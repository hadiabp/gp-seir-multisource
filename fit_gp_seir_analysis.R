#!/usr/bin/env Rscript
##############################################################################
# fit_gp_seir_analysis.R — GP-SEIR Complete Analysis
#
# Outputs:
#   CSV:  gp_seir_results.csv, ablation_results.csv, ablation_params.csv,
#         weight_optimization.csv, wdata_sensitivity.csv,
#         multistart_params.csv, wave_params.csv,
#         tikz_*.csv (8 files for pgfplots)
#   PDF:  gp_seir_plots.pdf (10 pages)
#
# Input:  gp_input_final_v2.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(fda); library(minpack.lm); library(zoo)
})

cat("============================================================\n")
cat("GP-SEIR: Complete Analysis\n")
cat("============================================================\n\n")

# ══════════════════════════════════════════════════════════════════
# SETTINGS
# ══════════════════════════════════════════════════════════════════

# ── RUN MODE ──────────────────────────────────────────────────────
# "full"   = run everything (Steps 1-7)
# "manual" = run single fit with weights below, then export
# "ablation_only" = Steps 1-2 only (main fit + ablation)
RUN_MODE <- "full"

# ── MANUAL WEIGHTS (only used when RUN_MODE = "manual") ──────────
W_CASES  <- 1.0
W_HOSP   <- 1.0
W_ICU    <- 1.0
W_RADAR  <- 1.0
# ──────────────────────────────────────────────────────────────────

# ── RANKING CRITERION ─────────────────────────────────────────────
# "prev_corr" = rank by prevalence correlation (fully independent)
# "rt_corr"   = rank by Rt correlation
SORT_CRITERION <- "prev_corr"
# ──────────────────────────────────────────────────────────────────

DATA_FILE  <- "gp_input_final_v2.csv"
DATE_START <- as.Date("2020-03-01")
DATE_END   <- as.Date("2021-03-01")

N_POP   <- 17400000
sigma   <- 1 / 5.5
gamma   <- 1 / 9.5

STATE_KNOT_DAYS <- 14;  STATE_NORDER <- 4
BETA_KNOT_DAYS  <- 28;  BETA_NORDER  <- 4

ETA        <- 1e6
KAPPA_BETA <- 1e3
KAPPA_BMAG <- 1e4
BETA_MAX   <- 0.6
KAPPA_MASS <- 1e6

N_STARTS   <- 5
MAX_INNER  <- 300
MAX_OUTER  <- 40
STEP0      <- 1.0

bounds <- list(
  rho    = c(0.05, 0.60),
  pH     = c(0.002, 0.05),
  pICU   = c(0.0003, 0.02),
  alphaR = c(0.10, 1.00)
)

# ══════════════════════════════════════════════════════════════════
# DATA
# ══════════════════════════════════════════════════════════════════

cat("--- Loading data ---\n")
df <- fread(DATA_FILE)[, date := as.Date(date)]
df <- df[date >= DATE_START & date <= DATE_END][order(date)]
n <- nrow(df)
df[, t := as.numeric(date - DATE_START)]
tv <- df$t;  rng <- range(tv);  dt <- 1;  dates <- df$date

cases_raw  <- as.numeric(df$cases_raw)
hosp_raw   <- as.numeric(df$hosp_nice)
icu_raw    <- as.numeric(df$icu_daily)
radar_frac <- as.numeric(df$radar_I_frac_sm7)
rt_rivm    <- as.numeric(df$rt_rivm)
prev_rivm  <- as.numeric(df$prev_json_avg)
tests_tot  <- as.numeric(df$tests_total)
tests_pos  <- as.numeric(df$tests_positive)
positivity <- as.numeric(df$positivity_rate)

y_cases <- log1p(pmax(0, cases_raw))
y_hosp  <- log1p(pmax(0, hosp_raw))
y_icu   <- log1p(pmax(0, icu_raw))
y_radar <- radar_frac

wsd <- function(y) {
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  1 / s^2
}
wb_cases <- wsd(y_cases)
wb_hosp  <- wsd(y_hosp)
wb_icu   <- wsd(y_icu)
wb_radar <- wsd(y_radar)

cat(sprintf("  %d days, ODE/cases ratio: %.0fx\n", n, sqrt(ETA * dt) / sqrt(wb_cases)))

# ══════════════════════════════════════════════════════════════════
# B-SPLINE BASES
# ══════════════════════════════════════════════════════════════════

mkb <- function(kd, no) {
  kn <- unique(c(rng[1], seq(rng[1], rng[2], by = kd), rng[2]))
  create.bspline.basis(rng, length(kn) + no - 2, no, kn)
}
Bs  <- mkb(STATE_KNOT_DAYS, STATE_NORDER)
Bb  <- mkb(BETA_KNOT_DAYS, BETA_NORDER)
nbs <- Bs$nbasis;  nbb <- Bb$nbasis
Ph  <- eval.basis(tv, Bs, 0)
dPh <- eval.basis(tv, Bs, 1)
Pb  <- eval.basis(tv, Bb, 0)
d2Pb <- eval.basis(tv, Bb, 2)

cat(sprintf("  State basis: %d, Beta basis: %d\n\n", nbs, nbb))

# ══════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════

to_u   <- function(x, lo, hi) qlogis((x - lo) / (hi - lo))
from_u <- function(u, lo, hi) lo + (hi - lo) * plogis(u)

seir_rhs <- function(Z, beta) {
  cbind(-beta * exp(Z[, 3]),
         beta * exp(Z[, 1] + Z[, 3] - Z[, 2]) - sigma,
         sigma * exp(Z[, 2] - Z[, 3]) - gamma,
         gamma * exp(Z[, 3] - Z[, 4]))
}

# ══════════════════════════════════════════════════════════════════
# INITIALIZATION
# ══════════════════════════════════════════════════════════════════

c0 <- ifelse(is.finite(cases_raw), cases_raw, 0)
E0 <- pmax(1e-10, c0 / (N_POP * 0.20 * sigma))
I0 <- as.numeric(rollmean(E0, 7, fill = median(E0), align = "right"))
R0 <- pmin(0.8, cumsum(pmax(c0, 0)) / N_POP)
S0 <- pmax(1e-6, 1 - E0 - I0 - R0)
m0 <- S0 + E0 + I0 + R0
S0 <- S0/m0; E0 <- E0/m0; I0 <- I0/m0; R0 <- R0/m0

Ci <- smooth.basis(tv, cbind(log(S0), log(E0), log(I0), log(R0)),
                   fdPar(Bs, int2Lfd(2), 1e-4))$fd$coefs
b0 <- rep(0.2, n); ok <- is.finite(rt_rivm)
if (any(ok)) b0[ok] <- pmax(0.05, pmin(0.5, rt_rivm[ok] * gamma))
ai <- as.numeric(smooth.basis(tv, log(b0),
                               fdPar(Bb, int2Lfd(2), 1e-2))$fd$coefs)

thi <- c(ai,
         to_u(0.20,   bounds$rho[1],    bounds$rho[2]),
         to_u(0.012,  bounds$pH[1],     bounds$pH[2]),
         to_u(0.0025, bounds$pICU[1],   bounds$pICU[2]),
         to_u(0.40,   bounds$alphaR[1], bounds$alphaR[2]))

# ══════════════════════════════════════════════════════════════════
# INNER SOLVE
# ══════════════════════════════════════════════════════════════════

solve_inner <- function(Cs, th, wc, wh, wi, wr) {
  av   <- th[1:nbb]
  rho  <- from_u(th[nbb + 1], bounds$rho[1],    bounds$rho[2])
  pH   <- from_u(th[nbb + 2], bounds$pH[1],     bounds$pH[2])
  pICU <- from_u(th[nbb + 3], bounds$pICU[1],   bounds$pICU[2])
  aR   <- from_u(th[nbb + 4], bounds$alphaR[1], bounds$alphaR[2])
  g    <- as.numeric(Pb %*% av)
  beta <- exp(g)

  fn <- function(cf) {
    C  <- matrix(cf, nrow = nbs, ncol = 4)
    Z  <- Ph %*% C;  dZ <- dPh %*% C
    S  <- exp(Z[, 1]); E <- exp(Z[, 2])
    I  <- exp(Z[, 3]); R <- exp(Z[, 4])

    mc <- log1p(pmax(0, N_POP * rho * sigma * E))
    mh <- log1p(pmax(0, N_POP * pH * gamma * I))
    mi <- log1p(pmax(0, N_POP * pICU * gamma * I))
    mr <- aR * I

    rc <- sqrt(wb_cases * wc) * ifelse(is.finite(y_cases), y_cases - mc, 0)
    rh <- sqrt(wb_hosp * wh)  * ifelse(is.finite(y_hosp),  y_hosp - mh, 0)
    ri <- sqrt(wb_icu * wi)   * ifelse(is.finite(y_icu),   y_icu - mi, 0)
    rr <- sqrt(wb_radar * wr) * ifelse(is.finite(y_radar), y_radar - mr, 0)

    rp <- sqrt(ETA * dt) * as.numeric(dZ - seir_rhs(Z, beta))
    rm <- sqrt(KAPPA_MASS) * (S + E + I + R - 1)
    rb <- sqrt(KAPPA_BMAG * dt) * pmax(beta - BETA_MAX, 0)
    rs <- sqrt(KAPPA_BETA * dt) * as.numeric(d2Pb %*% av)

    c(rc, rh, ri, rr, rp, rm, rb, rs)
  }

  r <- nls.lm(par = as.numeric(Cs), fn = fn,
              control = nls.lm.control(maxiter = MAX_INNER,
                                        ftol = 1e-10, ptol = 1e-10))
  list(C = matrix(r$par, nrow = nbs, ncol = 4), J = sum(r$fvec^2))
}

# ══════════════════════════════════════════════════════════════════
# CASCADING (OUTER PROBLEM)
# ══════════════════════════════════════════════════════════════════

cascade <- function(th0, C0, wc, wh, wi, wr) {
  th <- th0; Cc <- C0; bJ <- Inf
  for (oi in 1:MAX_OUTER) {
    inn <- solve_inner(Cc, th, wc, wh, wi, wr)
    Cc <- inn$C; Jc <- inn$J
    if (Jc < bJ) bJ <- Jc

    eps <- 1e-4
    gr  <- numeric(length(th))
    for (j in seq_along(th)) {
      tp    <- th; tp[j] <- tp[j] + eps
      gr[j] <- (solve_inner(Cc, tp, wc, wh, wi, wr)$J - Jc) / eps
    }
    gn <- sqrt(sum(gr^2))

    if (gn > 1e-12) {
      d  <- -gr / gn; st <- STEP0; found <- FALSE
      for (ls in 1:10) {
        tt <- th + st * d
        it <- solve_inner(Cc, tt, wc, wh, wi, wr)
        if (it$J < Jc - 1e-4 * st * gn) {
          th <- tt; Cc <- it$C; found <- TRUE; break
        }
        st <- st * 0.5
      }
      if (!found) th <- th + 1e-4 * d
    }
    if (gn < 1e-5 * max(1, sqrt(sum(th^2)))) break
  }
  list(theta = th, C = Cc, J = bJ)
}

# ══════════════════════════════════════════════════════════════════
# FULL DIAGNOSTICS (expanded: includes residuals, ODE per compartment)
# ══════════════════════════════════════════════════════════════════

extract_full <- function(res) {
  th <- res$theta; C <- res$C; av <- th[1:nbb]
  rho  <- from_u(th[nbb + 1], bounds$rho[1],    bounds$rho[2])
  pH   <- from_u(th[nbb + 2], bounds$pH[1],     bounds$pH[2])
  pICU <- from_u(th[nbb + 3], bounds$pICU[1],   bounds$pICU[2])
  aR   <- from_u(th[nbb + 4], bounds$alphaR[1], bounds$alphaR[2])

  Z  <- Ph %*% C; dZ <- dPh %*% C
  S  <- exp(Z[, 1]); E <- exp(Z[, 2])
  I  <- exp(Z[, 3]); R <- exp(Z[, 4])
  beta <- exp(as.numeric(Pb %*% av))
  Rt   <- beta * S / gamma
  mass <- S + E + I + R

  # ODE residuals per compartment
  ode_mat    <- dZ - seir_rhs(Z, beta)
  ode_rmse   <- sqrt(mean(ode_mat^2))
  ode_rmse_S <- sqrt(mean(ode_mat[, 1]^2))
  ode_rmse_E <- sqrt(mean(ode_mat[, 2]^2))
  ode_rmse_I <- sqrt(mean(ode_mat[, 3]^2))
  ode_rmse_R <- sqrt(mean(ode_mat[, 4]^2))

  # Model predictions on original scale
  mc <- N_POP * rho * sigma * E
  mh <- N_POP * pH * gamma * I
  mi <- N_POP * pICU * gamma * I
  mr <- aR * I

  # Data residuals on log(1+y) scale
  resid_cases <- ifelse(is.finite(y_cases), y_cases - log1p(pmax(0, mc)), NA)
  resid_hosp  <- ifelse(is.finite(y_hosp),  y_hosp  - log1p(pmax(0, mh)), NA)
  resid_icu   <- ifelse(is.finite(y_icu),   y_icu   - log1p(pmax(0, mi)), NA)
  resid_radar <- ifelse(is.finite(y_radar), y_radar  - mr, NA)

  # Data RMSE on original scale
  ok_c <- is.finite(cases_raw)
  ok_h <- is.finite(hosp_raw) & hosp_raw > 0
  ok_i <- is.finite(icu_raw)
  ok_r <- is.finite(radar_frac)
  cases_rmse <- sqrt(mean((cases_raw[ok_c] - mc[ok_c])^2))
  hosp_rmse  <- if (any(ok_h)) sqrt(mean((hosp_raw[ok_h] - mh[ok_h])^2)) else NA
  icu_rmse   <- if (any(ok_i)) sqrt(mean((icu_raw[ok_i] - mi[ok_i])^2)) else NA

  # Rt validation
  ok_rt <- is.finite(rt_rivm)
  rt_corr <- cor(rt_rivm[ok_rt], Rt[ok_rt])
  rt_spearman <- cor(rt_rivm[ok_rt], Rt[ok_rt], method = "spearman")
  rt_rmse <- sqrt(mean((rt_rivm[ok_rt] - Rt[ok_rt])^2))

  # Prevalence validation
  ok_prev <- is.finite(prev_rivm)
  I_count <- N_POP * I
  prev_corr <- if (any(ok_prev)) cor(prev_rivm[ok_prev], I_count[ok_prev]) else NA
  prev_spearman <- if (any(ok_prev)) cor(prev_rivm[ok_prev], I_count[ok_prev],
                                          method = "spearman") else NA

  # Effective time-varying detection rate: rho_eff(t) = cases / (N*sigma*E)
  # This shows how constant-rho assumption breaks down between waves
  denom <- N_POP * sigma * E
  rho_eff <- ifelse(is.finite(cases_raw) & denom > 1, cases_raw / denom, NA)

  # Per-period Rt
  cutoff <- as.Date("2020-06-12")
  bf <- dates < cutoff & ok_rt
  af <- dates >= cutoff & ok_rt
  rt_corr_pre  <- if (sum(bf) > 5) cor(rt_rivm[bf], Rt[bf]) else NA
  rt_corr_post <- if (sum(af) > 5) cor(rt_rivm[af], Rt[af]) else NA

  list(
    # Trajectories
    S = S, E = E, I = I, R = R, beta = beta, Rt = Rt, mass = mass,
    # Model predictions
    mu_cases = mc, mu_hosp = mh, mu_icu = mi, mu_radar = mr,
    # Parameters
    rho = rho, pH = pH, pICU = pICU, alphaR = aR,
    # ODE residuals (overall and per compartment)
    ode_mat = ode_mat, ode_rmse = ode_rmse,
    ode_rmse_S = ode_rmse_S, ode_rmse_E = ode_rmse_E,
    ode_rmse_I = ode_rmse_I, ode_rmse_R = ode_rmse_R,
    # Data residuals on log scale
    resid_cases = resid_cases, resid_hosp = resid_hosp,
    resid_icu = resid_icu, resid_radar = resid_radar,
    # Data RMSE on original scale
    cases_rmse = cases_rmse, hosp_rmse = hosp_rmse, icu_rmse = icu_rmse,
    # Validation metrics
    rt_corr = rt_corr, rt_rmse = rt_rmse, rt_spearman = rt_spearman,
    prev_corr = prev_corr, prev_spearman = prev_spearman,
    rt_corr_pre = rt_corr_pre, rt_corr_post = rt_corr_post,
    # Effective detection rate (time-varying)
    rho_eff = rho_eff,
    # Summary
    beta_range = c(min(beta), max(beta)),
    Rt_range   = c(min(Rt), max(Rt)),
    R_end = R[n],
    mass_err = max(abs(mass - 1))
  )
}

# ══════════════════════════════════════════════════════════════════
# MULTI-START WRAPPER (saves ALL starts, not just best)
# ══════════════════════════════════════════════════════════════════

run_fit_full <- function(wc, wh, wi, wr, label = "",
                         warm_th = NULL, warm_C = NULL) {
  cat(sprintf("  %-35s ", label))
  all_starts <- list()
  best <- NULL; bJ <- Inf
  th0 <- if (!is.null(warm_th)) warm_th else thi
  C0  <- if (!is.null(warm_C))  warm_C  else Ci

  for (s in 1:N_STARTS) {
    ts <- th0; Cs <- C0
    if (s > 1) {
      ts[1:nbb] <- ts[1:nbb] + rnorm(nbb, 0, 0.15)
      for (k in (nbb + 1):(nbb + 4))
        ts[k] <- ts[k] + rnorm(1, 0, 0.3)
    }
    r <- tryCatch(cascade(ts, Cs, wc, wh, wi, wr), error = function(e) NULL)
    if (!is.null(r)) {
      d <- extract_full(r)
      all_starts[[s]] <- list(
        J = r$J, rho = d$rho, pH = d$pH, pICU = d$pICU, alphaR = d$alphaR,
        rt_corr = d$rt_corr, rt_spearman = d$rt_spearman,
        prev_corr = d$prev_corr, prev_spearman = d$prev_spearman
      )
      if (r$J < bJ) { bJ <- r$J; best <- r }
    }
  }

  if (is.null(best)) { cat("FAILED\n"); return(NULL) }
  d <- extract_full(best)
  cat(sprintf("Rt=%.3f prev=%.3f cases=%.0f\n",
              d$rt_corr, d$prev_corr, d$cases_rmse))
  c(d, list(theta = best$theta, C_mat = best$C, J = best$J,
            all_starts = all_starts))
}

# Short wrapper (no start tracking, for ablation/optimization)
run_fit <- function(wc, wh, wi, wr, label = "",
                    warm_th = NULL, warm_C = NULL) {
  cat(sprintf("  %-35s ", label))
  best <- NULL; bJ <- Inf
  th0 <- if (!is.null(warm_th)) warm_th else thi
  C0  <- if (!is.null(warm_C))  warm_C  else Ci

  for (s in 1:N_STARTS) {
    ts <- th0; Cs <- C0
    if (s > 1) {
      ts[1:nbb] <- ts[1:nbb] + rnorm(nbb, 0, 0.15)
      for (k in (nbb + 1):(nbb + 4))
        ts[k] <- ts[k] + rnorm(1, 0, 0.3)
    }
    r <- tryCatch(cascade(ts, Cs, wc, wh, wi, wr), error = function(e) NULL)
    if (!is.null(r) && r$J < bJ) { bJ <- r$J; best <- r }
  }
  if (is.null(best)) { cat("FAILED\n"); return(NULL) }
  d <- extract_full(best)
  cat(sprintf("Rt=%.3f prev=%.3f cases=%.0f\n",
              d$rt_corr, d$prev_corr, d$cases_rmse))
  c(d, list(theta = best$theta, C_mat = best$C, J = best$J))
}


# ══════════════════════════════════════════════════════════════════
# STEP 1: MAIN FIT (with multi-start parameter tracking)
# ══════════════════════════════════════════════════════════════════

cat("\n=== STEP 1: Main fit (5 starts) ===\n")
set.seed(42)

if (RUN_MODE == "manual") {
  cat(sprintf("  Manual mode: w=(%s,%s,%s,%s)\n", W_CASES, W_HOSP, W_ICU, W_RADAR))
  main <- run_fit_full(W_CASES, W_HOSP, W_ICU, W_RADAR,
                       sprintf("manual (%.1f,%.1f,%.1f,%.1f)",
                               W_CASES, W_HOSP, W_ICU, W_RADAR))
} else {
  main <- run_fit_full(1, 1, 1, 1, "all sources (1,1,1,1)")
}

# [NEW 1] Save multi-start parameter variation
ms_df <- data.frame()
for (i in seq_along(main$all_starts)) {
  s <- main$all_starts[[i]]
  if (!is.null(s)) {
    ms_df <- rbind(ms_df, data.frame(
      start = i, J = s$J, rho = s$rho, pH = s$pH,
      pICU = s$pICU, alphaR = s$alphaR,
      rt_corr = s$rt_corr, prev_corr = s$prev_corr
    ))
  }
}
fwrite(ms_df, "multistart_params.csv")
cat(sprintf("\n  Parameters (best): rho=%.3f pH=%.4f pICU=%.4f alphaR=%.3f\n",
            main$rho, main$pH, main$pICU, main$alphaR))
cat(sprintf("  Validation: Rt=%.3f (pre=%.3f post=%.3f) Prev=%.3f\n",
            main$rt_corr, main$rt_corr_pre, main$rt_corr_post, main$prev_corr))
cat(sprintf("  Spearman:   Rt=%.3f  Prev=%.3f\n",
            main$rt_spearman, main$prev_spearman))
cat(sprintf("  ODE per compartment: S=%.4f E=%.4f I=%.4f R=%.4f\n",
            main$ode_rmse_S, main$ode_rmse_E, main$ode_rmse_I, main$ode_rmse_R))
# Detection rate summary
ok_rho <- is.finite(main$rho_eff)
if (any(ok_rho)) {
  w1 <- dates < as.Date("2020-06-01") & ok_rho
  w2 <- dates >= as.Date("2020-10-01") & ok_rho
  cat(sprintf("  Detection rate: wave1 median=%.3f, wave2 median=%.3f (fitted rho=%.3f)\n",
              median(main$rho_eff[w1], na.rm=TRUE),
              median(main$rho_eff[w2], na.rm=TRUE), main$rho))
}
if (nrow(ms_df) > 1) {
  cat(sprintf("  Multi-start rho range: [%.3f, %.3f]\n",
              min(ms_df$rho), max(ms_df$rho)))
  cat(sprintf("  Multi-start alphaR range: [%.3f, %.3f]\n",
              min(ms_df$alphaR), max(ms_df$alphaR)))
}

# ══════════════════════════════════════════════════════════════════
# STEP 2: SOURCE ABLATION (with parameter tracking per config)
# ══════════════════════════════════════════════════════════════════

abl <- data.frame()
abl_betas <- list()

if (RUN_MODE %in% c("full", "ablation_only")) {
cat("\n=== STEP 2: Source ablation ===\n")
set.seed(42)

cfgs <- list(
  list(1, 0, 0, 0, "cases only"),
  list(1, 1, 0, 0, "cases+hosp"),
  list(1, 1, 1, 0, "cases+hosp+icu"),
  list(1, 1, 1, 1, "cases+hosp+icu+radar"),
  list(0, 1, 1, 0, "hosp+icu"),
  list(0, 1, 1, 1, "hosp+icu+radar"),
  list(1, 0, 0, 1, "cases+radar"),
  list(0, 0, 0, 1, "radar only")
)

abl       <- data.frame()
abl_betas <- list()  # [NEW 5] store beta for each config

for (cfg in cfgs) {
  d <- run_fit(cfg[[1]], cfg[[2]], cfg[[3]], cfg[[4]], cfg[[5]],
               warm_th = main$theta, warm_C = main$C_mat)
  if (!is.null(d)) {
    abl <- rbind(abl, data.frame(
      sources = cfg[[5]],
      # [NEW 6] all estimated parameters
      rho = d$rho, pH = d$pH, pICU = d$pICU, alphaR = d$alphaR,
      # validation
      rt_corr = d$rt_corr, rt_spearman = d$rt_spearman, rt_rmse = d$rt_rmse,
      rt_corr_pre = d$rt_corr_pre, rt_corr_post = d$rt_corr_post,
      prev_corr = d$prev_corr, prev_spearman = d$prev_spearman,
      # fit quality
      cases_rmse = d$cases_rmse, hosp_rmse = d$hosp_rmse,
      ode_rmse = d$ode_rmse,
      # summary
      R_end = d$R_end, beta_min = d$beta_range[1], beta_max = d$beta_range[2]
    ))
    abl_betas[[cfg[[5]]]] <- d$beta  # [NEW 5]
  }
}
abl <- abl[order(-abl[[SORT_CRITERION]]), ]
fwrite(abl, "ablation_results.csv")
# [NEW 6] separate file with just parameters for easy comparison
fwrite(abl[, c("sources","rho","pH","pICU","alphaR","R_end")], "ablation_params.csv")

cat("\n  Ablation table:\n")
cat(sprintf("  %-25s %6s %6s %6s %6s %6s %6s\n",
            "Sources", "Rt", "Rt_pre", "Rt_pst", "Prev", "RtRMSE", "R_end"))
for (i in 1:nrow(abl)) {
  with(abl[i, ], cat(sprintf("  %-25s %6.3f %6.3f %6.3f %6.3f %6.3f %5.1f%%\n",
    sources, rt_corr, rt_corr_pre, rt_corr_post, prev_corr, rt_rmse, R_end * 100)))
}

} # end if ablation

# ══════════════════════════════════════════════════════════════════
# STEP 3: WEIGHT OPTIMIZATION
# ══════════════════════════════════════════════════════════════════

wopt <- data.frame()

if (RUN_MODE == "full") {
cat("\n=== STEP 3: Weight optimization (all 4 sources) ===\n")
set.seed(42)

# Full search over all 4 weights
# Cases and RADAR: 4 levels (most important based on ablation)
# Hospital and ICU: 2 levels (low marginal value, but still explored)
wopt <- data.frame()
for (wc in c(0.1, 0.5, 1.0, 2.0)) {
  for (wh in c(0, 1.0)) {
    for (wi in c(0, 1.0)) {
      for (wr in c(0, 0.5, 1.0, 2.0)) {
        # Skip: must have at least one data source
        if (wc == 0 & wh == 0 & wi == 0 & wr == 0) next
        label <- sprintf("c=%.1f h=%.0f i=%.0f r=%.1f", wc, wh, wi, wr)
        d <- run_fit(wc, wh, wi, wr, label,
                     warm_th = main$theta, warm_C = main$C_mat)
        if (!is.null(d)) {
          wopt <- rbind(wopt, data.frame(
            w_cases = wc, w_hosp = wh, w_icu = wi, w_radar = wr,
            rt_corr = d$rt_corr, rt_spearman = d$rt_spearman, rt_rmse = d$rt_rmse,
            prev_corr = d$prev_corr, prev_spearman = d$prev_spearman,
            cases_rmse = d$cases_rmse, R_end = d$R_end,
            rho = d$rho, pH = d$pH, pICU = d$pICU, alphaR = d$alphaR,
            ode_rmse = d$ode_rmse
          ))
        }
      }
    }
  }
}
wopt <- wopt[order(-wopt[[SORT_CRITERION]]), ]
fwrite(wopt, "weight_optimization.csv")
cat(sprintf("\n  %d weight combinations explored\n", nrow(wopt)))
cat("  Top 5 by ", SORT_CRITERION, ":\n", sep = "")
print(head(wopt[, c("w_cases","w_hosp","w_icu","w_radar",
                     "rt_corr","prev_corr","rho","alphaR")], 5))
} # end if weight optimization

# ══════════════════════════════════════════════════════════════════
# STEP 4: PENALTY WEIGHT SENSITIVITY
# ══════════════════════════════════════════════════════════════════

wsens <- data.frame()
wsens_fits <- list()

if (RUN_MODE == "full") {
cat("\n=== STEP 4: W_DATA sensitivity ===\n")
set.seed(42)

# Define day vector here (needed for tikz_rt_wdata export below)
day <- as.numeric(dates - DATE_START)

W_DATA_GRID <- c(0.01, 1, 100, 1e4, 1e6, 1e8)
wsens      <- data.frame()
wsens_fits <- list()

for (i in seq_along(W_DATA_GRID)) {
  wd <- W_DATA_GRID[i]
  eff <- sqrt(ETA * dt) / sqrt(wb_cases * wd)
  cat(sprintf("  W_DATA=%.0e (ratio=%.0fx) ... ", wd, eff))
  r <- tryCatch(cascade(main$theta, main$C_mat, wd, wd, wd, wd),
                error = function(e) NULL)
  if (!is.null(r)) {
    d <- extract_full(r)
    wsens_fits[[i]] <- d
    cat(sprintf("Rt=%.3f prev=%.3f\n", d$rt_corr, d$prev_corr))
    wsens <- rbind(wsens, data.frame(
      W_DATA = wd, eff_ratio = eff,
      rt_corr = d$rt_corr, rt_spearman = d$rt_spearman, rt_rmse = d$rt_rmse,
      prev_corr = d$prev_corr, prev_spearman = d$prev_spearman,
      cases_rmse = d$cases_rmse, ode_rmse = d$ode_rmse,
      beta_min = d$beta_range[1], beta_max = d$beta_range[2],
      R_end = d$R_end
    ))
  } else { cat("FAILED\n") }
}
fwrite(wsens, "wdata_sensitivity.csv")

# Save Rt trajectories per W_DATA value (for overlay plot in LaTeX)
rt_wd_df <- data.table(day = day, rt_rivm = rt_rivm)
for (i in seq_along(wsens_fits)) {
  if (!is.null(wsens_fits[[i]])) {
    col_name <- sprintf("Rt_wd%.0e", W_DATA_GRID[i])
    rt_wd_df[[col_name]] <- wsens_fits[[i]]$Rt
  }
}
fwrite(rt_wd_df, "tikz_rt_wdata.csv")
} # end if W_DATA sensitivity

# ══════════════════════════════════════════════════════════════════
# STEP 5: PER-WAVE REFITTING
#
# Fits independent GP-SEIR models for each wave with their own
# B-spline bases, giving separate (rho, pH, pICU, alphaR) per wave.
# Also computes sliced diagnostics from the main fit for comparison.
# ══════════════════════════════════════════════════════════════════

# Self-contained sub-period fitting function
# Builds new bases, slices data, runs full parameter cascading
fit_subperiod <- function(date_start, date_end, label,
                          main_fit, full_dates, full_df) {

  cat(sprintf("  %-30s ", label))

  # ── Slice data ──
  idx   <- full_dates >= as.Date(date_start) & full_dates <= as.Date(date_end)
  n_sub <- sum(idx)
  if (n_sub < 45) { cat("too few days\n"); return(NULL) }

  dates_sub <- full_dates[idx]
  tv_sub    <- as.numeric(dates_sub - as.Date(date_start))
  rng_sub   <- range(tv_sub)

  cases_sub <- cases_raw[idx]
  hosp_sub  <- hosp_raw[idx]
  icu_sub   <- icu_raw[idx]
  radar_sub <- radar_frac[idx]
  rt_sub    <- rt_rivm[idx]
  prev_sub  <- prev_rivm[idx]

  yc_sub <- log1p(pmax(0, cases_sub))
  yh_sub <- log1p(pmax(0, hosp_sub))
  yi_sub <- log1p(pmax(0, icu_sub))
  yr_sub <- radar_sub

  wsd_l <- function(y) { s <- sd(y, na.rm=TRUE); if(!is.finite(s)||s<=0) s<-1; 1/s^2 }
  wbc <- wsd_l(yc_sub); wbh <- wsd_l(yh_sub)
  wbi <- wsd_l(yi_sub); wbr <- wsd_l(yr_sub)

  # ── Build sub-period B-spline bases ──
  mkb_sub <- function(kd, no) {
    kn <- unique(c(rng_sub[1], seq(rng_sub[1], rng_sub[2], by=kd), rng_sub[2]))
    create.bspline.basis(rng_sub, length(kn)+no-2, no, kn)
  }
  Bs_sub  <- mkb_sub(STATE_KNOT_DAYS, STATE_NORDER)
  Bb_sub  <- mkb_sub(BETA_KNOT_DAYS, BETA_NORDER)
  nbs_sub <- Bs_sub$nbasis
  nbb_sub <- Bb_sub$nbasis
  Ph_sub  <- eval.basis(tv_sub, Bs_sub, 0)
  dPh_sub <- eval.basis(tv_sub, Bs_sub, 1)
  Pb_sub  <- eval.basis(tv_sub, Bb_sub, 0)
  d2Pb_sub <- eval.basis(tv_sub, Bb_sub, 2)

  cat(sprintf("(%d days, %d+%d bases) ", n_sub, nbs_sub, nbb_sub))

  # ── Initialize from main fit (interpolated to sub-period) ──
  S_init <- main_fit$S[idx]; E_init <- main_fit$E[idx]
  I_init <- main_fit$I[idx]; R_init <- main_fit$R[idx]
  b_init <- main_fit$beta[idx]

  Ci_sub <- tryCatch(
    smooth.basis(tv_sub, cbind(log(S_init), log(E_init), log(I_init), log(R_init)),
                 fdPar(Bs_sub, int2Lfd(2), 1e-4))$fd$coefs,
    error = function(e) matrix(0, nbs_sub, 4)
  )
  ai_sub <- tryCatch(
    as.numeric(smooth.basis(tv_sub, log(b_init),
                             fdPar(Bb_sub, int2Lfd(2), 1e-2))$fd$coefs),
    error = function(e) rep(-1.5, nbb_sub)
  )

  thi_sub <- c(ai_sub,
               to_u(main_fit$rho,    bounds$rho[1],    bounds$rho[2]),
               to_u(main_fit$pH,     bounds$pH[1],     bounds$pH[2]),
               to_u(main_fit$pICU,   bounds$pICU[1],   bounds$pICU[2]),
               to_u(main_fit$alphaR, bounds$alphaR[1], bounds$alphaR[2]))

  # ── Local inner solve (closure over sub-period variables) ──
  solve_inner_sub <- function(Cs, th, wc, wh, wi, wr) {
    av   <- th[1:nbb_sub]
    rho  <- from_u(th[nbb_sub+1], bounds$rho[1], bounds$rho[2])
    pH   <- from_u(th[nbb_sub+2], bounds$pH[1],  bounds$pH[2])
    pICU <- from_u(th[nbb_sub+3], bounds$pICU[1], bounds$pICU[2])
    aR   <- from_u(th[nbb_sub+4], bounds$alphaR[1], bounds$alphaR[2])
    g    <- as.numeric(Pb_sub %*% av); beta <- exp(g)

    fn <- function(cf) {
      C  <- matrix(cf, nrow=nbs_sub, ncol=4)
      Z  <- Ph_sub %*% C; dZ <- dPh_sub %*% C
      S  <- exp(Z[,1]); E <- exp(Z[,2]); I <- exp(Z[,3]); R <- exp(Z[,4])
      mc <- log1p(pmax(0, N_POP*rho*sigma*E))
      mh <- log1p(pmax(0, N_POP*pH*gamma*I))
      mi <- log1p(pmax(0, N_POP*pICU*gamma*I))
      mr <- aR*I
      rc <- sqrt(wbc*wc)*ifelse(is.finite(yc_sub), yc_sub-mc, 0)
      rh <- sqrt(wbh*wh)*ifelse(is.finite(yh_sub), yh_sub-mh, 0)
      ri <- sqrt(wbi*wi)*ifelse(is.finite(yi_sub), yi_sub-mi, 0)
      rr <- sqrt(wbr*wr)*ifelse(is.finite(yr_sub), yr_sub-mr, 0)
      rp <- sqrt(ETA*dt)*as.numeric(dZ - seir_rhs(Z, beta))
      rm <- sqrt(KAPPA_MASS)*(S+E+I+R-1)
      rb <- sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX, 0)
      rs <- sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb_sub %*% av)
      c(rc, rh, ri, rr, rp, rm, rb, rs)
    }
    r <- nls.lm(par=as.numeric(Cs), fn=fn,
                control=nls.lm.control(maxiter=MAX_INNER, ftol=1e-10, ptol=1e-10))
    list(C=matrix(r$par, nrow=nbs_sub, ncol=4), J=sum(r$fvec^2))
  }

  # ── Local cascading ──
  cascade_sub <- function(th0, C0, wc, wh, wi, wr) {
    th <- th0; Cc <- C0; bJ <- Inf
    for (oi in 1:MAX_OUTER) {
      inn <- solve_inner_sub(Cc, th, wc, wh, wi, wr)
      Cc <- inn$C; Jc <- inn$J; if(Jc < bJ) bJ <- Jc
      eps <- 1e-4; gr <- numeric(length(th))
      for (j in seq_along(th)) {
        tp <- th; tp[j] <- tp[j]+eps
        gr[j] <- (solve_inner_sub(Cc, tp, wc, wh, wi, wr)$J - Jc)/eps
      }
      gn <- sqrt(sum(gr^2))
      if (gn > 1e-12) {
        d <- -gr/gn; st <- STEP0; found <- FALSE
        for (ls in 1:10) {
          tt <- th+st*d; it <- solve_inner_sub(Cc, tt, wc, wh, wi, wr)
          if (it$J < Jc-1e-4*st*gn) {th<-tt; Cc<-it$C; found<-TRUE; break}
          st <- st*0.5
        }
        if (!found) th <- th+1e-4*d
      }
      if (gn < 1e-5*max(1, sqrt(sum(th^2)))) break
    }
    list(theta=th, C=Cc, J=bJ)
  }

  # ── Run multi-start ──
  best <- NULL; bJ <- Inf
  for (s in 1:N_STARTS) {
    ts <- thi_sub; Cs <- Ci_sub
    if (s > 1) {
      ts[1:nbb_sub] <- ts[1:nbb_sub] + rnorm(nbb_sub, 0, 0.15)
      for (k in (nbb_sub+1):(nbb_sub+4)) ts[k] <- ts[k] + rnorm(1, 0, 0.3)
    }
    r <- tryCatch(cascade_sub(ts, Cs, 1, 1, 1, 1), error=function(e) NULL)
    if (!is.null(r) && r$J < bJ) { bJ <- r$J; best <- r }
  }

  if (is.null(best)) { cat("FAILED\n"); return(NULL) }

  # ── Extract results ──
  th <- best$theta; C <- best$C; av <- th[1:nbb_sub]
  rho  <- from_u(th[nbb_sub+1], bounds$rho[1], bounds$rho[2])
  pH   <- from_u(th[nbb_sub+2], bounds$pH[1],  bounds$pH[2])
  pICU <- from_u(th[nbb_sub+3], bounds$pICU[1], bounds$pICU[2])
  aR   <- from_u(th[nbb_sub+4], bounds$alphaR[1], bounds$alphaR[2])

  Z <- Ph_sub %*% C; S <- exp(Z[,1]); E <- exp(Z[,2])
  I <- exp(Z[,3]); R <- exp(Z[,4])
  beta <- exp(as.numeric(Pb_sub %*% av))
  Rt   <- beta * S / gamma
  mc   <- N_POP * rho * sigma * E

  ok_rt <- is.finite(rt_sub)
  ok_prev <- is.finite(prev_sub)
  ok_c <- is.finite(cases_sub)
  rt_corr   <- if(sum(ok_rt)>5) cor(rt_sub[ok_rt], Rt[ok_rt]) else NA
  prev_corr <- if(sum(ok_prev)>5) cor(prev_sub[ok_prev], N_POP*I[ok_prev]) else NA
  cases_rmse <- if(sum(ok_c)>5) sqrt(mean((cases_sub[ok_c]-mc[ok_c])^2)) else NA

  detection <- if(sum(ok_c & ok_prev)>5)
    median(cases_sub[ok_c & ok_prev] / prev_sub[ok_c & ok_prev], na.rm=TRUE) else NA

  cat(sprintf("rho=%.3f pH=%.4f pICU=%.4f aR=%.3f Rt=%.3f prev=%.3f\n",
              rho, pH, pICU, aR, rt_corr, prev_corr))

  data.frame(
    period = label, n_days = n_sub,
    rho = rho, pH = pH, pICU = pICU, alphaR = aR,
    rt_corr = rt_corr, prev_corr = prev_corr,
    cases_rmse = cases_rmse, detection_rate = detection,
    beta_mean = mean(beta), beta_min = min(beta), beta_max = max(beta),
    I_mean = mean(I), R_end = R[n_sub],
    mass_err = max(abs(S+E+I+R-1))
  )
}

wave_params <- data.frame()
wave_refit  <- data.frame()

if (RUN_MODE %in% c("full", "ablation_only")) {
cat("\n=== STEP 5: Per-wave analysis ===\n")

periods <- list(
  list(start = "2020-03-01", end = "2020-06-01", label = "wave1 (Mar-Jun)"),
  list(start = "2020-06-01", end = "2020-10-01", label = "inter-wave (Jun-Oct)"),
  list(start = "2020-10-01", end = "2021-03-01", label = "wave2 (Oct-Mar)")
)

# ── Part A: Sliced diagnostics from main fit (fast, no refitting) ──
cat("  --- Sliced diagnostics from main fit ---\n")
for (p in periods) {
  idx <- dates >= as.Date(p$start) & dates <= as.Date(p$end)
  n_p <- sum(idx)
  if (n_p < 30) next

  cat(sprintf("  %-30s (%d days) ... ", p$label, n_p))

  ok_rt_p   <- is.finite(rt_rivm) & idx
  ok_prev_p <- is.finite(prev_rivm) & idx
  ok_c_p    <- is.finite(cases_raw) & idx

  rt_corr_p <- if(sum(ok_rt_p)>5) cor(rt_rivm[ok_rt_p], main$Rt[ok_rt_p]) else NA
  prev_corr_p <- if(sum(ok_prev_p)>5) cor(prev_rivm[ok_prev_p], N_POP*main$I[ok_prev_p]) else NA
  cases_rmse_p <- if(sum(ok_c_p)>5) sqrt(mean((cases_raw[ok_c_p]-main$mu_cases[ok_c_p])^2)) else NA
  detection_p <- if(sum(ok_c_p & ok_prev_p)>5)
    median(cases_raw[ok_c_p & ok_prev_p] / prev_rivm[ok_c_p & ok_prev_p], na.rm=TRUE) else NA

  cat(sprintf("Rt=%.3f prev=%.3f detect=%.3f\n", rt_corr_p, prev_corr_p, detection_p))

  wave_params <- rbind(wave_params, data.frame(
    period = p$label, n_days = n_p, type = "sliced",
    rho = main$rho, pH = main$pH, pICU = main$pICU, alphaR = main$alphaR,
    rt_corr = rt_corr_p, prev_corr = prev_corr_p,
    cases_rmse = cases_rmse_p, detection_rate = detection_p,
    beta_mean = mean(main$beta[idx]), I_mean = mean(main$I[idx]),
    R_end = main$R[max(which(idx))]
  ))
}

# ── Part B: Independent per-wave refitting (slow, new parameters) ──
cat("  --- Independent per-wave refitting ---\n")
set.seed(42)
for (p in periods) {
  result <- tryCatch(
    fit_subperiod(p$start, p$end, p$label, main, dates, df),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
  )
  if (!is.null(result)) {
    wave_refit <- rbind(wave_refit, result)
  }
}

fwrite(wave_params, "wave_params.csv")
if (nrow(wave_refit) > 0) {
  fwrite(wave_refit, "wave_refit_params.csv")
  cat("\n  Per-wave refitting results:\n")
  cat(sprintf("  %-25s %6s %6s %6s %6s %6s %6s\n",
              "Period", "rho", "pH", "pICU", "alphaR", "Rt", "Prev"))
  for (i in 1:nrow(wave_refit)) {
    with(wave_refit[i,], cat(sprintf("  %-25s %6.3f %6.4f %6.4f %6.3f %6.3f %6.3f\n",
      period, rho, pH, pICU, alphaR, rt_corr, prev_corr)))
  }
}
} # end if per-period

# ══════════════════════════════════════════════════════════════════
# STEP 6: EXPORT TIKZ DATA
# ══════════════════════════════════════════════════════════════════

cat("\n=== STEP 6: TikZ data ===\n")

# Ensure day is defined (Step 4 may have been skipped)
day <- as.numeric(dates - DATE_START)

# Main trajectories
fwrite(data.table(
  day = day, S = main$S, E = main$E, I = main$I, R = main$R,
  beta = main$beta, Rt_model = main$Rt, rt_rivm = rt_rivm
), "tikz_seir.csv")

# Fit panels
fwrite(data.table(day = day, cases_raw = cases_raw,
                  mu_cases = main$mu_cases), "tikz_cases.csv")
fwrite(data.table(day = day, hosp_raw = hosp_raw,
                  mu_hosp = main$mu_hosp), "tikz_hosp.csv")
fwrite(data.table(day = day, icu_raw = icu_raw,
                  mu_icu = main$mu_icu), "tikz_icu.csv")
fwrite(data.table(day = day, radar_frac = radar_frac,
                  mu_radar = main$mu_radar), "tikz_radar.csv")

# ACF of case residuals (pre-computed for pgfplots)
rc_clean <- main$resid_cases[is.finite(main$resid_cases)]
acf_obj  <- acf(rc_clean, lag.max = 30, plot = FALSE)
fwrite(data.table(lag = as.numeric(acf_obj$lag),
                  acf = as.numeric(acf_obj$acf)), "tikz_acf.csv")

# QQ data: sorted residuals vs theoretical normal quantiles
rc_sorted <- sort(rc_clean)
n_qq      <- length(rc_sorted)
qq_theor  <- qnorm(ppoints(n_qq))
fwrite(data.table(theoretical = qq_theor,
                  sample = rc_sorted), "tikz_qq.csv")

# Validation
fwrite(data.table(day = day, rt_rivm = rt_rivm,
                  rt_model = main$Rt), "tikz_rt.csv")
fwrite(data.table(day = day, prev_rivm = prev_rivm,
                  I_model = N_POP * main$I), "tikz_prev.csv")

# Ablation (only if run)
if (nrow(abl) > 0) {
  fwrite(abl[, c("sources", "rt_corr", "rt_rmse", "prev_corr")],
         "tikz_ablation.csv")
}

# [NEW 2] Residuals
fwrite(data.table(
  day = day,
  resid_cases = main$resid_cases, resid_hosp = main$resid_hosp,
  resid_icu = main$resid_icu, resid_radar = main$resid_radar,
  ode_S = main$ode_mat[, 1], ode_E = main$ode_mat[, 2],
  ode_I = main$ode_mat[, 3], ode_R = main$ode_mat[, 4]
), "tikz_residuals.csv")

# [NEW 5] Beta across ablation configs (only if run)
if (length(abl_betas) > 0) {
  beta_df <- data.table(day = day)
  for (nm in names(abl_betas)) {
    col <- gsub("[^a-zA-Z]", "_", nm)
    beta_df[[col]] <- abl_betas[[nm]]
  }
  fwrite(beta_df, "tikz_beta_ablation.csv")
}

# [NEW 7] Testing data
fwrite(data.table(
  day = day, tests_total = tests_tot, tests_positive = tests_pos,
  positivity = positivity
), "tikz_testing.csv")

# [NEW 8] Cumulative infections
fwrite(data.table(day = day, R_cumulative = main$R,
                  R_pct = main$R * 100), "tikz_cumulative.csv")

# Effective detection rate over time
fwrite(data.table(day = day, rho_eff = main$rho_eff,
                  rho_const = main$rho), "tikz_detection.csv")

# [NEW 9] Weight optimization heatmap
# [NEW 9] Weight optimization results (all 4 weights)
if (nrow(wopt) > 0) {
  fwrite(wopt[, c("w_cases", "w_hosp", "w_icu", "w_radar",
                   "rt_corr", "prev_corr", "rho", "alphaR")],
         "tikz_weight_heatmap.csv")
}

# Full results
fwrite(data.table(
  date = dates, day = day,
  S = main$S, E = main$E, I = main$I, R = main$R,
  beta = main$beta, Rt_model = main$Rt,
  mu_cases = main$mu_cases, mu_hosp = main$mu_hosp,
  mu_icu = main$mu_icu, mu_radar = main$mu_radar,
  cases_raw = cases_raw, hosp_raw = hosp_raw,
  icu_raw = icu_raw, radar_frac = radar_frac,
  rt_rivm = rt_rivm, prev_rivm = prev_rivm,
  tests_total = tests_tot, positivity = positivity,
  resid_cases = main$resid_cases, resid_hosp = main$resid_hosp,
  resid_icu = main$resid_icu, resid_radar = main$resid_radar,
  ode_S = main$ode_mat[, 1], ode_E = main$ode_mat[, 2],
  ode_I = main$ode_mat[, 3], ode_R = main$ode_mat[, 4],
  rho_eff = main$rho_eff
), "gp_seir_results.csv")

cat("  Saved all CSV files\n")

# ══════════════════════════════════════════════════════════════════
# STEP 7: PDF PLOTS (10 pages)
# ══════════════════════════════════════════════════════════════════

cat("\n=== STEP 7: Generating PDF ===\n")

lockdown1 <- as.Date("2020-03-15")
relax1    <- as.Date("2020-06-01")
lockdown2 <- as.Date("2020-10-14")
add_lock  <- function() {
  usr <- par("usr")
  rect(lockdown1, usr[3], relax1, usr[4], col = "#CC000020", border = NA)
  rect(lockdown2, usr[3], as.Date("2021-03-01"), usr[4],
       col = "#CC000020", border = NA)
}

pdf("gp_seir_plots.pdf", width = 14, height = 10)

# ── Page 1: Main fit (4 data streams) ──
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(dates, cases_raw, pch = 16, cex = 0.3, col = "gray60",
     xlab = "", ylab = "Daily cases",
     main = sprintf("A. Cases (rho=%.3f, RMSE=%.0f)", main$rho, main$cases_rmse))
add_lock(); lines(dates, main$mu_cases, col = "red", lwd = 2)
legend("topright", c("Data", "Model"), col = c("gray60", "red"),
       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), cex = 0.7)

plot(dates, hosp_raw, pch = 16, cex = 0.3, col = "gray60",
     xlab = "", ylab = "Admissions/day",
     main = sprintf("B. Hospital (pH=%.4f)", main$pH))
add_lock(); lines(dates, main$mu_hosp, col = "purple", lwd = 2)

ok_i <- is.finite(icu_raw)
plot(dates[ok_i], icu_raw[ok_i], pch = 16, cex = 0.3, col = "gray60",
     xlab = "", ylab = "Admissions/day",
     main = sprintf("C. ICU (pICU=%.4f)", main$pICU))
add_lock(); lines(dates, main$mu_icu, col = "darkred", lwd = 2)

ok_r <- is.finite(radar_frac)
plot(dates[ok_r], radar_frac[ok_r], pch = 16, cex = 0.3, col = "gray60",
     xlab = "", ylab = "I(t) fraction",
     main = sprintf("D. RADAR (alphaR=%.3f)", main$alphaR))
add_lock(); lines(dates, main$mu_radar, col = "darkgreen", lwd = 2)

# ── Page 2: Validation + dynamics ──
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

ok_rt <- is.finite(rt_rivm)
plot(dates, rt_rivm, type = "l", col = "gray40", lwd = 1,
     xlab = "", ylab = expression(R[t]),
     main = sprintf("E. Rt validation (corr=%.3f, RMSE=%.3f)",
                    main$rt_corr, main$rt_rmse),
     ylim = c(0, max(c(main$Rt, rt_rivm[ok_rt])) * 1.05))
add_lock(); lines(dates, main$Rt, col = "darkgreen", lwd = 2)
abline(h = 1, lty = 2, col = "gray")
legend("topright", c("RIVM Rt", "Model Rt"),
       col = c("gray40", "darkgreen"), lty = 1, lwd = c(1, 2), cex = 0.7)

ok_prev <- is.finite(prev_rivm)
plot(dates[ok_prev], prev_rivm[ok_prev] / 1000, type = "l",
     col = "gray40", lwd = 1, xlab = "", ylab = "Thousands",
     main = sprintf("F. Prevalence validation (corr=%.3f)", main$prev_corr),
     ylim = c(0, max(c(prev_rivm[ok_prev],
                        N_POP * main$I[ok_prev])) / 1000 * 1.05))
add_lock(); lines(dates, N_POP * main$I / 1000, col = "steelblue", lwd = 2)
legend("topright", c("RIVM prevalence", "Model N*I(t)"),
       col = c("gray40", "steelblue"), lty = 1, lwd = c(1, 2), cex = 0.7)

plot(dates, main$beta, type = "l", col = "purple", lwd = 2,
     xlab = "", ylab = expression(beta(t)),
     main = "G. Transmission rate")
add_lock(); abline(h = gamma, lty = 2, col = "gray60")

plot(dates, main$S, type = "l", col = "blue", lwd = 2,
     xlab = "", ylab = "Fraction", main = "H. SEIR compartments",
     ylim = c(0, 1))
lines(dates, main$E * 50, col = "orange", lwd = 2)
lines(dates, main$I * 50, col = "red", lwd = 2)
lines(dates, main$R, col = "darkgreen", lwd = 2)
legend("right", c("S", "E x 50", "I x 50", "R"),
       col = c("blue", "orange", "red", "darkgreen"),
       lty = 1, lwd = 2, cex = 0.7)

# ── Page 3: [NEW 2] Residual diagnostics ──
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(dates, main$resid_cases, type = "h", col = "gray40",
     xlab = "", ylab = "log(1+y) residual",
     main = "I. Case residuals")
abline(h = 0, col = "red")

plot(dates, main$resid_hosp, type = "h", col = "gray40",
     xlab = "", ylab = "log(1+y) residual",
     main = "J. Hospital residuals")
abline(h = 0, col = "red")

# Residual autocorrelation
acf_vals <- acf(main$resid_cases[is.finite(main$resid_cases)],
                lag.max = 30, plot = FALSE)
plot(acf_vals, main = "K. Case residual ACF", xlab = "Lag (days)")

# QQ plot of case residuals
qqnorm(main$resid_cases[is.finite(main$resid_cases)],
       main = "L. Case residual Q-Q",
       pch = 16, cex = 0.4, col = "gray40")
qqline(main$resid_cases[is.finite(main$resid_cases)], col = "red")

# ── Page 4: [NEW 4] ODE residuals per compartment ──
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(dates, main$ode_mat[, 1], type = "l", col = "blue", lwd = 1,
     xlab = "", ylab = "ODE residual",
     main = sprintf("M. ODE residual: S (RMSE=%.4f)", main$ode_rmse_S))
add_lock(); abline(h = 0, col = "red", lty = 2)

plot(dates, main$ode_mat[, 2], type = "l", col = "orange", lwd = 1,
     xlab = "", ylab = "ODE residual",
     main = sprintf("N. ODE residual: E (RMSE=%.4f)", main$ode_rmse_E))
add_lock(); abline(h = 0, col = "red", lty = 2)

plot(dates, main$ode_mat[, 3], type = "l", col = "red", lwd = 1,
     xlab = "", ylab = "ODE residual",
     main = sprintf("O. ODE residual: I (RMSE=%.4f)", main$ode_rmse_I))
add_lock(); abline(h = 0, col = "red", lty = 2)

plot(dates, main$ode_mat[, 4], type = "l", col = "darkgreen", lwd = 1,
     xlab = "", ylab = "ODE residual",
     main = sprintf("P. ODE residual: R (RMSE=%.4f)", main$ode_rmse_R))
add_lock(); abline(h = 0, col = "red", lty = 2)

# ── Page 5: Ablation bar chart ──
if (nrow(abl) > 1) {
par(mfrow = c(1, 1), mar = c(4, 10, 3, 2))
abl_p <- abl[order(abl[[SORT_CRITERION]]), ]
bp <- barplot(abl_p[[SORT_CRITERION]], horiz = TRUE, names.arg = abl_p$sources,
              las = 1, col = "steelblue", xlim = c(0, 1),
              main = paste0("Q. Source ablation: ", SORT_CRITERION),
              xlab = SORT_CRITERION)
abline(v = 0.8, lty = 2, col = "gray")
text(abl_p[[SORT_CRITERION]] - 0.03, bp,
     sprintf("%.3f", abl_p[[SORT_CRITERION]]), cex = 0.7)
}

# ── Page 6: [NEW 5] Beta across ablation configs ──
if (length(abl_betas) > 1) {
par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))
cols_abl <- c("red", "blue", "purple", "darkgreen",
              "orange", "cyan", "brown", "gray40")
plot(NULL, xlim = range(dates), ylim = c(0, 0.35),
     xlab = "", ylab = expression(beta(t)),
     main = "R. Transmission rate by source configuration")
add_lock()
nm_sorted <- abl$sources  # already sorted by SORT_CRITERION descending
for (i in seq_along(nm_sorted)) {
  nm <- nm_sorted[i]
  if (nm %in% names(abl_betas))
    lines(dates, abl_betas[[nm]], col = cols_abl[i], lwd = 1.5)
}
abline(h = gamma, lty = 2, col = "gray60")
legend("topright", nm_sorted, col = cols_abl[1:length(nm_sorted)],
       lty = 1, lwd = 1.5, cex = 0.55, bg = "white")
} # end if ablation betas

# ── Page 7: Testing, detection rate, and cumulative infections ──
par(mfrow = c(3, 1), mar = c(4, 4, 3, 4))
ok_t <- is.finite(tests_tot)
plot(dates[ok_t], tests_tot[ok_t], type = "l", col = "steelblue", lwd = 2,
     xlab = "", ylab = "Tests/day",
     main = "S. Daily testing volume and positivity rate")
add_lock()
par(new = TRUE)
ok_p <- is.finite(positivity)
plot(dates[ok_p], positivity[ok_p] * 100, type = "l", col = "red",
     lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(4, col = "red", col.axis = "red")
mtext("Positivity (%)", side = 4, line = 2.5, col = "red")
legend("topleft", c("Tests performed", "Positivity rate"),
       col = c("steelblue", "red"), lty = 1, lwd = 2, cex = 0.7)

# Effective detection rate
ok_rho <- is.finite(main$rho_eff)
if (any(ok_rho)) {
  rho_plot <- pmin(main$rho_eff, 1)  # cap at 1 for plotting
  plot(dates[ok_rho], rho_plot[ok_rho], type = "l", col = "darkblue", lwd = 2,
       xlab = "", ylab = expression(rho[eff](t)),
       main = sprintf("S2. Effective detection rate (fitted rho=%.3f)", main$rho),
       ylim = c(0, min(1, max(rho_plot[ok_rho], na.rm = TRUE) * 1.2)))
  add_lock()
  abline(h = main$rho, lty = 2, col = "red", lwd = 1.5)
  legend("topleft", c(expression(rho[eff](t)), expression(rho ~ "(constant, fitted)")),
         col = c("darkblue", "red"), lty = c(1, 2), lwd = c(2, 1.5), cex = 0.7)
}

# Cumulative infections
plot(dates, main$R * 100, type = "l", col = "darkgreen", lwd = 2,
     xlab = "", ylab = "Recovered (%)",
     main = sprintf("T. Cumulative infections (R_end=%.1f%%)", main$R_end * 100),
     ylim = c(0, max(main$R * 100) * 1.1))
add_lock()
abline(h = 15, lty = 2, col = "gray")
text(dates[50], 16, "Seroprevalence ~15%", cex = 0.7, col = "gray40")

# ── Page 8: W_DATA sensitivity panels ──
if (nrow(wsens) > 1) {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  x <- log10(wsens$W_DATA)

  plot(x, wsens$ode_rmse, type = "b", pch = 16, lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "ODE RMSE",
       main = "U. ODE compliance")
  plot(x, wsens$cases_rmse, type = "b", pch = 16, col = "red", lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "Cases RMSE",
       main = "V. Data fit (cases)")
  plot(x, wsens$rt_corr, type = "b", pch = 16, col = "darkgreen", lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "Rt corr",
       main = "W. Rt validation")
  abline(h = 0.8, lty = 2, col = "gray")
  plot(x, wsens$prev_corr, type = "b", pch = 16, col = "steelblue", lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "Prev corr",
       main = "X. Prevalence validation")
  abline(h = 0.8, lty = 2, col = "gray")
  plot(x, wsens$beta_max - wsens$beta_min, type = "b", pch = 16,
       col = "purple", lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "Beta range",
       main = "Y. Beta flexibility")
  plot(x, wsens$R_end * 100, type = "b", pch = 16, col = "brown", lwd = 2,
       xlab = expression(log[10](W[data])), ylab = "R(end) %",
       main = "Z. Cumulative infections")
  abline(h = 15, lty = 2, col = "gray")
}

# ── Page 9: Rt overlay for W_DATA ──
if (length(wsens_fits) > 1) {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  cols <- rainbow(length(wsens_fits))
  ylim_rt <- c(0, min(5, max(sapply(wsens_fits, function(x)
    if (!is.null(x)) max(x$Rt) else 2))))
  plot(dates, rt_rivm, type = "l", col = "black", lwd = 2.5,
       xlab = "", ylab = expression(R[t]),
       main = "AA. Rt trajectories across W_DATA", ylim = ylim_rt)
  for (i in seq_along(wsens_fits)) {
    if (!is.null(wsens_fits[[i]]))
      lines(dates, wsens_fits[[i]]$Rt, col = cols[i], lwd = 1.5, lty = 2)
  }
  abline(h = 1, lty = 3, col = "gray60")
  labs <- sprintf("W=%.0e (r=%.2f)", wsens$W_DATA, wsens$rt_corr)
  legend("topright", c("RIVM Rt", labs),
         col = c("black", cols[1:nrow(wsens)]),
         lty = c(1, rep(2, nrow(wsens))),
         lwd = c(2.5, rep(1.5, nrow(wsens))), cex = 0.55)
}

# ── Page 10: Weight optimization heatmap ──
if (nrow(wopt) > 1) {
par(mfrow = c(1, 1), mar = c(5, 5, 3, 5))
wc_vals <- sort(unique(wopt$w_cases))
wr_vals <- sort(unique(wopt$w_radar))
mat_val <- matrix(NA, nrow = length(wc_vals), ncol = length(wr_vals))
for (i in seq_along(wc_vals)) {
  for (j in seq_along(wr_vals)) {
    rows <- wopt[wopt$w_cases == wc_vals[i] & wopt$w_radar == wr_vals[j], ]
    if (nrow(rows) > 0)
      mat_val[i, j] <- max(rows[[SORT_CRITERION]], na.rm = TRUE)
  }
}
image(wc_vals, wr_vals, mat_val,
      xlab = expression(w[cases]), ylab = expression(w[radar]),
      main = paste0("AB. Best ", SORT_CRITERION, " (over w_hosp, w_icu)"),
      col = hcl.colors(20, "Blues3", rev = TRUE))
for (i in seq_along(wc_vals)) {
  for (j in seq_along(wr_vals)) {
    if (!is.na(mat_val[i, j]))
      text(wc_vals[i], wr_vals[j], sprintf("%.3f", mat_val[i, j]),
           cex = 0.7)
  }
}
} # end if weight heatmap

dev.off()
cat("  Saved: gp_seir_plots.pdf (10 pages)\n")

# ══════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("ALL DONE. Output files:\n\n")
cat("  CSV (results):\n")
cat("    gp_seir_results.csv       — main fit (all trajectories + residuals)\n")
cat("    ablation_results.csv      — source ablation (validation + fit quality)\n")
cat("    ablation_params.csv       — estimated parameters per ablation config\n")
cat("    weight_optimization.csv   — per-source weight search\n")
cat("    wdata_sensitivity.csv     — penalty weight sensitivity\n")
cat("    multistart_params.csv     — parameter variation across starts\n")
cat("    wave_params.csv           — per-period diagnostics (sliced from main fit)\n")
cat("    wave_refit_params.csv     — per-wave independent refit parameters\n")
cat("\n  CSV (TikZ/pgfplots):\n")
cat("    tikz_seir.csv             — SEIR + beta + Rt trajectories\n")
cat("    tikz_cases.csv            — cases data vs model\n")
cat("    tikz_hosp.csv             — hospital data vs model\n")
cat("    tikz_icu.csv              — ICU data vs model\n")
cat("    tikz_radar.csv            — RADAR data vs model\n")
cat("    tikz_rt.csv               — Rt validation\n")
cat("    tikz_prev.csv             — prevalence validation\n")
cat("    tikz_ablation.csv         — ablation bar chart data\n")
cat("    tikz_residuals.csv        — data + ODE residuals over time\n")
cat("    tikz_acf.csv              — case residual autocorrelation\n")
cat("    tikz_qq.csv               — case residual Q-Q plot data\n")
cat("    tikz_beta_ablation.csv    — beta(t) per ablation config\n")
cat("    tikz_testing.csv          — testing volume + positivity\n")
cat("    tikz_cumulative.csv       — cumulative infections R(t)\n")
cat("    tikz_detection.csv        — effective detection rate rho_eff(t)\n")
cat("    tikz_rt_wdata.csv         — Rt trajectories per W_DATA value\n")
cat("    tikz_weight_heatmap.csv   — weight optimization grid\n")
cat("\n  PDF:\n")
cat("    gp_seir_plots.pdf         — 10 pages of figures\n")
cat("============================================================\n")
