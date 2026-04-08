#!/usr/bin/env Rscript
##############################################################################
# gp_seir_wave_analysis.R
#
# Comprehensive per-wave GP-SEIR analysis for the paper:
#   "Calibration of Epidemic Models via Multi-Source Generalized Profiling:
#    Data Source Valuation Through Cooperative Game Theory"
#
# This script fits the GP-SEIR model independently to each pandemic wave
# and performs the full suite of analyses needed for the paper:
#   - Main fit (4-source) per wave
#   - Source ablation (all 15 subsets) per wave
#   - Shapley values and interaction indices per wave
#   - Leave-one-out analysis per wave
#   - Superadditivity assessment per wave
#   - Initialization sensitivity per wave
#   - Start-date sensitivity (wave 1)
#   - Profile likelihood per wave
#   - Sigma/gamma sensitivity per wave
#   - Residual diagnostics per wave
#   - Cross-correlations per wave
#   - ODE compliance per wave
#   - Derived quantities (flows, growth rates, detection rates)
#   - Cross-wave comparison
#   - TikZ-ready CSV outputs for all figures
#
# All outputs are prefixed with "wave{N}_" (N = 1, 2, 3 for inter-wave).
# A final cross-wave comparison file summarizes results across periods.
#
# Dependencies: data.table, fda, minpack.lm, zoo (all on CRAN)
#
# Usage:
#   Rscript gp_seir_wave_analysis.R
#
# Code repository: https://github.com/hadiabp/gp-seir-multisource
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fda)
  library(minpack.lm)
  library(zoo)
})

cat("================================================================\n")
cat("GP-SEIR Comprehensive Wave Analysis\n")
cat("================================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: GLOBAL SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

DATA_FILE <- "gp_input_final_v2.csv"

# Fixed epidemiological parameters
N_POP <- 17400000
SIGMA <- 1 / 5.5       # Latent rate (mean latent period 5.5 days)
GAMMA <- 1 / 9.5       # Removal rate (mean infectious period 9.5 days)
T_GEN <- 1/SIGMA + 1/GAMMA  # Generation time (~15 days)

# Spline settings
STATE_KNOT_DAYS <- 14
STATE_NORDER    <- 4   # Cubic B-splines
BETA_KNOT_DAYS  <- 28
BETA_NORDER     <- 4

# Penalty weights
ETA        <- 1e6      # ODE penalty
KAPPA_BETA <- 1e3      # Roughness of log-transmission
KAPPA_BMAG <- 1e4      # Soft upper bound on beta
BETA_MAX   <- 0.6      # Maximum plausible beta
KAPPA_MASS <- 1e6      # Mass conservation

# Optimization settings
N_WARM_STARTS  <- 5
N_COLD_STARTS  <- 3
MAX_INNER      <- 300    # Levenberg-Marquardt iterations
MAX_OUTER      <- 40     # Gradient descent iterations
STEP0          <- 1.0    # Initial step size
BETA_PERTURB   <- 0.30   # Perturbation SD for beta coefficients
PARAM_PERTURB  <- 0.50   # Perturbation SD for scaling parameters

# Delay convolution settings
USE_DELAY_KERNEL   <- TRUE
DELAY_CASE_MEAN    <- 3.0;  DELAY_CASE_SD <- 1.5
DELAY_HOSP_MEAN    <- 2.0;  DELAY_HOSP_SD <- 1.0
DELAY_ICU_MEAN     <- 3.5;  DELAY_ICU_SD  <- 1.5
DELAY_MAX_LAG      <- 14

# Parameter bounds (logistic transform)
BOUNDS <- list(
  rho    = c(0.05, 0.60),
  pH     = c(0.002, 0.05),
  pICU   = c(0.0003, 0.02),
  alphaR = c(0.10, 1.00)
)

# Wave definitions
# Note: fit_start/fit_end define the B-spline fitting window. For wave 1,
# the fit extends to June 30 for boundary stability, but the wave-1 period
# for reporting purposes is March-June 1. The validation metrics are computed
# over the full fitting window. For sliced metrics matching the full-period
# analysis, use the gp_seir_complete.R code instead.
WAVE_DEFS <- list(
  wave1 = list(
    id          = 1,
    label       = "Wave 1 (Mar-Jun 2020)",
    label_short = "wave1",
    fit_start   = as.Date("2020-03-01"),
    fit_end     = as.Date("2020-06-30")
  ),
  interwave = list(
    id          = 3,
    label       = "Inter-wave (Jun-Oct 2020)",
    label_short = "interwave",
    fit_start   = as.Date("2020-06-01"),
    fit_end     = as.Date("2020-10-01")
  ),
  wave2 = list(
    id          = 2,
    label       = "Wave 2 (Oct 2020-Mar 2021)",
    label_short = "wave2",
    fit_start   = as.Date("2020-10-01"),
    fit_end     = as.Date("2021-03-01")
  )
)

# Which waves get the full ablation + Shapley analysis (expensive!)
# Set to c("wave1", "wave2") to skip interwave ablation
ABLATION_WAVES <- c("wave1", "wave2", "interwave")

# Profile likelihood settings
PROFILE_NPTS <- 15

# Source names and all 15 non-empty subsets
SOURCE_NAMES <- c("cases", "hosp", "icu", "radar")
ALL_SUBSETS <- list()
for (k in 1:4) {
  combos <- combn(SOURCE_NAMES, k, simplify = FALSE)
  ALL_SUBSETS <- c(ALL_SUBSETS, combos)
}

cat(sprintf("  sigma=1/%.1f, gamma=1/%.1f, T_gen=%.1f days\n", 1/SIGMA, 1/GAMMA, T_GEN))
cat(sprintf("  USE_DELAY_KERNEL=%s\n", USE_DELAY_KERNEL))
cat(sprintf("  %d source subsets to evaluate\n\n", length(ALL_SUBSETS)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# --- Delay kernel construction ---
make_delay_kernel <- function(mean_d, sd_d, max_lag = DELAY_MAX_LAG) {
  if (mean_d <= 0) return(1)
  lags <- 0:max_lag
  mu_ln  <- log(mean_d^2 / sqrt(sd_d^2 + mean_d^2))
  sig_ln <- sqrt(log(1 + sd_d^2 / mean_d^2))
  w <- dlnorm(lags, mu_ln, sig_ln)
  w / sum(w)
}

# --- Discrete convolution with delay kernel ---
apply_delay <- function(y, kernel) {
  klen <- length(kernel)
  if (klen == 1) return(y)
  n_y <- length(y)
  y_conv <- numeric(n_y)
  for (t_idx in 1:n_y) {
    for (k_idx in 1:klen) {
      src <- t_idx - (k_idx - 1)
      if (src >= 1) y_conv[t_idx] <- y_conv[t_idx] + kernel[k_idx] * y[src]
    }
  }
  y_conv
}

# --- Logistic transforms for bounded parameters ---
to_u   <- function(x, lo, hi) qlogis((x - lo) / (hi - lo))
from_u <- function(u, lo, hi) lo + (hi - lo) * plogis(u)

# --- Weighted SD helper ---
wsd <- function(y) {
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  1 / s^2
}

# --- SEIR log-space right-hand side ---
seir_rhs <- function(Z, beta, sig = SIGMA, gam = GAMMA) {
  cbind(
    -beta * exp(Z[, 3]),
    beta * exp(Z[, 1] + Z[, 3] - Z[, 2]) - sig,
    sig * exp(Z[, 2] - Z[, 3]) - gam,
    gam * exp(Z[, 3] - Z[, 4])
  )
}

# --- Create B-spline basis ---
make_basis <- function(rng, knot_days, norder) {
  kn <- unique(c(rng[1], seq(rng[1], rng[2], by = knot_days), rng[2]))
  create.bspline.basis(rng, length(kn) + norder - 2, norder, kn)
}

# --- Subset indicator vector ---
subset_to_weights <- function(sources) {
  c(
    cases = as.numeric("cases" %in% sources),
    hosp  = as.numeric("hosp"  %in% sources),
    icu   = as.numeric("icu"   %in% sources),
    radar = as.numeric("radar" %in% sources)
  )
}

# --- Safe correlation ---
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}

# --- Safe RMSE ---
safe_rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: CORE FITTING ENGINE
# ══════════════════════════════════════════════════════════════════════════════

#' Fit GP-SEIR to a time window with specified source weights
#'
#' @param df        Data.table with the required columns
#' @param date_start Start date
#' @param date_end  End date
#' @param wc,wh,wi,wr Binary source weights (0 or 1)
#' @param init_method "hospital", "constant", or "random_smooth"
#' @param seed      Random seed
#' @param n_warm    Number of warm starts
#' @param n_cold    Number of cold starts
#' @param sig       Latent rate (default SIGMA)
#' @param gam       Removal rate (default GAMMA)
#' @param verbose   Print progress
#' @return List with fitted results
fit_gp_seir <- function(df, date_start, date_end,
                        wc = 1, wh = 1, wi = 1, wr = 1,
                        init_method = "hospital",
                        seed = 42,
                        n_warm = N_WARM_STARTS,
                        n_cold = N_COLD_STARTS,
                        sig = SIGMA,
                        gam = GAMMA,
                        bounds_override = NULL,
                        w_data_scale = 1.0,
                        verbose = TRUE) {

  # Use overridden bounds if provided, otherwise global
  bnds <- if (!is.null(bounds_override)) bounds_override else BOUNDS

  # --- Subset and prepare data ---
  dd <- copy(df)[date >= date_start & date <= date_end][order(date)]
  n  <- nrow(dd)
  if (n < 30) stop("Too few data points (", n, ") in selected window.")
  dd[, t := as.numeric(date - date_start)]
  tv    <- dd$t
  rng   <- range(tv)
  dt    <- 1
  dates <- dd$date

  cases_raw  <- as.numeric(dd$cases_raw)
  hosp_raw   <- as.numeric(dd$hosp_nice)
  icu_raw    <- as.numeric(dd$icu_daily)
  radar_frac <- as.numeric(dd$radar_I_frac_sm7)
  rt_rivm    <- as.numeric(dd$rt_rivm)
  prev_rivm  <- as.numeric(dd$prev_json_avg)

  # --- Transform to fitting scale ---
  y_cases <- log1p(pmax(0, cases_raw))
  y_hosp  <- log1p(pmax(0, hosp_raw))
  y_icu   <- log1p(pmax(0, icu_raw))
  y_radar <- radar_frac

  wb_cases <- wsd(y_cases) * w_data_scale
  wb_hosp  <- wsd(y_hosp)  * w_data_scale
  wb_icu   <- wsd(y_icu)   * w_data_scale
  wb_radar <- wsd(y_radar) * w_data_scale

  # --- Build spline bases ---
  Bs  <- make_basis(rng, STATE_KNOT_DAYS, STATE_NORDER)
  Bb  <- make_basis(rng, BETA_KNOT_DAYS, BETA_NORDER)
  nbs <- Bs$nbasis
  nbb <- Bb$nbasis

  Ph   <- eval.basis(tv, Bs, 0)
  dPh  <- eval.basis(tv, Bs, 1)
  Pb   <- eval.basis(tv, Bb, 0)
  d2Pb <- eval.basis(tv, Bb, 2)

  # --- Delay kernels ---
  if (USE_DELAY_KERNEL) {
    h_case <- make_delay_kernel(DELAY_CASE_MEAN, DELAY_CASE_SD)
    h_hosp <- make_delay_kernel(DELAY_HOSP_MEAN, DELAY_HOSP_SD)
    h_icu  <- make_delay_kernel(DELAY_ICU_MEAN,  DELAY_ICU_SD)
  } else {
    h_case <- 1; h_hosp <- 1; h_icu <- 1
  }

  # --- State initialization ---
  c0 <- ifelse(is.finite(cases_raw), cases_raw, 0)
  E0 <- pmax(1e-10, c0 / (N_POP * 0.20 * sig))
  win_roll <- min(7, max(3, n %/% 10))
  I0 <- as.numeric(rollmean(E0, win_roll, fill = median(E0), align = "right"))
  R0 <- pmin(0.8, cumsum(pmax(c0, 0)) / N_POP)
  S0 <- pmax(1e-6, 1 - E0 - I0 - R0)
  m0 <- S0 + E0 + I0 + R0
  S0 <- S0 / m0; E0 <- E0 / m0; I0 <- I0 / m0; R0 <- R0 / m0

  Ci <- smooth.basis(tv, cbind(log(S0), log(E0), log(I0), log(R0)),
                     fdPar(Bs, int2Lfd(2), 1e-4))$fd$coefs

  # --- Beta initialization ---
  if (init_method == "hospital") {
    h0_safe <- ifelse(is.finite(hosp_raw), pmax(0.5, hosp_raw), 0.5)
    win_h <- min(28, max(7, n %/% 5))
    hosp_sm <- as.numeric(rollmean(h0_safe, win_h, fill = NA, align = "center"))
    hosp_sm[is.na(hosp_sm)] <- hosp_sm[which(!is.na(hosp_sm))[1]]
    log_hosp <- log(pmax(1, hosp_sm))
    dlog_h <- c(diff(log_hosp), 0)
    win_d <- min(21, max(5, n %/% 6))
    dlog_sm <- as.numeric(rollmean(dlog_h, win_d, fill = 0, align = "center"))
    dlog_sm[is.na(dlog_sm)] <- 0
    Tc <- 1/sig + 1/gam
    b0 <- gam * (1 + dlog_sm * Tc)
    b0 <- pmax(0.05, pmin(0.30, b0))
    win_b <- min(14, max(3, n %/% 10))
    b0 <- as.numeric(rollmean(b0, win_b, fill = NA, align = "center"))
    b0[is.na(b0)] <- median(b0, na.rm = TRUE)
  } else if (init_method == "constant") {
    b0 <- rep(gam, n)
  } else if (init_method == "random_smooth") {
    set.seed(seed + 1000)
    b0 <- runif(n, min = 0.05, max = 0.25)
    win_b <- min(14, max(3, n %/% 10))
    b0 <- as.numeric(rollmean(b0, win_b, fill = NA, align = "center"))
    b0[is.na(b0)] <- median(b0, na.rm = TRUE)
  } else {
    stop("Unknown init_method: ", init_method)
  }

  ai <- as.numeric(smooth.basis(tv, log(pmax(b0, 1e-6)),
                                fdPar(Bb, int2Lfd(2), 1e-2))$fd$coefs)

  # Clamp initial scaling parameters to lie within bounds (essential when
  # bounds are narrowed for profile likelihood)
  clamp <- function(x, lo, hi) max(lo + 1e-8 * (hi - lo), min(hi - 1e-8 * (hi - lo), x))

  thi <- c(ai,
           to_u(clamp(0.20,   bnds$rho[1],    bnds$rho[2]),    bnds$rho[1],    bnds$rho[2]),
           to_u(clamp(0.012,  bnds$pH[1],     bnds$pH[2]),     bnds$pH[1],     bnds$pH[2]),
           to_u(clamp(0.0025, bnds$pICU[1],   bnds$pICU[2]),   bnds$pICU[1],   bnds$pICU[2]),
           to_u(clamp(0.40,   bnds$alphaR[1], bnds$alphaR[2]), bnds$alphaR[1], bnds$alphaR[2]))
  beta_init_curve <- exp(as.numeric(Pb %*% ai))

  # --- Inner solver (Levenberg-Marquardt) ---
  solve_inner <- function(Cs, th, w_c, w_h, w_i, w_r) {
    av   <- th[1:nbb]
    rho  <- from_u(th[nbb + 1], bnds$rho[1],    bnds$rho[2])
    pH   <- from_u(th[nbb + 2], bnds$pH[1],     bnds$pH[2])
    pICU <- from_u(th[nbb + 3], bnds$pICU[1],   bnds$pICU[2])
    aR   <- from_u(th[nbb + 4], bnds$alphaR[1], bnds$alphaR[2])
    g    <- as.numeric(Pb %*% av)
    beta <- exp(g)

    fn <- function(cf) {
      C  <- matrix(cf, nrow = nbs, ncol = 4)
      Z  <- Ph %*% C;  dZ <- dPh %*% C
      Sv <- exp(Z[, 1]); Ev <- exp(Z[, 2])
      Iv <- exp(Z[, 3]); Rv <- exp(Z[, 4])

      fc <- N_POP * rho  * sig * Ev
      fh <- N_POP * pH   * gam * Iv
      fi <- N_POP * pICU * gam * Iv
      mc <- log1p(pmax(0, apply_delay(fc, h_case)))
      mh <- log1p(pmax(0, apply_delay(fh, h_hosp)))
      mi <- log1p(pmax(0, apply_delay(fi, h_icu)))
      mr <- aR * Iv

      rc <- sqrt(wb_cases * w_c) * ifelse(is.finite(y_cases), y_cases - mc, 0)
      rh <- sqrt(wb_hosp  * w_h) * ifelse(is.finite(y_hosp),  y_hosp  - mh, 0)
      ri <- sqrt(wb_icu   * w_i) * ifelse(is.finite(y_icu),   y_icu   - mi, 0)
      rr <- sqrt(wb_radar * w_r) * ifelse(is.finite(y_radar), y_radar - mr, 0)

      rp <- sqrt(ETA * dt) * as.numeric(dZ - seir_rhs(Z, beta, sig, gam))
      rm <- sqrt(KAPPA_MASS) * (Sv + Ev + Iv + Rv - 1)
      rb <- sqrt(KAPPA_BMAG * dt) * pmax(beta - BETA_MAX, 0)
      rs <- sqrt(KAPPA_BETA * dt) * as.numeric(d2Pb %*% av)

      c(rc, rh, ri, rr, rp, rm, rb, rs)
    }

    r <- nls.lm(par = as.numeric(Cs), fn = fn,
                control = nls.lm.control(maxiter = MAX_INNER,
                                         ftol = 1e-10, ptol = 1e-10))
    list(C = matrix(r$par, nrow = nbs, ncol = 4), J = sum(r$fvec^2))
  }

  # --- Parameter cascading (outer loop) ---
  cascade <- function(th0, C0, w_c, w_h, w_i, w_r) {
    th <- th0; Cc <- C0; bJ <- Inf
    for (oi in 1:MAX_OUTER) {
      inn <- solve_inner(Cc, th, w_c, w_h, w_i, w_r)
      Cc <- inn$C; Jc <- inn$J
      if (Jc < bJ) bJ <- Jc

      eps <- 1e-4; gr <- numeric(length(th))
      for (j in seq_along(th)) {
        tp <- th; tp[j] <- tp[j] + eps
        gr[j] <- (solve_inner(Cc, tp, w_c, w_h, w_i, w_r)$J - Jc) / eps
      }
      gn <- sqrt(sum(gr^2))
      if (gn > 1e-12) {
        d <- -gr / gn; st <- STEP0; found <- FALSE
        for (ls in 1:10) {
          tt <- th + st * d
          it <- solve_inner(Cc, tt, w_c, w_h, w_i, w_r)
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

  # --- Extract all quantities from a fitted result ---
  extract_results <- function(res) {
    th <- res$theta; C <- res$C; av <- th[1:nbb]
    rho  <- from_u(th[nbb + 1], bnds$rho[1],    bnds$rho[2])
    pH   <- from_u(th[nbb + 2], bnds$pH[1],     bnds$pH[2])
    pICU <- from_u(th[nbb + 3], bnds$pICU[1],   bnds$pICU[2])
    aR   <- from_u(th[nbb + 4], bnds$alphaR[1], bnds$alphaR[2])

    Z  <- Ph %*% C;  dZ <- dPh %*% C
    Sv <- exp(Z[, 1]); Ev <- exp(Z[, 2])
    Iv <- exp(Z[, 3]); Rv <- exp(Z[, 4])
    beta <- exp(as.numeric(Pb %*% av))
    Rt   <- beta * Sv / gam
    mass <- Sv + Ev + Iv + Rv

    ode_mat  <- dZ - seir_rhs(Z, beta, sig, gam)
    ode_rmse <- sqrt(mean(ode_mat^2))

    fc <- N_POP * rho  * sig * Ev
    fh <- N_POP * pH   * gam * Iv
    fi <- N_POP * pICU * gam * Iv
    mc <- apply_delay(fc, h_case)
    mh <- apply_delay(fh, h_hosp)
    mi <- apply_delay(fi, h_icu)
    mr <- aR * Iv

    I_count <- N_POP * Iv

    rt_corr_p  <- safe_cor(rt_rivm, Rt, "pearson")
    rt_corr_s  <- safe_cor(rt_rivm, Rt, "spearman")
    rt_rmse_v  <- safe_rmse(rt_rivm, Rt)
    prev_corr_p <- safe_cor(prev_rivm, I_count, "pearson")
    prev_corr_s <- safe_cor(prev_rivm, I_count, "spearman")
    prev_rmse_v <- safe_rmse(prev_rivm, I_count)

    ok_rt <- is.finite(rt_rivm) & is.finite(Rt)
    rt_bias <- if (sum(ok_rt) > 3) mean(Rt[ok_rt] - rt_rivm[ok_rt]) else NA_real_
    threshold_agree <- if (sum(ok_rt) > 3) {
      mean((Rt[ok_rt] > 1) == (rt_rivm[ok_rt] > 1))
    } else NA_real_

    denom <- N_POP * sig * Ev
    rho_eff <- ifelse(is.finite(cases_raw) & denom > 1, cases_raw / denom, NA)

    # Flows
    new_infections <- N_POP * beta * Sv * Iv
    onset_flow     <- N_POP * sig * Ev
    removal_flow   <- N_POP * gam * Iv

    # Growth rate
    growth_rate <- c(diff(log(pmax(Iv, 1e-15))), NA)

    # Residuals
    resid_cases <- ifelse(is.finite(y_cases), y_cases - log1p(pmax(0, mc)), NA)
    resid_hosp  <- ifelse(is.finite(y_hosp),  y_hosp  - log1p(pmax(0, mh)), NA)
    resid_icu   <- ifelse(is.finite(y_icu),   y_icu   - log1p(pmax(0, mi)), NA)
    resid_radar <- ifelse(is.finite(y_radar), y_radar - mr, NA)

    # Cases RMSE on raw scale
    cases_rmse <- safe_rmse(cases_raw, mc)
    hosp_rmse  <- safe_rmse(hosp_raw, mh)

    list(
      # States
      S = Sv, E = Ev, I = Iv, R = Rv, beta = beta, Rt = Rt, mass = mass,
      # Model predictions
      mu_cases = mc, mu_hosp = mh, mu_icu = mi, mu_radar = mr,
      # Parameters
      rho = rho, pH = pH, pICU = pICU, alphaR = aR,
      # Validation
      rt_corr = rt_corr_p, rt_spearman = rt_corr_s,
      rt_rmse = rt_rmse_v, rt_bias = rt_bias,
      prev_corr = prev_corr_p, prev_spearman = prev_corr_s,
      prev_rmse = prev_rmse_v,
      threshold_agree = threshold_agree,
      # Quality
      ode_rmse = ode_rmse, cases_rmse = cases_rmse, hosp_rmse = hosp_rmse,
      mass_err = max(abs(mass - 1)),
      # Ranges
      beta_range = c(min(beta), max(beta)),
      Rt_range   = c(min(Rt), max(Rt)),
      R_end = Rv[n],
      # Flows
      new_infections = new_infections, onset_flow = onset_flow,
      removal_flow = removal_flow, growth_rate = growth_rate,
      # Residuals
      resid_cases = resid_cases, resid_hosp = resid_hosp,
      resid_icu = resid_icu, resid_radar = resid_radar,
      # Detection
      rho_eff = rho_eff,
      # ODE residuals per compartment
      ode_S = ode_mat[, 1], ode_E = ode_mat[, 2],
      ode_I = ode_mat[, 3], ode_R = ode_mat[, 4]
    )
  }

  # --- Multi-start optimization ---
  set.seed(seed)
  best <- NULL; bJ <- Inf; all_results <- list(); all_J <- numeric(0)

  if (verbose) cat(sprintf("    %d warm + %d cold starts...", n_warm, n_cold))

  for (s in 1:n_warm) {
    ts <- thi; Cs <- Ci
    if (s > 1) {
      ts[1:nbb] <- ts[1:nbb] + rnorm(nbb, 0, BETA_PERTURB)
      for (k in (nbb + 1):(nbb + 4)) ts[k] <- ts[k] + rnorm(1, 0, PARAM_PERTURB)
    }
    r <- tryCatch(cascade(ts, Cs, wc, wh, wi, wr), error = function(e) NULL)
    if (!is.null(r)) {
      all_J <- c(all_J, r$J)
      all_results[[length(all_results) + 1]] <- r
      if (r$J < bJ) { bJ <- r$J; best <- r }
    }
  }

  for (s in 1:n_cold) {
    ts <- thi; Cs <- Ci
    ts[1:nbb] <- ts[1:nbb] + rnorm(nbb, 0, 0.50)
    for (k in (nbb + 1):(nbb + 4)) ts[k] <- ts[k] + rnorm(1, 0, 0.80)
    r <- tryCatch(cascade(ts, Cs, wc, wh, wi, wr), error = function(e) NULL)
    if (!is.null(r)) {
      all_J <- c(all_J, r$J)
      all_results[[length(all_results) + 1]] <- r
      if (r$J < bJ) { bJ <- r$J; best <- r }
    }
  }

  if (is.null(best)) stop("All starts failed.")

  main <- extract_results(best)

  # Extract parameters from all converged starts for multi-start tracking
  multistart_data <- list()
  for (k in seq_along(all_results)) {
    mk <- extract_results(all_results[[k]])
    multistart_data[[k]] <- data.table(
      start = k, J = all_J[k],
      rho = mk$rho, pH = mk$pH, pICU = mk$pICU, alphaR = mk$alphaR,
      rt_corr = mk$rt_corr, prev_corr = mk$prev_corr
    )
  }
  multistart_dt <- if (length(multistart_data) > 0) rbindlist(multistart_data) else NULL

  if (verbose) {
    cat(sprintf(" J=%.1f | Rt=%.3f | Prev=%.3f\n",
                bJ, ifelse(is.na(main$rt_corr), NA, main$rt_corr),
                ifelse(is.na(main$prev_corr), NA, main$prev_corr)))
  }

  list(
    dates       = dates,
    tv          = tv,
    n           = n,
    J           = bJ,
    main        = main,
    all_J       = all_J,
    n_converged = sum(all_J <= bJ * 1.10),
    beta_init   = beta_init_curve,
    multistart  = multistart_dt,
    raw = list(
      cases_raw = cases_raw, hosp_raw = hosp_raw,
      icu_raw = icu_raw, radar_frac = radar_frac,
      rt_rivm = rt_rivm, prev_rivm = prev_rivm,
      y_cases = y_cases, y_hosp = y_hosp,
      y_icu = y_icu, y_radar = y_radar
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Loading data ---\n")
if (!file.exists(DATA_FILE)) {
  stop("Data file not found: ", DATA_FILE,
       "\nPlease place gp_input_final_v2.csv in the working directory.")
}
FULL_DF <- fread(DATA_FILE)[, date := as.Date(date)][order(date)]
cat(sprintf("  %d rows (%s to %s)\n\n", nrow(FULL_DF),
            min(FULL_DF$date), max(FULL_DF$date)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: MAIN FIT FOR EACH WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("STEP 1: Main fit (all 4 sources) per wave\n")
cat("================================================================\n")

wave_fits <- list()
for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  fit <- fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                     wc = 1, wh = 1, wi = 1, wr = 1,
                     init_method = "hospital", seed = 42 + wdef$id)
  wave_fits[[wname]] <- fit

  # Save main results CSV
  out_prefix <- wdef$label_short
  beta_abs_diff    <- fit$main$beta - fit$beta_init
  beta_rel_diff_pct <- 100 * beta_abs_diff / pmax(fit$beta_init, 1e-8)
  beta_init_fit_corr <- cor(fit$beta_init, fit$main$beta)

  results_df <- data.table(
    date           = fit$dates,
    day            = fit$tv,
    S              = fit$main$S,
    E              = fit$main$E,
    I              = fit$main$I,
    R              = fit$main$R,
    beta_init      = fit$beta_init,
    beta           = fit$main$beta,
    beta_abs_diff  = beta_abs_diff,
    beta_rel_diff_pct = beta_rel_diff_pct,
    Rt_model       = fit$main$Rt,
    mu_cases       = fit$main$mu_cases,
    mu_hosp        = fit$main$mu_hosp,
    mu_icu         = fit$main$mu_icu,
    mu_radar       = fit$main$mu_radar,
    cases_raw      = fit$raw$cases_raw,
    hosp_raw       = fit$raw$hosp_raw,
    icu_raw        = fit$raw$icu_raw,
    radar_frac     = fit$raw$radar_frac,
    rt_rivm        = fit$raw$rt_rivm,
    prev_rivm      = fit$raw$prev_rivm,
    prev_model     = N_POP * fit$main$I,
    rho_eff        = fit$main$rho_eff,
    new_infections = fit$main$new_infections,
    onset_flow     = fit$main$onset_flow,
    removal_flow   = fit$main$removal_flow,
    growth_rate    = fit$main$growth_rate,
    resid_cases    = fit$main$resid_cases,
    resid_hosp     = fit$main$resid_hosp,
    resid_icu      = fit$main$resid_icu,
    resid_radar    = fit$main$resid_radar,
    ode_S          = fit$main$ode_S,
    ode_E          = fit$main$ode_E,
    ode_I          = fit$main$ode_I,
    ode_R          = fit$main$ode_R,
    mass           = fit$main$S + fit$main$E + fit$main$I + fit$main$R
  )
  fwrite(results_df, sprintf("%s_results.csv", out_prefix))

  # Save TikZ-ready CSVs
  fwrite(results_df[, .(day, cases_raw, mu_cases)], sprintf("%s_tikz_cases.csv", out_prefix))
  fwrite(results_df[, .(day, hosp_raw, mu_hosp)], sprintf("%s_tikz_hosp.csv", out_prefix))
  fwrite(results_df[, .(day, icu_raw, mu_icu)], sprintf("%s_tikz_icu.csv", out_prefix))
  fwrite(results_df[, .(day, radar_frac, mu_radar)], sprintf("%s_tikz_radar.csv", out_prefix))
  fwrite(results_df[, .(day, rt_rivm, rt_model = Rt_model)], sprintf("%s_tikz_rt.csv", out_prefix))
  fwrite(results_df[, .(day, prev_rivm, I_model = prev_model)], sprintf("%s_tikz_prev.csv", out_prefix))
  fwrite(results_df[, .(day, S, E, I, R, beta)], sprintf("%s_tikz_seir.csv", out_prefix))

  # Rt scatter plot (model vs RIVM)
  ok_rt <- is.finite(results_df$rt_rivm) & is.finite(results_df$Rt_model)
  if (any(ok_rt)) {
    fwrite(results_df[ok_rt, .(rt_rivm, rt_model = Rt_model)],
           sprintf("%s_tikz_rt_scatter.csv", out_prefix))
  }

  # Bland-Altman plot data
  if (any(ok_rt)) {
    ba_dt <- results_df[ok_rt, .(
      mean_rt = (Rt_model + rt_rivm) / 2,
      diff_rt = Rt_model - rt_rivm
    )]
    ba_bias <- mean(ba_dt$diff_rt)
    ba_sd   <- sd(ba_dt$diff_rt)
    fwrite(cbind(ba_dt, bias = ba_bias, loa_lo = ba_bias - 1.96*ba_sd,
                 loa_hi = ba_bias + 1.96*ba_sd),
           sprintf("%s_tikz_bland_altman.csv", out_prefix))
  }

  # Detection rate and cumulative infections
  fwrite(results_df[, .(day, rho_eff, R_pct = R * 100)],
         sprintf("%s_tikz_detection.csv", out_prefix))

  # Initial vs fitted beta(t)
  fwrite(results_df[, .(day, beta_init, beta, beta_abs_diff, beta_rel_diff_pct)],
         sprintf("%s_tikz_beta_init.csv", out_prefix))

  cat(sprintf("  Saved: %s_results.csv + TikZ CSVs\n", out_prefix))
  cat(sprintf("    beta init: corr=%.3f, mean|diff|=%.4f, median|diff%%|=%.1f%%\n",
              beta_init_fit_corr,
              mean(abs(beta_abs_diff)),
              median(abs(beta_rel_diff_pct), na.rm = TRUE)))

  # Save multi-start parameter tracking
  if (!is.null(fit$multistart) && nrow(fit$multistart) >= 3) {
    fwrite(fit$multistart, sprintf("%s_multistart_params.csv", out_prefix))
    # Correlation matrix across starts
    ms_num <- fit$multistart[, .(J, rho, pH, pICU, alphaR, rt_corr, prev_corr)]
    ms_cor <- cor(ms_num, use = "pairwise.complete.obs")
    fwrite(as.data.table(ms_cor, keep.rownames = "variable"),
           sprintf("%s_multistart_correlation.csv", out_prefix))
  } else if (!is.null(fit$multistart)) {
    fwrite(fit$multistart, sprintf("%s_multistart_params.csv", out_prefix))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: SOURCE ABLATION (ALL 15 SUBSETS) PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 2: Source ablation (15 subsets) per wave\n")
cat("================================================================\n")

ablation_all <- list()
for (wname in names(WAVE_DEFS)) {
  if (!(wname %in% ABLATION_WAVES)) {
    cat(sprintf("\n  === %s === SKIPPED (not in ABLATION_WAVES)\n", WAVE_DEFS[[wname]]$label))
    next
  }
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  abl_rows <- list()
  abl_beta_trajs <- list()  # For TikZ beta comparison figure
  for (si in seq_along(ALL_SUBSETS)) {
    sources <- ALL_SUBSETS[[si]]
    src_label <- paste(sources, collapse = "+")
    w <- subset_to_weights(sources)

    cat(sprintf("    [%2d/15] %s ...", si, src_label))
    fit <- tryCatch(
      fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                  wc = w["cases"], wh = w["hosp"],
                  wi = w["icu"], wr = w["radar"],
                  seed = 100 * wdef$id + si, verbose = FALSE),
      error = function(e) { cat(" FAILED\n"); NULL }
    )

    if (!is.null(fit)) {
      m <- fit$main
      abl_rows[[length(abl_rows) + 1]] <- data.table(
        wave          = wdef$label_short,
        sources       = src_label,
        w_cases       = w["cases"],
        w_hosp        = w["hosp"],
        w_icu         = w["icu"],
        w_radar       = w["radar"],
        J             = fit$J,
        n_converged   = fit$n_converged,
        rho           = m$rho,
        pH            = m$pH,
        pICU          = m$pICU,
        alphaR        = m$alphaR,
        rt_corr       = m$rt_corr,
        rt_spearman   = m$rt_spearman,
        rt_rmse       = m$rt_rmse,
        prev_corr     = m$prev_corr,
        prev_spearman = m$prev_spearman,
        cases_rmse    = m$cases_rmse,
        ode_rmse      = m$ode_rmse,
        R_end         = m$R_end,
        beta_min      = m$beta_range[1],
        beta_max      = m$beta_range[2],
        rt_min        = m$Rt_range[1],
        rt_max        = m$Rt_range[2],
        threshold_agree = m$threshold_agree
      )

      # Save beta trajectory for single-source and all-4 configs
      if (length(sources) == 1 || length(sources) == 4) {
        abl_beta_trajs[[length(abl_beta_trajs) + 1]] <- data.table(
          day     = fit$tv,
          sources = src_label,
          beta    = m$beta,
          Rt      = m$Rt
        )
      }

      cat(sprintf(" J=%.1f Rt=%.3f Prev=%.3f\n",
                  fit$J,
                  ifelse(is.na(m$rt_corr), NA, m$rt_corr),
                  ifelse(is.na(m$prev_corr), NA, m$prev_corr)))
    }
  }

  abl_dt <- rbindlist(abl_rows)
  ablation_all[[wname]] <- abl_dt
  fwrite(abl_dt, sprintf("%s_ablation_results.csv", wdef$label_short))

  # Save ablation beta trajectories for TikZ figure
  if (length(abl_beta_trajs) > 0) {
    fwrite(rbindlist(abl_beta_trajs),
           sprintf("%s_tikz_beta_ablation.csv", wdef$label_short))
  }

  cat(sprintf("  Saved: %s_ablation_results.csv + beta trajectories\n", wdef$label_short))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: SHAPLEY VALUES AND INTERACTION INDICES PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 3: Shapley values and interaction indices per wave\n")
cat("================================================================\n")

compute_shapley <- function(abl_dt, metric_col) {
  n_players <- 4
  sources <- SOURCE_NAMES
  phi <- numeric(n_players)
  names(phi) <- sources

  # Build lookup: source_set -> metric value
  get_v <- function(src_set) {
    if (length(src_set) == 0) return(0)
    key <- paste(sort(src_set), collapse = "+")
    row <- abl_dt[sources == key]
    if (nrow(row) == 0) return(NA_real_)
    row[[metric_col]][1]
  }

  for (i in seq_along(sources)) {
    player <- sources[i]
    others <- sources[-i]
    total <- 0

    for (s_size in 0:(n_players - 1)) {
      if (s_size == 0) {
        coalitions <- list(character(0))
      } else {
        coalitions <- combn(others, s_size, simplify = FALSE)
      }
      weight <- factorial(s_size) * factorial(n_players - s_size - 1) /
                factorial(n_players)

      for (S in coalitions) {
        v_with    <- get_v(c(S, player))
        v_without <- get_v(S)
        if (is.finite(v_with) && is.finite(v_without)) {
          total <- total + weight * (v_with - v_without)
        }
      }
    }
    phi[i] <- total
  }
  phi
}

compute_interactions <- function(abl_dt, metric_col) {
  n_players <- 4
  sources <- SOURCE_NAMES
  pairs <- combn(sources, 2, simplify = FALSE)
  result <- data.table(
    source_i = character(0), source_j = character(0), interaction = numeric(0)
  )

  get_v <- function(src_set) {
    if (length(src_set) == 0) return(0)
    key <- paste(sort(src_set), collapse = "+")
    row <- abl_dt[sources == key]
    if (nrow(row) == 0) return(NA_real_)
    row[[metric_col]][1]
  }

  for (pair in pairs) {
    i <- pair[1]; j <- pair[2]
    others <- setdiff(sources, pair)
    total <- 0

    for (s_size in 0:(n_players - 2)) {
      if (s_size == 0) {
        coalitions <- list(character(0))
      } else {
        coalitions <- combn(others, s_size, simplify = FALSE)
      }
      weight <- factorial(s_size) * factorial(n_players - s_size - 2) /
                factorial(n_players - 1)

      for (S in coalitions) {
        v_ij <- get_v(c(S, i, j))
        v_i  <- get_v(c(S, i))
        v_j  <- get_v(c(S, j))
        v_0  <- get_v(S)
        if (all(is.finite(c(v_ij, v_i, v_j, v_0)))) {
          total <- total + weight * (v_ij - v_i - v_j + v_0)
        }
      }
    }
    result <- rbind(result, data.table(source_i = i, source_j = j,
                                       interaction = total))
  }
  result
}

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  abl_dt <- ablation_all[[wname]]
  if (is.null(abl_dt) || nrow(abl_dt) == 0) next

  cat(sprintf("\n  === %s ===\n", wdef$label))

  phi_rt   <- compute_shapley(abl_dt, "rt_corr")
  phi_prev <- compute_shapley(abl_dt, "prev_corr")

  shap_dt <- data.table(
    source     = SOURCE_NAMES,
    phi_rt     = phi_rt,
    phi_prev   = phi_prev
  )
  fwrite(shap_dt, sprintf("%s_shapley_values.csv", wdef$label_short))

  int_rt   <- compute_interactions(abl_dt, "rt_corr")
  int_prev <- compute_interactions(abl_dt, "prev_corr")
  int_dt <- merge(int_rt, int_prev, by = c("source_i", "source_j"),
                  suffixes = c("_rt", "_prev"))
  fwrite(int_dt, sprintf("%s_shapley_interactions.csv", wdef$label_short))

  cat("  Shapley Rt:  ", paste(sprintf("%s=%.3f", names(phi_rt), phi_rt), collapse = ", "), "\n")
  cat("  Shapley Prev:", paste(sprintf("%s=%.3f", names(phi_prev), phi_prev), collapse = ", "), "\n")

  # --- Leave-one-out ---
  all4 <- abl_dt[sources == paste(SOURCE_NAMES, collapse = "+")]
  if (nrow(all4) > 0) {
    loo_rows <- list()
    for (src in SOURCE_NAMES) {
      remaining <- setdiff(SOURCE_NAMES, src)
      key <- paste(sort(remaining), collapse = "+")
      row3 <- abl_dt[sources == key]
      if (nrow(row3) > 0) {
        loo_rows[[length(loo_rows) + 1]] <- data.table(
          dropped     = src,
          remaining   = key,
          rt_all      = all4$rt_corr[1],
          rt_without  = row3$rt_corr[1],
          rt_drop     = all4$rt_corr[1] - row3$rt_corr[1],
          prev_all    = all4$prev_corr[1],
          prev_without = row3$prev_corr[1],
          prev_drop   = all4$prev_corr[1] - row3$prev_corr[1]
        )
      }
    }
    loo_dt <- rbindlist(loo_rows)
    fwrite(loo_dt, sprintf("%s_leave_one_out.csv", wdef$label_short))
  }

  # --- Superadditivity ---
  singles_rt   <- sapply(SOURCE_NAMES, function(s) {
    r <- abl_dt[sources == s]; if (nrow(r) > 0) r$rt_corr[1] else NA
  })
  singles_prev <- sapply(SOURCE_NAMES, function(s) {
    r <- abl_dt[sources == s]; if (nrow(r) > 0) r$prev_corr[1] else NA
  })

  if (nrow(all4) > 0) {
    super_dt <- data.table(
      metric     = c("rt_corr", "prev_corr"),
      v_all      = c(all4$rt_corr[1], all4$prev_corr[1]),
      v_sum      = c(sum(singles_rt, na.rm = TRUE), sum(singles_prev, na.rm = TRUE)),
      ratio      = c(all4$rt_corr[1] / sum(singles_rt, na.rm = TRUE),
                     all4$prev_corr[1] / sum(singles_prev, na.rm = TRUE)),
      best_single = c(max(singles_rt, na.rm = TRUE), max(singles_prev, na.rm = TRUE)),
      redundancy = c(1 - all4$rt_corr[1] / sum(singles_rt, na.rm = TRUE),
                     1 - all4$prev_corr[1] / sum(singles_prev, na.rm = TRUE))
    )
    fwrite(super_dt, sprintf("%s_superadditivity.csv", wdef$label_short))
  }

  cat(sprintf("  Saved: %s_shapley_values/interactions/leave_one_out/superadditivity.csv\n",
              wdef$label_short))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: INITIALIZATION SENSITIVITY PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 4: Initialization sensitivity per wave\n")
cat("================================================================\n")

init_methods <- c("hospital", "constant", "random_smooth")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  init_rows <- list()
  init_trajs <- list()

  for (mi in seq_along(init_methods)) {
    meth <- init_methods[mi]
    cat(sprintf("    init=%s ...", meth))

    if (meth == "hospital" && !is.null(wave_fits[[wname]])) {
      fit <- wave_fits[[wname]]
    } else {
      fit <- tryCatch(
        fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                    init_method = meth,
                    seed = 200 * wdef$id + mi,
                    verbose = FALSE),
        error = function(e) { cat(" FAILED\n"); NULL }
      )
    }

    if (!is.null(fit)) {
      m <- fit$main
      # Beta initialization deviation analysis
      b_abs_diff <- fit$main$beta - fit$beta_init
      b_rel_diff <- 100 * b_abs_diff / pmax(fit$beta_init, 1e-8)
      b_corr     <- cor(fit$beta_init, fit$main$beta)

      init_rows[[length(init_rows) + 1]] <- data.table(
        wave        = wdef$label_short,
        init_method = meth,
        J           = fit$J,
        rho         = m$rho,
        pH          = m$pH,
        pICU        = m$pICU,
        alphaR      = m$alphaR,
        rt_corr     = m$rt_corr,
        prev_corr   = m$prev_corr,
        beta_min    = m$beta_range[1],
        beta_max    = m$beta_range[2],
        R_end       = m$R_end,
        beta_init_fit_corr   = b_corr,
        beta_mean_abs_diff   = mean(abs(b_abs_diff)),
        beta_median_rel_pct  = median(abs(b_rel_diff), na.rm = TRUE)
      )
      init_trajs[[length(init_trajs) + 1]] <- data.table(
        date        = fit$dates,
        init_method = meth,
        beta_init   = fit$beta_init,
        beta        = fit$main$beta,
        beta_diff   = b_abs_diff,
        Rt          = fit$main$Rt
      )
      cat(sprintf(" J=%.1f Rt=%.3f (beta corr=%.3f, median|Δ%%|=%.1f%%)\n",
                  fit$J,
                  ifelse(is.na(m$rt_corr), NA, m$rt_corr),
                  b_corr,
                  median(abs(b_rel_diff), na.rm = TRUE)))
    }
  }

  if (length(init_rows) > 0) {
    fwrite(rbindlist(init_rows),
           sprintf("%s_init_sensitivity.csv", wdef$label_short))
    fwrite(rbindlist(init_trajs),
           sprintf("%s_init_sensitivity_trajectories.csv", wdef$label_short))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: START-DATE SENSITIVITY (WAVE 1 ONLY)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 5: Start-date sensitivity (wave 1)\n")
cat("================================================================\n")

start_offsets <- c(0, 7, 14)
start_rows <- list()
start_trajs <- list()

for (offset in start_offsets) {
  ds <- WAVE_DEFS$wave1$fit_start + offset
  de <- WAVE_DEFS$wave1$fit_end
  lbl <- sprintf("start+%dd", offset)
  cat(sprintf("  %s ...", lbl))

  if (offset == 0 && !is.null(wave_fits$wave1)) {
    fit <- wave_fits$wave1
  } else {
    fit <- tryCatch(
      fit_gp_seir(FULL_DF, ds, de, seed = 300 + offset,
                  n_cold = max(N_COLD_STARTS, 5), verbose = FALSE),
      error = function(e) { cat(" FAILED\n"); NULL }
    )
  }

  if (!is.null(fit)) {
    m <- fit$main
    start_rows[[length(start_rows) + 1]] <- data.table(
      start_label = lbl, date_start = as.character(ds),
      n_days = fit$n, J = fit$J,
      rt_corr = m$rt_corr, prev_corr = m$prev_corr,
      rt_min = m$Rt_range[1], rt_max = m$Rt_range[2],
      n_near = fit$n_converged, n_total = length(fit$all_J)
    )
    start_trajs[[length(start_trajs) + 1]] <- data.table(
      date = fit$dates, start_label = lbl,
      beta = fit$main$beta, Rt = fit$main$Rt,
      Rt_rivm = fit$raw$rt_rivm
    )
    cat(sprintf(" J=%.1f Rt=%.3f\n", fit$J,
                ifelse(is.na(m$rt_corr), NA, m$rt_corr)))
  }
}

if (length(start_rows) > 0) {
  fwrite(rbindlist(start_rows), "wave1_startdate_sensitivity.csv")
  fwrite(rbindlist(start_trajs), "wave1_startdate_sensitivity_trajectories.csv")
  cat("  Saved: wave1_startdate_sensitivity*.csv\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: PROFILE LIKELIHOOD PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 6: Profile likelihood per wave\n")
cat("================================================================\n")

for (wname in c("wave1", "wave2")) {
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  fit0 <- wave_fits[[wname]]
  if (is.null(fit0)) next

  m0 <- fit0$main
  profile_rows <- list()

  # Profile each scaling parameter
  param_specs <- list(
    list(name = "rho",    lo = BOUNDS$rho[1],    hi = BOUNDS$rho[2]),
    list(name = "pH",     lo = BOUNDS$pH[1],     hi = BOUNDS$pH[2]),
    list(name = "pICU",   lo = BOUNDS$pICU[1],   hi = BOUNDS$pICU[2]),
    list(name = "alphaR", lo = BOUNDS$alphaR[1], hi = BOUNDS$alphaR[2])
  )

  for (ps in param_specs) {
    cat(sprintf("    Profiling %s (%d points)...", ps$name, PROFILE_NPTS))

    lo_margin <- ps$lo + 0.05 * (ps$hi - ps$lo)
    hi_margin <- ps$hi - 0.05 * (ps$hi - ps$lo)
    grid_vals <- seq(lo_margin, hi_margin, length.out = PROFILE_NPTS)

    for (gv in grid_vals) {
      # Run a single-start fit with this parameter pinned
      fit_p <- tryCatch({
        # Pin the parameter by setting tight bounds
        bounds_mod <- BOUNDS
        bounds_mod[[ps$name]] <- c(gv - 1e-6, gv + 1e-6)

        # Refit with pinned parameter (using 2 warm + 1 cold)
        fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                    seed = 400 * wdef$id + round(gv * 1000),
                    n_warm = 2, n_cold = 1, verbose = FALSE,
                    bounds_override = bounds_mod)
      }, error = function(e) NULL)

      if (!is.null(fit_p)) {
        profile_rows[[length(profile_rows) + 1]] <- data.table(
          wave       = wdef$label_short,
          parameter  = ps$name,
          fixed_value = gv,
          J          = fit_p$J,
          rt_corr    = fit_p$main$rt_corr,
          prev_corr  = fit_p$main$prev_corr,
          rho        = fit_p$main$rho,
          pH         = fit_p$main$pH,
          pICU       = fit_p$main$pICU,
          alphaR     = fit_p$main$alphaR,
          R_end      = fit_p$main$R_end
        )
      }
    }
    cat(" done\n")
  }

  if (length(profile_rows) > 0) {
    profile_dt <- rbindlist(profile_rows)
    fwrite(profile_dt, sprintf("%s_profile_likelihood.csv", wdef$label_short))

    # Split into per-parameter TikZ CSVs (paper's pgfplots reference these)
    for (pname in unique(profile_dt$parameter)) {
      pdt <- profile_dt[parameter == pname]
      fwrite(pdt, sprintf("%s_tikz_profile_%s.csv", wdef$label_short, pname))
    }
    cat(sprintf("  Saved: %s_profile_likelihood.csv + per-parameter TikZ CSVs\n",
                wdef$label_short))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11: SIGMA/GAMMA SENSITIVITY PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 7: Sigma/gamma sensitivity per wave\n")
cat("================================================================\n")

sg_configs <- list(
  list(label = "default (sigma=1/5.5, gamma=1/9.5)", sig = 1/5.5, gam = 1/9.5),
  list(label = "sigma=1/4 (shorter latent)",         sig = 1/4,   gam = 1/9.5),
  list(label = "sigma=1/7 (longer latent)",          sig = 1/7,   gam = 1/9.5),
  list(label = "gamma=1/7 (shorter infectious)",     sig = 1/5.5, gam = 1/7),
  list(label = "gamma=1/12 (longer infectious)",     sig = 1/5.5, gam = 1/12)
)

for (wname in c("wave1", "wave2")) {
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  sg_rows <- list()
  for (cfg in sg_configs) {
    cat(sprintf("    %s ...", cfg$label))
    fit_sg <- tryCatch(
      fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                  sig = cfg$sig, gam = cfg$gam,
                  seed = 500 * wdef$id, n_warm = 3, n_cold = 2,
                  verbose = FALSE),
      error = function(e) { cat(" FAILED\n"); NULL }
    )
    if (!is.null(fit_sg)) {
      m <- fit_sg$main
      sg_rows[[length(sg_rows) + 1]] <- data.table(
        wave     = wdef$label_short,
        label    = cfg$label,
        sigma    = cfg$sig,
        gamma    = cfg$gam,
        J        = fit_sg$J,
        rho      = m$rho,
        pH       = m$pH,
        pICU     = m$pICU,
        alphaR   = m$alphaR,
        rt_corr  = m$rt_corr,
        prev_corr = m$prev_corr,
        R_end    = m$R_end
      )
      cat(sprintf(" J=%.1f Rt=%.3f\n", fit_sg$J,
                  ifelse(is.na(m$rt_corr), NA, m$rt_corr)))
    }
  }

  if (length(sg_rows) > 0) {
    fwrite(rbindlist(sg_rows),
           sprintf("%s_sigma_gamma_sensitivity.csv", wdef$label_short))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11b: W_DATA SENSITIVITY PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 7b: W_DATA (penalty weight) sensitivity per wave\n")
cat("================================================================\n")

WDATA_VALUES <- c(0.01, 1, 100, 10000, 1e6, 1e8)

for (wname in c("wave1", "wave2")) {
  wdef <- WAVE_DEFS[[wname]]
  cat(sprintf("\n  === %s ===\n", wdef$label))

  wd_rows <- list()
  for (wd in WDATA_VALUES) {
    # Compute effective ODE-to-data ratio
    eff_ratio <- sqrt(ETA / (wd * wsd(log1p(pmax(0, as.numeric(
      copy(FULL_DF)[date >= wdef$fit_start & date <= wdef$fit_end]$cases_raw
    ))))))

    cat(sprintf("    W_DATA=%.0e (ratio=%.0fx) ...", wd, eff_ratio))
    fit_wd <- tryCatch(
      fit_gp_seir(FULL_DF, wdef$fit_start, wdef$fit_end,
                  w_data_scale = wd,
                  seed = 600 * wdef$id + log10(max(wd, 1e-2)),
                  n_warm = 3, n_cold = 2, verbose = FALSE),
      error = function(e) { cat(" FAILED\n"); NULL }
    )

    if (!is.null(fit_wd)) {
      m <- fit_wd$main
      wd_rows[[length(wd_rows) + 1]] <- data.table(
        wave       = wdef$label_short,
        W_DATA     = wd,
        eff_ratio  = eff_ratio,
        J          = fit_wd$J,
        rt_corr    = m$rt_corr,
        rt_spearman = m$rt_spearman,
        rt_rmse    = m$rt_rmse,
        prev_corr  = m$prev_corr,
        prev_spearman = m$prev_spearman,
        cases_rmse = m$cases_rmse,
        ode_rmse   = m$ode_rmse,
        beta_min   = m$beta_range[1],
        beta_max   = m$beta_range[2],
        R_end      = m$R_end
      )
      cat(sprintf(" J=%.1f Rt=%.3f Prev=%.3f\n",
                  fit_wd$J,
                  ifelse(is.na(m$rt_corr), NA, m$rt_corr),
                  ifelse(is.na(m$prev_corr), NA, m$prev_corr)))
    }
  }

  if (length(wd_rows) > 0) {
    fwrite(rbindlist(wd_rows),
           sprintf("%s_wdata_sensitivity.csv", wdef$label_short))
    cat(sprintf("  Saved: %s_wdata_sensitivity.csv\n", wdef$label_short))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12: RESIDUAL DIAGNOSTICS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 8: Residual diagnostics per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  diag_rows <- list()
  streams <- list(
    list(name = "cases", resid = m$resid_cases),
    list(name = "hosp",  resid = m$resid_hosp),
    list(name = "icu",   resid = m$resid_icu),
    list(name = "radar", resid = m$resid_radar)
  )

  for (s in streams) {
    r <- s$resid[is.finite(s$resid)]
    if (length(r) < 10) next

    sw_test <- tryCatch(shapiro.test(r[1:min(5000, length(r))]),
                        error = function(e) list(statistic = NA, p.value = NA))
    lb_test <- tryCatch(Box.test(r, lag = 10, type = "Ljung-Box"),
                        error = function(e) list(statistic = NA, p.value = NA))
    acf_vals <- tryCatch(acf(r, lag.max = 5, plot = FALSE)$acf[-1],
                         error = function(e) rep(NA, 5))

    skewness_val <- (mean((r - mean(r))^3)) / (sd(r)^3)
    kurtosis_val <- (mean((r - mean(r))^4)) / (sd(r)^4) - 3

    diag_rows[[length(diag_rows) + 1]] <- data.table(
      stream      = s$name,
      n           = length(r),
      mean        = mean(r),
      sd          = sd(r),
      skewness    = skewness_val,
      kurtosis    = kurtosis_val,
      acf_lag1    = acf_vals[1],
      acf_lag2    = acf_vals[2],
      acf_lag3    = acf_vals[3],
      shapiro_W   = as.numeric(sw_test$statistic),
      shapiro_p   = sw_test$p.value,
      ljungbox_stat = as.numeric(lb_test$statistic),
      ljungbox_p  = lb_test$p.value
    )
  }

  if (length(diag_rows) > 0) {
    fwrite(rbindlist(diag_rows),
           sprintf("%s_residual_diagnostics.csv", wdef$label_short))
  }

  # TikZ residual CSVs
  resid_dt <- data.table(
    day = fit$tv,
    resid_cases = m$resid_cases,
    resid_hosp  = m$resid_hosp,
    resid_icu   = m$resid_icu,
    resid_radar = m$resid_radar
  )
  fwrite(resid_dt, sprintf("%s_tikz_residuals.csv", wdef$label_short))

  # ACF data for cases and hospital (for paper figures)
  for (sname in c("cases", "hosp")) {
    r <- if (sname == "cases") m$resid_cases else m$resid_hosp
    r <- r[is.finite(r)]
    if (length(r) > 20) {
      acf_obj <- tryCatch(acf(r, lag.max = 30, plot = FALSE),
                          error = function(e) NULL)
      if (!is.null(acf_obj)) {
        fwrite(data.table(lag = acf_obj$lag[-1], acf = acf_obj$acf[-1]),
               sprintf("%s_tikz_acf_%s.csv", wdef$label_short, sname))
      }
      # Q-Q data
      qq <- qqnorm(r, plot.it = FALSE)
      fwrite(data.table(theoretical = qq$x, sample = qq$y),
             sprintf("%s_tikz_qq_%s.csv", wdef$label_short, sname))
    }
  }
}

cat("  Saved: residual diagnostics per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13: CROSS-CORRELATIONS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 9: Cross-correlations per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next

  streams <- list(
    cases = fit$raw$cases_raw,
    hosp  = fit$raw$hosp_raw,
    icu   = fit$raw$icu_raw,
    radar = fit$raw$radar_frac
  )

  cc_rows <- list()
  pairs <- combn(names(streams), 2, simplify = FALSE)
  for (pair in pairs) {
    x <- streams[[pair[1]]]
    y <- streams[[pair[2]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 20) next

    lag0 <- cor(x[ok], y[ok])

    # Find best lag via CCF
    ccf_obj <- tryCatch(
      ccf(x[ok], y[ok], lag.max = 14, plot = FALSE),
      error = function(e) NULL
    )
    if (!is.null(ccf_obj)) {
      best_idx <- which.max(abs(ccf_obj$acf))
      cc_rows[[length(cc_rows) + 1]] <- data.table(
        stream_i   = pair[1],
        stream_j   = pair[2],
        pearson_lag0 = lag0,
        best_lag   = ccf_obj$lag[best_idx],
        best_ccf   = ccf_obj$acf[best_idx]
      )
    }
  }

  if (length(cc_rows) > 0) {
    fwrite(rbindlist(cc_rows),
           sprintf("%s_cross_correlations.csv", wdef$label_short))
  }
}

cat("  Saved: cross-correlations per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 14: ODE COMPLIANCE PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 10: ODE compliance per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  ode_dt <- data.table(
    compartment = c("S", "E", "I", "R"),
    rmse        = c(sqrt(mean(m$ode_S^2)), sqrt(mean(m$ode_E^2)),
                    sqrt(mean(m$ode_I^2)), sqrt(mean(m$ode_R^2))),
    max_abs     = c(max(abs(m$ode_S)), max(abs(m$ode_E)),
                    max(abs(m$ode_I)), max(abs(m$ode_R))),
    mean        = c(mean(m$ode_S), mean(m$ode_E),
                    mean(m$ode_I), mean(m$ode_R))
  )
  fwrite(ode_dt, sprintf("%s_ode_compliance.csv", wdef$label_short))
}

cat("  Saved: ODE compliance per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 15: STREAM FIT METRICS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 11: Stream fit metrics per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  sfm_rows <- list()
  stream_specs <- list(
    list(name = "cases", obs = fit$raw$cases_raw, pred = m$mu_cases, scale = "raw"),
    list(name = "hosp",  obs = fit$raw$hosp_raw,  pred = m$mu_hosp,  scale = "raw"),
    list(name = "icu",   obs = fit$raw$icu_raw,   pred = m$mu_icu,   scale = "raw"),
    list(name = "radar", obs = fit$raw$radar_frac, pred = m$mu_radar, scale = "fraction")
  )

  for (ss in stream_specs) {
    ok <- is.finite(ss$obs) & is.finite(ss$pred)
    if (sum(ok) < 5) next
    x <- ss$obs[ok]; y <- ss$pred[ok]
    sfm_rows[[length(sfm_rows) + 1]] <- data.table(
      stream   = ss$name,
      scale    = ss$scale,
      pearson  = cor(x, y),
      spearman = cor(x, y, method = "spearman"),
      RMSE     = sqrt(mean((x - y)^2)),
      MAE      = mean(abs(x - y)),
      bias     = mean(y - x),
      n_obs    = sum(ok)
    )
  }

  if (length(sfm_rows) > 0) {
    fwrite(rbindlist(sfm_rows),
           sprintf("%s_stream_fit_metrics.csv", wdef$label_short))
  }
}

cat("  Saved: stream fit metrics per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 16: Rt LAG ANALYSIS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 12: Rt lag analysis per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next

  rt_model <- fit$main$Rt
  rt_rivm  <- fit$raw$rt_rivm
  n <- length(rt_model)

  lag_rows <- list()
  for (lag in -14:14) {
    if (lag >= 0) {
      idx_m <- 1:(n - lag)
      idx_r <- (1 + lag):n
    } else {
      idx_m <- (1 - lag):n
      idx_r <- 1:(n + lag)
    }
    ok <- is.finite(rt_model[idx_m]) & is.finite(rt_rivm[idx_r])
    if (sum(ok) > 10) {
      lag_rows[[length(lag_rows) + 1]] <- data.table(
        lag     = lag,
        pearson  = cor(rt_model[idx_m][ok], rt_rivm[idx_r][ok]),
        spearman = cor(rt_model[idx_m][ok], rt_rivm[idx_r][ok], method = "spearman"),
        rmse    = sqrt(mean((rt_model[idx_m][ok] - rt_rivm[idx_r][ok])^2))
      )
    }
  }

  if (length(lag_rows) > 0) {
    fwrite(rbindlist(lag_rows),
           sprintf("%s_rt_lag_analysis.csv", wdef$label_short))
  }
}

cat("  Saved: Rt lag analysis per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 17: Rt THRESHOLD CROSSINGS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 13: Rt threshold crossings per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next

  Rt <- fit$main$Rt
  n <- length(Rt)
  cross_rows <- list()

  for (i in 2:n) {
    if ((Rt[i-1] >= 1 && Rt[i] < 1) || (Rt[i-1] < 1 && Rt[i] >= 1)) {
      cross_rows[[length(cross_rows) + 1]] <- data.table(
        day       = fit$tv[i],
        date      = fit$dates[i],
        direction = if (Rt[i] < 1) "falling" else "rising",
        Rt_before = Rt[i-1],
        Rt_after  = Rt[i]
      )
    }
  }

  if (length(cross_rows) > 0) {
    fwrite(rbindlist(cross_rows),
           sprintf("%s_rt_crossings.csv", wdef$label_short))
  }
}

cat("  Saved: Rt crossings per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 18: COMPREHENSIVE VALIDATION METRICS PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 14: Comprehensive validation metrics per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  ok_rt <- is.finite(fit$raw$rt_rivm) & is.finite(m$Rt)
  ok_prev <- is.finite(fit$raw$prev_rivm) & is.finite(N_POP * m$I)

  rt_comp <- data.table(
    metric = c("pearson", "spearman", "R2", "RMSE", "bias",
               "threshold_agree", "n_days"),
    value  = c(m$rt_corr, m$rt_spearman,
               ifelse(sum(ok_rt) > 3,
                      1 - sum((m$Rt[ok_rt] - fit$raw$rt_rivm[ok_rt])^2) /
                        sum((fit$raw$rt_rivm[ok_rt] - mean(fit$raw$rt_rivm[ok_rt]))^2),
                      NA),
               m$rt_rmse, m$rt_bias, m$threshold_agree, sum(ok_rt))
  )
  fwrite(rt_comp, sprintf("%s_rt_comprehensive.csv", wdef$label_short))

  prev_comp <- data.table(
    metric = c("pearson", "spearman", "R2", "RMSE"),
    value  = c(m$prev_corr, m$prev_spearman,
               ifelse(sum(ok_prev) > 3, {
                 I_count <- N_POP * m$I
                 1 - sum((I_count[ok_prev] - fit$raw$prev_rivm[ok_prev])^2) /
                   sum((fit$raw$prev_rivm[ok_prev] -
                          mean(fit$raw$prev_rivm[ok_prev]))^2)
               }, NA),
               m$prev_rmse)
  )
  fwrite(prev_comp, sprintf("%s_prev_comprehensive.csv", wdef$label_short))
}

cat("  Saved: comprehensive validation metrics per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 19: CROSS-WAVE COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 15: Cross-wave comparison\n")
cat("================================================================\n")

comparison_rows <- list()
for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  comparison_rows[[length(comparison_rows) + 1]] <- data.table(
    wave           = wdef$label_short,
    wave_label     = wdef$label,
    n_days         = fit$n,
    J              = fit$J,
    rho            = m$rho,
    pH             = m$pH,
    pICU           = m$pICU,
    alphaR         = m$alphaR,
    rt_corr        = m$rt_corr,
    rt_spearman    = m$rt_spearman,
    rt_rmse        = m$rt_rmse,
    rt_bias        = m$rt_bias,
    prev_corr      = m$prev_corr,
    prev_spearman  = m$prev_spearman,
    beta_min       = m$beta_range[1],
    beta_max       = m$beta_range[2],
    rt_min         = m$Rt_range[1],
    rt_max         = m$Rt_range[2],
    R_end          = m$R_end,
    ode_rmse       = m$ode_rmse,
    mass_err       = m$mass_err,
    cases_rmse     = m$cases_rmse,
    hosp_rmse      = m$hosp_rmse,
    threshold_agree = m$threshold_agree,
    attack_rate    = m$R_end,
    peak_I_count   = max(N_POP * m$I),
    peak_I_date    = as.character(fit$dates[which.max(m$I)]),
    mean_rho_eff   = mean(m$rho_eff, na.rm = TRUE),
    pct_above_Rt1  = 100 * mean(m$Rt > 1),
    beta_init_corr      = cor(fit$beta_init, m$beta),
    beta_mean_abs_diff  = mean(abs(m$beta - fit$beta_init)),
    beta_median_rel_pct = median(abs(100 * (m$beta - fit$beta_init) /
                                       pmax(fit$beta_init, 1e-8)), na.rm = TRUE)
  )
}

if (length(comparison_rows) > 0) {
  comp_dt <- rbindlist(comparison_rows)
  fwrite(comp_dt, "cross_wave_comparison.csv")
  cat("  Saved: cross_wave_comparison.csv\n")

  cat("\n  === Cross-wave parameter comparison ===\n")
  for (i in 1:nrow(comp_dt)) {
    r <- comp_dt[i]
    cat(sprintf("  %s: rho=%.3f pH=%.4f pICU=%.4f alphaR=%.3f | Rt_r=%.3f Prev_r=%.3f\n",
                r$wave, r$rho, r$pH, r$pICU, r$alphaR, r$rt_corr, r$prev_corr))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 20: SUMMARY CSV PER WAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 16: Summary files per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main
  abl <- ablation_all[[wname]]

  summary_rows <- list()
  add_row <- function(cat, met, val) {
    summary_rows[[length(summary_rows) + 1]] <<- data.table(
      category = cat, metric = met, value = val
    )
  }

  add_row("main", "J", fit$J)
  add_row("main", "rho", m$rho)
  add_row("main", "pH", m$pH)
  add_row("main", "pICU", m$pICU)
  add_row("main", "alphaR", m$alphaR)
  add_row("main", "rt_pearson", m$rt_corr)
  add_row("main", "rt_spearman", m$rt_spearman)
  add_row("main", "rt_rmse", m$rt_rmse)
  add_row("main", "prev_pearson", m$prev_corr)
  add_row("main", "prev_spearman", m$prev_spearman)
  add_row("main", "R_end", m$R_end)
  add_row("main", "beta_min", m$beta_range[1])
  add_row("main", "beta_max", m$beta_range[2])
  add_row("main", "ode_rmse", m$ode_rmse)
  add_row("main", "cases_rmse", m$cases_rmse)
  add_row("main", "threshold_agreement", m$threshold_agree)
  add_row("main", "delay_kernel", as.numeric(USE_DELAY_KERNEL))
  add_row("main", "n_days", fit$n)
  add_row("main", "peak_I_count", max(N_POP * m$I))
  add_row("main", "peak_I_date", as.character(fit$dates[which.max(m$I)]))
  add_row("main", "mean_rho_eff", mean(m$rho_eff, na.rm = TRUE))
  add_row("main", "pct_above_Rt1", 100 * mean(m$Rt > 1))
  add_row("beta_init", "corr_init_fitted", cor(fit$beta_init, m$beta))
  add_row("beta_init", "mean_abs_diff", mean(abs(m$beta - fit$beta_init)))
  add_row("beta_init", "median_rel_pct", median(abs(100 * (m$beta - fit$beta_init) /
                                                       pmax(fit$beta_init, 1e-8)), na.rm = TRUE))

  # Shapley summary
  shap_file <- sprintf("%s_shapley_values.csv", wdef$label_short)
  if (file.exists(shap_file)) {
    shap <- fread(shap_file)
    for (i in 1:nrow(shap)) {
      add_row("shapley_rt", shap$source[i], shap$phi_rt[i])
      add_row("shapley_prev", shap$source[i], shap$phi_prev[i])
    }
  }

  # Ablation summary
  if (!is.null(abl) && nrow(abl) > 0) {
    singles <- abl[nchar(gsub("[^+]", "", sources)) == 0]
    if (nrow(singles) > 0) {
      add_row("ablation", "best_single_rt_source",
              singles$sources[which.max(singles$rt_corr)])
      add_row("ablation", "best_single_rt_value",
              max(singles$rt_corr, na.rm = TRUE))
      add_row("ablation", "best_single_prev_source",
              singles$sources[which.max(singles$prev_corr)])
      add_row("ablation", "best_single_prev_value",
              max(singles$prev_corr, na.rm = TRUE))
    }
  }

  fwrite(rbindlist(summary_rows),
         sprintf("%s_summary.csv", wdef$label_short))
}

cat("  Saved: summary per wave\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 21: ABLATION TikZ CSVs (for bar charts)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 17: Ablation TikZ CSVs\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  abl <- ablation_all[[wname]]
  if (is.null(abl) || nrow(abl) == 0) next

  # Sort by prevalence correlation for the bar chart
  abl_sorted <- abl[order(prev_corr)]
  fwrite(abl_sorted[, .(sources, rt_corr, prev_corr, J)],
         sprintf("%s_tikz_ablation.csv", wdef$label_short))

  # Ablation Rt chart (sorted by Rt)
  abl_rt_sorted <- abl[order(rt_corr)]
  fwrite(abl_rt_sorted[, .(sources, rt_corr, prev_corr, J)],
         sprintf("%s_tikz_ablation_rt.csv", wdef$label_short))
}

cat("  Saved: ablation TikZ CSVs\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 22: MAIN FIT PLOTS (PDF)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("STEP 18: Generating PDF plots per wave\n")
cat("================================================================\n")

for (wname in names(WAVE_DEFS)) {
  wdef <- WAVE_DEFS[[wname]]
  fit <- wave_fits[[wname]]
  if (is.null(fit)) next
  m <- fit$main

  pdf_file <- sprintf("%s_fit.pdf", wdef$label_short)
  pdf(pdf_file, width = 13, height = 13)
  par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))

  # SEIR compartments
  plot(fit$dates, m$S, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
       xlab = "Date", ylab = "Population fraction",
       main = paste0("SEIR Compartments - ", wdef$label))
  lines(fit$dates, m$E, col = "orange", lwd = 2)
  lines(fit$dates, m$I, col = "red", lwd = 2)
  lines(fit$dates, m$R, col = "darkgreen", lwd = 2)
  legend("right", c("S", "E", "I", "R"),
         col = c("blue", "orange", "red", "darkgreen"), lwd = 2, cex = 0.8)

  # Rt
  ylim_rt <- c(0, max(3, max(c(m$Rt, fit$raw$rt_rivm), na.rm = TRUE) * 1.1))
  plot(fit$dates, m$Rt, type = "l", col = "darkred", lwd = 2, ylim = ylim_rt,
       xlab = "Date", ylab = expression(R[t]),
       main = sprintf("Reproduction number (r=%.3f)",
                      ifelse(is.na(m$rt_corr), NA, m$rt_corr)))
  ok <- is.finite(fit$raw$rt_rivm)
  if (any(ok)) points(fit$dates[ok], fit$raw$rt_rivm[ok], col = "gray50", pch = 16, cex = 0.5)
  abline(h = 1, lty = 2, col = "gray40")
  legend("topright", c("Model", "RIVM"), col = c("darkred", "gray50"),
         lwd = c(2, NA), pch = c(NA, 16), cex = 0.8)

  # Cases
  plot(fit$dates, fit$raw$cases_raw, pch = 16, cex = 0.3, col = "gray60",
       xlab = "Date", ylab = "Daily cases", main = "Cases")
  lines(fit$dates, m$mu_cases, col = "steelblue", lwd = 2)

  # Hospital
  plot(fit$dates, fit$raw$hosp_raw, pch = 16, cex = 0.3, col = "gray60",
       xlab = "Date", ylab = "Admissions/day", main = "Hospital admissions")
  lines(fit$dates, m$mu_hosp, col = "darkorange", lwd = 2)

  # ICU
  plot(fit$dates, fit$raw$icu_raw, pch = 16, cex = 0.3, col = "gray60",
       xlab = "Date", ylab = "Admissions/day", main = "ICU admissions")
  lines(fit$dates, m$mu_icu, col = "purple", lwd = 2)

  # RADAR
  ok_r <- is.finite(fit$raw$radar_frac)
  if (any(ok_r)) {
    plot(fit$dates[ok_r], fit$raw$radar_frac[ok_r], pch = 16, cex = 0.3, col = "gray60",
         xlab = "Date", ylab = "I(t) fraction", main = "COVID RADAR")
    lines(fit$dates, m$mu_radar, col = "forestgreen", lwd = 2)
  } else {
    plot.new(); title("RADAR: no data")
  }

  # Beta
  plot(fit$dates, m$beta, type = "l", col = "purple", lwd = 2,
       xlab = "Date", ylab = expression(beta(t)), main = "Transmission rate")
  lines(fit$dates, fit$beta_init, col = "gray50", lwd = 1, lty = 2)
  legend("topright", c("Fitted", "Initial"), col = c("purple", "gray50"),
         lwd = c(2, 1), lty = c(1, 2), cex = 0.8)

  # Prevalence validation
  ok_p <- is.finite(fit$raw$prev_rivm)
  if (any(ok_p)) {
    plot(fit$dates[ok_p], fit$raw$prev_rivm[ok_p] / 1000, pch = 16, cex = 0.4, col = "gray50",
         xlab = "Date", ylab = "Prevalence (thousands)",
         main = sprintf("Prevalence (r=%.3f)",
                        ifelse(is.na(m$prev_corr), NA, m$prev_corr)))
    lines(fit$dates, N_POP * m$I / 1000, col = "blue", lwd = 2)
    legend("topright", c("RIVM", "Model"), col = c("gray50", "blue"),
           pch = c(16, NA), lwd = c(NA, 2), cex = 0.8)
  } else {
    plot.new(); title("Prevalence: no validation data")
  }

  dev.off()
  cat(sprintf("  Saved: %s\n", pdf_file))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 23: EPIDEMIOLOGICAL PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

epi_dt <- data.table(
  parameter = c("sigma", "gamma", "gen_time", "latent_period",
                "infectious_period", "N_pop",
                "delay_case_mean", "delay_hosp_mean", "delay_icu_mean",
                "use_delay_kernel"),
  value     = c(SIGMA, GAMMA, T_GEN, 1/SIGMA, 1/GAMMA, N_POP,
                DELAY_CASE_MEAN, DELAY_HOSP_MEAN, DELAY_ICU_MEAN,
                as.numeric(USE_DELAY_KERNEL))
)
fwrite(epi_dt, "epi_parameters.csv")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 24: SAVE WORKSPACE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("Saving workspace...\n")
save.image("gp_seir_wave_workspace.RData")
cat("Saved: gp_seir_wave_workspace.RData\n")

cat("\n================================================================\n")
cat("COMPLETE. Output files:\n")
cat("================================================================\n")
output_files <- list.files(pattern = "^(wave|cross|epi).*\\.(csv|pdf|RData)$")
for (f in sort(output_files)) cat(sprintf("  %s\n", f))
cat("================================================================\n")
