#!/usr/bin/env Rscript
##############################################################################
# extend_rho_tv_crossval.R
#
# HIGH-IMPACT EXTENSIONS — run after fit_gp_seir_analysis.R
#
# 1. Time-varying rho(t): adds ~17 B-spline coefficients for detection rate
# 2. Cross-validation: fit on days 1-300, predict days 301-366
# 3. EpiEstim comparison: Cori et al. (2013) Rt estimation
#
# Reads: gp_input_final_v2.csv, gp_seir_results.csv (for comparison)
# Produces: rho_tv_results.csv, crossval_results.csv, epiestim_rt.csv,
#           tikz_rho_tv.csv, tikz_crossval.csv, tikz_epiestim.csv,
#           gp_seir_plots_extensions.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(fda); library(minpack.lm); library(zoo)
})

cat("============================================================\n")
cat("GP-SEIR Extensions\n")
cat("============================================================\n\n")

# Load settings (same as main script)
DATA_FILE <- "gp_input_final_v2.csv"
DATE_START <- as.Date("2020-03-01"); DATE_END <- as.Date("2021-03-01")
N_POP <- 17400000; sigma <- 1/5.5; gamma <- 1/9.5
STATE_KNOT_DAYS <- 14; STATE_NORDER <- 4
BETA_KNOT_DAYS <- 28; BETA_NORDER <- 4
RHO_KNOT_DAYS <- 28; RHO_NORDER <- 4  # NEW: knots for rho(t)
ETA <- 1e6; KAPPA_BETA <- 1e3; KAPPA_BMAG <- 1e4; BETA_MAX <- 0.6; KAPPA_MASS <- 1e6
KAPPA_RHO <- 1e2  # roughness penalty on rho(t)
MAX_INNER <- 300; MAX_OUTER <- 40; STEP0 <- 1.0
bounds <- list(pH=c(0.002,0.05), pICU=c(0.0003,0.02), alphaR=c(0.10,1.00))

# Load data
df <- fread(DATA_FILE)[, date := as.Date(date)]
df <- df[date >= DATE_START & date <= DATE_END][order(date)]
n <- nrow(df); df[, t := as.numeric(date - DATE_START)]
tv <- df$t; rng <- range(tv); dt <- 1; dates <- df$date
cases_raw <- as.numeric(df$cases_raw); hosp_raw <- as.numeric(df$hosp_nice)
icu_raw <- as.numeric(df$icu_daily); radar_frac <- as.numeric(df$radar_I_frac_sm7)
rt_rivm <- as.numeric(df$rt_rivm); prev_rivm <- as.numeric(df$prev_json_avg)
y_cases <- log1p(pmax(0, cases_raw)); y_hosp <- log1p(pmax(0, hosp_raw))
y_icu <- log1p(pmax(0, icu_raw)); y_radar <- radar_frac
wsd <- function(y){s<-sd(y,na.rm=TRUE);if(!is.finite(s)||s<=0)s<-1;1/s^2}
wb_cases<-wsd(y_cases);wb_hosp<-wsd(y_hosp);wb_icu<-wsd(y_icu);wb_radar<-wsd(y_radar)

# Bases
mkb <- function(kd,no){kn<-unique(c(rng[1],seq(rng[1],rng[2],by=kd),rng[2]));create.bspline.basis(rng,length(kn)+no-2,no,kn)}
Bs <- mkb(STATE_KNOT_DAYS, STATE_NORDER)
Bb <- mkb(BETA_KNOT_DAYS, BETA_NORDER)
Br <- mkb(RHO_KNOT_DAYS, RHO_NORDER)  # NEW: rho basis
nbs <- Bs$nbasis; nbb <- Bb$nbasis; nbr <- Br$nbasis
Ph <- eval.basis(tv,Bs,0); dPh <- eval.basis(tv,Bs,1)
Pb <- eval.basis(tv,Bb,0); d2Pb <- eval.basis(tv,Bb,2)
Pr <- eval.basis(tv,Br,0); d2Pr <- eval.basis(tv,Br,2)  # NEW

cat(sprintf("  State: %d, Beta: %d, Rho: %d basis functions\n", nbs, nbb, nbr))

to_u <- function(x,lo,hi) qlogis((x-lo)/(hi-lo))
from_u <- function(u,lo,hi) lo+(hi-lo)*plogis(u)
seir_rhs <- function(Z,beta) cbind(-beta*exp(Z[,3]),beta*exp(Z[,1]+Z[,3]-Z[,2])-sigma,sigma*exp(Z[,2]-Z[,3])-gamma,gamma*exp(Z[,3]-Z[,4]))

# ══════════════════════════════════════════════════════════════════
# EXTENSION 1: TIME-VARYING rho(t)
#
# Instead of a single constant rho, we model:
#   rho(t) = 0.60 * logistic(Pr %*% a_rho)
# This gives rho(t) ∈ (0, 0.60) by construction.
# The outer parameters become: (a_beta, a_rho, pH, pICU, alphaR)
# ══════════════════════════════════════════════════════════════════

cat("\n=== Extension 1: Time-varying rho(t) ===\n")

# Initialize from constant-rho fit
old <- fread("gp_seir_results.csv")
Ci_tv <- smooth.basis(tv, cbind(log(old$S), log(old$E), log(old$I), log(old$R)),
                       fdPar(Bs, int2Lfd(2), 1e-4))$fd$coefs
ai_tv <- as.numeric(smooth.basis(tv, log(old$beta),
                                  fdPar(Bb, int2Lfd(2), 1e-2))$fd$coefs)
# Initialize rho spline at logit(0.20/0.60) ≈ -1.1
ri_tv <- rep(qlogis(0.20/0.60), nbr)

# Outer params: a_beta (nbb), a_rho (nbr), pH, pICU, alphaR
thi_tv <- c(ai_tv, ri_tv,
            to_u(0.012, bounds$pH[1], bounds$pH[2]),
            to_u(0.0025, bounds$pICU[1], bounds$pICU[2]),
            to_u(0.406, bounds$alphaR[1], bounds$alphaR[2]))

solve_inner_tv <- function(Cs, th) {
  av <- th[1:nbb]; ar <- th[(nbb+1):(nbb+nbr)]
  pH   <- from_u(th[nbb+nbr+1], bounds$pH[1], bounds$pH[2])
  pICU <- from_u(th[nbb+nbr+2], bounds$pICU[1], bounds$pICU[2])
  aR   <- from_u(th[nbb+nbr+3], bounds$alphaR[1], bounds$alphaR[2])
  beta <- exp(as.numeric(Pb %*% av))
  rho_t <- 0.60 * plogis(as.numeric(Pr %*% ar))  # time-varying rho

  fn <- function(cf) {
    C <- matrix(cf, nrow=nbs, ncol=4); Z <- Ph%*%C; dZ <- dPh%*%C
    S <- exp(Z[,1]); E <- exp(Z[,2]); I <- exp(Z[,3]); R <- exp(Z[,4])
    mc <- log1p(pmax(0, N_POP * rho_t * sigma * E))  # rho_t instead of rho
    mh <- log1p(pmax(0, N_POP*pH*gamma*I))
    mi <- log1p(pmax(0, N_POP*pICU*gamma*I))
    mr <- aR*I
    rc <- sqrt(wb_cases)*ifelse(is.finite(y_cases), y_cases-mc, 0)
    rh <- sqrt(wb_hosp)*ifelse(is.finite(y_hosp), y_hosp-mh, 0)
    ri <- sqrt(wb_icu)*ifelse(is.finite(y_icu), y_icu-mi, 0)
    rr <- sqrt(wb_radar)*ifelse(is.finite(y_radar), y_radar-mr, 0)
    rp <- sqrt(ETA*dt)*as.numeric(dZ - seir_rhs(Z, beta))
    rm <- sqrt(KAPPA_MASS)*(S+E+I+R-1)
    rb <- sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX, 0)
    rs_beta <- sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb %*% av)
    rs_rho  <- sqrt(KAPPA_RHO*dt)*as.numeric(d2Pr %*% ar)  # rho roughness
    c(rc, rh, ri, rr, rp, rm, rb, rs_beta, rs_rho)
  }
  r <- nls.lm(par=as.numeric(Cs), fn=fn,
              control=nls.lm.control(maxiter=MAX_INNER, ftol=1e-10, ptol=1e-10))
  list(C=matrix(r$par, nrow=nbs, ncol=4), J=sum(r$fvec^2))
}

cascade_tv <- function(th0, C0) {
  th <- th0; Cc <- C0; bJ <- Inf
  for (oi in 1:MAX_OUTER) {
    inn <- solve_inner_tv(Cc, th); Cc <- inn$C; Jc <- inn$J
    if (Jc < bJ) bJ <- Jc
    eps <- 1e-4; gr <- numeric(length(th))
    for (j in seq_along(th)) {
      tp <- th; tp[j] <- tp[j]+eps
      gr[j] <- (solve_inner_tv(Cc, tp)$J - Jc)/eps
    }
    gn <- sqrt(sum(gr^2))
    if (gn > 1e-12) {
      d <- -gr/gn; st <- STEP0; found <- FALSE
      for (ls in 1:10) {
        tt <- th+st*d; it <- solve_inner_tv(Cc, tt)
        if (it$J < Jc-1e-4*st*gn) {th<-tt; Cc<-it$C; found<-TRUE; break}
        st <- st*0.5
      }
      if (!found) th <- th+1e-4*d
    }
    cat(sprintf("  Outer %d: J=%.2f gn=%.4f\n", oi, bJ, gn))
    if (gn < 1e-5*max(1, sqrt(sum(th^2)))) break
  }
  list(theta=th, C=Cc, J=bJ)
}

set.seed(42)
cat("  Running time-varying rho(t) fit...\n")
res_tv <- tryCatch(cascade_tv(thi_tv, Ci_tv), error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(res_tv)) {
  th <- res_tv$theta; C <- res_tv$C
  av <- th[1:nbb]; ar <- th[(nbb+1):(nbb+nbr)]
  pH   <- from_u(th[nbb+nbr+1], bounds$pH[1], bounds$pH[2])
  pICU <- from_u(th[nbb+nbr+2], bounds$pICU[1], bounds$pICU[2])
  aR   <- from_u(th[nbb+nbr+3], bounds$alphaR[1], bounds$alphaR[2])
  beta <- exp(as.numeric(Pb %*% av))
  rho_t <- 0.60 * plogis(as.numeric(Pr %*% ar))
  Z <- Ph %*% C; S <- exp(Z[,1]); E <- exp(Z[,2]); I <- exp(Z[,3]); R <- exp(Z[,4])
  Rt <- beta * S / gamma

  ok_rt <- is.finite(rt_rivm); ok_prev <- is.finite(prev_rivm)
  rt_corr_tv <- cor(rt_rivm[ok_rt], Rt[ok_rt])
  prev_corr_tv <- if(any(ok_prev)) cor(prev_rivm[ok_prev], N_POP*I[ok_prev]) else NA

  # Case residual ACF
  mc_tv <- N_POP * rho_t * sigma * E
  resid_tv <- ifelse(is.finite(y_cases), y_cases - log1p(pmax(0, mc_tv)), NA)
  acf_tv <- acf(resid_tv[is.finite(resid_tv)], lag.max = 7, plot = FALSE)$acf[2]

  cat(sprintf("  rho(t) range: [%.3f, %.3f]\n", min(rho_t), max(rho_t)))
  cat(sprintf("  Rt corr: %.3f (constant-rho: 0.824)\n", rt_corr_tv))
  cat(sprintf("  Prev corr: %.3f (constant-rho: 0.847)\n", prev_corr_tv))
  cat(sprintf("  ACF(1): %.3f (constant-rho: 0.941)\n", acf_tv))

  fwrite(data.table(day = tv, rho_t = rho_t, beta = beta,
                    S = S, E = E, I = I, R = R, Rt = Rt,
                    mu_cases_tv = mc_tv, resid_cases_tv = resid_tv),
         "rho_tv_results.csv")
  fwrite(data.table(day = tv, rho_t = rho_t, rho_const = 0.201), "tikz_rho_tv.csv")
  cat("  Saved rho_tv_results.csv and tikz_rho_tv.csv\n")
}

# ══════════════════════════════════════════════════════════════════
# EXTENSION 2: CROSS-VALIDATION (fit 1-300, predict 301-366)
# ══════════════════════════════════════════════════════════════════

cat("\n=== Extension 2: Cross-validation ===\n")

# Build bases for training period only (days 0-299)
train_end <- 299
train_idx <- tv <= train_end
test_idx  <- tv > train_end

rng_train <- c(0, train_end)
mkb_tr <- function(kd, no) {
  kn <- unique(c(rng_train[1], seq(rng_train[1], rng_train[2], by=kd), rng_train[2]))
  create.bspline.basis(rng_train, length(kn)+no-2, no, kn)
}
Bs_tr <- mkb_tr(STATE_KNOT_DAYS, STATE_NORDER)
Bb_tr <- mkb_tr(BETA_KNOT_DAYS, BETA_NORDER)
nbs_tr <- Bs_tr$nbasis; nbb_tr <- Bb_tr$nbasis
tv_tr <- tv[train_idx]
Ph_tr <- eval.basis(tv_tr, Bs_tr, 0); dPh_tr <- eval.basis(tv_tr, Bs_tr, 1)
Pb_tr <- eval.basis(tv_tr, Bb_tr, 0); d2Pb_tr <- eval.basis(tv_tr, Bb_tr, 2)

cat(sprintf("  Training: days 0-%d (%d days), Test: days %d-%d (%d days)\n",
            train_end, sum(train_idx), train_end+1, max(tv), sum(test_idx)))

# Slice training data
yc_tr <- y_cases[train_idx]; yh_tr <- y_hosp[train_idx]
yi_tr <- y_icu[train_idx]; yr_tr <- y_radar[train_idx]
wbc_tr <- wsd(yc_tr); wbh_tr <- wsd(yh_tr)
wbi_tr <- wsd(yi_tr); wbr_tr <- wsd(yr_tr)

# Initialize from main fit (sliced to training period)
S_tr <- old$S[train_idx]; E_tr <- old$E[train_idx]
I_tr <- old$I[train_idx]; R_tr <- old$R[train_idx]
b_tr <- old$beta[train_idx]

bounds_cv <- list(rho=c(0.05,0.60), pH=c(0.002,0.05), pICU=c(0.0003,0.02), alphaR=c(0.10,1.00))

Ci_tr <- smooth.basis(tv_tr, cbind(log(S_tr),log(E_tr),log(I_tr),log(R_tr)),
                      fdPar(Bs_tr, int2Lfd(2), 1e-4))$fd$coefs
ai_tr <- as.numeric(smooth.basis(tv_tr, log(b_tr),
                                  fdPar(Bb_tr, int2Lfd(2), 1e-2))$fd$coefs)
thi_tr <- c(ai_tr,
            to_u(0.20, bounds_cv$rho[1], bounds_cv$rho[2]),
            to_u(0.012, bounds_cv$pH[1], bounds_cv$pH[2]),
            to_u(0.0025, bounds_cv$pICU[1], bounds_cv$pICU[2]),
            to_u(0.40, bounds_cv$alphaR[1], bounds_cv$alphaR[2]))

solve_inner_cv <- function(Cs, th) {
  av <- th[1:nbb_tr]
  rho  <- from_u(th[nbb_tr+1], bounds_cv$rho[1], bounds_cv$rho[2])
  pH   <- from_u(th[nbb_tr+2], bounds_cv$pH[1], bounds_cv$pH[2])
  pICU <- from_u(th[nbb_tr+3], bounds_cv$pICU[1], bounds_cv$pICU[2])
  aR   <- from_u(th[nbb_tr+4], bounds_cv$alphaR[1], bounds_cv$alphaR[2])
  beta <- exp(as.numeric(Pb_tr %*% av))

  fn <- function(cf) {
    C <- matrix(cf, nrow=nbs_tr, ncol=4); Z <- Ph_tr%*%C; dZ <- dPh_tr%*%C
    S <- exp(Z[,1]); E <- exp(Z[,2]); I <- exp(Z[,3]); R <- exp(Z[,4])
    mc <- log1p(pmax(0, N_POP*rho*sigma*E))
    mh <- log1p(pmax(0, N_POP*pH*gamma*I))
    mi <- log1p(pmax(0, N_POP*pICU*gamma*I))
    mr <- aR*I
    rc <- sqrt(wbc_tr)*ifelse(is.finite(yc_tr), yc_tr-mc, 0)
    rh <- sqrt(wbh_tr)*ifelse(is.finite(yh_tr), yh_tr-mh, 0)
    ri <- sqrt(wbi_tr)*ifelse(is.finite(yi_tr), yi_tr-mi, 0)
    rr <- sqrt(wbr_tr)*ifelse(is.finite(yr_tr), yr_tr-mr, 0)
    rp <- sqrt(ETA*dt)*as.numeric(dZ - seir_rhs(Z, beta))
    rm <- sqrt(KAPPA_MASS)*(S+E+I+R-1)
    rb <- sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX, 0)
    rs <- sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb_tr %*% av)
    c(rc, rh, ri, rr, rp, rm, rb, rs)
  }
  r <- nls.lm(par=as.numeric(Cs), fn=fn,
              control=nls.lm.control(maxiter=MAX_INNER, ftol=1e-10, ptol=1e-10))
  list(C=matrix(r$par, nrow=nbs_tr, ncol=4), J=sum(r$fvec^2))
}

cascade_cv <- function(th0, C0) {
  th <- th0; Cc <- C0; bJ <- Inf
  for (oi in 1:MAX_OUTER) {
    inn <- solve_inner_cv(Cc, th); Cc <- inn$C; Jc <- inn$J
    if (Jc < bJ) bJ <- Jc
    eps <- 1e-4; gr <- numeric(length(th))
    for (j in seq_along(th)) {
      tp <- th; tp[j] <- tp[j]+eps
      gr[j] <- (solve_inner_cv(Cc, tp)$J - Jc)/eps
    }
    gn <- sqrt(sum(gr^2))
    if (gn > 1e-12) {
      d <- -gr/gn; st <- STEP0; found <- FALSE
      for (ls in 1:10) {
        tt <- th+st*d; it <- solve_inner_cv(Cc, tt)
        if (it$J < Jc-1e-4*st*gn) {th<-tt; Cc<-it$C; found<-TRUE; break}
        st <- st*0.5
      }
      if (!found) th <- th+1e-4*d
    }
    if (gn < 1e-5*max(1, sqrt(sum(th^2)))) break
  }
  list(theta=th, C=Cc, J=bJ)
}

set.seed(42)
cat("  Running cross-validation fit (days 0-299)...\n")
res_cv <- tryCatch(cascade_cv(thi_tr, Ci_tr), error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(res_cv)) {
  th <- res_cv$theta; C <- res_cv$C
  av <- th[1:nbb_tr]
  rho  <- from_u(th[nbb_tr+1], bounds_cv$rho[1], bounds_cv$rho[2])
  pH   <- from_u(th[nbb_tr+2], bounds_cv$pH[1], bounds_cv$pH[2])
  pICU <- from_u(th[nbb_tr+3], bounds_cv$pICU[1], bounds_cv$pICU[2])
  aR   <- from_u(th[nbb_tr+4], bounds_cv$alphaR[1], bounds_cv$alphaR[2])

  # Training metrics
  Z_tr <- Ph_tr %*% C; S_tr <- exp(Z_tr[,1]); I_tr_fit <- exp(Z_tr[,3])
  beta_tr <- exp(as.numeric(Pb_tr %*% av)); Rt_tr <- beta_tr * S_tr / gamma

  ok_rt_tr <- is.finite(rt_rivm[train_idx])
  rt_corr_train <- cor(rt_rivm[train_idx][ok_rt_tr], Rt_tr[ok_rt_tr])

  # For test period: extrapolate the last fitted beta and SEIR state
  # (simple: use the last fitted values as constant forecast)
  last_beta <- beta_tr[length(beta_tr)]
  last_S <- S_tr[length(S_tr)]
  last_Rt <- last_beta * last_S / gamma
  # Predicted Rt for test period
  test_days <- sum(test_idx)
  Rt_pred <- rep(last_Rt, test_days)  # naive constant extrapolation

  ok_rt_te <- is.finite(rt_rivm[test_idx])
  rt_corr_test <- cor(rt_rivm[test_idx][ok_rt_te], Rt_pred[ok_rt_te])
  rt_rmse_test <- sqrt(mean((rt_rivm[test_idx][ok_rt_te] - Rt_pred[ok_rt_te])^2))

  cat(sprintf("  Train Rt corr: %.3f\n", rt_corr_train))
  cat(sprintf("  Test Rt RMSE: %.3f (naive constant forecast)\n", rt_rmse_test))
  cat(sprintf("  Test Rt corr: %.3f\n", rt_corr_test))
  cat(sprintf("  Fitted rho: %.3f, pH: %.4f\n", rho, pH))

  fwrite(data.table(
    metric = c("rt_corr_train", "rt_corr_test", "rt_rmse_test",
               "rho", "pH", "pICU", "alphaR", "last_Rt"),
    value = c(rt_corr_train, rt_corr_test, rt_rmse_test,
              rho, pH, pICU, aR, last_Rt)
  ), "crossval_results.csv")

  # TikZ: training fit + test forecast
  cv_df <- data.table(day = tv, rt_rivm = rt_rivm, type = "data")
  cv_df[train_idx, rt_model := Rt_tr]
  cv_df[test_idx, rt_forecast := Rt_pred]
  fwrite(cv_df, "tikz_crossval.csv")
  cat("  Saved crossval_results.csv and tikz_crossval.csv\n")
}

# ══════════════════════════════════════════════════════════════════
# EXTENSION 3: EpiEstim Rt COMPARISON
# ══════════════════════════════════════════════════════════════════

cat("\n=== Extension 3: EpiEstim Rt ===\n")

# Check if EpiEstim is available
epiestim_available <- requireNamespace("EpiEstim", quietly = TRUE)

if (epiestim_available) {
  library(EpiEstim)

  # Estimate Rt using Cori method
  # Serial interval: mean 6.5 days, sd 3.5 days (COVID-19 consensus)
  incid <- pmax(1, cases_raw)  # EpiEstim needs positive integers
  si_config <- make_config(
    mean_si = 6.5, std_si = 3.5,
    t_start = seq(2, n - 6), t_end = seq(8, n)
  )

  rt_est <- estimate_R(incid, method = "parametric_si", config = si_config)

  # Extract mean Rt and CI
  rt_cori <- data.table(
    day = rt_est$R$t_end,
    rt_cori_mean = rt_est$R$`Mean(R)`,
    rt_cori_q025 = rt_est$R$`Quantile.0.025(R)`,
    rt_cori_q975 = rt_est$R$`Quantile.0.975(R)`
  )

  # Merge with model Rt
  all_rt <- merge(
    data.table(day = tv, rt_rivm = rt_rivm, rt_model = old$Rt_model),
    rt_cori, by = "day", all.x = TRUE
  )

  fwrite(all_rt, "tikz_epiestim.csv")

  # Correlations
  ok <- complete.cases(all_rt[, .(rt_rivm, rt_model, rt_cori_mean)])
  cat(sprintf("  Cori Rt vs RIVM: r=%.3f\n", cor(all_rt$rt_rivm[ok], all_rt$rt_cori_mean[ok])))
  cat(sprintf("  Model Rt vs RIVM: r=%.3f\n", cor(all_rt$rt_rivm[ok], all_rt$rt_model[ok])))
  cat(sprintf("  Model Rt vs Cori: r=%.3f\n", cor(all_rt$rt_model[ok], all_rt$rt_cori_mean[ok])))
  cat("  Saved tikz_epiestim.csv\n")
} else {
  cat("  EpiEstim not installed. Install with: install.packages('EpiEstim')\n")
  cat("  Skipping EpiEstim comparison.\n")
}

cat("\n============================================================\n")
cat("Extensions done.\n")
cat("============================================================\n")
