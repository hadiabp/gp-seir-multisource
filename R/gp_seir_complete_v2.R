#!/usr/bin/env Rscript
##############################################################################
# gp_seir_complete_v2.R — GP-SEIR Full-Period Analysis (Version 2)
#
# Full analysis on the 366-day window (March 2020 – March 2021).
#
# Output prefix: FALSE_ or TRUE_ depending on USE_DELAY_KERNEL.
#
# Dependencies: data.table, fda, minpack.lm, zoo
#               EpiEstim (optional – graceful skip if absent)
#
# Usage:
#   Rscript gp_seir_complete_v2.R
#   Rscript gp_seir_complete_v2.R --delay  # enable delay kernels
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(fda); library(minpack.lm); library(zoo)
})

HAVE_EPIESTIM <- requireNamespace("EpiEstim", quietly = TRUE)
if (!HAVE_EPIESTIM)
  cat("NOTE: EpiEstim not installed. Cori comparison will be skipped.\n",
      "      Install with: install.packages('EpiEstim')\n\n")

cat("============================================================\n")
cat("GP-SEIR v2: Full-Period Analysis with Banzhaf + Cori\n")
cat("============================================================\n\n")

# ─────────────────────────────────────────────────────────────────
# §1  SETTINGS
# ─────────────────────────────────────────────────────────────────

DATA_FILE  <- "gp_input_final_v2.csv"
DATE_START <- as.Date("2020-03-01")
DATE_END   <- as.Date("2021-03-01")

N_POP   <- 17400000
sigma   <- 1 / 5.5     # working assumption; mean latent period 5.5 d
gamma   <- 1 / 9.5     # working assumption; mean infectious period 9.5 d

STATE_KNOT_DAYS <- 14;  STATE_NORDER <- 4
BETA_KNOT_DAYS  <- 28;  BETA_NORDER  <- 4

ETA        <- 1e6      # ODE penalty
KAPPA_BETA <- 1e3      # roughness of log-beta
KAPPA_BMAG <- 1e4      # soft cap on beta magnitude
BETA_MAX   <- 0.6
KAPPA_MASS <- 1e6      # mass conservation

SEEDS <- list(main=42, ablation=42, wgrid=42, init=42, profile=42,
              sg_inner=123, sg_sweep=42, waves=42)
N_STARTS      <- 5     # warm starts per config
N_COLD_STARTS <- 3     # cold starts per ablation config
MAX_INNER     <- 300
MAX_OUTER     <- 40
USE_RIVM_INIT <- FALSE  # MUST be FALSE – prevents validation leakage
STEP0         <- 1.0
BETA_PERTURB  <- 0.30
PARAM_PERTURB <- 0.50

# Delay convolution (log-normal kernels)
args_cmd <- commandArgs(trailingOnly = TRUE)
USE_DELAY_KERNEL <- "--delay" %in% args_cmd
DELAY_CASE_MEAN <- 3.0;  DELAY_CASE_SD <- 1.5
DELAY_HOSP_MEAN <- 2.0;  DELAY_HOSP_SD <- 1.0
DELAY_ICU_MEAN  <- 3.5;  DELAY_ICU_SD  <- 1.5
DELAY_MAX_LAG   <- 14

OUT_PREFIX <- ifelse(USE_DELAY_KERNEL, "TRUE_", "FALSE_")
out <- function(f) paste0(OUT_PREFIX, f)

bounds <- list(
  rho    = c(0.05, 0.60),
  pH     = c(0.002, 0.05),
  pICU   = c(0.0003, 0.02),
  alphaR = c(0.10, 1.00)
)

# Cori / EpiEstim settings (COVID-19 serial interval, Bi et al. 2020)
CORI_SI_MEAN <- 5.1;  CORI_SI_SD <- 2.8;  CORI_TAU <- 7L

# Game-theory settings
SOURCE_NAMES  <- c("cases", "hosp", "icu", "radar")
N_PLAYERS     <- 4L
ALL_SUBSETS   <- unlist(lapply(1:4, function(k)
                   combn(SOURCE_NAMES, k, simplify = FALSE)), recursive = FALSE)
# Banzhaf stability: flag window if |sum(beta_i)| < threshold OR sign mismatch
BANZHAF_STABLE_THRESHOLD <- 0.05

cat(sprintf("  sigma=1/%.1f, gamma=1/%.1f, gen_time=%.1f d\n",
            1/sigma, 1/gamma, 1/sigma + 1/gamma))
cat(sprintf("  USE_DELAY_KERNEL=%s, output prefix: %s\n\n",
            USE_DELAY_KERNEL, OUT_PREFIX))


# ─────────────────────────────────────────────────────────────────
# §2  DATA LOADING
# ─────────────────────────────────────────────────────────────────

cat("--- Loading data ---\n")
if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)
df <- fread(DATA_FILE)[, date := as.Date(date)]
df <- df[date >= DATE_START & date <= DATE_END][order(date)]
n  <- nrow(df)
df[, t := as.numeric(date - DATE_START)]
tv <- df$t;  rng <- range(tv);  dt <- 1;  dates <- df$date

cases_raw  <- as.numeric(df$cases_raw)
hosp_raw   <- as.numeric(df$hosp_nice)
icu_raw    <- as.numeric(df$icu_daily)
radar_frac <- as.numeric(df$radar_I_frac_sm7)

# VALIDATION DATA — loaded here, never passed to solve_inner / cascade
rt_rivm   <- as.numeric(df$rt_rivm)
prev_rivm <- as.numeric(df$prev_json_avg)
tests_tot <- as.numeric(df$tests_total)
tests_pos <- as.numeric(df$tests_positive)
positivity <- as.numeric(df$positivity_rate)

y_cases <- log1p(pmax(0, cases_raw))
y_hosp  <- log1p(pmax(0, hosp_raw))
y_icu   <- log1p(pmax(0, icu_raw))
y_radar <- radar_frac

wsd <- function(y) { s <- sd(y, na.rm=TRUE); if (!is.finite(s)||s<=0) s<-1; 1/s^2 }
wb_cases <- wsd(y_cases); wb_hosp <- wsd(y_hosp)
wb_icu   <- wsd(y_icu);   wb_radar <- wsd(y_radar)

cat(sprintf("  %d days (%s to %s)\n", n, DATE_START, DATE_END))
cat(sprintf("  ODE/cases ratio: %.0fx\n\n", sqrt(ETA*dt)/sqrt(wb_cases)))


# ─────────────────────────────────────────────────────────────────
# §3  B-SPLINE BASES
# ─────────────────────────────────────────────────────────────────

mkb <- function(kd, no) {
  kn <- unique(c(rng[1], seq(rng[1], rng[2], by=kd), rng[2]))
  create.bspline.basis(rng, length(kn)+no-2, no, kn)
}
Bs   <- mkb(STATE_KNOT_DAYS, STATE_NORDER)
Bb   <- mkb(BETA_KNOT_DAYS,  BETA_NORDER)
nbs  <- Bs$nbasis;  nbb  <- Bb$nbasis
Ph   <- eval.basis(tv, Bs, 0);  dPh  <- eval.basis(tv, Bs, 1)
Pb   <- eval.basis(tv, Bb, 0);  d2Pb <- eval.basis(tv, Bb, 2)
cat(sprintf("  State basis: %d, Beta basis: %d\n", nbs, nbb))


# ─────────────────────────────────────────────────────────────────
# §4  DELAY KERNELS
# ─────────────────────────────────────────────────────────────────

make_delay_kernel <- function(mean_d, sd_d, max_lag = DELAY_MAX_LAG) {
  if (mean_d <= 0) return(1)
  lags <- 0:max_lag
  mu_ln  <- log(mean_d^2 / sqrt(sd_d^2 + mean_d^2))
  sig_ln <- sqrt(log(1 + sd_d^2 / mean_d^2))
  w <- dlnorm(lags, mu_ln, sig_ln); w / sum(w)
}
apply_delay <- function(y, kernel) {
  if (length(kernel)==1) return(y)
  klen <- length(kernel); n_y <- length(y); yc <- numeric(n_y)
  for (t in 1:n_y)
    for (k in 1:klen) { s <- t-(k-1); if (s>=1) yc[t] <- yc[t]+kernel[k]*y[s] }
  yc
}

if (USE_DELAY_KERNEL) {
  h_case <- make_delay_kernel(DELAY_CASE_MEAN, DELAY_CASE_SD)
  h_hosp <- make_delay_kernel(DELAY_HOSP_MEAN, DELAY_HOSP_SD)
  h_icu  <- make_delay_kernel(DELAY_ICU_MEAN,  DELAY_ICU_SD)
  cat(sprintf("  Delay kernels: case %.1fd, hosp %.1fd, icu %.1fd\n",
              DELAY_CASE_MEAN, DELAY_HOSP_MEAN, DELAY_ICU_MEAN))
} else {
  h_case <- 1; h_hosp <- 1; h_icu <- 1
  cat("  Delay kernels: disabled\n")
}


# ─────────────────────────────────────────────────────────────────
# §5  HELPERS
# ─────────────────────────────────────────────────────────────────

to_u   <- function(x, lo, hi) qlogis((x-lo)/(hi-lo))
from_u <- function(u, lo, hi) lo+(hi-lo)*plogis(u)

seir_rhs <- function(Z, beta)
  cbind(-beta*exp(Z[,3]),
         beta*exp(Z[,1]+Z[,3]-Z[,2])-sigma,
         sigma*exp(Z[,2]-Z[,3])-gamma,
         gamma*exp(Z[,3]-Z[,4]))

safe_cor  <- function(x,y,m="pearson"){ok<-is.finite(x)&is.finite(y);if(sum(ok)<5)return(NA_real_);cor(x[ok],y[ok],method=m)}
safe_rmse <- function(x,y){ok<-is.finite(x)&is.finite(y);if(sum(ok)<3)return(NA_real_);sqrt(mean((x[ok]-y[ok])^2))}
safe_skill <- function(x,y){ok<-is.finite(x)&is.finite(y);if(sum(ok)<5)return(NA_real_);mse<-mean((x[ok]-y[ok])^2);vref<-var(y[ok]);if(!is.finite(vref)||vref<1e-12)return(NA_real_);1-mse/vref}


# ─────────────────────────────────────────────────────────────────
# §6  GAME-THEORY FUNCTIONS (Shapley + Banzhaf + interactions)
# ─────────────────────────────────────────────────────────────────

# v(S) lookup from ablation data frame
v_lookup <- function(abl_df, src_set, metric_col) {
  if (length(src_set)==0) return(0)
  key <- paste(sort(src_set), collapse="+")
  idx <- which(abl_df$sources==key)
  if (length(idx)==0) return(0)
  val <- abl_df[[metric_col]][idx[1]]
  ifelse(is.finite(val), val, 0)
}

# Exact Shapley values for n=4 players
compute_shapley <- function(abl_df, metric_col) {
  phi <- numeric(N_PLAYERS); names(phi) <- SOURCE_NAMES
  for (i in seq_along(SOURCE_NAMES)) {
    others <- SOURCE_NAMES[-i]; total <- 0
    for (sz in 0:(N_PLAYERS-1)) {
      coals <- if (sz==0) list(character(0)) else combn(others, sz, simplify=FALSE)
      weight <- factorial(sz)*factorial(N_PLAYERS-sz-1)/factorial(N_PLAYERS)
      for (S in coals) {
        vw  <- v_lookup(abl_df, c(S, SOURCE_NAMES[i]), metric_col)
        vnw <- v_lookup(abl_df, S, metric_col)
        if (is.finite(vw) && is.finite(vnw)) total <- total+weight*(vw-vnw)
      }
    }
    phi[i] <- total
  }
  phi
}

# Banzhaf values: equal weight 1/2^(n-1) across all coalitions
compute_banzhaf <- function(abl_df, metric_col) {
  norm <- 1/(2^(N_PLAYERS-1))
  beta_v <- numeric(N_PLAYERS); names(beta_v) <- SOURCE_NAMES
  for (i in seq_along(SOURCE_NAMES)) {
    others <- SOURCE_NAMES[-i]; total <- 0
    for (sz in 0:(N_PLAYERS-1)) {
      coals <- if (sz==0) list(character(0)) else combn(others, sz, simplify=FALSE)
      for (S in coals) {
        vw  <- v_lookup(abl_df, c(S, SOURCE_NAMES[i]), metric_col)
        vnw <- v_lookup(abl_df, S, metric_col)
        if (is.finite(vw) && is.finite(vnw)) total <- total+(vw-vnw)
      }
    }
    beta_v[i] <- norm*total
  }
  beta_v
}

# Normalised Banzhaf: sum rescaled to v(N)
normalise_banzhaf <- function(beta_v, vN) {
  s <- sum(beta_v)
  sign_mismatch <- abs(s)>1e-10 && abs(vN)>0.01 && s*vN < 0
  if (abs(s) < BANZHAF_STABLE_THRESHOLD || sign_mismatch)
    return(list(values=setNames(rep(NA_real_,N_PLAYERS),SOURCE_NAMES), stable=FALSE))
  list(values=beta_v*vN/s, stable=TRUE)
}

# Shapley interaction indices I_ij
compute_interactions <- function(abl_df, metric_col) {
  pairs <- combn(SOURCE_NAMES, 2, simplify=FALSE)
  result <- list()
  for (pair in pairs) {
    i <- pair[1]; j <- pair[2]; others <- setdiff(SOURCE_NAMES, pair); total <- 0
    for (sz in 0:(N_PLAYERS-2)) {
      coals <- if (sz==0) list(character(0)) else combn(others, sz, simplify=FALSE)
      weight <- factorial(sz)*factorial(N_PLAYERS-sz-2)/factorial(N_PLAYERS-1)
      for (S in coals) {
        v_ij  <- v_lookup(abl_df, c(S,i,j), metric_col)
        v_i   <- v_lookup(abl_df, c(S,i),   metric_col)
        v_j   <- v_lookup(abl_df, c(S,j),   metric_col)
        v_0   <- v_lookup(abl_df, S,         metric_col)
        if (all(is.finite(c(v_ij,v_i,v_j,v_0))))
          total <- total+weight*(v_ij-v_i-v_j+v_0)
      }
    }
    result[[length(result)+1]] <- data.table(source_i=i,source_j=j,interaction=total)
  }
  rbindlist(result)
}

# Full game-theory output for one metric: phi, beta, normBeta, v_grand
game_theory_for <- function(abl_df, metric_col, metric_name) {
  vN  <- v_lookup(abl_df, SOURCE_NAMES, metric_col)
  phi <- compute_shapley(abl_df, metric_col)
  bv  <- compute_banzhaf(abl_df, metric_col)
  nb  <- normalise_banzhaf(bv, vN)
  ints <- compute_interactions(abl_df, metric_col)

  df_vals <- data.table(
    metric=metric_col, metric_name=metric_name, source=SOURCE_NAMES,
    shapley=phi, banzhaf=bv, norm_banzhaf=nb$values,
    banzhaf_stable=nb$stable,
    singleton=sapply(SOURCE_NAMES, function(s) v_lookup(abl_df,s,metric_col)),
    v_grand=vN
  )
  list(values=df_vals, interactions=ints, vN=vN, phi=phi, banzhaf=bv,
       norm_banzhaf=nb$values, stable=nb$stable)
}


# ─────────────────────────────────────────────────────────────────
# §7  CORI METHOD (EpiEstim wrapper)
# ─────────────────────────────────────────────────────────────────

run_cori <- function(case_vec, date_vec,
                     si_mean=CORI_SI_MEAN, si_sd=CORI_SI_SD, tau=CORI_TAU) {
  if (!HAVE_EPIESTIM) return(NULL)
  inc <- pmax(0, round(ifelse(is.finite(case_vec), case_vec, 0)))
  cfg <- tryCatch(
    EpiEstim::make_config(method="parametric_si", mean_si=si_mean, std_si=si_sd,
      t_start=2:(length(inc)-tau+1), t_end=(tau+1):length(inc)),
    error=function(e) NULL)
  if (is.null(cfg)) return(NULL)
  res <- tryCatch(EpiEstim::estimate_R(inc, method="parametric_si", config=cfg),
                  error=function(e){ cat("  EpiEstim error:", e$message, "\n"); NULL })
  if (is.null(res)) return(NULL)
  r_dt <- as.data.table(res$R)
  r_dt[, date := date_vec[t_end]]
  r_dt[, .(date, Rt_cori=`Mean(R)`,
            Rt_cori_lo=`Quantile.0.025(R)`,
            Rt_cori_hi=`Quantile.0.975(R)`)]
}


# ─────────────────────────────────────────────────────────────────
# §8  INITIALIZATION (hospital log-derivative — no validation leakage)
# ─────────────────────────────────────────────────────────────────

if (USE_RIVM_INIT) stop("USE_RIVM_INIT=TRUE is disabled.")

c0 <- ifelse(is.finite(cases_raw), cases_raw, 0)
E0 <- pmax(1e-10, c0/(N_POP*0.20*sigma))
I0 <- as.numeric(rollmean(E0, 7, fill=median(E0), align="right"))
R0 <- pmin(0.8, cumsum(pmax(c0,0))/N_POP)
S0 <- pmax(1e-6, 1-E0-I0-R0); m0<-S0+E0+I0+R0
S0<-S0/m0; E0<-E0/m0; I0<-I0/m0; R0<-R0/m0
Ci <- smooth.basis(tv, cbind(log(S0),log(E0),log(I0),log(R0)),
                   fdPar(Bs, int2Lfd(2), 1e-4))$fd$coefs

h0_safe <- ifelse(is.finite(hosp_raw), pmax(0.5,hosp_raw), 0.5)
hosp_sm <- as.numeric(rollmean(h0_safe, 28, fill=NA, align="center"))
hosp_sm[is.na(hosp_sm)] <- hosp_sm[which(!is.na(hosp_sm))[1]]
dlog_h  <- c(diff(log(pmax(1,hosp_sm))), 0)
dlog_sm <- as.numeric(rollmean(dlog_h, 21, fill=0, align="center"))
dlog_sm[is.na(dlog_sm)] <- 0
b0 <- gamma*(1+dlog_sm*(1/sigma+1/gamma))
b0 <- pmax(0.05, pmin(0.35, b0))
b0 <- as.numeric(rollmean(b0, 14, fill=NA, align="center"))
b0[is.na(b0)] <- 0.15
ai  <- as.numeric(smooth.basis(tv, log(b0), fdPar(Bb, int2Lfd(2), 1e-2))$fd$coefs)
thi <- c(ai,
         to_u(0.20,   bounds$rho[1],    bounds$rho[2]),
         to_u(0.012,  bounds$pH[1],     bounds$pH[2]),
         to_u(0.0025, bounds$pICU[1],   bounds$pICU[2]),
         to_u(0.40,   bounds$alphaR[1], bounds$alphaR[2]))
beta_init_curve <- exp(as.numeric(Pb %*% ai))
cat(sprintf("  Beta init: [%.3f, %.3f]\n\n", min(beta_init_curve), max(beta_init_curve)))


# ─────────────────────────────────────────────────────────────────
# §9  INNER SOLVER AND CASCADE (identical to v1)
# ─────────────────────────────────────────────────────────────────

solve_inner <- function(Cs, th, wc, wh, wi, wr) {
  av <- th[1:nbb]
  rho  <- from_u(th[nbb+1], bounds$rho[1],    bounds$rho[2])
  pH   <- from_u(th[nbb+2], bounds$pH[1],     bounds$pH[2])
  pICU <- from_u(th[nbb+3], bounds$pICU[1],   bounds$pICU[2])
  aR   <- from_u(th[nbb+4], bounds$alphaR[1], bounds$alphaR[2])
  beta <- exp(as.numeric(Pb %*% av))
  fn <- function(cf) {
    C <- matrix(cf,nrow=nbs,ncol=4); Z<-Ph%*%C; dZ<-dPh%*%C
    S<-exp(Z[,1]); E<-exp(Z[,2]); I<-exp(Z[,3]); R<-exp(Z[,4])
    fc<-N_POP*rho*sigma*E; fh<-N_POP*pH*gamma*I; fi<-N_POP*pICU*gamma*I
    mc<-log1p(pmax(0,apply_delay(fc,h_case))); mh<-log1p(pmax(0,apply_delay(fh,h_hosp)))
    mi<-log1p(pmax(0,apply_delay(fi,h_icu)));  mr<-aR*I
    rc<-sqrt(wb_cases*wc)*ifelse(is.finite(y_cases),y_cases-mc,0)
    rh<-sqrt(wb_hosp*wh)*ifelse(is.finite(y_hosp),y_hosp-mh,0)
    ri<-sqrt(wb_icu*wi)*ifelse(is.finite(y_icu),y_icu-mi,0)
    rr<-sqrt(wb_radar*wr)*ifelse(is.finite(y_radar),y_radar-mr,0)
    rp<-sqrt(ETA*dt)*as.numeric(dZ-seir_rhs(Z,beta))
    rm<-sqrt(KAPPA_MASS)*(S+E+I+R-1)
    rb<-sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX,0)
    rs<-sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb%*%av)
    c(rc,rh,ri,rr,rp,rm,rb,rs)
  }
  r <- nls.lm(par=as.numeric(Cs), fn=fn,
              control=nls.lm.control(maxiter=MAX_INNER,ftol=1e-10,ptol=1e-10))
  list(C=matrix(r$par,nrow=nbs,ncol=4), J=sum(r$fvec^2))
}

cascade <- function(th0, C0, wc, wh, wi, wr) {
  th <- th0; Cc <- C0; bJ <- Inf
  for (oi in 1:MAX_OUTER) {
    inn <- solve_inner(Cc,th,wc,wh,wi,wr); Cc<-inn$C; Jc<-inn$J; if(Jc<bJ) bJ<-Jc
    eps <- 1e-4; gr <- numeric(length(th))
    for (j in seq_along(th)) { tp<-th; tp[j]<-tp[j]+eps; gr[j]<-(solve_inner(Cc,tp,wc,wh,wi,wr)$J-Jc)/eps }
    gn <- sqrt(sum(gr^2))
    if (gn>1e-12) {
      d <- -gr/gn; st<-STEP0; found<-FALSE
      for (ls in 1:10) { tt<-th+st*d; it<-solve_inner(Cc,tt,wc,wh,wi,wr)
        if (it$J<Jc-1e-4*st*gn){th<-tt;Cc<-it$C;found<-TRUE;break}; st<-st*0.5 }
      if (!found) th<-th+1e-4*d
    }
    if (gn<1e-5*max(1,sqrt(sum(th^2)))) break
  }
  list(theta=th, C=Cc, J=bJ)
}


# ─────────────────────────────────────────────────────────────────
# §10  EXTRACT_FULL — post-hoc diagnostics (validation data used here)
# ─────────────────────────────────────────────────────────────────

extract_full <- function(res) {
  th <- res$theta; C <- res$C; av <- th[1:nbb]
  rho  <- from_u(th[nbb+1],bounds$rho[1],   bounds$rho[2])
  pH   <- from_u(th[nbb+2],bounds$pH[1],    bounds$pH[2])
  pICU <- from_u(th[nbb+3],bounds$pICU[1],  bounds$pICU[2])
  aR   <- from_u(th[nbb+4],bounds$alphaR[1],bounds$alphaR[2])
  Z  <- Ph%*%C; dZ <- dPh%*%C
  S  <- exp(Z[,1]); E <- exp(Z[,2]); I <- exp(Z[,3]); R <- exp(Z[,4])
  beta <- exp(as.numeric(Pb%*%av)); Rt <- beta*S/gamma; mass <- S+E+I+R
  ode_mat  <- dZ - seir_rhs(Z,beta); ode_rmse <- sqrt(mean(ode_mat^2))
  fc<-N_POP*rho*sigma*E; fh<-N_POP*pH*gamma*I; fi<-N_POP*pICU*gamma*I
  mc<-apply_delay(fc,h_case); mh<-apply_delay(fh,h_hosp); mi<-apply_delay(fi,h_icu); mr<-aR*I
  I_count <- N_POP*I
  resid_cases <- ifelse(is.finite(y_cases),y_cases-log1p(pmax(0,mc)),NA)
  resid_hosp  <- ifelse(is.finite(y_hosp), y_hosp-log1p(pmax(0,mh)),NA)
  resid_icu   <- ifelse(is.finite(y_icu),  y_icu-log1p(pmax(0,mi)),NA)
  resid_radar <- ifelse(is.finite(y_radar),y_radar-mr,NA)
  ok_c  <- is.finite(cases_raw); ok_h <- is.finite(hosp_raw)&hosp_raw>0

  # Full-period validation metrics
  ok_rt   <- is.finite(rt_rivm)
  ok_prev <- is.finite(prev_rivm)
  rt_corr     <- safe_cor(rt_rivm,Rt)
  rt_spearman <- safe_cor(rt_rivm,Rt,"spearman")
  rt_rmse_v   <- safe_rmse(rt_rivm,Rt)
  rt_bias_v   <- if(sum(ok_rt)>3) mean(Rt[ok_rt]-rt_rivm[ok_rt]) else NA_real_
  rt_skill_v  <- safe_skill(Rt,rt_rivm)
  prev_corr     <- safe_cor(prev_rivm,I_count)
  prev_spearman <- safe_cor(prev_rivm,I_count,"spearman")
  prev_skill_v  <- safe_skill(I_count,prev_rivm)
  threshold_v <- if(sum(ok_rt)>3) mean((Rt[ok_rt]>1)==(rt_rivm[ok_rt]>1)) else NA_real_

  # Pre/post June 12 split (RIVM methodology switch)
  cut <- as.Date("2020-06-12")
  pre  <- ok_rt & dates<cut; post <- ok_rt & dates>=cut
  rt_corr_pre  <- safe_cor(rt_rivm[pre],  Rt[pre])
  rt_corr_post <- safe_cor(rt_rivm[post], Rt[post])
  rt_bias_pre  <- if(sum(pre)>3)  mean(Rt[pre]-rt_rivm[pre])   else NA_real_
  rt_bias_post <- if(sum(post)>3) mean(Rt[post]-rt_rivm[post]) else NA_real_

  # Per-wave sliced metrics
  w1  <- ok_rt  & dates<"2020-06-01"
  iw  <- ok_rt  & dates>="2020-06-01" & dates<"2020-10-01"
  w2  <- ok_rt  & dates>="2020-10-01"
  pw1 <- ok_prev & dates<"2020-06-01"
  piw <- ok_prev & dates>="2020-06-01" & dates<"2020-10-01"
  pw2 <- ok_prev & dates>="2020-10-01"
  rt_wave1   <- safe_cor(rt_rivm[w1],   Rt[w1])
  rt_inter   <- safe_cor(rt_rivm[iw],   Rt[iw])
  rt_wave2   <- safe_cor(rt_rivm[w2],   Rt[w2])
  rt_sp_wave1 <- safe_cor(rt_rivm[w1],  Rt[w1],"spearman")
  rt_sp_inter <- safe_cor(rt_rivm[iw],  Rt[iw],"spearman")
  rt_sp_wave2 <- safe_cor(rt_rivm[w2],  Rt[w2],"spearman")
  rt_skill_w1 <- safe_skill(Rt[w1],  rt_rivm[w1])
  rt_skill_iw <- safe_skill(Rt[iw],  rt_rivm[iw])
  rt_skill_w2 <- safe_skill(Rt[w2],  rt_rivm[w2])
  prev_wave1 <- safe_cor(prev_rivm[pw1],I_count[pw1])
  prev_inter <- safe_cor(prev_rivm[piw],I_count[piw])
  prev_wave2 <- safe_cor(prev_rivm[pw2],I_count[pw2])
  prev_sp_w1 <- safe_cor(prev_rivm[pw1],I_count[pw1],"spearman")
  prev_sp_iw <- safe_cor(prev_rivm[piw],I_count[piw],"spearman")
  prev_sp_w2 <- safe_cor(prev_rivm[pw2],I_count[pw2],"spearman")
  prev_sk_w1 <- safe_skill(I_count[pw1],prev_rivm[pw1])
  prev_sk_iw <- safe_skill(I_count[piw],prev_rivm[piw])
  prev_sk_w2 <- safe_skill(I_count[pw2],prev_rivm[pw2])

  denom <- N_POP*sigma*E
  rho_eff <- ifelse(is.finite(cases_raw)&denom>1, cases_raw/denom, NA)

  list(
    S=S,E=E,I=I,R=R,beta=beta,Rt=Rt,mass=mass,
    mu_cases=mc,mu_hosp=mh,mu_icu=mi,mu_radar=mr,
    rho=rho,pH=pH,pICU=pICU,alphaR=aR,
    ode_mat=ode_mat,ode_rmse=ode_rmse,
    resid_cases=resid_cases,resid_hosp=resid_hosp,resid_icu=resid_icu,resid_radar=resid_radar,
    cases_rmse=safe_rmse(cases_raw,mc),hosp_rmse=safe_rmse(hosp_raw,mh),mass_err=max(abs(mass-1)),
    # Validation
    rt_corr=rt_corr,rt_spearman=rt_spearman,rt_rmse=rt_rmse_v,rt_bias=rt_bias_v,rt_skill=rt_skill_v,
    prev_corr=prev_corr,prev_spearman=prev_spearman,prev_skill=prev_skill_v,
    threshold_agree=threshold_v,
    rt_corr_pre=rt_corr_pre,rt_corr_post=rt_corr_post,
    rt_bias_pre=rt_bias_pre,rt_bias_post=rt_bias_post,
    rt_wave1=rt_wave1,rt_inter=rt_inter,rt_wave2=rt_wave2,
    rt_sp_wave1=rt_sp_wave1,rt_sp_inter=rt_sp_inter,rt_sp_wave2=rt_sp_wave2,
    rt_skill_w1=rt_skill_w1,rt_skill_iw=rt_skill_iw,rt_skill_w2=rt_skill_w2,
    prev_wave1=prev_wave1,prev_inter=prev_inter,prev_wave2=prev_wave2,
    prev_sp_w1=prev_sp_w1,prev_sp_iw=prev_sp_iw,prev_sp_w2=prev_sp_w2,
    prev_sk_w1=prev_sk_w1,prev_sk_iw=prev_sk_iw,prev_sk_w2=prev_sk_w2,
    rho_eff=rho_eff,
    beta_range=c(min(beta),max(beta)),Rt_range=c(min(Rt),max(Rt)),
    R_end=R[n]
  )
}


# ─────────────────────────────────────────────────────────────────
# §11  MULTI-START WRAPPERS
# ─────────────────────────────────────────────────────────────────

run_fit_full <- function(wc,wh,wi,wr,label="",warm_th=NULL,warm_C=NULL) {
  cat(sprintf("  %-35s ", label))
  best<-NULL; bJ<-Inf; all_starts<-list()
  th0<-if(!is.null(warm_th)) warm_th else thi
  C0 <-if(!is.null(warm_C))  warm_C  else Ci
  for (s in 1:N_STARTS) {
    ts<-th0; Cs<-C0
    if (s>1) { ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,BETA_PERTURB)
               for(k in(nbb+1):(nbb+4)) ts[k]<-ts[k]+rnorm(1,0,PARAM_PERTURB) }
    r<-tryCatch(cascade(ts,Cs,wc,wh,wi,wr),error=function(e)NULL)
    if (!is.null(r)) {
      d<-extract_full(r)
      all_starts[[s]]<-list(J=r$J,rho=d$rho,pH=d$pH,pICU=d$pICU,alphaR=d$alphaR,
                              rt_corr=d$rt_corr,prev_corr=d$prev_corr)
      if (r$J<bJ){bJ<-r$J;best<-r}
    }
  }
  if (is.null(best)){cat("FAILED\n");return(NULL)}
  d<-extract_full(best)
  cat(sprintf("Rt=%.3f prev=%.3f J=%.1f\n",d$rt_corr,d$prev_corr,bJ))
  c(d,list(theta=best$theta,C_mat=best$C,J=best$J,all_starts=all_starts))
}

run_fit <- function(wc,wh,wi,wr,label="",warm_th=NULL,warm_C=NULL) {
  cat(sprintf("  %-35s ", label))
  best<-NULL; bJ<-Inf; all_J<-numeric(0)
  if (!is.null(warm_th)) {
    th0<-warm_th; C0<-if(!is.null(warm_C)) warm_C else Ci
    for (s in 1:N_STARTS) {
      ts<-th0; Cs<-C0
      if (s>1) { ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,BETA_PERTURB)
                 for(k in(nbb+1):(nbb+4)) ts[k]<-ts[k]+rnorm(1,0,PARAM_PERTURB) }
      r<-tryCatch(cascade(ts,Cs,wc,wh,wi,wr),error=function(e)NULL)
      if (!is.null(r)){all_J<-c(all_J,r$J);if(r$J<bJ){bJ<-r$J;best<-r}}
    }
  }
  for (s in 1:max(N_COLD_STARTS,if(is.null(warm_th)) N_STARTS else 0)) {
    ts<-thi; Cs<-Ci
    if (s>1) { ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,0.40)
               for(k in(nbb+1):(nbb+4)) ts[k]<-ts[k]+rnorm(1,0,0.6) }
    r<-tryCatch(cascade(ts,Cs,wc,wh,wi,wr),error=function(e)NULL)
    if (!is.null(r)){all_J<-c(all_J,r$J);if(r$J<bJ){bJ<-r$J;best<-r}}
  }
  if (is.null(best)){cat("FAILED\n");return(NULL)}
  d<-extract_full(best); n_near<-sum(all_J<=bJ*1.10)
  cat(sprintf("Rt=%.3f prev=%.3f J=%.0f skill=%.3f (%d/%d near)\n",
              d$rt_corr,d$prev_corr,bJ,ifelse(is.na(d$rt_skill),NA,d$rt_skill),
              n_near,length(all_J)))
  c(d,list(theta=best$theta,C_mat=best$C,J=best$J,
            n_converged=n_near,n_total=length(all_J)))
}


# ─────────────────────────────────────────────────────────────────
# STEP 1: MAIN FIT (all 4 sources)
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 1: Main fit ===\n")
set.seed(SEEDS$main)
main <- run_fit_full(1,1,1,1,"all sources (1,1,1,1)")

ms_df <- rbindlist(Filter(Negate(is.null), lapply(seq_along(main$all_starts), function(i) {
  s<-main$all_starts[[i]]; if(is.null(s)) return(NULL)
  as.data.table(c(list(start=i),s))
})))
if (nrow(ms_df)>0) fwrite(ms_df, out("multistart_params.csv"))

cat(sprintf("\n  rho=%.3f pH=%.4f pICU=%.4f alphaR=%.3f\n",
            main$rho,main$pH,main$pICU,main$alphaR))
cat(sprintf("  Rt: r=%.3f, rs=%.3f, RMSE=%.3f, bias=%.3f, skill=%.3f\n",
            main$rt_corr,main$rt_spearman,main$rt_rmse,main$rt_bias,main$rt_skill))
cat(sprintf("  Prev: r=%.3f, rs=%.3f, skill=%.3f\n",
            main$prev_corr,main$prev_spearman,main$prev_skill))
cat(sprintf("  Rt range: [%.2f, %.2f], R(end)=%.1f%%\n",
            main$Rt_range[1],main$Rt_range[2],main$R_end*100))
save.image(out("gp_seir_checkpoint_main.RData"))


# ─────────────────────────────────────────────────────────────────
# STEP 2: SOURCE ABLATION (all 15 subsets)
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 2: Full source ablation (15 subsets) ===\n")
set.seed(SEEDS$ablation)
abl_rows <- list(); abl_betas <- list()

for (si in seq_along(ALL_SUBSETS)) {
  srcs <- ALL_SUBSETS[[si]]; src_key <- paste(sort(srcs),collapse="+")
  w <- setNames(as.numeric(SOURCE_NAMES %in% srcs), SOURCE_NAMES)
  d <- run_fit(w["cases"],w["hosp"],w["icu"],w["radar"],src_key,
               warm_th=main$theta, warm_C=main$C_mat)
  if (!is.null(d)) {
    abl_rows[[src_key]] <- data.table(
      sources=src_key,
      w_cases=w["cases"],w_hosp=w["hosp"],w_icu=w["icu"],w_radar=w["radar"],
      J=d$J, n_converged=d$n_converged, n_total=d$n_total,
      rho=d$rho,pH=d$pH,pICU=d$pICU,alphaR=d$alphaR,
      # Full-period Pearson, Spearman, skill, RMSE, bias
      rt_corr=d$rt_corr,     rt_spearman=d$rt_spearman,
      rt_rmse=d$rt_rmse,     rt_bias=d$rt_bias,   rt_skill=d$rt_skill,
      prev_corr=d$prev_corr, prev_spearman=d$prev_spearman, prev_skill=d$prev_skill,
      # Pre/post June split
      rt_corr_pre=d$rt_corr_pre, rt_corr_post=d$rt_corr_post,
      # Per-wave sliced Pearson and skill
      rt_wave1=d$rt_wave1, rt_inter=d$rt_inter, rt_wave2=d$rt_wave2,
      rt_sp_w1=d$rt_sp_wave1, rt_sp_iw=d$rt_sp_inter, rt_sp_w2=d$rt_sp_wave2,
      rt_sk_w1=d$rt_skill_w1, rt_sk_iw=d$rt_skill_iw, rt_sk_w2=d$rt_skill_w2,
      prev_wave1=d$prev_wave1, prev_inter=d$prev_inter, prev_wave2=d$prev_wave2,
      prev_sp_w1=d$prev_sp_w1,prev_sp_iw=d$prev_sp_iw,prev_sp_w2=d$prev_sp_w2,
      prev_sk_w1=d$prev_sk_w1,prev_sk_iw=d$prev_sk_iw,prev_sk_w2=d$prev_sk_w2,
      # Other
      cases_rmse=d$cases_rmse,hosp_rmse=d$hosp_rmse,ode_rmse=d$ode_rmse,
      R_end=d$R_end,beta_min=d$beta_range[1],beta_max=d$beta_range[2],
      neg_rt_rmse=-d$rt_rmse,neg_J=-d$J
    )
    abl_betas[[src_key]] <- d$beta
  } else {
    # Record failure (v(S)=0 for Shapley)
    abl_rows[[src_key]] <- data.table(
      sources=src_key,w_cases=w["cases"],w_hosp=w["hosp"],w_icu=w["icu"],w_radar=w["radar"],
      J=NA,n_converged=0L,n_total=0L,rho=NA,pH=NA,pICU=NA,alphaR=NA,
      rt_corr=0,rt_spearman=0,rt_rmse=NA,rt_bias=NA,rt_skill=0,
      prev_corr=0,prev_spearman=0,prev_skill=0,
      rt_corr_pre=NA,rt_corr_post=NA,
      rt_wave1=0,rt_inter=0,rt_wave2=0,rt_sp_w1=NA,rt_sp_iw=NA,rt_sp_w2=NA,
      rt_sk_w1=0,rt_sk_iw=0,rt_sk_w2=0,
      prev_wave1=0,prev_inter=0,prev_wave2=0,prev_sp_w1=NA,prev_sp_iw=NA,prev_sp_w2=NA,
      prev_sk_w1=0,prev_sk_iw=0,prev_sk_w2=0,
      cases_rmse=NA,hosp_rmse=NA,ode_rmse=NA,R_end=NA,beta_min=NA,beta_max=NA,
      neg_rt_rmse=NA,neg_J=NA
    )
  }
}
abl <- rbindlist(abl_rows)
fwrite(abl, out("ablation_results.csv"))
save.image(out("gp_seir_checkpoint_ablation.RData"))
cat(sprintf("\n  %d/%d subsets converged\n", sum(is.finite(abl$rho)), nrow(abl)))


# ─────────────────────────────────────────────────────────────────
# STEP 3: SHAPLEY + BANZHAF VALUES (all metrics)
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 3: Shapley and Banzhaf values ===\n")

METRICS <- list(
  list(col="rt_corr",      name="Full-period Rt Pearson"),
  list(col="rt_spearman",  name="Full-period Rt Spearman"),
  list(col="rt_skill",     name="Full-period Rt skill score"),
  list(col="prev_corr",    name="Full-period Prev Pearson"),
  list(col="prev_spearman",name="Full-period Prev Spearman"),
  list(col="prev_skill",   name="Full-period Prev skill score"),
  list(col="rt_wave1",     name="Wave-1 Rt Pearson"),
  list(col="rt_inter",     name="Inter-wave Rt Pearson"),
  list(col="rt_wave2",     name="Wave-2 Rt Pearson"),
  list(col="prev_wave1",   name="Wave-1 Prev Pearson"),
  list(col="prev_inter",   name="Inter-wave Prev Pearson"),
  list(col="prev_wave2",   name="Wave-2 Prev Pearson"),
  list(col="rt_sk_w1",     name="Wave-1 Rt skill score"),
  list(col="rt_sk_w2",     name="Wave-2 Rt skill score"),
  list(col="rt_corr_pre",  name="Pre-June Rt Pearson"),
  list(col="rt_corr_post", name="Post-June Rt Pearson")
)

all_gt_rows <- list(); all_int_rows <- list()

for (m in METRICS) {
  gt <- game_theory_for(abl, m$col, m$name)
  all_gt_rows[[m$col]] <- gt$values
  all_int_rows[[m$col]] <- cbind(gt$interactions, data.table(metric=m$col, metric_name=m$name))

  # Efficiency check
  phi_sum <- sum(gt$phi, na.rm=TRUE)
  vN      <- gt$vN
  err     <- abs(phi_sum - vN)
  cat(sprintf("  %-38s v(N)=%+.4f sum(phi)=%+.4f err=%.1e\n", m$name, vN, phi_sum, err))
  if (err > 1e-3) cat(sprintf("    WARNING: efficiency violation! err=%.4f\n", err))
}

gt_all  <- rbindlist(all_gt_rows)
int_all <- rbindlist(all_int_rows)

fwrite(gt_all,  out("shapley_banzhaf_values.csv"))
fwrite(int_all, out("shapley_interactions_all.csv"))

# Backward-compatible Shapley-only table (wide format)
shap_wide <- dcast(gt_all[metric %in% c("rt_corr","rt_spearman","prev_corr",
                                         "rt_wave1","rt_wave2","prev_wave1","prev_wave2")],
                   source ~ metric, value.var="shapley")
fwrite(shap_wide, out("shapley_values.csv"))

# Shapley/Banzhaf comparison + rank agreement (primary metrics only)
cat("\n  Shapley vs Banzhaf rank-1 agreement:\n")
for (m in c("rt_corr","prev_corr","rt_skill","prev_skill")) {
  sub <- gt_all[metric==m]
  if (nrow(sub)==0 || !any(is.finite(sub$shapley))) next
  r1_phi <- sub$source[which.max(sub$shapley)]
  nb_ok  <- sub[is.finite(norm_banzhaf) & banzhaf_stable==TRUE]
  r1_nb  <- if (nrow(nb_ok)>0) nb_ok$source[which.max(nb_ok$norm_banzhaf)] else "unstable"
  cat(sprintf("    %-35s Shapley=%s  Banzhaf=%s  %s\n",
              m, r1_phi, r1_nb, if(r1_phi==r1_nb)"✓" else "✗"))
}

# Leave-one-out (derived from 3-source subsets)
loo_rows <- list()
for (m in c("rt_corr","prev_corr")) {
  vN <- v_lookup(abl,SOURCE_NAMES,m)
  for (src in SOURCE_NAMES) {
    rem <- setdiff(SOURCE_NAMES,src)
    loo_rows[[length(loo_rows)+1]] <- data.table(
      metric=m, dropped=src,
      all4=vN, without=v_lookup(abl,rem,m),
      drop_value=vN-v_lookup(abl,rem,m))
  }
}
fwrite(rbindlist(loo_rows), out("leave_one_out.csv"))

# Superadditivity
super_rows <- list()
for (m in c("rt_corr","prev_corr","rt_skill","prev_skill")) {
  vN    <- v_lookup(abl,SOURCE_NAMES,m)
  v_sing <- sapply(SOURCE_NAMES, function(s) v_lookup(abl,s,m))
  best_pair <- max(sapply(combn(4,2,simplify=FALSE), function(idx)
    v_lookup(abl,SOURCE_NAMES[idx],m)))
  super_rows[[length(super_rows)+1]] <- data.table(
    metric=m, v_all=vN, sum_singles=sum(v_sing),
    best_single=max(v_sing), best_single_source=SOURCE_NAMES[which.max(v_sing)],
    best_pair=best_pair,
    ratio=vN/sum(v_sing), redundancy_pct=(1-vN/sum(v_sing))*100)
}
fwrite(rbindlist(super_rows), out("superadditivity.csv"))

cat(sprintf("\n  Full-period Rt Pearson Shapley: %s\n",
  paste(sprintf("%s=%+.3f",SOURCE_NAMES,gt_all[metric=="rt_corr"]$shapley),collapse=" ")))
cat(sprintf("  Full-period Rt skill Shapley:   %s\n",
  paste(sprintf("%s=%+.3f",SOURCE_NAMES,gt_all[metric=="rt_skill"]$shapley),collapse=" ")))


# ─────────────────────────────────────────────────────────────────
# STEP 4: PENALTY WEIGHT SENSITIVITY
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 4: W_DATA sensitivity ===\n")
set.seed(SEEDS$wgrid)
W_DATA_GRID <- c(0.01,1,100,1e4,1e6,1e8)
wsens_rows <- list(); wsens_fits <- list()
for (wd in W_DATA_GRID) {
  cat(sprintf("  W_DATA=%.0e ... ", wd))
  r <- tryCatch(cascade(main$theta,main$C_mat,wd,wd,wd,wd),error=function(e)NULL)
  if (!is.null(r)) {
    d <- extract_full(r)
    wsens_fits[[as.character(wd)]] <- d
    cat(sprintf("Rt=%.3f prev=%.3f skill=%.3f\n",d$rt_corr,d$prev_corr,d$rt_skill))
    wsens_rows[[length(wsens_rows)+1]] <- data.table(
      W_DATA=wd,eff_ratio=sqrt(ETA*dt)/sqrt(wb_cases*wd),
      rt_corr=d$rt_corr,rt_spearman=d$rt_spearman,rt_rmse=d$rt_rmse,rt_skill=d$rt_skill,
      prev_corr=d$prev_corr,prev_spearman=d$prev_spearman,prev_skill=d$prev_skill,
      cases_rmse=d$cases_rmse,ode_rmse=d$ode_rmse,
      beta_min=d$beta_range[1],beta_max=d$beta_range[2],R_end=d$R_end)
  } else cat("FAILED\n")
}
wsens <- rbindlist(wsens_rows)
fwrite(wsens, out("wdata_sensitivity.csv"))
rt_wd_dt <- data.table(day=as.integer(tv),rt_rivm=rt_rivm)
for (nm in names(wsens_fits)) {
  d <- wsens_fits[[nm]]
  if (!is.null(d)) rt_wd_dt[[paste0("Rt_wd",nm)]] <- d$Rt
}
fwrite(rt_wd_dt, out("tikz_rt_wdata.csv"))
save.image(out("gp_seir_checkpoint_sensitivity.RData"))


# ─────────────────────────────────────────────────────────────────
# STEP 4b: BETA INITIALIZATION SENSITIVITY
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 4b: Beta initialization sensitivity ===\n")
set.seed(SEEDS$init)
init_rows <- list()

# (1) Hospital log-derivative — main fit
init_rows[[1]] <- data.table(init_method="hosp_logderiv",J=main$J,
  rho=main$rho,alphaR=main$alphaR,
  rt_corr=main$rt_corr,prev_corr=main$prev_corr,rt_skill=main$rt_skill,
  beta_min=main$beta_range[1],beta_max=main$beta_range[2],R_end=main$R_end)
cat(sprintf("  (1) hosp_logderiv (main): Rt=%.3f prev=%.3f skill=%.3f\n",
            main$rt_corr,main$prev_corr,main$rt_skill))

# (2) Constant beta = gamma
cat("  (2) constant_gamma ... ")
ai_c <- as.numeric(smooth.basis(tv,rep(log(gamma),n),fdPar(Bb,int2Lfd(2),1e-2))$fd$coefs)
thi_c <- c(ai_c,thi[(nbb+1):(nbb+4)])
best_c<-NULL; bJ_c<-Inf
for (s in 1:N_STARTS) {
  ts<-thi_c; Cs<-Ci
  if(s>1){ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,0.30);for(k in(nbb+1):(nbb+4))ts[k]<-ts[k]+rnorm(1,0,0.5)}
  r<-tryCatch(cascade(ts,Cs,1,1,1,1),error=function(e)NULL)
  if(!is.null(r)&&r$J<bJ_c){bJ_c<-r$J;best_c<-r}
}
if(!is.null(best_c)){d<-extract_full(best_c)
  cat(sprintf("Rt=%.3f prev=%.3f skill=%.3f\n",d$rt_corr,d$prev_corr,d$rt_skill))
  init_rows[[2]]<-data.table(init_method="constant_gamma",J=best_c$J,
    rho=d$rho,alphaR=d$alphaR,rt_corr=d$rt_corr,prev_corr=d$prev_corr,rt_skill=d$rt_skill,
    beta_min=d$beta_range[1],beta_max=d$beta_range[2],R_end=d$R_end)
}else cat("FAILED\n")

# (3) Random uniform
cat("  (3) random_uniform ... ")
set.seed(SEEDS$sg_inner)
b0_r<-runif(n,0.05,0.25)
ai_r<-as.numeric(smooth.basis(tv,log(b0_r),fdPar(Bb,int2Lfd(2),1e-2))$fd$coefs)
thi_r<-c(ai_r,thi[(nbb+1):(nbb+4)])
best_r<-NULL; bJ_r<-Inf
for(s in 1:N_STARTS){ts<-thi_r;Cs<-Ci
  if(s>1){ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,0.40);for(k in(nbb+1):(nbb+4))ts[k]<-ts[k]+rnorm(1,0,0.5)}
  r<-tryCatch(cascade(ts,Cs,1,1,1,1),error=function(e)NULL)
  if(!is.null(r)&&r$J<bJ_r){bJ_r<-r$J;best_r<-r}}
if(!is.null(best_r)){d<-extract_full(best_r)
  cat(sprintf("Rt=%.3f prev=%.3f skill=%.3f\n",d$rt_corr,d$prev_corr,d$rt_skill))
  init_rows[[3]]<-data.table(init_method="random_uniform",J=best_r$J,
    rho=d$rho,alphaR=d$alphaR,rt_corr=d$rt_corr,prev_corr=d$prev_corr,rt_skill=d$rt_skill,
    beta_min=d$beta_range[1],beta_max=d$beta_range[2],R_end=d$R_end)
}else cat("FAILED\n")

beta_init_dt <- rbindlist(Filter(Negate(is.null),init_rows))
fwrite(beta_init_dt, out("beta_init_sensitivity.csv"))


# ─────────────────────────────────────────────────────────────────
# STEP 4c: PROFILE LIKELIHOOD
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 4c: Profile likelihood ===\n")
set.seed(SEEDS$profile)
param_grids <- list(
  rho    = seq(0.08, 0.45, length.out=8),
  pH     = seq(0.004, 0.025, length.out=7),
  pICU   = seq(0.001, 0.008, length.out=6),
  alphaR = seq(0.15, 0.80, length.out=8)
)
prof_rows <- list()
for (pname in names(param_grids)) {
  cat(sprintf("  Profiling %s (%d pts) ... ", pname, length(param_grids[[pname]])))
  pidx <- which(c("rho","pH","pICU","alphaR")==pname)
  for (pval in param_grids[[pname]]) {
    th_p <- main$theta; th_p[nbb+pidx] <- to_u(pval,bounds[[pname]][1],bounds[[pname]][2])
    th0<-th_p; Cc<-main$C_mat; bJ<-Inf
    for (oi in 1:min(20,MAX_OUTER)) {
      inn<-solve_inner(Cc,th0,1,1,1,1); Cc<-inn$C; Jc<-inn$J; if(Jc<bJ) bJ<-Jc
      eps<-1e-4; gr<-numeric(length(th0))
      for(j in seq_along(th0)){if(j==nbb+pidx){gr[j]<-0;next};tp<-th0;tp[j]<-tp[j]+eps;gr[j]<-(solve_inner(Cc,tp,1,1,1,1)$J-Jc)/eps}
      gn<-sqrt(sum(gr^2)); if(gn<1e-12) break
      d_dir<- -gr/gn; st<-STEP0; found<-FALSE
      for(ls in 1:8){tt<-th0+st*d_dir;tt[nbb+pidx]<-th_p[nbb+pidx]
        it<-solve_inner(Cc,tt,1,1,1,1)
        if(it$J<Jc-1e-4*st*gn){th0<-tt;Cc<-it$C;found<-TRUE;break};st<-st*0.5}
      if(!found){th0<-th0+1e-4*d_dir;th0[nbb+pidx]<-th_p[nbb+pidx]}
      if(gn<1e-4*max(1,sqrt(sum(th0^2)))) break
    }
    dp<-extract_full(list(theta=th0,C=Cc,J=bJ))
    prof_rows[[length(prof_rows)+1]]<-data.table(
      parameter=pname,fixed_value=pval,J=bJ,
      rt_corr=dp$rt_corr,prev_corr=dp$prev_corr,
      rt_skill=dp$rt_skill,prev_skill=dp$prev_skill,
      rho=dp$rho,pH=dp$pH,pICU=dp$pICU,alphaR=dp$alphaR,R_end=dp$R_end)
  }
  cat("done\n")
}
prof_dt <- rbindlist(prof_rows)
fwrite(prof_dt, out("profile_likelihood.csv"))
fwrite(prof_dt[,.(parameter,fixed_value,J,rt_corr,prev_corr,rt_skill)], out("tikz_profile.csv"))
save.image(out("gp_seir_checkpoint_profile.RData"))


# ─────────────────────────────────────────────────────────────────
# STEP 4d: σ/γ SENSITIVITY
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 4d: sigma/gamma sensitivity ===\n")
set.seed(SEEDS$sg_sweep)
sg_configs <- list(
  list(s=1/5.5,g=1/9.5,lab="default (1/5.5, 1/9.5)"),
  list(s=1/4,  g=1/9.5,lab="sigma=1/4"),
  list(s=1/7,  g=1/9.5,lab="sigma=1/7"),
  list(s=1/5.5,g=1/7,  lab="gamma=1/7"),
  list(s=1/5.5,g=1/12, lab="gamma=1/12")
)
sg_rows <- list()
sg_rows[[1]] <- data.table(sigma_val=sigma,gamma_val=gamma,label="default (1/5.5, 1/9.5)",
  J=main$J,rho=main$rho,pH=main$pH,pICU=main$pICU,alphaR=main$alphaR,
  rt_corr=main$rt_corr,prev_corr=main$prev_corr,rt_skill=main$rt_skill,R_end=main$R_end)

for (cfg in sg_configs[-1]) {
  cat(sprintf("  %-35s ", cfg$lab))
  # Temporarily substitute sigma/gamma using local environment
  sv <- cfg$s; gv <- cfg$g
  Tc_tmp <- 1/sv+1/gv
  b_tmp <- gv*(1+dlog_sm*Tc_tmp); b_tmp<-pmax(0.05,pmin(0.35,b_tmp))
  b_tmp<-as.numeric(rollmean(b_tmp,14,fill=NA,align="center")); b_tmp[is.na(b_tmp)]<-0.15
  ai_sg<-as.numeric(smooth.basis(tv,log(b_tmp),fdPar(Bb,int2Lfd(2),1e-2))$fd$coefs)
  thi_sg<-c(ai_sg,to_u(0.20,bounds$rho[1],bounds$rho[2]),to_u(0.012,bounds$pH[1],bounds$pH[2]),
            to_u(0.0025,bounds$pICU[1],bounds$pICU[2]),to_u(0.40,bounds$alphaR[1],bounds$alphaR[2]))
  # Custom solve_inner that uses local sv/gv
  solve_inner_sg <- function(Cs,th,wc=1,wh=1,wi=1,wr=1){
    av<-th[1:nbb];rho<-from_u(th[nbb+1],bounds$rho[1],bounds$rho[2]);pH<-from_u(th[nbb+2],bounds$pH[1],bounds$pH[2]);pICU<-from_u(th[nbb+3],bounds$pICU[1],bounds$pICU[2]);aR<-from_u(th[nbb+4],bounds$alphaR[1],bounds$alphaR[2]);beta<-exp(as.numeric(Pb%*%av))
    fn<-function(cf){C<-matrix(cf,nrow=nbs,ncol=4);Z<-Ph%*%C;dZ<-dPh%*%C;S<-exp(Z[,1]);E<-exp(Z[,2]);I<-exp(Z[,3]);R<-exp(Z[,4]);fc<-N_POP*rho*sv*E;fh<-N_POP*pH*gv*I;fi<-N_POP*pICU*gv*I;mc<-log1p(pmax(0,apply_delay(fc,h_case)));mh<-log1p(pmax(0,apply_delay(fh,h_hosp)));mi<-log1p(pmax(0,apply_delay(fi,h_icu)));mr<-aR*I;rc<-sqrt(wb_cases*wc)*ifelse(is.finite(y_cases),y_cases-mc,0);rh<-sqrt(wb_hosp*wh)*ifelse(is.finite(y_hosp),y_hosp-mh,0);ri<-sqrt(wb_icu*wi)*ifelse(is.finite(y_icu),y_icu-mi,0);rr<-sqrt(wb_radar*wr)*ifelse(is.finite(y_radar),y_radar-mr,0);seir_rhs_loc<-function(Z2,b)cbind(-b*exp(Z2[,3]),b*exp(Z2[,1]+Z2[,3]-Z2[,2])-sv,sv*exp(Z2[,2]-Z2[,3])-gv,gv*exp(Z2[,3]-Z2[,4]));rp<-sqrt(ETA*dt)*as.numeric(dZ-seir_rhs_loc(Z,beta));rm<-sqrt(KAPPA_MASS)*(S+E+I+R-1);rb<-sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX,0);rs<-sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb%*%av);c(rc,rh,ri,rr,rp,rm,rb,rs)}
    r<-nls.lm(par=as.numeric(Cs),fn=fn,control=nls.lm.control(maxiter=MAX_INNER,ftol=1e-10,ptol=1e-10));list(C=matrix(r$par,nrow=nbs,ncol=4),J=sum(r$fvec^2))
  }
  cascade_sg <- function(th0,C0){th<-th0;Cc<-C0;bJ<-Inf;for(oi in 1:MAX_OUTER){inn<-solve_inner_sg(Cc,th);Cc<-inn$C;Jc<-inn$J;if(Jc<bJ)bJ<-Jc;eps<-1e-4;gr<-numeric(length(th));for(j in seq_along(th)){tp<-th;tp[j]<-tp[j]+eps;gr[j]<-(solve_inner_sg(Cc,tp)$J-Jc)/eps};gn<-sqrt(sum(gr^2));if(gn>1e-12){d<- -gr/gn;st<-STEP0;found<-FALSE;for(ls in 1:10){tt<-th+st*d;it<-solve_inner_sg(Cc,tt);if(it$J<Jc-1e-4*st*gn){th<-tt;Cc<-it$C;found<-TRUE;break};st<-st*0.5};if(!found)th<-th+1e-4*d};if(gn<1e-5*max(1,sqrt(sum(th^2))))break};list(theta=th,C=Cc,J=bJ)}
  best_sg<-NULL; bJ_sg<-Inf
  for(s in 1:N_STARTS){ts<-thi_sg;Cs<-Ci;if(s>1){ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,BETA_PERTURB);for(k in(nbb+1):(nbb+4))ts[k]<-ts[k]+rnorm(1,0,PARAM_PERTURB)};r<-tryCatch(cascade_sg(ts,Cs),error=function(e)NULL);if(!is.null(r)&&r$J<bJ_sg){bJ_sg<-r$J;best_sg<-r}}
  if(!is.null(best_sg)){
    # Compute Rt with local gamma for this result
    th<-best_sg$theta;C<-best_sg$C;av<-th[1:nbb];Z<-Ph%*%C;S<-exp(Z[,1]);Rt_sg<-exp(as.numeric(Pb%*%av))*S/gv
    ok_rt2<-is.finite(rt_rivm);rt_c<-safe_cor(rt_rivm,Rt_sg);rt_sk<-safe_skill(Rt_sg,rt_rivm)
    rho_sg<-from_u(th[nbb+1],bounds$rho[1],bounds$rho[2])
    cat(sprintf("Rt=%.3f skill=%.3f J=%.0f\n",rt_c,rt_sk,best_sg$J))
    sg_rows[[length(sg_rows)+1]]<-data.table(sigma_val=sv,gamma_val=gv,label=cfg$lab,
      J=best_sg$J,rho=rho_sg,pH=NA,pICU=NA,alphaR=NA,
      rt_corr=rt_c,prev_corr=NA,rt_skill=rt_sk,R_end=NA)
  } else cat("FAILED\n")
}
sg_dt <- rbindlist(sg_rows, fill=TRUE)
fwrite(sg_dt, out("sigma_gamma_sensitivity.csv"))
fwrite(sg_dt, out("tikz_sigma_gamma.csv"))
cat(sprintf("  Verified: sigma=%.4f, gamma=%.4f (unchanged)\n",sigma,gamma))


# ─────────────────────────────────────────────────────────────────
# STEP 5: EPIESTIM / CORI COMPARISON
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 5: Cori-method Rt comparison ===\n")
cori_dt <- run_cori(cases_raw, dates, si_mean=CORI_SI_MEAN, si_sd=CORI_SI_SD, tau=CORI_TAU)

if (!is.null(cori_dt) && nrow(cori_dt)>0) {
  # Three-way comparison: GP-SEIR vs Cori vs RIVM
  gp_dt <- data.table(date=dates, Rt_gp=main$Rt, Rt_rivm=rt_rivm)
  merged3 <- merge(gp_dt, cori_dt, by="date", all.x=TRUE)
  fwrite(merged3, out("tikz_rt_threeway.csv"))

  ok_gc <- is.finite(merged3$Rt_gp)  & is.finite(merged3$Rt_cori)
  ok_gr <- is.finite(merged3$Rt_gp)  & is.finite(merged3$Rt_rivm)
  ok_cr <- is.finite(merged3$Rt_cori)& is.finite(merged3$Rt_rivm)

  comp3 <- rbindlist(list(
    data.table(comparison="GP_vs_Cori",
      pearson=safe_cor(merged3$Rt_gp,merged3$Rt_cori),
      spearman=safe_cor(merged3$Rt_gp,merged3$Rt_cori,"spearman"),
      RMSE=safe_rmse(merged3$Rt_gp,merged3$Rt_cori),
      bias=if(sum(ok_gc)>3) mean(merged3$Rt_gp[ok_gc]-merged3$Rt_cori[ok_gc]) else NA,
      skill=safe_skill(merged3$Rt_gp[ok_gc],merged3$Rt_cori[ok_gc]),
      threshold_agree=if(sum(ok_gc)>3) mean((merged3$Rt_gp[ok_gc]>1)==(merged3$Rt_cori[ok_gc]>1)) else NA,
      n_days=sum(ok_gc)),
    data.table(comparison="Cori_vs_RIVM",
      pearson=safe_cor(merged3$Rt_cori,merged3$Rt_rivm),
      spearman=safe_cor(merged3$Rt_cori,merged3$Rt_rivm,"spearman"),
      RMSE=safe_rmse(merged3$Rt_cori,merged3$Rt_rivm),
      bias=if(sum(ok_cr)>3) mean(merged3$Rt_cori[ok_cr]-merged3$Rt_rivm[ok_cr]) else NA,
      skill=safe_skill(merged3$Rt_cori[ok_cr],merged3$Rt_rivm[ok_cr]),
      threshold_agree=if(sum(ok_cr)>3) mean((merged3$Rt_cori[ok_cr]>1)==(merged3$Rt_rivm[ok_cr]>1)) else NA,
      n_days=sum(ok_cr)),
    data.table(comparison="GP_vs_RIVM",
      pearson=safe_cor(merged3$Rt_gp,merged3$Rt_rivm),
      spearman=safe_cor(merged3$Rt_gp,merged3$Rt_rivm,"spearman"),
      RMSE=safe_rmse(merged3$Rt_gp,merged3$Rt_rivm),
      bias=if(sum(ok_gr)>3) mean(merged3$Rt_gp[ok_gr]-merged3$Rt_rivm[ok_gr]) else NA,
      skill=safe_skill(merged3$Rt_gp[ok_gr],merged3$Rt_rivm[ok_gr]),
      threshold_agree=if(sum(ok_gr)>3) mean((merged3$Rt_gp[ok_gr]>1)==(merged3$Rt_rivm[ok_gr]>1)) else NA,
      n_days=sum(ok_gr))
  ))
  fwrite(comp3, out("rt_method_comparison.csv"))

  # Lag analysis: GP vs Cori
  lag_rows <- list()
  for (lag in -14:14) {
    n2 <- length(merged3$Rt_gp)
    if (lag>=0){idx_g<-1:(n2-lag);idx_c<-(1+lag):n2}else{idx_g<-(1-lag):n2;idx_c<-1:(n2+lag)}
    ok2 <- is.finite(merged3$Rt_gp[idx_g]) & is.finite(merged3$Rt_cori[idx_c])
    if (sum(ok2)>20)
      lag_rows[[length(lag_rows)+1]] <- data.table(
        lag=lag, pearson=safe_cor(merged3$Rt_gp[idx_g][ok2],merged3$Rt_cori[idx_c][ok2]),
        rmse=safe_rmse(merged3$Rt_gp[idx_g][ok2],merged3$Rt_cori[idx_c][ok2]))
  }
  if (length(lag_rows)>0) fwrite(rbindlist(lag_rows), out("tikz_rt_lag_comparison.csv"))

  for (i in 1:nrow(comp3))
    cat(sprintf("  %-15s r=%.3f RMSE=%.3f bias=%.3f threshold=%.1f%%\n",
                comp3$comparison[i],comp3$pearson[i],comp3$RMSE[i],
                comp3$bias[i],comp3$threshold_agree[i]*100))

  # Also save full Cori series (for EpiEstim figure)
  rt_epi <- merge(data.table(day=as.integer(tv),rt_rivm=rt_rivm,rt_model=main$Rt,date=dates),
                  cori_dt, by="date", all.x=TRUE)
  fwrite(rt_epi, out("tikz_epiestim.csv"))
} else {
  cat("  Cori estimation skipped or failed.\n")
  fwrite(data.table(day=as.integer(tv),rt_rivm=rt_rivm,rt_model=main$Rt,
                    Rt_cori=NA,Rt_cori_lo=NA,Rt_cori_hi=NA), out("tikz_epiestim.csv"))
}


# ─────────────────────────────────────────────────────────────────
# STEP 6: COMPREHENSIVE Rt AND PREVALENCE METRICS
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 6: Comprehensive validation metrics ===\n")

ok_rt2  <- is.finite(rt_rivm) & is.finite(main$Rt)
ok_prev2 <- is.finite(prev_rivm) & is.finite(N_POP*main$I)
r_r <- rt_rivm[ok_rt2]; r_m <- main$Rt[ok_rt2]
I_ct <- N_POP*main$I

rt_met <- rbindlist(list(
  data.table(metric="pearson",       value=safe_cor(rt_rivm,main$Rt)),
  data.table(metric="spearman",      value=safe_cor(rt_rivm,main$Rt,"spearman")),
  data.table(metric="RMSE",          value=safe_rmse(rt_rivm,main$Rt)),
  data.table(metric="bias",          value=if(sum(ok_rt2)>3) mean(r_m-r_r) else NA),
  data.table(metric="skill_score",   value=main$rt_skill),
  data.table(metric="R2",            value=if(sum(ok_rt2)>3){ss_r<-sum((r_r-r_m)^2);ss_t<-sum((r_r-mean(r_r))^2);1-ss_r/ss_t}else NA),
  data.table(metric="MAE",           value=if(sum(ok_rt2)>3) mean(abs(r_r-r_m)) else NA),
  data.table(metric="threshold_agree",value=main$threshold_agree),
  data.table(metric="pearson_pre_jun",value=main$rt_corr_pre),
  data.table(metric="pearson_post_jun",value=main$rt_corr_post),
  data.table(metric="bias_pre_jun",  value=main$rt_bias_pre),
  data.table(metric="bias_post_jun", value=main$rt_bias_post),
  data.table(metric="pearson_wave1", value=main$rt_wave1),
  data.table(metric="pearson_interwave",value=main$rt_inter),
  data.table(metric="pearson_wave2", value=main$rt_wave2),
  data.table(metric="spearman_wave1",value=main$rt_sp_wave1),
  data.table(metric="spearman_wave2",value=main$rt_sp_wave2),
  data.table(metric="skill_wave1",   value=main$rt_skill_w1),
  data.table(metric="skill_interwave",value=main$rt_skill_iw),
  data.table(metric="skill_wave2",   value=main$rt_skill_w2),
  data.table(metric="RMSE_wave1",    value=safe_rmse(rt_rivm[w1idx<-dates<"2020-06-01"&ok_rt2],main$Rt[w1idx])),
  data.table(metric="RMSE_wave2",    value=safe_rmse(rt_rivm[w2idx<-dates>="2020-10-01"&ok_rt2],main$Rt[w2idx]))
))
fwrite(rt_met, out("rt_comprehensive.csv"))

prev_met <- rbindlist(list(
  data.table(metric="pearson",      value=main$prev_corr),
  data.table(metric="spearman",     value=main$prev_spearman),
  data.table(metric="skill_score",  value=main$prev_skill),
  data.table(metric="R2",           value=if(sum(ok_prev2)>3){ss_p<-sum((prev_rivm[ok_prev2]-I_ct[ok_prev2])^2);ss_tp<-sum((prev_rivm[ok_prev2]-mean(prev_rivm[ok_prev2]))^2);1-ss_p/ss_tp}else NA),
  data.table(metric="RMSE",         value=safe_rmse(prev_rivm,I_ct)),
  data.table(metric="bias",         value=if(sum(ok_prev2)>3) mean(I_ct[ok_prev2]-prev_rivm[ok_prev2]) else NA),
  data.table(metric="scale_ratio",  value=if(sum(ok_prev2)>3) mean(I_ct[ok_prev2])/mean(prev_rivm[ok_prev2]) else NA),
  data.table(metric="pearson_wave1",value=main$prev_wave1),
  data.table(metric="pearson_interwave",value=main$prev_inter),
  data.table(metric="pearson_wave2",value=main$prev_wave2),
  data.table(metric="skill_wave1",  value=main$prev_sk_w1),
  data.table(metric="skill_wave2",  value=main$prev_sk_w2)
))
fwrite(prev_met, out("prev_comprehensive.csv"))

cat(sprintf("  Rt:   r=%.3f, skill=%.3f, RMSE=%.3f, bias=%.3f, threshold=%.1f%%\n",
            main$rt_corr,main$rt_skill,main$rt_rmse,main$rt_bias,main$threshold_agree*100))
cat(sprintf("  Prev: r=%.3f, skill=%.3f\n",main$prev_corr,main$prev_skill))


# ─────────────────────────────────────────────────────────────────
# STEP 7: PER-WAVE REFITTING (independent sub-period fits)
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 7: Per-wave sliced and refitted diagnostics ===\n")
wave_periods <- list(
  list(start="2020-03-01",end="2020-06-01",label="wave1 (Mar-Jun)"),
  list(start="2020-06-01",end="2020-10-01",label="inter-wave (Jun-Oct)"),
  list(start="2020-10-01",end="2021-03-01",label="wave2 (Oct-Mar)")
)

# Sliced from main fit
wave_sliced_rows <- list()
for (p in wave_periods) {
  idx <- dates>=as.Date(p$start) & dates<as.Date(p$end)
  ok_r2<-is.finite(rt_rivm)&idx; ok_p2<-is.finite(prev_rivm)&idx
  wave_sliced_rows[[p$label]] <- data.table(
    period=p$label, n_days=sum(idx), type="sliced",
    rho=main$rho,pH=main$pH,pICU=main$pICU,alphaR=main$alphaR,
    rt_corr=safe_cor(rt_rivm[ok_r2],main$Rt[ok_r2]),
    rt_skill=safe_skill(main$Rt[ok_r2],rt_rivm[ok_r2]),
    prev_corr=safe_cor(prev_rivm[ok_p2],I_ct[ok_p2]),
    prev_skill=safe_skill(I_ct[ok_p2],prev_rivm[ok_p2]),
    beta_mean=mean(main$beta[idx]),R_end=main$R[max(which(idx))])
}
fwrite(rbindlist(wave_sliced_rows), out("wave_params.csv"))

# Independent subperiod refit function
fit_subperiod <- function(date_start, date_end, label) {
  cat(sprintf("  %-30s ", label))
  idx   <- dates >= as.Date(date_start) & dates < as.Date(date_end)
  n_sub <- sum(idx); if (n_sub < 40) { cat("too few days\n"); return(NULL) }
  tv_s  <- as.numeric(dates[idx]-as.Date(date_start)); rng_s <- range(tv_s)
  mkb_s <- function(kd,no){kn<-unique(c(rng_s[1],seq(rng_s[1],rng_s[2],by=kd),rng_s[2]));create.bspline.basis(rng_s,length(kn)+no-2,no,kn)}
  Bs_s<-mkb_s(STATE_KNOT_DAYS,STATE_NORDER);Bb_s<-mkb_s(BETA_KNOT_DAYS,BETA_NORDER)
  nbs_s<-Bs_s$nbasis;nbb_s<-Bb_s$nbasis
  Ph_s<-eval.basis(tv_s,Bs_s,0);dPh_s<-eval.basis(tv_s,Bs_s,1)
  Pb_s<-eval.basis(tv_s,Bb_s,0);d2Pb_s<-eval.basis(tv_s,Bb_s,2)
  yc_s<-y_cases[idx];yh_s<-y_hosp[idx];yi_s<-y_icu[idx];yr_s<-y_radar[idx]
  rc_s<-cases_raw[idx];rh_s<-hosp_raw[idx];rt_s<-rt_rivm[idx];pv_s<-prev_rivm[idx]
  wbc_s<-wsd(yc_s);wbh_s<-wsd(yh_s);wbi_s<-wsd(yi_s);wbr_s<-wsd(yr_s)
  S_i<-main$S[idx];E_i<-main$E[idx];I_i<-main$I[idx];R_i<-main$R[idx];b_i<-main$beta[idx]
  Ci_s<-tryCatch(smooth.basis(tv_s,cbind(log(S_i),log(E_i),log(I_i),log(R_i)),fdPar(Bs_s,int2Lfd(2),1e-4))$fd$coefs,error=function(e)matrix(0,nbs_s,4))
  ai_s<-tryCatch(as.numeric(smooth.basis(tv_s,log(b_i),fdPar(Bb_s,int2Lfd(2),1e-2))$fd$coefs),error=function(e)rep(-1.5,nbb_s))
  thi_s<-c(ai_s,to_u(main$rho,bounds$rho[1],bounds$rho[2]),to_u(main$pH,bounds$pH[1],bounds$pH[2]),to_u(main$pICU,bounds$pICU[1],bounds$pICU[2]),to_u(main$alphaR,bounds$alphaR[1],bounds$alphaR[2]))
  si_fn<-function(Cs,th){av<-th[1:nbb_s];rho<-from_u(th[nbb_s+1],bounds$rho[1],bounds$rho[2]);pH<-from_u(th[nbb_s+2],bounds$pH[1],bounds$pH[2]);pICU<-from_u(th[nbb_s+3],bounds$pICU[1],bounds$pICU[2]);aR<-from_u(th[nbb_s+4],bounds$alphaR[1],bounds$alphaR[2]);beta<-exp(as.numeric(Pb_s%*%av));fn<-function(cf){C<-matrix(cf,nrow=nbs_s,ncol=4);Z<-Ph_s%*%C;dZ<-dPh_s%*%C;S<-exp(Z[,1]);E<-exp(Z[,2]);I<-exp(Z[,3]);R<-exp(Z[,4]);fc<-N_POP*rho*sigma*E;fh<-N_POP*pH*gamma*I;fi<-N_POP*pICU*gamma*I;mc<-log1p(pmax(0,apply_delay(fc,h_case)));mh<-log1p(pmax(0,apply_delay(fh,h_hosp)));mi<-log1p(pmax(0,apply_delay(fi,h_icu)));mr<-aR*I;rc<-sqrt(wbc_s)*ifelse(is.finite(yc_s),yc_s-mc,0);rh<-sqrt(wbh_s)*ifelse(is.finite(yh_s),yh_s-mh,0);ri<-sqrt(wbi_s)*ifelse(is.finite(yi_s),yi_s-mi,0);rr<-sqrt(wbr_s)*ifelse(is.finite(yr_s),yr_s-mr,0);rp<-sqrt(ETA*dt)*as.numeric(dZ-seir_rhs(Z,beta));rm<-sqrt(KAPPA_MASS)*(S+E+I+R-1);rb<-sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX,0);rs<-sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb_s%*%av);c(rc,rh,ri,rr,rp,rm,rb,rs)};r<-nls.lm(par=as.numeric(Cs),fn=fn,control=nls.lm.control(maxiter=MAX_INNER,ftol=1e-10,ptol=1e-10));list(C=matrix(r$par,nrow=nbs_s,ncol=4),J=sum(r$fvec^2))}
  casc_s<-function(th0,C0){th<-th0;Cc<-C0;bJ<-Inf;for(oi in 1:MAX_OUTER){inn<-si_fn(Cc,th);Cc<-inn$C;Jc<-inn$J;if(Jc<bJ)bJ<-Jc;eps<-1e-4;gr<-numeric(length(th));for(j in seq_along(th)){tp<-th;tp[j]<-tp[j]+eps;gr[j]<-(si_fn(Cc,tp)$J-Jc)/eps};gn<-sqrt(sum(gr^2));if(gn>1e-12){d<- -gr/gn;st<-STEP0;found<-FALSE;for(ls in 1:10){tt<-th+st*d;it<-si_fn(Cc,tt);if(it$J<Jc-1e-4*st*gn){th<-tt;Cc<-it$C;found<-TRUE;break};st<-st*0.5};if(!found)th<-th+1e-4*d};if(gn<1e-5*max(1,sqrt(sum(th^2))))break};list(theta=th,C=Cc,J=bJ)}
  best_s<-NULL;bJ_s<-Inf
  for(s in 1:(N_STARTS+N_COLD_STARTS)){ts<-thi_s;Cs<-Ci_s;if(s>1){ts[1:nbb_s]<-ts[1:nbb_s]+rnorm(nbb_s,0,if(s<=N_STARTS)BETA_PERTURB else 0.5);for(k in(nbb_s+1):(nbb_s+4))ts[k]<-ts[k]+rnorm(1,0,if(s<=N_STARTS)PARAM_PERTURB else 0.8)};r<-tryCatch(casc_s(ts,Cs),error=function(e)NULL);if(!is.null(r)&&r$J<bJ_s){bJ_s<-r$J;best_s<-r}}
  if(is.null(best_s)){cat("FAILED\n");return(NULL)}
  th<-best_s$theta;C<-best_s$C;av<-th[1:nbb_s]
  rho<-from_u(th[nbb_s+1],bounds$rho[1],bounds$rho[2]);pH<-from_u(th[nbb_s+2],bounds$pH[1],bounds$pH[2]);pICU<-from_u(th[nbb_s+3],bounds$pICU[1],bounds$pICU[2]);aR<-from_u(th[nbb_s+4],bounds$alphaR[1],bounds$alphaR[2])
  Z<-Ph_s%*%C;S<-exp(Z[,1]);E<-exp(Z[,2]);I<-exp(Z[,3]);R<-exp(Z[,4]);beta<-exp(as.numeric(Pb_s%*%av));Rt_s<-beta*S/gamma
  ok_rt3<-is.finite(rt_s);ok_pv3<-is.finite(pv_s)
  rt_c3<-safe_cor(rt_s[ok_rt3],Rt_s[ok_rt3]);rt_sk3<-safe_skill(Rt_s[ok_rt3],rt_s[ok_rt3])
  pv_c3<-safe_cor(pv_s[ok_pv3],N_POP*I[ok_pv3]);pv_sk3<-safe_skill(N_POP*I[ok_pv3],pv_s[ok_pv3])
  cat(sprintf("rho=%.3f pH=%.4f alphaR=%.3f Rt=%.3f skill=%.3f\n",rho,pH,aR,rt_c3,rt_sk3))
  data.table(period=label,n_days=n_sub,type="refit",rho=rho,pH=pH,pICU=pICU,alphaR=aR,
    rt_corr=rt_c3,rt_skill=rt_sk3,prev_corr=pv_c3,prev_skill=pv_sk3,
    beta_mean=mean(beta),beta_min=min(beta),beta_max=max(beta),R_end=R[n_sub],
    mass_err=max(abs(S+E+I+R-1)))
}

set.seed(SEEDS$waves)
wave_refit_rows <- lapply(wave_periods, function(p)
  tryCatch(fit_subperiod(p$start,p$end,p$label), error=function(e){cat("  ERROR:",e$message,"\n");NULL}))
wave_refit_dt <- rbindlist(Filter(Negate(is.null),wave_refit_rows))
if (nrow(wave_refit_dt)>0) fwrite(wave_refit_dt, out("wave_refit_params.csv"))
save.image(out("gp_seir_checkpoint_waves.RData"))


# ─────────────────────────────────────────────────────────────────
# STEP 8: ALL TikZ EXPORTS
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 8: TikZ exports ===\n")
day <- as.integer(tv)
wave_label <- ifelse(dates<"2020-06-01","wave1",ifelse(dates<"2020-10-01","inter","wave2"))

fwrite(data.table(day=day,date=dates,S=main$S,E=main$E,I=main$I,R=main$R,beta=main$beta,Rt_model=main$Rt,rt_rivm=rt_rivm),out("tikz_seir.csv"))
fwrite(data.table(day=day,cases_raw=cases_raw,mu_cases=main$mu_cases),out("tikz_cases.csv"))
fwrite(data.table(day=day,hosp_raw=hosp_raw,mu_hosp=main$mu_hosp),out("tikz_hosp.csv"))
fwrite(data.table(day=day,icu_raw=icu_raw,mu_icu=main$mu_icu),out("tikz_icu.csv"))
fwrite(data.table(day=day,radar_frac=radar_frac,mu_radar=main$mu_radar),out("tikz_radar.csv"))
fwrite(data.table(day=day,rt_rivm=rt_rivm,rt_model=main$Rt),out("tikz_rt.csv"))
fwrite(data.table(day=day,prev_rivm=prev_rivm,I_model=N_POP*main$I),out("tikz_prev.csv"))
fwrite(data.table(day=day,prev_rivm=prev_rivm,I_model=N_POP*main$I,wave=wave_label),out("tikz_prev_comparison.csv"))
# Residuals and ODE
fwrite(data.table(day=day,resid_cases=main$resid_cases,resid_hosp=main$resid_hosp,resid_icu=main$resid_icu,resid_radar=main$resid_radar,ode_S=main$ode_mat[,1],ode_E=main$ode_mat[,2],ode_I=main$ode_mat[,3],ode_R=main$ode_mat[,4]),out("tikz_residuals.csv"))
rc_clean <- main$resid_cases[is.finite(main$resid_cases)]
if(length(rc_clean)>20){acf_obj<-acf(rc_clean,lag.max=30,plot=FALSE);fwrite(data.table(lag=as.numeric(acf_obj$lag),acf=as.numeric(acf_obj$acf)),out("tikz_acf.csv"));n_qq<-length(rc_clean);fwrite(data.table(theoretical=qnorm(ppoints(n_qq)),sample=sort(rc_clean)),out("tikz_qq.csv"))}
# Validation scatter
fwrite(data.table(rt_rivm=rt_rivm,rt_model=main$Rt,wave=wave_label,day=day),out("tikz_rt_scatter.csv"))
ok_ba<-is.finite(rt_rivm)&is.finite(main$Rt)
fwrite(data.table(day=day,mean_rt=(rt_rivm+main$Rt)/2,diff_rt=main$Rt-rt_rivm,wave=wave_label,ok=as.integer(ok_ba)),out("tikz_bland_altman.csv"))
# Detection
fwrite(data.table(day=day,rho_eff=main$rho_eff,rho_const=main$rho),out("tikz_detection.csv"))
fwrite(data.table(day=day,R_cumulative=main$R,R_pct=main$R*100),out("tikz_cumulative.csv"))
# Ablation
fwrite(abl[order(-abl$prev_corr),.(sources,rt_corr,prev_corr,rt_skill,prev_skill,J)],out("tikz_ablation.csv"))
fwrite(abl[order(-abl$rt_corr),.(sources,rt_corr,prev_corr,rt_skill,J)],out("tikz_ablation_rt.csv"))
# Beta ablation
if(length(abl_betas)>0){beta_dt<-data.table(day=day);for(nm in names(abl_betas)) beta_dt[[gsub("[^a-zA-Z]","_",nm)]]<-abl_betas[[nm]];fwrite(beta_dt,out("tikz_beta_ablation.csv"))}
# Beta init
fwrite(data.table(day=day,beta_init=beta_init_curve,beta_fitted=main$beta),out("tikz_beta_init.csv"))
# Testing
fwrite(data.table(day=day,tests_total=tests_tot,tests_positive=tests_pos,positivity=positivity),out("tikz_testing.csv"))
# Mass conservation
fwrite(data.table(day=day,date=dates,mass_error=main$S+main$E+main$I+main$R-1),out("tikz_mass_conservation.csv"))
# Shapley + Banzhaf bar chart (primary metrics, wide format)
primary_gt <- gt_all[metric %in% c("rt_corr","rt_skill","prev_corr","prev_skill")]
fwrite(primary_gt, out("tikz_shapley_banzhaf.csv"))
fwrite(gt_all[metric=="rt_corr"],  out("tikz_shapley.csv"))  # backward compat
# Interactions
int_primary <- int_all[metric %in% c("rt_corr","prev_corr")]
fwrite(int_primary, out("tikz_shapley_interactions.csv"))
# Rt Pearson vs skill tradeoff
rt_tradeoff <- abl[,.(sources,rt_corr,rt_skill,prev_corr,prev_skill)]
fwrite(rt_tradeoff, out("tikz_rt_prev_tradeoff.csv"))
cat("  All TikZ CSVs saved.\n")


# ─────────────────────────────────────────────────────────────────
# STEP 9: DERIVED QUANTITIES (comprehensive)
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 9: Derived quantities ===\n")

# Flows
flow_new_inf <- N_POP*main$beta*main$S*main$I
flow_onset   <- N_POP*sigma*main$E
flow_removal <- N_POP*gamma*main$I
growth_rate  <- c(diff(log(pmax(main$I,1e-15))),NA)
doubling_t   <- ifelse(growth_rate>0,log(2)/growth_rate,NA)
halving_t    <- ifelse(growth_rate<0,-log(2)/growth_rate,NA)

# Peaks
w1_idx<-which(dates<"2020-06-01"); iw_idx<-which(dates>="2020-06-01"&dates<"2020-10-01"); w2_idx<-which(dates>="2020-10-01")
peak_I_val  <- N_POP*max(main$I); peak_I_date <- dates[which.max(main$I)]
peak_I_w1   <- dates[w1_idx[which.max(main$I[w1_idx])]]
peak_I_w2   <- dates[w2_idx[which.max(main$I[w2_idx])]]
attack_w1   <- main$R[max(w1_idx)]-main$R[min(w1_idx)]
attack_iw   <- main$R[max(iw_idx)]-main$R[min(iw_idx)]
attack_w2   <- main$R[max(w2_idx)]-main$R[min(w2_idx)]
attack_tot  <- main$R[n]-main$R[1]
rt_crossings <- which(diff(main$Rt>1)!=0)
R0_range    <- main$beta/gamma

# Derived time series
fwrite(data.table(date=dates,day=day,new_infections=flow_new_inf,onset_flow=flow_onset,removal_flow=flow_removal,growth_rate=growth_rate,doubling_time=doubling_t,halving_time=halving_t,mass_error=main$S+main$E+main$I+main$R-1,rho_eff=main$rho_eff),out("derived_timeseries.csv"))

# Stream fit metrics
stream_met_rows <- list()
for (sname in c("cases","hosp","icu","radar")) {
  obs<-switch(sname,cases=cases_raw,hosp=hosp_raw,icu=icu_raw,radar=radar_frac)
  pred<-switch(sname,cases=main$mu_cases,hosp=main$mu_hosp,icu=main$mu_icu,radar=main$mu_radar)
  ok2<-is.finite(obs)&is.finite(pred)&(if(sname=="hosp")obs>0 else TRUE)
  if(sum(ok2)>5) stream_met_rows[[sname]]<-data.table(stream=sname,pearson=safe_cor(obs[ok2],pred[ok2]),spearman=safe_cor(obs[ok2],pred[ok2],"spearman"),RMSE=safe_rmse(obs[ok2],pred[ok2]),bias=mean(pred[ok2]-obs[ok2]),n_obs=sum(ok2))
}
fwrite(rbindlist(stream_met_rows), out("stream_fit_metrics.csv"))

# Cross-correlations
ccorr_rows<-list()
streams_l<-list(cases=y_cases,hosp=y_hosp,icu=y_icu,radar=y_radar)
for(i in 1:3) for(j in(i+1):4){x<-streams_l[[i]];y<-streams_l[[j]];ok2<-is.finite(x)&is.finite(y)
  if(sum(ok2)>30){ccf_obj<-ccf(x[ok2],y[ok2],lag.max=21,plot=FALSE);bi<-which.max(abs(ccf_obj$acf));ccorr_rows[[length(ccorr_rows)+1]]<-data.table(stream_i=names(streams_l)[i],stream_j=names(streams_l)[j],pearson_lag0=cor(x[ok2],y[ok2]),best_lag=ccf_obj$lag[bi],best_ccf=ccf_obj$acf[bi])}}
fwrite(rbindlist(ccorr_rows), out("cross_correlations.csv"))

# Residual diagnostics
resid_diag_rows<-list()
for(sname in c("cases","hosp","icu","radar")){rv<-switch(sname,cases=main$resid_cases,hosp=main$resid_hosp,icu=main$resid_icu,radar=main$resid_radar);rv2<-rv[is.finite(rv)];if(length(rv2)>20){sw<-tryCatch(shapiro.test(rv2[1:min(5000,length(rv2))]),error=function(e)list(statistic=NA,p.value=NA));acf1<-acf(rv2,lag.max=1,plot=FALSE)$acf[2];lb<-tryCatch(Box.test(rv2,lag=10,type="Ljung-Box"),error=function(e)list(statistic=NA,p.value=NA));resid_diag_rows[[sname]]<-data.table(stream=sname,n=length(rv2),mean=mean(rv2),sd=sd(rv2),skewness=mean((rv2-mean(rv2))^3)/sd(rv2)^3,kurtosis=mean((rv2-mean(rv2))^4)/sd(rv2)^4-3,acf_lag1=acf1,shapiro_W=as.numeric(sw$statistic),shapiro_p=sw$p.value,lb_stat=as.numeric(lb$statistic),lb_p=lb$p.value)}}
fwrite(rbindlist(resid_diag_rows), out("residual_diagnostics.csv"))

# Rt lag analysis
rt_lag_rows<-list()
for(shift in -10:20){if(shift>=0){r_sl<-rt_rivm[(shift+1):n];m_sl<-main$Rt[1:(n-shift)]}else{r_sl<-rt_rivm[1:(n+shift)];m_sl<-main$Rt[(1-shift):n]};ok2<-is.finite(r_sl)&is.finite(m_sl);if(sum(ok2)>50)rt_lag_rows[[length(rt_lag_rows)+1]]<-data.table(lag=shift,pearson=safe_cor(r_sl[ok2],m_sl[ok2]),spearman=safe_cor(r_sl[ok2],m_sl[ok2],"spearman"),rmse=safe_rmse(r_sl[ok2],m_sl[ok2]))}
fwrite(rbindlist(rt_lag_rows), out("rt_lag_analysis.csv"))

# ODE compliance
fwrite(data.table(compartment=c("S","E","I","R"),rmse=apply(main$ode_mat,2,function(x)sqrt(mean(x^2))),max_abs=apply(main$ode_mat,2,function(x)max(abs(x)))),out("ode_compliance.csv"))

# Wave summary
fwrite(data.table(
  period=c("wave1","inter-wave","wave2","total"),
  start=c("2020-03-01","2020-06-01","2020-10-01","2020-03-01"),
  end=c("2020-06-01","2020-10-01","2021-03-01","2021-03-01"),
  attack_rate=c(attack_w1,attack_iw,attack_w2,attack_tot),
  attack_rate_pct=c(attack_w1,attack_iw,attack_w2,attack_tot)*100,
  peak_I_date=c(as.character(peak_I_w1),NA,as.character(peak_I_w2),as.character(peak_I_date)),
  peak_I=c(N_POP*max(main$I[w1_idx]),N_POP*max(main$I[iw_idx]),N_POP*max(main$I[w2_idx]),peak_I_val),
  mean_Rt=c(mean(main$Rt[w1_idx]),mean(main$Rt[iw_idx]),mean(main$Rt[w2_idx]),mean(main$Rt)),
  max_Rt=c(max(main$Rt[w1_idx]),max(main$Rt[iw_idx]),max(main$Rt[w2_idx]),max(main$Rt)),
  pct_above_Rt1=c(mean(main$Rt[w1_idx]>1),mean(main$Rt[iw_idx]>1),mean(main$Rt[w2_idx]>1),mean(main$Rt>1))*100,
  total_cases=c(sum(cases_raw[w1_idx],na.rm=TRUE),sum(cases_raw[iw_idx],na.rm=TRUE),sum(cases_raw[w2_idx],na.rm=TRUE),sum(cases_raw,na.rm=TRUE)),
  total_hosp=c(sum(hosp_raw[w1_idx],na.rm=TRUE),sum(hosp_raw[iw_idx],na.rm=TRUE),sum(hosp_raw[w2_idx],na.rm=TRUE),sum(hosp_raw,na.rm=TRUE)),
  R_end=c(main$R[max(w1_idx)],main$R[max(iw_idx)],main$R[max(w2_idx)],main$R[n])
),out("wave_summary.csv"))

# Rt crossings
if(length(rt_crossings)>0) fwrite(data.table(day=rt_crossings,date=dates[rt_crossings],direction=ifelse(diff(main$Rt>1)[rt_crossings]>0,"rising","falling"),Rt_before=main$Rt[rt_crossings],Rt_after=main$Rt[pmin(rt_crossings+1,n)]),out("rt_crossings.csv"))

# Epi parameters
fwrite(data.table(parameter=c("sigma","gamma","gen_time","latent_period","infectious_period","N_pop","CORI_si_mean","CORI_si_sd","CORI_tau","use_delay_kernel","banzhaf_stable_threshold"),value=c(sigma,gamma,1/sigma+1/gamma,1/sigma,1/gamma,N_POP,CORI_SI_MEAN,CORI_SI_SD,CORI_TAU,as.numeric(USE_DELAY_KERNEL),BANZHAF_STABLE_THRESHOLD)),out("epi_parameters.csv"))

# Multi-start correlation
if(nrow(ms_df)>=3){ms_c<-cor(ms_df[,.(J,rho,pH,pICU,alphaR,rt_corr,prev_corr)],use="pairwise.complete.obs");fwrite(as.data.table(as.data.frame(ms_c),keep.rownames="variable"),out("multistart_correlation.csv"))}

# Prevalence
fwrite(data.table(day=day,date=dates,prev_rivm=prev_rivm,I_model=N_POP*main$I,wave=wave_label),out("tikz_prev_comparison.csv"))

cat(sprintf("  Peaks: I at %s (%.0f k), attack: w1=%.1f%% iw=%.1f%% w2=%.1f%%\n",
            peak_I_date,peak_I_val/1000,attack_w1*100,attack_iw*100,attack_w2*100))
cat(sprintf("  Rt crossings: %d, R0 range [%.2f, %.2f]\n",length(rt_crossings),min(R0_range),max(R0_range)))


# ─────────────────────────────────────────────────────────────────
# STEP 10: FULL RESULTS CSV + SUMMARY
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 10: Full results CSV and summary ===\n")

fwrite(data.table(
  date=dates,day=day,wave=wave_label,
  S=main$S,E=main$E,I=main$I,R=main$R,
  beta=main$beta,Rt_model=main$Rt,
  mu_cases=main$mu_cases,mu_hosp=main$mu_hosp,mu_icu=main$mu_icu,mu_radar=main$mu_radar,
  cases_raw=cases_raw,hosp_raw=hosp_raw,icu_raw=icu_raw,radar_frac=radar_frac,
  rt_rivm=rt_rivm,prev_rivm=prev_rivm,tests_total=tests_tot,positivity=positivity,
  resid_cases=main$resid_cases,resid_hosp=main$resid_hosp,resid_icu=main$resid_icu,resid_radar=main$resid_radar,
  ode_S=main$ode_mat[,1],ode_E=main$ode_mat[,2],ode_I=main$ode_mat[,3],ode_R=main$ode_mat[,4],
  rho_eff=main$rho_eff,
  new_infections=flow_new_inf,onset_flow=flow_onset,removal_flow=flow_removal,
  growth_rate=growth_rate,doubling_time=doubling_t,
  I_count=N_POP*main$I,E_count=N_POP*main$E,R_pct=main$R*100
),out("gp_seir_results.csv"))

# Comprehensive summary (all key numbers for paper writing)
sum_rows <- list(
  data.table(category="main",metric="J",value=main$J),
  data.table(category="main",metric="rho",value=main$rho),
  data.table(category="main",metric="pH",value=main$pH),
  data.table(category="main",metric="pICU",value=main$pICU),
  data.table(category="main",metric="alphaR",value=main$alphaR),
  data.table(category="main",metric="rt_pearson",value=main$rt_corr),
  data.table(category="main",metric="rt_spearman",value=main$rt_spearman),
  data.table(category="main",metric="rt_rmse",value=main$rt_rmse),
  data.table(category="main",metric="rt_bias",value=main$rt_bias),
  data.table(category="main",metric="rt_skill",value=main$rt_skill),
  data.table(category="main",metric="prev_pearson",value=main$prev_corr),
  data.table(category="main",metric="prev_spearman",value=main$prev_spearman),
  data.table(category="main",metric="prev_skill",value=main$prev_skill),
  data.table(category="main",metric="R_end",value=main$R_end),
  data.table(category="main",metric="beta_min",value=main$beta_range[1]),
  data.table(category="main",metric="beta_max",value=main$beta_range[2]),
  data.table(category="main",metric="ode_rmse",value=main$ode_rmse),
  data.table(category="main",metric="threshold_agree",value=main$threshold_agree),
  data.table(category="main",metric="delay_kernel",value=as.numeric(USE_DELAY_KERNEL))
)
# Shapley and Banzhaf for primary metrics
for (m in c("rt_corr","prev_corr","rt_skill","prev_skill")) {
  sub <- gt_all[metric==m]
  for(i in seq_len(nrow(sub))) {
    sum_rows[[length(sum_rows)+1]]<-data.table(category=paste0("shapley_",m),metric=sub$source[i],value=sub$shapley[i])
    sum_rows[[length(sum_rows)+1]]<-data.table(category=paste0("banzhaf_",m),metric=sub$source[i],value=sub$norm_banzhaf[i])
  }
}
# Superadditivity
for(i in seq_len(nrow(rbindlist(super_rows)))){r<-rbindlist(super_rows)[i];sum_rows[[length(sum_rows)+1]]<-data.table(category="superadditivity",metric=paste0(r$metric,"_redundancy_pct"),value=r$redundancy_pct)}
# Derived
for(nm_v in list(list("gen_time",1/sigma+1/gamma),list("peak_I_count",peak_I_val),list("attack_total",attack_tot),list("attack_wave1",attack_w1),list("attack_wave2",attack_w2),list("n_rt_crossings",length(rt_crossings)),list("pct_time_above_Rt1",mean(main$Rt>1)*100),list("R0_max",max(R0_range))))
  sum_rows[[length(sum_rows)+1]]<-data.table(category="derived",metric=nm_v[[1]],value=nm_v[[2]])
sum_df <- rbindlist(sum_rows, fill=TRUE)
fwrite(sum_df, out("summary.csv"))
cat("  summary.csv saved.\n")


# ─────────────────────────────────────────────────────────────────
# STEP 11: PDF PLOTS
# ─────────────────────────────────────────────────────────────────

cat("\n=== STEP 11: PDF plots ===\n")

lockdown1 <- as.Date("2020-03-15"); relax1 <- as.Date("2020-06-01")
lockdown2 <- as.Date("2020-10-14")
add_lock <- function(){usr<-par("usr");rect(lockdown1,usr[3],relax1,usr[4],col="#CC000020",border=NA);rect(lockdown2,usr[3],as.Date("2021-03-01"),usr[4],col="#CC000020",border=NA)}

pdf(out("gp_seir_plots.pdf"), width=14, height=10)

# Page 1: Main fit
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(dates,cases_raw,pch=16,cex=0.3,col="gray60",xlab="",ylab="Daily cases",main=sprintf("A. Cases (rho=%.3f)",main$rho));add_lock();lines(dates,main$mu_cases,col="red",lwd=2)
plot(dates,hosp_raw,pch=16,cex=0.3,col="gray60",xlab="",ylab="Admissions/d",main=sprintf("B. Hospital (pH=%.4f)",main$pH));add_lock();lines(dates,main$mu_hosp,col="purple",lwd=2)
ok_i<-is.finite(icu_raw);plot(dates[ok_i],icu_raw[ok_i],pch=16,cex=0.3,col="gray60",xlab="",ylab="Admissions/d",main=sprintf("C. ICU (pICU=%.4f)",main$pICU));add_lock();lines(dates,main$mu_icu,col="darkred",lwd=2)
ok_r<-is.finite(radar_frac);plot(dates[ok_r],radar_frac[ok_r],pch=16,cex=0.3,col="gray60",xlab="",ylab="I(t) fraction",main=sprintf("D. RADAR (alphaR=%.3f)",main$alphaR));add_lock();lines(dates,main$mu_radar,col="darkgreen",lwd=2)

# Page 2: Validation + SEIR
par(mfrow=c(2,2),mar=c(4,4,3,1))
ok_rt3<-is.finite(rt_rivm);plot(dates,rt_rivm,type="l",col="gray40",lwd=1,xlab="",ylab=expression(R[t]),main=sprintf("E. Rt validation (r=%.3f, skill=%.3f)",main$rt_corr,main$rt_skill),ylim=c(0,max(c(main$Rt,rt_rivm[ok_rt3])*1.05,na.rm=TRUE)));add_lock();lines(dates,main$Rt,col="darkgreen",lwd=2);abline(h=1,lty=2,col="gray")
ok_pr3<-is.finite(prev_rivm);if(any(ok_pr3)){plot(dates[ok_pr3],prev_rivm[ok_pr3]/1000,type="l",col="gray40",lwd=1,xlab="",ylab="Thousands",main=sprintf("F. Prevalence (r=%.3f, skill=%.3f)",main$prev_corr,main$prev_skill));add_lock();lines(dates,N_POP*main$I/1000,col="steelblue",lwd=2)}
plot(dates,main$beta,type="l",col="purple",lwd=2,xlab="",ylab=expression(beta(t)),main="G. Transmission rate");add_lock()
plot(dates,main$S,type="l",col="blue",lwd=2,xlab="",ylab="Fraction",main="H. SEIR",ylim=c(0,1));lines(dates,main$E*50,col="orange",lwd=2);lines(dates,main$I*50,col="red",lwd=2);lines(dates,main$R,col="darkgreen",lwd=2);legend("right",c("S","E×50","I×50","R"),col=c("blue","orange","red","darkgreen"),lty=1,lwd=2,cex=0.7)

# Page 3: Cori three-way Rt
if (file.exists(out("tikz_rt_threeway.csv"))) {
  tw <- fread(out("tikz_rt_threeway.csv"))
  par(mfrow=c(1,1),mar=c(4,4,3,1))
  ok3 <- is.finite(tw$Rt_rivm)
  plot(tw$date[ok3],tw$Rt_rivm[ok3],type="l",col="gray30",lwd=1.5,
       xlab="",ylab=expression(R[t]),
       main=sprintf("Rt comparison: GP (r=%.3f) vs Cori vs RIVM",main$rt_corr),
       ylim=c(0,max(tw$Rt_gp,tw$Rt_rivm[ok3],na.rm=TRUE)*1.1))
  add_lock(); abline(h=1,lty=2,col="gray70")
  lines(tw$date,tw$Rt_gp,col="darkgreen",lwd=2)
  ok_c2<-is.finite(tw$Rt_cori)
  if(any(ok_c2)) lines(tw$date[ok_c2],tw$Rt_cori[ok_c2],col="steelblue",lwd=2,lty=2)
  legend("topright",c("RIVM","GP-SEIR model","EpiEstim Cori"),col=c("gray30","darkgreen","steelblue"),lty=c(1,1,2),lwd=c(1.5,2,2),cex=0.8)
}

# Page 4: Shapley vs Banzhaf bar charts
par(mfrow=c(2,2),mar=c(5,5,3,2))
for(m_name in c("rt_corr","prev_corr","rt_skill","prev_skill")) {
  sub <- gt_all[metric==m_name]
  if(nrow(sub)<4) next
  phi_v  <- setNames(sub$shapley, sub$source)
  nb_v   <- setNames(sub$norm_banzhaf, sub$source)
  ylim_r <- range(c(phi_v,nb_v,0),na.rm=TRUE)*c(1.2,1.2)
  xp <- barplot(rbind(phi_v,nb_v[names(phi_v)]),beside=TRUE,
                col=c("steelblue","coral"),
                names.arg=names(phi_v),las=2,
                main=sprintf("%s\nShapley(blue) vs normBanzhaf(coral)",m_name),
                ylab="Value",ylim=ylim_r)
  abline(h=0,col="gray40")
  if(all(!is.na(nb_v))) legend("topright",c("Shapley","Norm.Banzhaf"),fill=c("steelblue","coral"),cex=0.7)
}

# Page 5: Residuals
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(dates,main$resid_cases,type="h",col="gray40",xlab="",ylab="Residual",main="I. Case residuals");abline(h=0,col="red")
plot(dates,main$resid_hosp,type="h",col="gray40",xlab="",ylab="Residual",main="J. Hospital residuals");abline(h=0,col="red")
if(length(rc_clean)>20){acf(rc_clean,lag.max=30,main="K. Case residual ACF");qqnorm(rc_clean,main="L. Case Q-Q",pch=16,cex=0.4,col="gray40");qqline(rc_clean,col="red")}

# Page 6: Profile likelihood
if(nrow(prof_dt)>0){
  par(mfrow=c(2,4),mar=c(4,4,3,1))
  for(pname in c("rho","pH","pICU","alphaR")) {
    for(mname in c("J","rt_skill")) {
      sub2<-prof_dt[parameter==pname]
      yv<-if(mname=="J") sub2$J-min(sub2$J) else sub2$rt_skill
      ylab2<-if(mname=="J") expression(Delta*J) else "Rt skill"
      plot(sub2$fixed_value,yv,type="b",pch=16,lwd=2,col="darkblue",xlab=pname,ylab=ylab2,main=sprintf("Profile %s: %s",pname,mname))
      abline(h=if(mname=="J")0 else main$rt_skill,col="red",lty=2)
    }
  }
}

# Page 7: W_DATA sensitivity
if(nrow(wsens)>1){par(mfrow=c(2,3),mar=c(4,4,3,1));x<-log10(wsens$W_DATA);plot(x,wsens$ode_rmse,type="b",pch=16,lwd=2,xlab=expression(log[10](W[data])),ylab="ODE RMSE",main="ODE compliance");plot(x,wsens$cases_rmse,type="b",pch=16,col="red",lwd=2,xlab=expression(log[10](W[data])),ylab="Cases RMSE",main="Data fit");plot(x,wsens$rt_corr,type="b",pch=16,col="darkgreen",lwd=2,xlab=expression(log[10](W[data])),ylab="Rt corr",main="Rt Pearson");plot(x,wsens$rt_skill,type="b",pch=16,col="purple",lwd=2,xlab=expression(log[10](W[data])),ylab="Rt skill",main="Rt skill score");plot(x,wsens$prev_corr,type="b",pch=16,col="steelblue",lwd=2,xlab=expression(log[10](W[data])),ylab="Prev corr",main="Prevalence");plot(x,wsens$R_end*100,type="b",pch=16,col="brown",lwd=2,xlab=expression(log[10](W[data])),ylab="R(end) %",main="Cumulative")}

# Page 8: Rt scatter + Bland-Altman
par(mfrow=c(1,2),mar=c(4,4,3,1))
cols_w<-ifelse(wave_label=="wave1","red",ifelse(wave_label=="inter","blue","darkgreen"))
ok2<-is.finite(rt_rivm)&is.finite(main$Rt)
plot(rt_rivm[ok2],main$Rt[ok2],pch=16,cex=0.5,col=cols_w[ok2],xlab=expression(R[t]^RIVM),ylab=expression(R[t]^model),main=sprintf("Rt scatter (r=%.3f, skill=%.3f)",main$rt_corr,main$rt_skill),xlim=c(0.3,2.3),ylim=c(0.3,2.3));abline(0,1,col="gray40",lty=2);abline(h=1,v=1,col="gray80",lty=3);legend("topleft",c("Wave1","Inter","Wave2"),col=c("red","blue","darkgreen"),pch=16,cex=0.7)
mr2<-(rt_rivm[ok2]+main$Rt[ok2])/2;dr2<-main$Rt[ok2]-rt_rivm[ok2];plot(mr2,dr2,pch=16,cex=0.4,col=cols_w[ok2],xlab="Mean Rt",ylab="Model−RIVM",main=sprintf("Bland-Altman (bias=%.3f)",mean(dr2)));abline(h=mean(dr2),col="red");abline(h=mean(dr2)+c(-1.96,1.96)*sd(dr2),col="red",lty=2)

# Page 9: Detection + cumulative
par(mfrow=c(2,1),mar=c(4,4,3,1))
ok_rho<-is.finite(main$rho_eff);if(any(ok_rho)){plot(dates[ok_rho],pmin(main$rho_eff[ok_rho],1),type="l",col="darkblue",lwd=2,xlab="",ylab=expression(rho[eff](t)),main=sprintf("Detection rate (rho=%.3f)",main$rho),ylim=c(0,1));add_lock();abline(h=main$rho,lty=2,col="red",lwd=1.5)}
plot(dates,main$R*100,type="l",col="darkgreen",lwd=2,xlab="",ylab="Recovered (%)",main=sprintf("Cumulative infections (R_end=%.1f%%)",main$R_end*100));add_lock()

dev.off()
cat(sprintf("  Saved: %s\n", out("gp_seir_plots.pdf")))


# ─────────────────────────────────────────────────────────────────
# FINAL WORKSPACE
# ─────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("ALL DONE.\n\n")
cat(sprintf("  Rt:   r=%.3f  rs=%.3f  RMSE=%.3f  bias=%.3f  skill=%.3f\n",
            main$rt_corr,main$rt_spearman,main$rt_rmse,main$rt_bias,main$rt_skill))
cat(sprintf("  Prev: r=%.3f  rs=%.3f  skill=%.3f\n",
            main$prev_corr,main$prev_spearman,main$prev_skill))
cat(sprintf("  Params: rho=%.3f pH=%.4f pICU=%.4f alphaR=%.3f\n",
            main$rho,main$pH,main$pICU,main$alphaR))
cat(sprintf("  R(end)=%.1f%%  peak I=%s (%.0f k)\n",
            main$R_end*100,peak_I_date,peak_I_val/1000))
cat("\n  Shapley + Banzhaf (Rt Pearson):\n")
for(i in seq_along(SOURCE_NAMES)){sub<-gt_all[metric=="rt_corr"&source==SOURCE_NAMES[i]];cat(sprintf("    %-8s phi=%+.4f  normBeta=%+.4f  stable=%s\n",SOURCE_NAMES[i],sub$shapley,sub$norm_banzhaf,sub$banzhaf_stable))}
cat("\n  Shapley (Rt skill score):\n")
for(i in seq_along(SOURCE_NAMES)){sub<-gt_all[metric=="rt_skill"&source==SOURCE_NAMES[i]];cat(sprintf("    %-8s phi=%+.4f\n",SOURCE_NAMES[i],sub$shapley))}
cat(sprintf("\n  Ablation: %d/%d subsets converged\n",sum(is.finite(abl$rho)),nrow(abl)))
cat(sprintf("  Output prefix: %s\n",OUT_PREFIX))
if(USE_DELAY_KERNEL) cat(sprintf("  Delay kernels: case=%.1fd, hosp=%.1fd, icu=%.1fd\n",DELAY_CASE_MEAN,DELAY_HOSP_MEAN,DELAY_ICU_MEAN))
cat("\n  Output files:\n")
cat("    gp_seir_results.csv        main time series + residuals\n")
cat("    ablation_results.csv       all 15 subsets × all metrics\n")
cat("    shapley_banzhaf_values.csv Shapley + Banzhaf for 16 metrics\n")
cat("    shapley_values.csv         Shapley only (backward-compat)\n")
cat("    shapley_interactions_all.csv  pairwise interaction indices\n")
cat("    leave_one_out.csv          drop-one-source analysis\n")
cat("    superadditivity.csv        redundancy analysis\n")
cat("    rt_comprehensive.csv       Rt Pearson/Spearman/skill/RMSE/bias\n")
cat("    prev_comprehensive.csv     prevalence metrics\n")
cat("    rt_method_comparison.csv   GP vs Cori vs RIVM\n")
cat("    wave_params.csv / _refit   per-wave parameters\n")
cat("    wave_summary.csv           attack rates, peaks, etc.\n")
cat("    wdata_sensitivity.csv      ODE weight sensitivity\n")
cat("    beta_init_sensitivity.csv  initialization sensitivity\n")
cat("    profile_likelihood.csv     approximate profile likelihood\n")
cat("    sigma_gamma_sensitivity.csv sigma/gamma sensitivity\n")
cat("    residual_diagnostics.csv   ACF, Shapiro-Wilk, Ljung-Box\n")
cat("    stream_fit_metrics.csv     per-stream Pearson/RMSE/bias\n")
cat("    cross_correlations.csv     pairwise stream CCF\n")
cat("    rt_lag_analysis.csv        model–RIVM lag at ±10..20d\n")
cat("    derived_timeseries.csv     flows, growth rates, doubling time\n")
cat("    epi_parameters.csv         all fixed parameters\n")
cat("    tikz_*.csv                 26 TikZ-ready data files\n")
cat("    gp_seir_plots.pdf          9-page diagnostic PDF\n")
cat("============================================================\n")

save.image(out("gp_seir_workspace.RData"))
cat(sprintf("  Workspace saved: %s\n",out("gp_seir_workspace.RData")))
