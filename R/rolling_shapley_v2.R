#!/usr/bin/env Rscript
##############################################################################
# rolling_shapley_v2.R
#
# Rolling-window Shapley AND Banzhaf analysis for GP-SEIR multi-source model.
#
#
# USAGE
#   Rscript rolling_shapley_v2.R
#   Rscript rolling_shapley_v2.R --window=5 --step=1 --init=constant
#   Rscript rolling_shapley_v2.R --resume
#
# OUTPUT FILES (prefix = rolling_shapley_w{W}_s{S}_{init})
#   {prefix}_subset_results.csv        all metrics per window per subset
#   {prefix}_shapley_all.csv           Shapley for all metrics (long)
#   {prefix}_banzhaf_all.csv           Banzhaf (raw + normalised) for all metrics (long)
#   {prefix}_shapley_rt.csv / _prev    Shapley for Rt and prevalence
#   {prefix}_banzhaf_rt.csv / _prev    normBanzhaf for Rt and prevalence
#   {prefix}_tikz_rt.csv               Shapley wide format for figure
#   {prefix}_tikz_prev.csv             Shapley wide format for figure
#   {prefix}_tikz_banzhaf_rt.csv       normBanzhaf wide format for figure
#   {prefix}_tikz_banzhaf_prev.csv     normBanzhaf wide format for figure
#   {prefix}_windows.csv               window-level all-4 summary
#   {prefix}_window_quality.csv        completeness + rank agreement per window
#   {prefix}_interactions.csv          pairwise Shapley interaction indices
#   {prefix}_rankings.csv              source rankings (Shapley and Banzhaf)
#   {prefix}_redundancy.csv            redundancy over time
#   {prefix}_singletons.csv            single-source performance
#   {prefix}_leave_one_out.csv         drop-one-source analysis
#   {prefix}_best_subset.csv           best source subset per window
#   {prefix}_cori_comparison.csv       GP vs Cori vs RIVM Rt per window
#   {prefix}_all4_trajectories.csv     full trajectories for all-4 windows
#   {prefix}_selected_trajectories.csv singleton + all-4 trajectories
#
# Dependencies: data.table, fda, minpack.lm, zoo, EpiEstim (optional)
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(fda); library(minpack.lm); library(zoo)
})
HAVE_EPIESTIM <- requireNamespace("EpiEstim", quietly=TRUE)
if (!HAVE_EPIESTIM) cat("NOTE: EpiEstim not available. Rolling Cori Rt will be skipped.\n\n")

cat("================================================================\n")
cat("Rolling-Window Shapley + Banzhaf + Cori Analysis v2\n")
cat("================================================================\n\n")

# ══════════════════════════════════════════════════════════════
# SECTION 1: PARAMETERS
# ══════════════════════════════════════════════════════════════

DATA_FILE    <- "gp_input_final_v2.csv"
WINDOW_WEEKS <- 5L; STEP_WEEKS <- 1L; INIT_METHOD <- "constant"
GLOBAL_START <- as.Date("2020-03-01"); GLOBAL_END <- as.Date("2021-03-01")
N_WARM <- 4L; N_COLD <- 3L; N_LEVEL <- 3L; RESUME <- FALSE

args <- commandArgs(trailingOnly=TRUE)
for (a in args) {
  if (grepl("^--window=",a))  WINDOW_WEEKS <- as.integer(sub("^--window=","",a))
  if (grepl("^--step=",a))    STEP_WEEKS   <- as.integer(sub("^--step=","",a))
  if (grepl("^--init=",a))    INIT_METHOD  <- sub("^--init=","",a)
  if (grepl("^--start=",a))   GLOBAL_START <- as.Date(sub("^--start=","",a))
  if (grepl("^--end=",a))     GLOBAL_END   <- as.Date(sub("^--end=","",a))
  if (grepl("^--nwarm=",a))   N_WARM  <- as.integer(sub("^--nwarm=","",a))
  if (grepl("^--ncold=",a))   N_COLD  <- as.integer(sub("^--ncold=","",a))
  if (grepl("^--nlevel=",a))  N_LEVEL <- as.integer(sub("^--nlevel=","",a))
  if (a=="--resume")          RESUME <- TRUE
}
WINDOW_DAYS <- WINDOW_WEEKS*7L; STEP_DAYS <- STEP_WEEKS*7L
stopifnot(INIT_METHOD %in% c("hospital","constant","random_smooth"))
OUT_PREFIX <- sprintf("rolling_shapley_w%d_s%d_%s", WINDOW_WEEKS, STEP_WEEKS, INIT_METHOD)
CHECKPOINT_FILE <- normalizePath(
  file.path(getwd(), sprintf("%s_checkpoint.rds", OUT_PREFIX)),
  winslash = "/",
  mustWork = FALSE
)

cat(sprintf("  Checkpoint file: %s\n", CHECKPOINT_FILE))

# Fixed epidemiological parameters (shared with wave analysis)
N_POP  <- 17400000; SIGMA  <- 1/5.5; GAMMA  <- 1/9.5
STATE_KNOT_DAYS <- 14; STATE_NORDER <- 4; BETA_KNOT_DAYS <- 28; BETA_NORDER <- 4
ETA <- 1e6; KAPPA_BETA <- 1e3; KAPPA_BMAG <- 1e4; BETA_MAX <- 0.6; KAPPA_MASS <- 1e6
USE_DELAY_KERNEL <- TRUE
DELAY_CASE_MEAN <- 3.0; DELAY_CASE_SD <- 1.5
DELAY_HOSP_MEAN <- 2.0; DELAY_HOSP_SD <- 1.0
DELAY_ICU_MEAN  <- 3.5; DELAY_ICU_SD  <- 1.5; DELAY_MAX_LAG <- 14
BOUNDS <- list(rho=c(0.05,0.60), pH=c(0.002,0.05), pICU=c(0.0003,0.02), alphaR=c(0.10,1.00))

SOURCE_NAMES <- c("cases","hosp","icu","radar")
ALL_SUBSETS  <- unlist(lapply(1:4, function(k) combn(SOURCE_NAMES,k,simplify=FALSE)), recursive=FALSE)
SUBSET_KEYS  <- sapply(ALL_SUBSETS, function(s) paste(sort(s), collapse="+"))
N_SUBSETS    <- length(ALL_SUBSETS)

# Cori settings
CORI_SI_MEAN <- 5.1; CORI_SI_SD <- 2.8; CORI_TAU <- 7L

# Banzhaf stability threshold: flag window if |sum(beta_i)| < this
BANZHAF_STABLE_THRESHOLD <- 0.05

cat(sprintf("  Window=%d wks (%d days), step=%d wk, init=%s\n",
            WINDOW_WEEKS, WINDOW_DAYS, STEP_WEEKS, INIT_METHOD))
cat(sprintf("  Output prefix: %s\n\n", OUT_PREFIX))


# ══════════════════════════════════════════════════════════════
# SECTION 2: HELPERS (same as wave analysis v2)
# ══════════════════════════════════════════════════════════════

make_delay_kernel <- function(mean_d, sd_d, max_lag=DELAY_MAX_LAG) {
  if (mean_d<=0) return(1); lags <- 0:max_lag
  mu_ln <- log(mean_d^2/sqrt(sd_d^2+mean_d^2)); sig_ln <- sqrt(log(1+sd_d^2/mean_d^2))
  w <- dlnorm(lags,mu_ln,sig_ln); w/sum(w)
}
apply_delay <- function(y, kernel) {
  if (length(kernel)==1) return(y); n_y <- length(y); y_conv <- numeric(n_y)
  for (t in 1:n_y) for (k in seq_along(kernel)) { s <- t-(k-1); if(s>=1) y_conv[t]<-y_conv[t]+kernel[k]*y[s] }
  y_conv
}
to_u   <- function(x,lo,hi) qlogis((x-lo)/(hi-lo))
from_u <- function(u,lo,hi) lo+(hi-lo)*plogis(u)
wsd    <- function(y){ s<-sd(y,na.rm=TRUE); if(!is.finite(s)||s<=0) s<-1; 1/s^2 }
seir_rhs <- function(Z,beta,sig=SIGMA,gam=GAMMA)
  cbind(-beta*exp(Z[,3]), beta*exp(Z[,1]+Z[,3]-Z[,2])-sig,
        sig*exp(Z[,2]-Z[,3])-gam, gam*exp(Z[,3]-Z[,4]))
make_basis <- function(rng,knot_days,norder) {
  kn <- unique(c(rng[1],seq(rng[1],rng[2],by=knot_days),rng[2]))
  create.bspline.basis(rng,length(kn)+norder-2,norder,kn)
}
subset_to_weights <- function(sources)
  c(cases=as.numeric("cases"%in%sources), hosp=as.numeric("hosp"%in%sources),
    icu=as.numeric("icu"%in%sources),   radar=as.numeric("radar"%in%sources))
safe_cor  <- function(x,y,method="pearson"){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<5) return(NA_real_); cor(x[ok],y[ok],method=method) }
safe_rmse <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<3) return(NA_real_); sqrt(mean((x[ok]-y[ok])^2)) }


# ══════════════════════════════════════════════════════════════
# SECTION 3: GAME THEORY FUNCTIONS
# ══════════════════════════════════════════════════════════════

compute_shapley_lookup <- function(results_list, metric) {
  n_players <- 4L; phi <- numeric(n_players); names(phi) <- SOURCE_NAMES
  get_v <- function(src_set) {
    if (length(src_set)==0) return(0)
    key <- paste(sort(src_set),collapse="+"); r <- results_list[[key]]
    if (is.null(r)||!is.finite(r[[metric]])) return(NA_real_); r[[metric]]
  }
  for (i in seq_along(SOURCE_NAMES)) {
    player <- SOURCE_NAMES[i]; others <- SOURCE_NAMES[-i]; total <- 0
    for (s_size in 0:(n_players-1)) {
      coalitions <- if(s_size==0) list(character(0)) else combn(others,s_size,simplify=FALSE)
      weight <- factorial(s_size)*factorial(n_players-s_size-1)/factorial(n_players)
      for (S in coalitions) {
        vw <- get_v(c(S,player)); vnw <- get_v(S)
        if(is.finite(vw)&&is.finite(vnw)) total <- total+weight*(vw-vnw)
      }
    }
    phi[i] <- total
  }
  phi
}

compute_banzhaf_lookup <- function(results_list, metric) {
  n_players <- 4L; norm <- 1/(2^(n_players-1))
  beta <- numeric(n_players); names(beta) <- SOURCE_NAMES
  get_v <- function(src_set) {
    if (length(src_set)==0) return(0)
    key <- paste(sort(src_set),collapse="+"); r <- results_list[[key]]
    if (is.null(r)||!is.finite(r[[metric]])) return(NA_real_); r[[metric]]
  }
  for (i in seq_along(SOURCE_NAMES)) {
    player <- SOURCE_NAMES[i]; others <- SOURCE_NAMES[-i]; total <- 0
    for (s_size in 0:(n_players-1)) {
      coalitions <- if(s_size==0) list(character(0)) else combn(others,s_size,simplify=FALSE)
      for (S in coalitions) {
        vw <- get_v(c(S,player)); vnw <- get_v(S)
        if(is.finite(vw)&&is.finite(vnw)) total <- total+(vw-vnw)
      }
    }
    beta[i] <- norm*total
  }
  beta
}

normalise_banzhaf <- function(beta_vec, vN) {
  s <- sum(beta_vec)
  if (abs(s) < 1e-12) return(setNames(rep(0,length(beta_vec)),names(beta_vec)))
  beta_vec*vN/s
}

compute_pair_interactions <- function(results_list, metric) {
  n_players <- 4L; pairs <- combn(SOURCE_NAMES,2,simplify=FALSE); result <- list()
  get_v <- function(src_set) {
    if(length(src_set)==0) return(0)
    key <- paste(sort(src_set),collapse="+"); r <- results_list[[key]]
    if(is.null(r)||!is.finite(r[[metric]])) return(NA_real_); r[[metric]]
  }
  for (pair in pairs) {
    i <- pair[1]; j <- pair[2]; others <- setdiff(SOURCE_NAMES,pair); total <- 0
    for (s_size in 0:(n_players-2)) {
      coalitions <- if(s_size==0) list(character(0)) else combn(others,s_size,simplify=FALSE)
      weight <- factorial(s_size)*factorial(n_players-s_size-2)/factorial(n_players-1)
      for (S in coalitions) {
        v_ij<-get_v(c(S,i,j)); v_i<-get_v(c(S,i)); v_j<-get_v(c(S,j)); v_0<-get_v(S)
        if(all(is.finite(c(v_ij,v_i,v_j,v_0)))) total <- total+weight*(v_ij-v_i-v_j+v_0)
      }
    }
    result[[length(result)+1]] <- data.table(source_i=i,source_j=j,interaction=total)
  }
  rbindlist(result)
}

make_results_lookup <- function(win_rows) {
  out <- list()
  for (r in seq_len(nrow(win_rows))) { row <- win_rows[r]; out[[row$sources]] <- as.list(row) }
  out
}


# ══════════════════════════════════════════════════════════════
# SECTION 4: WINDOW-LEVEL GP-SEIR FIT
# ══════════════════════════════════════════════════════════════

fit_subset_window <- function(df, date_start, date_end, wc, wh, wi, wr,
                              init_method=INIT_METHOD, seed=42) {
  bnds <- BOUNDS
  dd <- copy(df)[date>=date_start & date<=date_end][order(date)]
  n <- nrow(dd); if(n<14) return(NULL)
  dd[, t:=as.numeric(date-date_start)]
  tv <- dd$t; rng <- range(tv); dt <- 1; dates <- dd$date

  cases_raw  <- as.numeric(dd$cases_raw); hosp_raw <- as.numeric(dd$hosp_nice)
  icu_raw    <- as.numeric(dd$icu_daily); radar_frac <- as.numeric(dd$radar_I_frac_sm7)
  rt_rivm    <- as.numeric(dd$rt_rivm);  prev_rivm  <- as.numeric(dd$prev_json_avg)

  y_cases <- log1p(pmax(0,cases_raw)); y_hosp <- log1p(pmax(0,hosp_raw))
  y_icu   <- log1p(pmax(0,icu_raw));   y_radar <- radar_frac
  wb_c <- wsd(y_cases); wb_h <- wsd(y_hosp); wb_i <- wsd(y_icu); wb_r <- wsd(y_radar)

  Bs <- make_basis(rng,STATE_KNOT_DAYS,STATE_NORDER)
  Bb <- make_basis(rng,BETA_KNOT_DAYS,BETA_NORDER)
  nbs <- Bs$nbasis; nbb <- Bb$nbasis
  Ph <- eval.basis(tv,Bs,0); dPh <- eval.basis(tv,Bs,1)
  Pb <- eval.basis(tv,Bb,0); d2Pb <- eval.basis(tv,Bb,2)

  h_case <- if(USE_DELAY_KERNEL) make_delay_kernel(DELAY_CASE_MEAN,DELAY_CASE_SD) else 1
  h_hosp <- if(USE_DELAY_KERNEL) make_delay_kernel(DELAY_HOSP_MEAN,DELAY_HOSP_SD) else 1
  h_icu  <- if(USE_DELAY_KERNEL) make_delay_kernel(DELAY_ICU_MEAN,DELAY_ICU_SD)  else 1

  c0 <- ifelse(is.finite(cases_raw),cases_raw,0)
  E0 <- pmax(1e-10,c0/(N_POP*0.20*SIGMA))
  I0 <- as.numeric(rollmean(E0,min(5,n),fill=median(E0),align="right"))
  R0 <- pmin(0.8,cumsum(pmax(c0,0))/N_POP)
  S0 <- pmax(1e-6,1-E0-I0-R0); m0<-S0+E0+I0+R0
  S0<-S0/m0; E0<-E0/m0; I0<-I0/m0; R0<-R0/m0
  Ci <- tryCatch(smooth.basis(tv,cbind(log(S0),log(E0),log(I0),log(R0)),
                              fdPar(Bs,int2Lfd(2),1e-4))$fd$coefs,
                 error=function(e) matrix(0,nbs,4))

  if (init_method=="hospital") {
    h0 <- ifelse(is.finite(hosp_raw),pmax(0.5,hosp_raw),0.5)
    hs <- as.numeric(rollmean(h0,min(14,n),fill=NA,align="center"))
    hs[is.na(hs)] <- median(hs,na.rm=TRUE)
    dl <- c(diff(log(pmax(1,hs))),0)
    ds <- as.numeric(rollmean(dl,min(7,n),fill=0,align="center")); ds[is.na(ds)]<-0
    b0 <- pmax(0.05,pmin(0.30,GAMMA*(1+ds*(1/SIGMA+1/GAMMA))))
    b0 <- as.numeric(rollmean(b0,min(7,n),fill=NA,align="center")); b0[is.na(b0)]<-median(b0,na.rm=TRUE)
  } else if (init_method=="constant") {
    b0 <- rep(GAMMA,n)
  } else {
    set.seed(seed+1000); b0 <- runif(n,0.05,0.25)
    b0 <- as.numeric(rollmean(b0,min(7,n),fill=NA,align="center")); b0[is.na(b0)]<-median(b0,na.rm=TRUE)
  }
  ai <- tryCatch(as.numeric(smooth.basis(tv,log(pmax(b0,1e-6)),fdPar(Bb,int2Lfd(2),1e-2))$fd$coefs),
                 error=function(e) rep(log(GAMMA),nbb))
  clamp <- function(x,lo,hi) max(lo+1e-8*(hi-lo),min(hi-1e-8*(hi-lo),x))
  thi <- c(ai,
           to_u(clamp(0.20,bnds$rho[1],bnds$rho[2]),bnds$rho[1],bnds$rho[2]),
           to_u(clamp(0.012,bnds$pH[1],bnds$pH[2]),bnds$pH[1],bnds$pH[2]),
           to_u(clamp(0.0025,bnds$pICU[1],bnds$pICU[2]),bnds$pICU[1],bnds$pICU[2]),
           to_u(clamp(0.40,bnds$alphaR[1],bnds$alphaR[2]),bnds$alphaR[1],bnds$alphaR[2]))

  solve_inner <- function(Cs, th) {
    av<-th[1:nbb]; rho<-from_u(th[nbb+1],bnds$rho[1],bnds$rho[2])
    pH<-from_u(th[nbb+2],bnds$pH[1],bnds$pH[2]); pICU<-from_u(th[nbb+3],bnds$pICU[1],bnds$pICU[2])
    aR<-from_u(th[nbb+4],bnds$alphaR[1],bnds$alphaR[2])
    g<-as.numeric(Pb%*%av); beta<-exp(g)
    fn <- function(cf) {
      C<-matrix(cf,nrow=nbs,ncol=4); Z<-Ph%*%C; dZ<-dPh%*%C
      Sv<-exp(Z[,1]); Ev<-exp(Z[,2]); Iv<-exp(Z[,3]); Rv<-exp(Z[,4])
      fc<-N_POP*rho*SIGMA*Ev; fh<-N_POP*pH*GAMMA*Iv; fi<-N_POP*pICU*GAMMA*Iv
      mc<-log1p(pmax(0,apply_delay(fc,h_case))); mh<-log1p(pmax(0,apply_delay(fh,h_hosp)))
      mi<-log1p(pmax(0,apply_delay(fi,h_icu))); mr<-aR*Iv
      rc<-sqrt(wb_c*wc)*ifelse(is.finite(y_cases),y_cases-mc,0)
      rh<-sqrt(wb_h*wh)*ifelse(is.finite(y_hosp),y_hosp-mh,0)
      ri<-sqrt(wb_i*wi)*ifelse(is.finite(y_icu),y_icu-mi,0)
      rr<-sqrt(wb_r*wr)*ifelse(is.finite(y_radar),y_radar-mr,0)
      rp<-sqrt(ETA*dt)*as.numeric(dZ-seir_rhs(Z,beta))
      rm<-sqrt(KAPPA_MASS)*(Sv+Ev+Iv+Rv-1)
      rb<-sqrt(KAPPA_BMAG*dt)*pmax(beta-BETA_MAX,0)
      rs<-sqrt(KAPPA_BETA*dt)*as.numeric(d2Pb%*%av)
      c(rc,rh,ri,rr,rp,rm,rb,rs)
    }
    r <- nls.lm(par=as.numeric(Cs),fn=fn,
                control=nls.lm.control(maxiter=200,ftol=1e-8,ptol=1e-8))
    list(C=matrix(r$par,nrow=nbs,ncol=4),J=sum(r$fvec^2))
  }
  stopifnot(!any(c("rt_rivm","prev_rivm") %in% all.vars(body(solve_inner))))

  cascade <- function(th0,C0) {
    th<-th0; Cc<-C0; bJ<-Inf
    for(oi in 1:30){
      inn<-solve_inner(Cc,th); Cc<-inn$C; Jc<-inn$J; if(Jc<bJ) bJ<-Jc
      eps<-1e-4; gr<-numeric(length(th))
      for(j in seq_along(th)){ tp<-th; tp[j]<-tp[j]+eps; gr[j]<-(solve_inner(Cc,tp)$J-Jc)/eps }
      gn<-sqrt(sum(gr^2)); if(gn<1e-12) break
      d<--gr/gn; st<-1.0; found<-FALSE
      for(ls in 1:8){ tt<-th+st*d; it<-solve_inner(Cc,tt); if(it$J<Jc-1e-4*st*gn){th<-tt;Cc<-it$C;found<-TRUE;break}; st<-st*0.5 }
      if(!found) th<-th+1e-4*d
      if(gn<1e-5*max(1,sqrt(sum(th^2)))) break
    }
    list(theta=th,C=Cc,J=bJ)
  }

  set.seed(seed); best<-NULL; bJ<-Inf; all_J<-numeric(0)
  starts <- list(list(th=thi,C=Ci,type="warm"))
  for(s in 2:N_WARM){ ts<-thi; ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,0.30)
    for(k in (nbb+1):(nbb+4)) ts[k]<-ts[k]+rnorm(1,0,0.50)
    starts[[length(starts)+1]] <- list(th=ts,C=Ci,type="warm") }
  for(s in 1:N_COLD){ ts<-thi; ts[1:nbb]<-ts[1:nbb]+rnorm(nbb,0,0.50)
    for(k in (nbb+1):(nbb+4)) ts[k]<-ts[k]+rnorm(1,0,0.80)
    starts[[length(starts)+1]] <- list(th=ts,C=Ci,type="cold") }
  for(s in 1:N_LEVEL){ ts<-thi; shift<-rnorm(1,0,0.3); ts[1:nbb]<-ts[1:nbb]+shift
    starts[[length(starts)+1]] <- list(th=ts,C=Ci,type="level") }

  for(st in starts) {
    r <- tryCatch(cascade(st$th,st$C), error=function(e) NULL)
    if(!is.null(r)){ all_J<-c(all_J,r$J); if(r$J<bJ){bJ<-r$J;best<-r} }
  }
  if(is.null(best)) return(NULL)

  th<-best$theta; C<-best$C; av<-th[1:nbb]
  rho<-from_u(th[nbb+1],bnds$rho[1],bnds$rho[2])
  pH<-from_u(th[nbb+2],bnds$pH[1],bnds$pH[2])
  pICU<-from_u(th[nbb+3],bnds$pICU[1],bnds$pICU[2])
  aR<-from_u(th[nbb+4],bnds$alphaR[1],bnds$alphaR[2])
  Z<-Ph%*%C; Sv<-exp(Z[,1]); Ev<-exp(Z[,2]); Iv<-exp(Z[,3]); Rv<-exp(Z[,4])
  beta<-exp(as.numeric(Pb%*%av)); Rt<-beta*Sv/GAMMA; I_count<-N_POP*Iv
  ode_mat<-dPh%*%C-seir_rhs(Z,beta); ode_rmse<-sqrt(mean(ode_mat^2))
  fc<-N_POP*rho*SIGMA*Ev; fh<-N_POP*pH*GAMMA*Iv; fi<-N_POP*pICU*GAMMA*Iv
  mc<-apply_delay(fc,h_case); mh<-apply_delay(fh,h_hosp); mi<-apply_delay(fi,h_icu)
  ok_rt<-is.finite(rt_rivm)&is.finite(Rt)

  traj_dt <- data.table(
    day=tv, date=dates, S=Sv, E=Ev, I=Iv, R=Rv, beta=beta, Rt=Rt,
    mu_cases=mc, mu_hosp=mh, mu_icu=mi, mu_radar=aR*Iv,
    rt_rivm=rt_rivm, prev_rivm=prev_rivm
  )

  list(
    J=bJ, all_J=all_J,
    rt_pearson=safe_cor(rt_rivm,Rt), rt_spearman=safe_cor(rt_rivm,Rt,"spearman"),
    rt_rmse=safe_rmse(rt_rivm,Rt), neg_rt_rmse=-safe_rmse(rt_rivm,Rt),
    rt_bias=if(sum(ok_rt)>3) mean(Rt[ok_rt]-rt_rivm[ok_rt]) else NA_real_,
    rt_threshold=if(sum(ok_rt)>3) mean((Rt[ok_rt]>1)==(rt_rivm[ok_rt]>1)) else NA_real_,
    prev_pearson=safe_cor(prev_rivm,I_count), prev_spearman=safe_cor(prev_rivm,I_count,"spearman"),
    prev_rmse=safe_rmse(prev_rivm,I_count),
    cases_rmse=safe_rmse(cases_raw,mc), hosp_rmse=safe_rmse(hosp_raw,mh),
    ode_rmse=ode_rmse, mass_err=max(abs(Sv+Ev+Iv+Rv-1)),
    rho=rho, pH=pH, pICU=pICU, alphaR=aR,
    beta_min=min(beta), beta_max=max(beta), beta_mean=mean(beta),
    rt_min=min(Rt), rt_max=max(Rt), rt_mean=mean(Rt),
    R_end=Rv[n], attack_rate=Rv[n],
    peak_I=max(I_count), peak_I_date=as.character(dates[which.max(Iv)]),
    peak_Rt=max(Rt), peak_Rt_date=as.character(dates[which.max(Rt)]),
    pct_Rt_above1=100*mean(Rt>1),
    n_Rt_crossings=sum(diff(Rt>1)!=0),
    first_Rt_crossing=tryCatch(as.character(dates[which(diff(Rt>1)!=0)[1]+1]), error=function(e) NA),
    mean_rho_eff=mean(ifelse(is.finite(cases_raw)&N_POP*SIGMA*Ev>1,cases_raw/(N_POP*SIGMA*Ev),NA),na.rm=TRUE),
    neg_J=-bJ,
    # Skill score: 1 - MSE/Var(Rt_RIVM) -- natural zero, avoids Pearson sign problems
    rt_skill = {
      ok <- is.finite(rt_rivm) & is.finite(Rt)
      if (sum(ok)>5) { mse<-mean((Rt[ok]-rt_rivm[ok])^2); vref<-var(rt_rivm[ok])
        if (is.finite(vref) && vref>1e-12) 1-mse/vref else NA_real_ } else NA_real_
    },
    prev_skill = {
      I_v <- N_POP*Iv; ok <- is.finite(prev_rivm) & is.finite(I_v)
      if (sum(ok)>5) { mse<-mean((I_v[ok]-prev_rivm[ok])^2); vref<-var(prev_rivm[ok])
        if (is.finite(vref) && vref>1e-12) 1-mse/vref else NA_real_ } else NA_real_
    },
    trajectories=list(traj_dt),
    cases_raw=cases_raw, rt_rivm_vec=rt_rivm
  )
}


# ══════════════════════════════════════════════════════════════
# SECTION 5: LOAD DATA
# ══════════════════════════════════════════════════════════════

cat("Loading data...\n")
if(!file.exists(DATA_FILE)) stop("Data file not found: ",DATA_FILE)
FULL_DF <- fread(DATA_FILE)[,date:=as.Date(date)][order(date)]
cat(sprintf("  %d rows (%s to %s)\n\n",nrow(FULL_DF),min(FULL_DF$date),max(FULL_DF$date)))

window_starts <- seq(GLOBAL_START, GLOBAL_END-WINDOW_DAYS+1L, by=STEP_DAYS)
n_windows <- length(window_starts)
cat(sprintf("  %d windows | W=%d wks, S=%d wk, init=%s\n", n_windows, WINDOW_WEEKS, STEP_WEEKS, INIT_METHOD))
cat(sprintf("  Total subset-fits: %d\n\n", n_windows*N_SUBSETS))


# ══════════════════════════════════════════════════════════════
# SECTION 6: RESUME FROM CHECKPOINT
# ══════════════════════════════════════════════════════════════

all_subset_rows <- list(); completed_windows <- integer(0)
all4_traj_rows <- list(); selected_traj_rows <- list()

if (RESUME) {
  if (!file.exists(CHECKPOINT_FILE)) {
    stop(sprintf(
      paste0(
        "Checkpoint not found:\n  %s\n\n",
        "Resume must use the same --window, --step, and --init values ",
        "as the run that created the checkpoint."
      ),
      CHECKPOINT_FILE
    ))
  }
  chk <- readRDS(CHECKPOINT_FILE)
  all_subset_rows <- chk$all_subset_rows
  completed_windows <- chk$completed_windows
  if(!is.null(chk$all4_traj_rows)) all4_traj_rows <- chk$all4_traj_rows
  if(!is.null(chk$selected_traj_rows)) selected_traj_rows <- chk$selected_traj_rows
  cat(sprintf("Resumed: %d/%d windows completed\n\n", length(completed_windows), n_windows))
}


# ══════════════════════════════════════════════════════════════
# SECTION 7: MAIN LOOP
# ══════════════════════════════════════════════════════════════

t_start <- Sys.time()
for(wi in seq_along(window_starts)) {
  if(wi %in% completed_windows) next
  w_start <- window_starts[wi]; w_end <- w_start+WINDOW_DAYS-1L
  w_mid <- w_start+(WINDOW_DAYS-1L)%/%2L
  cat(sprintf("\n[%d/%d] %s to %s (mid=%s)\n", wi, n_windows, w_start, w_end, w_mid))
  subset_cache <- list()

  for(si in seq_len(N_SUBSETS)) {
    src_key <- SUBSET_KEYS[si]; w <- subset_to_weights(ALL_SUBSETS[[si]])
    if (src_key %in% names(subset_cache)) {
      fit <- subset_cache[[src_key]]$fit
    } else {
      fit <- tryCatch(
        fit_subset_window(FULL_DF,w_start,w_end,
                          wc=w["cases"],wh=w["hosp"],wi=w["icu"],wr=w["radar"],
                          init_method=INIT_METHOD, seed=wi*1000L+si),
        error=function(e) NULL)
      subset_cache[[src_key]] <- list(fit=fit)
    }
    if(!is.null(fit)) {
      all_subset_rows[[length(all_subset_rows)+1L]] <- data.table(
        window_idx=wi, window_start=w_start, window_end=w_end, mid_date=w_mid,
        sources=src_key, n_sources=length(ALL_SUBSETS[[si]]),
        J=fit$J, rt_pearson=fit$rt_pearson, rt_spearman=fit$rt_spearman,
        rt_rmse=fit$rt_rmse, rt_bias=fit$rt_bias, rt_threshold=fit$rt_threshold,
        prev_pearson=fit$prev_pearson, prev_spearman=fit$prev_spearman, prev_rmse=fit$prev_rmse,
        cases_rmse=fit$cases_rmse, hosp_rmse=fit$hosp_rmse,
        ode_rmse=fit$ode_rmse, mass_err=fit$mass_err,
        rho=fit$rho, pH=fit$pH, pICU=fit$pICU, alphaR=fit$alphaR,
        beta_min=fit$beta_min, beta_max=fit$beta_max, beta_mean=fit$beta_mean,
        rt_min=fit$rt_min, rt_max=fit$rt_max, rt_mean=fit$rt_mean,
        R_end=fit$R_end, attack_rate=fit$attack_rate,
        peak_I=fit$peak_I, peak_I_date=fit$peak_I_date,
        peak_Rt=fit$peak_Rt, peak_Rt_date=fit$peak_Rt_date,
        pct_Rt_above1=fit$pct_Rt_above1, n_Rt_crossings=fit$n_Rt_crossings,
        first_Rt_crossing=fit$first_Rt_crossing,
        mean_rho_eff=fit$mean_rho_eff, neg_rt_rmse=fit$neg_rt_rmse, neg_J=fit$neg_J,
        rt_skill=fit$rt_skill, prev_skill=fit$prev_skill)

      if(identical(src_key,paste(sort(SOURCE_NAMES),collapse="+"))) {
        tr <- copy(fit$trajectories[[1]])
        tr[,`:=`(window_idx=wi,window_start=w_start,window_end=w_end,mid_date=w_mid,sources=src_key)]
        all4_traj_rows[[length(all4_traj_rows)+1L]] <- tr
      }
      if(length(ALL_SUBSETS[[si]])==1L || length(ALL_SUBSETS[[si]])==4L) {
        tr <- copy(fit$trajectories[[1]])
        tr[,`:=`(window_idx=wi,window_start=w_start,window_end=w_end,mid_date=w_mid,sources=src_key)]
        selected_traj_rows[[length(selected_traj_rows)+1L]] <- tr
      }
    }
  }
  completed_windows <- c(completed_windows,wi)
  saveRDS(list(all_subset_rows=all_subset_rows,completed_windows=completed_windows,
               all4_traj_rows=all4_traj_rows,selected_traj_rows=selected_traj_rows),
          CHECKPOINT_FILE)
  elapsed <- as.numeric(difftime(Sys.time(),t_start,units="mins"))
  rate <- elapsed/length(completed_windows); n_left <- n_windows-length(completed_windows)
  cat(sprintf("  ETA: %.0f min\n", rate*n_left))
}


# ══════════════════════════════════════════════════════════════
# SECTION 8: COMPUTE SHAPLEY AND BANZHAF VALUES
# ══════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("Computing Shapley and Banzhaf values...\n")
cat("================================================================\n")

if(length(all_subset_rows)==0){cat("ERROR: No results.\n"); quit(status=1)}
subset_dt <- rbindlist(all_subset_rows)
fwrite(subset_dt, sprintf("%s_subset_results.csv",OUT_PREFIX))
cat(sprintf("  Saved: %s_subset_results.csv (%d rows)\n",OUT_PREFIX,nrow(subset_dt)))

PRIMARY_METRICS <- c("rt_pearson","prev_pearson")
ALL_METRICS <- c("rt_pearson","rt_spearman","rt_threshold","prev_pearson","prev_spearman",
                  "neg_rt_rmse","neg_J","rt_skill","prev_skill")

all_shapley <- list(); all_banzhaf <- list()
window_quality_rows <- list(); interaction_rows <- list()
ranking_rows <- list(); best_subset_rows <- list()
leave_one_out_rows <- list(); redundancy_rows <- list(); singleton_rows <- list()
cori_comp_rows <- list()

for(widx in sort(unique(subset_dt$window_idx))) {
  win_rows <- subset_dt[window_idx==widx]
  w_mid <- win_rows$mid_date[1]; w_start <- win_rows$window_start[1]; w_end <- win_rows$window_end[1]
  present_keys <- sort(unique(win_rows$sources))
  complete_subsets <- identical(present_keys, sort(SUBSET_KEYS))
  results_list <- make_results_lookup(win_rows)
  all4_key <- paste(sort(SOURCE_NAMES),collapse="+")
  all4_row  <- win_rows[sources==all4_key]
  vN_rt   <- if(nrow(all4_row)>0) all4_row$rt_pearson[1]   else NA_real_
  vN_prev <- if(nrow(all4_row)>0) all4_row$prev_pearson[1] else NA_real_

  for(metric in ALL_METRICS) {
    metric_complete <- complete_subsets && all(is.finite(win_rows[[metric]]))
    window_quality_rows[[length(window_quality_rows)+1L]] <- data.table(
      window_idx=widx, window_start=w_start, window_end=w_end, mid_date=w_mid,
      metric=metric, complete=metric_complete)

    if (!metric_complete) {
      for (src in SOURCE_NAMES) {
        all_shapley[[length(all_shapley)+1L]] <- data.table(
          window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,source=src,shapley=NA_real_)
        all_banzhaf[[length(all_banzhaf)+1L]] <- data.table(
          window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,source=src,banzhaf=NA_real_,norm_banzhaf=NA_real_,
          banzhaf_stable=FALSE)
      }
      next
    }

    phi  <- compute_shapley_lookup(results_list, metric)
    beta <- compute_banzhaf_lookup(results_list, metric)
    vN_m <- if(metric=="rt_pearson") vN_rt else if(metric=="prev_pearson") vN_prev else {
      r<-all4_row[[metric]]; if(length(r)>0&&is.finite(r[1])) r[1] else NA_real_
    }
    b_sum <- sum(beta)
    # A window is unstable if: (1) |sum| < threshold, OR (2) sign(sum) != sign(v(N))
    # which causes the normalisation to invert all values.
    sign_mismatch <- !is.na(vN_m) && abs(b_sum) > 1e-10 && abs(vN_m) > 0.01 &&
                     (b_sum * vN_m < 0)
    b_stable <- abs(b_sum) >= BANZHAF_STABLE_THRESHOLD && !sign_mismatch
    nb <- if(!is.na(vN_m) && b_stable) {
      normalise_banzhaf(beta, vN_m)
    } else {
      setNames(rep(NA_real_, length(SOURCE_NAMES)), SOURCE_NAMES)
    }

    for(src in SOURCE_NAMES) {
      all_shapley[[length(all_shapley)+1L]] <- data.table(
        window_start=w_start,window_end=w_end,mid_date=w_mid,
        metric=metric,source=src,shapley=phi[src])
      all_banzhaf[[length(all_banzhaf)+1L]] <- data.table(
        window_start=w_start,window_end=w_end,mid_date=w_mid,
        metric=metric,source=src,banzhaf=beta[src],norm_banzhaf=nb[src],
        banzhaf_stable=b_stable)
    }
  }

  if(complete_subsets) {
    # Rankings under both indices
    for(metric in PRIMARY_METRICS) {
      metric_complete <- all(is.finite(win_rows[[metric]]))
      if(!metric_complete) next
      phi  <- compute_shapley_lookup(results_list,metric)
      beta <- compute_banzhaf_lookup(results_list,metric)
      vN_m <- if(metric=="rt_pearson") vN_rt else vN_prev
      b_sum <- sum(beta)
      sign_mismatch <- !is.na(vN_m) && abs(b_sum) > 1e-10 && abs(vN_m) > 0.01 &&
                       (b_sum * vN_m < 0)
      b_stable <- abs(b_sum) >= BANZHAF_STABLE_THRESHOLD && !sign_mismatch
      nb <- if(!is.na(vN_m) && b_stable) {
        normalise_banzhaf(beta, vN_m)
      } else {
        setNames(rep(NA_real_, length(SOURCE_NAMES)), SOURCE_NAMES)
      }
      ord_phi <- order(-phi)
      ord_nb <- if (all(is.na(nb))) seq_along(SOURCE_NAMES) else order(-nb, na.last=TRUE)
      for(rank_idx in seq_along(ord_phi)) {
        src_phi <- names(phi)[ord_phi][rank_idx]
        src_nb  <- if (rank_idx <= sum(is.finite(nb))) names(nb)[ord_nb][rank_idx] else NA_character_
        ranking_rows[[length(ranking_rows)+1L]] <- data.table(
          window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,rank=rank_idx,
          source_shapley=src_phi, shapley=phi[src_phi],
          source_banzhaf=src_nb,
          norm_banzhaf=if(!is.na(src_nb)) nb[src_nb] else NA_real_,
          banzhaf_stable=b_stable,
          all4_value=vN_m,
          signal_ok=abs(vN_m) >= 0.15,
          rank_agrees=if(!is.na(src_nb)) (src_phi==src_nb) else NA)
      }
    }

    # Pairwise interactions (Shapley)
    for(metric in PRIMARY_METRICS) {
      if(!all(is.finite(win_rows[[metric]]))) next
      ints <- compute_pair_interactions(results_list,metric)
      ints[,`:=`(window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,metric=metric)]
      interaction_rows[[length(interaction_rows)+1L]] <- ints
    }

    # Best subset, LOO, redundancy, singletons
    for(metric in PRIMARY_METRICS) {
      if(!all(is.finite(win_rows[[metric]]))) next
      all4_val <- if(nrow(all4_row)>0) all4_row[[metric]][1] else NA_real_
      tmp <- win_rows[is.finite(get(metric))]
      if(nrow(tmp)>0) {
        bi <- which.max(tmp[[metric]])[1]
        best_subset_rows[[length(best_subset_rows)+1L]] <- data.table(
          window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,best_sources=tmp$sources[bi],best_n=tmp$n_sources[bi],
          best_value=tmp[[metric]][bi],all4_value=all4_val)
      }
      for(src in SOURCE_NAMES) {
        rem_key <- paste(sort(setdiff(SOURCE_NAMES,src)),collapse="+")
        rem_row <- win_rows[sources==rem_key]; sing_row <- win_rows[sources==src]
        leave_one_out_rows[[length(leave_one_out_rows)+1L]] <- data.table(
          window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,dropped=src,all4=all4_val,
          without=if(nrow(rem_row)>0) rem_row[[metric]][1] else NA_real_,
          drop_value=if(nrow(rem_row)>0) all4_val-rem_row[[metric]][1] else NA_real_)
        singleton_rows[[length(singleton_rows)+1L]] <- data.table(
          window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,source=src,
          singleton_value=if(nrow(sing_row)>0) sing_row[[metric]][1] else NA_real_,
          all4_value=all4_val)
      }
      sing_vals <- win_rows[n_sources==1][[metric]]
      if(length(sing_vals)==4 && all(is.finite(sing_vals)) && is.finite(all4_val) && sum(sing_vals)!=0)
        redundancy_rows[[length(redundancy_rows)+1L]] <- data.table(
          window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
          metric=metric,all4_value=all4_val,singleton_sum=sum(sing_vals),
          best_singleton=max(sing_vals),
          ratio_all4_to_sum=all4_val/sum(sing_vals),
          redundancy=1-all4_val/sum(sing_vals))
    }

    # Cori comparison per window
    if(HAVE_EPIESTIM && nrow(all4_row)>0) {
      dd_w <- FULL_DF[date>=w_start & date<=w_end][order(date)]
      inc <- pmax(0,round(ifelse(is.finite(dd_w$cases_raw),dd_w$cases_raw,0)))
      if(sum(inc)>10 && length(inc)>CORI_TAU+2) {
        cfg <- tryCatch(EpiEstim::make_config(method="parametric_si",
          mean_si=CORI_SI_MEAN, std_si=CORI_SI_SD,
          t_start=2:(length(inc)-CORI_TAU+1), t_end=(CORI_TAU+1):length(inc)),
          error=function(e) NULL)
        if(!is.null(cfg)) {
          cori_res <- tryCatch(EpiEstim::estimate_R(inc,method="parametric_si",config=cfg),
                               error=function(e) NULL)
          if(!is.null(cori_res)) {
            r_cori <- as.data.table(cori_res$R)
            r_cori[, date:=dd_w$date[t_end]]
            setnames(r_cori, c("Mean(R)","Quantile.0.025(R)","Quantile.0.975(R)"),
                     c("Rt_cori","Rt_cori_lo","Rt_cori_hi"))
            # Get all-4 trajectory for this window from saved trajectories
            traj_dt_w <- rbindlist(Filter(function(x) x$window_idx[1]==widx && x$sources[1]==all4_key,
                                          all4_traj_rows))
            if(nrow(traj_dt_w)>0) {
              merged_w <- merge(traj_dt_w[,.(date,Rt_gp=Rt,Rt_rivm=rt_rivm)],
                                r_cori[,.(date,Rt_cori)], by="date",all.x=TRUE)
              ok_gc <- is.finite(merged_w$Rt_gp) & is.finite(merged_w$Rt_cori)
              ok_cr <- is.finite(merged_w$Rt_cori) & is.finite(merged_w$Rt_rivm)
              cori_comp_rows[[length(cori_comp_rows)+1L]] <- data.table(
                window_idx=widx,window_start=w_start,window_end=w_end,mid_date=w_mid,
                gp_vs_cori_pearson=safe_cor(merged_w$Rt_gp,merged_w$Rt_cori),
                gp_vs_cori_rmse=safe_rmse(merged_w$Rt_gp,merged_w$Rt_cori),
                gp_vs_cori_bias=if(sum(ok_gc)>3) mean(merged_w$Rt_gp[ok_gc]-merged_w$Rt_cori[ok_gc]) else NA,
                cori_vs_rivm_pearson=safe_cor(merged_w$Rt_cori,merged_w$Rt_rivm),
                cori_vs_rivm_rmse=safe_rmse(merged_w$Rt_cori,merged_w$Rt_rivm),
                cori_vs_rivm_bias=if(sum(ok_cr)>3) mean(merged_w$Rt_cori[ok_cr]-merged_w$Rt_rivm[ok_cr]) else NA,
                gp_rt_pearson=all4_row$rt_pearson[1])
            }
          }
        }
      }
    }
  }
}


# ══════════════════════════════════════════════════════════════
# SECTION 9: SAVE ALL OUTPUTS
# ══════════════════════════════════════════════════════════════

cat("\nSaving outputs...\n")

shapley_dt <- rbindlist(all_shapley)
banzhaf_dt <- rbindlist(all_banzhaf)
fwrite(shapley_dt, sprintf("%s_shapley_all.csv",OUT_PREFIX))
fwrite(banzhaf_dt, sprintf("%s_banzhaf_all.csv",OUT_PREFIX))

# Shapley and Banzhaf TikZ wide tables
for(met in c("rt_pearson","prev_pearson")) {
  short <- sub("_pearson","",met)
  s_dt  <- shapley_dt[metric==met]
  b_dt  <- banzhaf_dt[metric==met & banzhaf_stable==TRUE]
  if(nrow(s_dt)>0) {
    fwrite(s_dt, sprintf("%s_shapley_%s.csv",OUT_PREFIX,short))
    wide <- dcast(s_dt, mid_date~source, value.var="shapley")
    fwrite(wide, sprintf("%s_tikz_%s.csv",OUT_PREFIX,short))
  }
  if(nrow(b_dt)>0) {
    fwrite(b_dt, sprintf("%s_banzhaf_%s.csv",OUT_PREFIX,short))
    # Wide format for TikZ figure: unstable windows show as nan (pgfplots skips them = gaps)
    wide_b <- dcast(
      banzhaf_dt[metric==met],  # use ALL rows, not just stable, to keep all dates
      mid_date~source, value.var="norm_banzhaf")
    # NA -> nan for pgfplots compatibility
    for(col in SOURCE_NAMES) set(wide_b, which(is.na(wide_b[[col]])), col, NaN)
    fwrite(wide_b, sprintf("%s_tikz_banzhaf_%s.csv",OUT_PREFIX,short))
  }
  cat(sprintf("  Shapley + Banzhaf TikZ for %s\n", short))
}

# All other outputs
save_if_nonempty <- function(lst, fname) {
  if(length(lst)>0){ dt <- rbindlist(lst); fwrite(dt,fname); cat(sprintf("  %s (%d rows)\n",fname,nrow(dt))) }
}
save_if_nonempty(window_quality_rows, sprintf("%s_window_quality.csv",OUT_PREFIX))
save_if_nonempty(interaction_rows,    sprintf("%s_interactions.csv",OUT_PREFIX))
save_if_nonempty(ranking_rows,        sprintf("%s_rankings.csv",OUT_PREFIX))
save_if_nonempty(best_subset_rows,    sprintf("%s_best_subset.csv",OUT_PREFIX))
save_if_nonempty(leave_one_out_rows,  sprintf("%s_leave_one_out.csv",OUT_PREFIX))
save_if_nonempty(redundancy_rows,     sprintf("%s_redundancy.csv",OUT_PREFIX))
save_if_nonempty(singleton_rows,      sprintf("%s_singletons.csv",OUT_PREFIX))
save_if_nonempty(cori_comp_rows,      sprintf("%s_cori_comparison.csv",OUT_PREFIX))
save_if_nonempty(all4_traj_rows,      sprintf("%s_all4_trajectories.csv",OUT_PREFIX))
save_if_nonempty(selected_traj_rows,  sprintf("%s_selected_trajectories.csv",OUT_PREFIX))

# Window quality with rank-1 agreement statistic
if(length(ranking_rows)>0) {
  rank_dt <- rbindlist(ranking_rows)
  tryCatch({
    agree_summary <- rank_dt[rank==1 & is.finite(rank_agrees),
      .(n_agree=sum(rank_agrees),
        n_total=.N,
        pct_agree=100*mean(rank_agrees)),
      by=.(metric)][, comparison := "all_rank1_windows"]

    agree_summary_stable <- rank_dt[rank==1 & is.finite(rank_agrees) & banzhaf_stable,
      .(n_agree=sum(rank_agrees),
        n_total=.N,
        pct_agree=100*mean(rank_agrees)),
      by=.(metric)][, comparison := "banzhaf_stable_rank1_windows"]

    agree_summary_signal <- rank_dt[rank==1 & is.finite(rank_agrees) & banzhaf_stable & signal_ok,
      .(n_agree=sum(rank_agrees),
        n_total=.N,
        pct_agree=100*mean(rank_agrees)),
      by=.(metric)][, comparison := "banzhaf_stable_signal_windows"]

    agree_summary <- rbindlist(
      list(agree_summary, agree_summary_stable, agree_summary_signal),
      use.names=TRUE, fill=TRUE
    )
    fwrite(agree_summary, sprintf("%s_rank1_agreement.csv",OUT_PREFIX))

    rt_row <- agree_summary[metric=="rt_pearson" & comparison=="banzhaf_stable_signal_windows"]
    prev_row <- agree_summary[metric=="prev_pearson" & comparison=="banzhaf_stable_signal_windows"]

    if (nrow(rt_row) > 0)
      cat(sprintf("  Rank-1 agreement (Rt, stable signal windows): %.0f%%\n", rt_row$pct_agree[1]))
    if (nrow(prev_row) > 0)
      cat(sprintf("  Rank-1 agreement (Prev, stable signal windows): %.0f%%\n", prev_row$pct_agree[1]))
  }, error = function(e) {
    cat(sprintf("WARNING: agree_summary failed: %s\n", conditionMessage(e)))
  })
}

# All-4 window summary
win_summary <- subset_dt[sources==paste(sort(SOURCE_NAMES),collapse="+"),
  .(window_start,window_end,mid_date,J,rt_pearson,rt_spearman,rt_rmse,rt_bias,
    prev_pearson,prev_spearman,prev_rmse,ode_rmse,mass_err,
    rho,pH,pICU,alphaR,beta_mean,rt_mean,R_end,attack_rate,
    peak_I,peak_I_date,peak_Rt,peak_Rt_date,pct_Rt_above1,
    n_Rt_crossings,mean_rho_eff)]
fwrite(win_summary, sprintf("%s_windows.csv",OUT_PREFIX))
cat(sprintf("  %s_windows.csv (%d windows)\n",OUT_PREFIX,nrow(win_summary)))

if(file.exists(CHECKPOINT_FILE)) file.remove(CHECKPOINT_FILE)

total_time <- as.numeric(difftime(Sys.time(),t_start,units="mins"))
cat(sprintf("\n================================================================\n"))
cat(sprintf("COMPLETE in %.1f minutes (%d windows, %d subset-fits)\n",
            total_time,length(unique(subset_dt$window_idx)),nrow(subset_dt)))
cat(sprintf("================================================================\n"))
