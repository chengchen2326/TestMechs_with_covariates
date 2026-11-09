# Is cross-fitting worth it, and is the optimistic standard error trustworthy?
#
# Three estimators on the same draws:
#   split    one split, stage 1 on 30%
#   cf_cons  cross-fitted, Cauchy-Schwarz standard error (always valid)
#   cf_pool  cross-fitted, zero-covariance standard error (the usual practice)
#
# The question the Monte Carlo answers is whether cf_pool keeps coverage. If it
# does, cross-fitting buys the full sqrt(2), which is most of what we lose to
# KR's test. If it does not, cf_cons is still a free improvement.
suppressMessages({library(purrr); library(Rglpk)})
SRC <- Sys.getenv("TM_SRC","../src")
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_simulate_data_binaryM.R","src_simulate_data_binaryM_x.R",
            "src_lb_frac_affected_x.R","src_lb_frac_affected_x_cf.R")) source(file.path(SRC,f))
source("violation_scale.R")

calib <- list(
  list(p_y_m0d0=c(0.65,0.35), p_y_m0d1=c(0.50,0.50),
       p_y_m1d0=c(0.50,0.50), p_y_m1d1=c(0.50,0.50), p_m_0=0.6, p_m_1=0.6),
  list(p_y_m0d0=c(0.50,0.50), p_y_m0d1=c(0.65,0.35),
       p_y_m1d0=c(0.50,0.50), p_y_m1d1=c(0.50,0.50), p_m_0=0.6, p_m_1=0.6))

REPS <- as.integer(Sys.getenv("TM_REPS","600"))
cat(sprintf("%4s %6s %-9s %9s %9s %9s %8s %8s\n",
            "t","n","estimator","bound","mean_eta","mean_se","power","coverage"))
for (tt in c(0, 0.5, 1)) {
  cells <- scale_violation_binary(calib, tt)
  pb <- population_bound_binaryM_x(cells, c(0.5,0.5), at_group=0)
  bd <- pb$eta_conditional
  for (n in c(1000, 4000)) {
    set.seed(20261109)
    M <- matrix(NA_real_, REPS, 6)
    for (r in seq_len(REPS)) {
      df <- simulate_data_binaryM_x(cells, p_x=c(0.5,0.5), n=n, p_d=0.5)
      s  <- try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=0,
                                   method="split", split_frac=0.3), silent=TRUE)
      cf <- try(lb_frac_affected_x_cf(df,"d","m","y",x="x",at_group=0), silent=TRUE)
      if (!inherits(s,"try-error")) M[r,1:2] <- c(s$eta_lcb, s$eta_se)
      if (!inherits(cf,"try-error")) M[r,3:6] <- c(cf$eta_lcb, cf$eta_se,
                                                   cf$eta_lcb_pooled, cf$eta_se_pooled)
      if (!inherits(cf,"try-error")) M[r,2] <- M[r,2]   # keep
      if (!inherits(cf,"try-error")) attr(M,"eta") <- cf$eta
    }
    # recompute point estimates for reporting
    set.seed(20261109); ests <- matrix(NA_real_, REPS, 2)
    for (r in seq_len(REPS)) {
      df <- simulate_data_binaryM_x(cells, p_x=c(0.5,0.5), n=n, p_d=0.5)
      s  <- try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=0,
                                   method="split", split_frac=0.3), silent=TRUE)
      cf <- try(lb_frac_affected_x_cf(df,"d","m","y",x="x",at_group=0), silent=TRUE)
      if (!inherits(s,"try-error"))  ests[r,1] <- s$eta
      if (!inherits(cf,"try-error")) ests[r,2] <- cf$eta
    }
    rows <- list(
      split   = list(lcb=M[,1], se=M[,2], est=ests[,1]),
      cf_cons = list(lcb=M[,3], se=M[,4], est=ests[,2]),
      cf_pool = list(lcb=M[,5], se=M[,6], est=ests[,2]))
    for (nm in names(rows)) { z <- rows[[nm]]
      cat(sprintf("%4.1f %6d %-9s %9.5f %9.5f %9.5f %8.3f %8.3f\n", tt, n, nm, bd,
          mean(z$est, na.rm=TRUE), mean(z$se, na.rm=TRUE),
          mean(z$lcb > 0, na.rm=TRUE), mean(z$lcb <= bd + 1e-12, na.rm=TRUE))) }
    # how much of the conservative SE is slack?
    cat(sprintf("     empirical sd of the cross-fitted estimate: %.5f  (cons SE %.5f, pooled SE %.5f)\n",
        sd(ests[,2], na.rm=TRUE), mean(M[,4], na.rm=TRUE), mean(M[,6], na.rm=TRUE)))
  }
}
