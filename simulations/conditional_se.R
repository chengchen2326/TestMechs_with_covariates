# Is the influence-function SE the CONDITIONAL standard deviation given Ahat?
# Fix fold 1 (hence Ahat) once, then redraw fold 2 many times. If the IF SE
# matches the dispersion of eta_hat around its own conditional mean, the
# estimator is right and the earlier gap was the estimand eta^{Ahat} moving
# with Ahat -- which sample splitting does not have to pay for, since any fixed
# A gives a valid lower bound on eta_k.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_lb_frac_affected_x.R","src_simulate_data_multiM_x.R")) source(file.path("../src",f))
source("calib_C.R"); source("violation_scale.R"); source("design_C.R")

B <- build_config("homog","maximal",5,5,t=0)
pb3 <- population_bound_multiM_x(B$cells, B$p_x, at_group=3)
cat(sprintf("population bound at k=3: %.6f\n\n", pb3$eta_conditional))
cat(sprintf("%2s %6s %12s %12s %9s %12s\n","k","n2","IF SE (mean)","cond. emp sd","ratio","cond cover"))
for (k in c(3,5)) for (n in c(4000, 16000)) {
  n1 <- floor(0.3*n); n2 <- n - n1
  # one fold-1 sample, fixed for the whole experiment
  set.seed(77)
  d1 <- simulate_data_multiM_x(B$cells, B$p_x, n=n1, p_d=0.5); d1$const <- 1L
  est <- se <- numeric(400)
  for (r in 1:400) {
    d2 <- simulate_data_multiM_x(B$cells, B$p_x, n=n2, p_d=0.5); d2$const <- 1L
    df <- rbind(d1, d2); fv <- c(rep(1L, nrow(d1)), rep(2L, nrow(d2)))
    z <- try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=k,method="split",folds=fv),
             silent=TRUE)
    if (inherits(z,"try-error")) { est[r] <- NA; se[r] <- NA } else { est[r] <- z$eta; se[r] <- z$eta_se }
  }
  pbk <- population_bound_multiM_x(B$cells, B$p_x, at_group=k)$eta_conditional
  lcb <- pmax(est - qnorm(.95)*se, 0)
  cat(sprintf("%2d %6d %12.6f %12.6f %9.3f %12.3f\n", k, n2,
      mean(se,na.rm=TRUE), sd(est,na.rm=TRUE), mean(se,na.rm=TRUE)/sd(est,na.rm=TRUE),
      mean(lcb <= pbk + 1e-12, na.rm=TRUE)))
}
