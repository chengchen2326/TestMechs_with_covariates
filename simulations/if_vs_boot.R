# Referee the influence function against a nonparametric bootstrap, on the one
# configuration where the audit says the reported SE deteriorates with n:
# homog (compliers present) at an interior mediator value.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_lb_frac_affected_x.R","src_simulate_data_multiM_x.R")) source(file.path("../src",f))
source("calib_C.R"); source("violation_scale.R"); source("design_C.R")

B <- build_config("homog","maximal",5,5,t=0)
set.seed(11)
cat(sprintf("%2s %6s %11s %11s %11s %9s\n","k","n","IF SE","bootstrap","emp sd (MC)","IF/boot"))
for (k in c(3,5)) for (n in c(4000, 16000)) {
  # one dataset -> IF SE and a bootstrap SE on the SAME fold split
  df <- simulate_data_multiM_x(B$cells, B$p_x, n=n, p_d=0.5); df$const <- 1L
  set.seed(1); fv <- rep(2L,n); fv[sample.int(n, floor(0.3*n))] <- 1L
  base <- lb_frac_affected_x(df,"d","m","y",x="x",at_group=k,method="split",folds=fv)
  bs <- replicate(300, {
    i <- sample.int(n, n, replace=TRUE)
    r <- try(lb_frac_affected_x(df[i,],"d","m","y",x="x",at_group=k,
                                method="split",folds=fv[i]), silent=TRUE)
    if (inherits(r,"try-error")) NA_real_ else r$eta })
  # Monte Carlo dispersion over fresh datasets, for reference
  mc <- replicate(300, {
    d2 <- simulate_data_multiM_x(B$cells, B$p_x, n=n, p_d=0.5); d2$const <- 1L
    r <- try(lb_frac_affected_x(d2,"d","m","y",x="x",at_group=k,method="split"), silent=TRUE)
    if (inherits(r,"try-error")) NA_real_ else r$eta })
  cat(sprintf("%2d %6d %11.6f %11.6f %11.6f %9.3f%s\n", k, n, base$eta_se,
      sd(bs,na.rm=TRUE), sd(mc,na.rm=TRUE), base$eta_se/sd(bs,na.rm=TRUE),
      if (base$eta_se/sd(bs,na.rm=TRUE) < 0.95) "  <-- IF below bootstrap" else ""))
}
