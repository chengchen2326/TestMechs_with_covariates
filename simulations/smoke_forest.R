suppressMessages({library(purrr);library(Rglpk);library(ranger)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_lb_frac_affected_x.R","src_lb_frac_affected_x_cts.R")) source(file.path("../src",f))
source("design_B.R")   # for dgp_B / truth_B / cut_cols
tr <- truth_B()
cat(sprintf("truth E[eta_0(X)] = %.4f, pooled = 0\n\n", tr$eta_cond))
set.seed(11); REPS <- 60
cat(sprintf("%-7s %-6s %-16s %8s %8s %9s\n","signal","n","arm","power","cover","mean_lcb"))
for (sg in c("x1x2")) for (n in c(2000, 8000)) {
  acc <- list()
  for (r in seq_len(REPS)) {
    df <- dgp_B(n, 5, signal=sg); xs <- paste0("x",1:5)
    z <- list(
      cts_glm    = try(lb_frac_affected_x_cts(df,"d","m","y",x=xs,at_group=0,pi_true=0.5,learner="glm"),silent=TRUE),
      cts_forest = try(lb_frac_affected_x_cts(df,"d","m","y",x=xs,at_group=0,pi_true=0.5,learner="forest"),silent=TRUE))
      da <- cut_cols(df, xs, 3)
      z$binned_all <- try(lb_frac_affected_x(da,"d","m","y",x=paste0("b_",xs),at_group=0,method="split"),silent=TRUE)
      do2 <- cut_cols(df, c("x1","x2"), 5)
      z$binned_oracle <- try(lb_frac_affected_x(do2,"d","m","y",x=c("b_x1","b_x2"),at_group=0,method="split"),silent=TRUE)
    for (nm in names(z)) { if (inherits(z[[nm]],"try-error")) next
      acc[[nm]] <- rbind(acc[[nm]], data.frame(lcb=z[[nm]]$eta_lcb)) }
  }
  for (nm in names(acc)) { A <- acc[[nm]]
    cat(sprintf("%-7s %-6d %-16s %8.3f %8.3f %9.4f\n", sg, n, nm,
        mean(A$lcb>0), mean(A$lcb <= tr$eta_cond + 1e-12), mean(A$lcb))) }
}
