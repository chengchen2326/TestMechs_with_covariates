suppressMessages({library(purrr); library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_simulate_data_binaryM.R","src_simulate_data_binaryM_x.R",
            "src_lb_frac_affected_x.R")) source(file.path("../src", f))

load("/home/claude/repo/data/burstzyn_data.rda")
bz <- get(ls(pattern="burstzyn|bursztyn")[1])
bz <- bz[!is.na(bz$index) & !is.na(bz$college_deg), ]
bz$const <- 1L
cat(sprintf("Bursztyn analysis sample: n = %d, treated = %d\n\n", nrow(bz), sum(bz$condition2)))

kr <- as.numeric(lb_frac_affected(bz, d="condition2", m="signed_up_number",
                                  y="applied_out_fl", at_group=0))
cat(sprintf("KR published / lb_frac_affected point estimate : nu_0 = %.4f\n\n", kr))

# Distribution over 400 random splits, since n = 284 makes the split material.
set.seed(20261108)
res <- list()
for (arm in c("pooled_split","split")) {
  xv <- if (arm == "pooled_split") "const" else "college_deg"
  M <- t(replicate(400, {
    r <- lb_frac_affected_x(bz, "condition2","signed_up_number","applied_out_fl",
                            x=xv, at_group=0, method=arm, alpha=0.05)
    c(nu=r$nu, se=r$nu_se, lcb=r$nu_lcb)
  }))
  res[[arm]] <- M
  cat(sprintf("%-13s nu_hat  median %.4f  [%.4f, %.4f] across splits\n",
              arm, median(M[,"nu"]), quantile(M[,"nu"],.1), quantile(M[,"nu"],.9)))
  cat(sprintf("%-13s 95%% LCB median %.4f  [%.4f, %.4f];  share of splits with LCB > 0: %.2f\n\n",
              "", median(M[,"lcb"]), quantile(M[,"lcb"],.1), quantile(M[,"lcb"],.9),
              mean(M[,"lcb"] > 0)))
}

# Plug-in, full sample, by cell -- what a naive covariate analysis would report
pl <- lb_frac_affected_x(bz,"condition2","signed_up_number","applied_out_fl",
                         x="college_deg", at_group=0, method="plugin")
cat(sprintf("plugin by college_deg, full sample: nu = %.4f  (eta = %.4f, no valid CI)\n", pl$nu, pl$eta))
for (v in c(0,1)) {
  s <- bz[bz$college_deg==v,]
  r <- as.numeric(lb_frac_affected(s, d="condition2", m="signed_up_number", y="applied_out_fl", at_group=0))
  cat(sprintf("   college_deg = %d: n = %3d, nu_0 = %.4f\n", v, nrow(s), r))
}
