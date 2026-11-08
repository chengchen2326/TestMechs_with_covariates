suppressMessages({library(haven);library(purrr);library(Rglpk);library(dplyr);library(quadprog);library(magrittr)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_discretize_y.R","src_nonparametric_bootstrap.R","src_bayes_bootstrap.R","src_cox_shi_nonuisance.R","src_test_sharp_null_coxandshi_binary_m.R",
            "src_lb_frac_affected_x.R")) source(file.path("../src", f))
args <- commandArgs(TRUE)
show <- function(lbl, df, d, m, y) {
  r <- tryCatch(test_sharp_null_coxandshi_binary_m(df=df, d=d, m=m, y=y),
                error=function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-34s ", lbl)); print(unlist(r))
}
load("/home/claude/repo/data/burstzyn_data.rda")
bz <- get(ls(pattern="burstzyn|bursztyn")[1]); bz <- bz[!is.na(bz$index),]
show("Bursztyn pooled (KR Cox-Shi)", bz, "condition2","signed_up_number","applied_out_fl")

load("/home/claude/repo/data/baranov_data.rda")
br <- mother_data
br <- br[!is.na(br$dep_sample) & br$dep_sample==1, ]
b2 <- br[!is.na(br$treat) & !is.na(br$grandmother) & !is.na(br$motherfinancial), ]
b2 <- haven::zap_labels(b2)
for (v in c("treat","grandmother","motherfinancial","grandmother_baseline"))
  if (v %in% names(b2)) b2[[v]] <- as.numeric(unclass(b2[[v]]))
b2 <- b2[!is.na(b2$grandmother_baseline), ]
b2$const <- 1L
cat(sprintf("\nBaranov n = %d (grandmother mediator)\n", nrow(b2)))
cat("  lb_frac_affected nu_0 = ")
cat(sprintf("%.4f\n", as.numeric(lb_frac_affected(b2, d="treat", m="grandmother",
        y="motherfinancial", at_group=0, num_Ybins=5))))
for (arm in c("pooled_split","split")) {
  xv <- if (arm=="pooled_split") "const" else "grandmother_baseline"
  if (xv != "const" && !(xv %in% names(b2))) { cat("  (no ", xv, ")\n"); next }
  b2$const <- 1L
  set.seed(1); M <- t(replicate(200, {
    r <- lb_frac_affected_x(b2,"treat","grandmother","motherfinancial", x=xv,
                            at_group=0, method=arm); c(r$nu, r$nu_lcb)}))
  cat(sprintf("  %-13s nu median %.4f | 95%% LCB median %.4f | share LCB>0 %.2f\n",
              arm, median(M[,1]), median(M[,2]), mean(M[,2]>0)))
}

cat("\n--- Baranov, KR Cox-Shi test (binary mediator = grandmother) ---\n")
b2$yb <- as.integer(b2$motherfinancial > median(b2$motherfinancial, na.rm=TRUE))
show("Baranov pooled (KR Cox-Shi, Y binarised)", b2, "treat","grandmother","yb")
cat("\n--- Baranov with Y binarised at the median, our procedure ---\n")
for (arm in c("pooled_split","split")) {
  xv <- if (arm=="pooled_split") "const" else "grandmother_baseline"
  set.seed(1); M <- t(replicate(300, {
    r <- lb_frac_affected_x(b2,"treat","grandmother","yb", x=xv,
                            at_group=0, method=arm); c(r$nu, r$nu_lcb)}))
  cat(sprintf("  %-13s nu median %.4f | 95%% LCB median %.4f | share LCB>0 %.2f\n",
              arm, median(M[,1]), median(M[,2]), mean(M[,2]>0)))
}
pl <- lb_frac_affected_x(b2,"treat","grandmother","yb", x="grandmother_baseline",
                         at_group=0, method="plugin")
cat(sprintf("  plugin by cell (full sample): nu = %.4f   [no valid CI]\n", pl$nu))
for (v in c(0,1)) { s <- b2[b2$grandmother_baseline==v,]
  cat(sprintf("     grandmother_baseline = %d: n = %3d, nu_0 = %.4f\n", v, nrow(s),
      as.numeric(lb_frac_affected(s, d="treat", m="grandmother", y="yb", at_group=0)))) }
