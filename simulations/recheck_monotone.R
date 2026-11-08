suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_lb_frac_affected_x.R")) source(file.path("../src",f))

cat("=== 1. Bursztyn regression check (must still be 4e-17) ===\n")
load("/home/claude/repo/data/burstzyn_data.rda")
bz <- get(ls(pattern="burstzyn|bursztyn")[1]); bz <- bz[!is.na(bz$index),]; bz$const <- 1L
a <- as.numeric(lb_frac_affected(bz,d="condition2",m="signed_up_number",y="applied_out_fl",at_group=0))
b <- lb_frac_affected_x(bz,"condition2","signed_up_number","applied_out_fl",x="const",at_group=0,method="plugin")
cat(sprintf("  %.10f vs %.10f   diff %.2e   nonmono cells %d\n\n", a, b$nu, abs(a-b$nu), b$n_cells_nonmonotone))

cat("=== 2. Baranov binary (grandmother): did the guard change anything? ===\n")
load("/home/claude/repo/data/baranov_data.rda")
br <- haven::zap_labels(mother_data); br <- br[!is.na(br$dep_sample) & br$dep_sample==1,]
br <- br[!is.na(br$treat)&!is.na(br$grandmother)&!is.na(br$motherfinancial)&!is.na(br$grandmother_baseline),]
for (v in c("treat","grandmother","motherfinancial","grandmother_baseline")) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial)); br$const <- 1L
cat(sprintf("  n = %d;  P(M=1|D=0) = %.3f,  P(M=1|D=1) = %.3f  -> monotone: %s\n",
    nrow(br), mean(br$grandmother[br$treat==0]), mean(br$grandmother[br$treat==1]),
    mean(br$grandmother[br$treat==1]) >= mean(br$grandmother[br$treat==0])))
for (xv in c("const","grandmother_baseline")) {
  set.seed(3)
  M <- t(vapply(seq_len(400), function(i) {
    r <- lb_frac_affected_x(br,"treat","grandmother","yb", x=xv, at_group=0,
                            method=if(xv=="const") "pooled_split" else "split")
    c(r$nu_lcb, r$n_cells_nonmonotone, r$mass_nonmonotone)}, numeric(3)))
  cat(sprintf("  %-22s LCB>0 %.2f | med LCB %.4f | nonmono cells/draw %.2f | mass %.3f\n",
      xv, mean(M[,1]>0), median(M[,1]), mean(M[,2]), mean(M[,3])))
}

cat("\n=== 3. finer Bursztyn splits: how much mass is lost to non-monotone cells? ===\n")
bz2 <- bz[!is.na(bz$college_deg)&!is.na(bz$employed_wife)&!is.na(bz$vjobs_self)&!is.na(bz$employed_now),]
specs <- list(c("college_deg"), c("college_deg","employed_wife"),
              c("college_deg","employed_wife","vjobs_self"),
              c("college_deg","employed_wife","vjobs_self","employed_now"))
for (sp in specs) {
  set.seed(5)
  M <- t(vapply(seq_len(200), function(i) {
    r <- lb_frac_affected_x(bz2,"condition2","signed_up_number","applied_out_fl",
                            x=sp, at_group=0, method="split")
    c(r$eta_lcb, r$n_cells, r$n_cells_degenerate, r$n_cells_nonmonotone, r$mass_degenerate, r$max_cell_defier_share)}, numeric(6)))
  cat(sprintf("  %d cov: cells %.1f | degenerate %.1f | lost mass %.3f | cells w/ defiers %.1f | max cell defier %.3f | LCB>0 %.2f\n",
      length(sp), mean(M[,2]), mean(M[,3]), mean(M[,5]), mean(M[,4]), mean(M[,6]), mean(M[,1]>0)))
}
