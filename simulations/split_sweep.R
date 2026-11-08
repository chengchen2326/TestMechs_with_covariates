# How much of the sample should stage 1 get?  Any A is dual-feasible, so giving
# stage 1 less data costs tightness, never validity.  Stage 2 is where the CLT
# happens, so the optimum should be well below 1/2.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_simulate_data_binaryM.R","src_simulate_data_binaryM_x.R",
            "src_lb_frac_affected_x.R")) source(file.path("../src", f))

cells <- list(
  list(p_y_m0d0=c(0.65,0.35), p_y_m0d1=c(0.50,0.50),
       p_y_m1d0=c(0.50,0.50), p_y_m1d1=c(0.50,0.50), p_m_0=0.6, p_m_1=0.6),
  list(p_y_m0d0=c(0.50,0.50), p_y_m0d1=c(0.65,0.35),
       p_y_m1d0=c(0.50,0.50), p_y_m1d1=c(0.50,0.50), p_m_0=0.6, p_m_1=0.6))
truth <- population_bound_binaryM_x(cells, c(0.5,0.5), 0)

fracs <- c(0.10, 0.20, 0.30, 0.50, 0.70)
REPS <- as.integer(Sys.getenv("TM_REPS","600"))
cat("=== sign calibration (truth E[eta(X)] = 0.06, pooled = 0) ===\n")
cat(sprintf("%6s %6s %8s %8s %8s %9s\n","n","frac","power","mean_lcb","cover","mean_eta"))
set.seed(7)
for (n in c(284, 500, 1000, 2000)) {
  for (fr in fracs) {
    M <- t(replicate(REPS, {
      df <- simulate_data_binaryM_x(cells, p_x=c(0.5,0.5), n=n, p_d=0.5)
      r <- lb_frac_affected_x(df,"d","m","y",x="x",at_group=0,
                              method="split", split_frac=fr)
      c(r$eta_lcb, r$eta)}))
    cat(sprintf("%6d %6.2f %8.3f %8.4f %8.3f %9.4f\n", n, fr,
                mean(M[,1]>0), mean(M[,1]),
                mean(M[,1] <= truth$eta_conditional + 1e-12), mean(M[,2])))
  }
}

cat("\n=== real data: share of splits with 95% LCB > 0 ===\n")
load("/home/claude/repo/data/burstzyn_data.rda")
bz <- get(ls(pattern="burstzyn|bursztyn")[1]); bz <- bz[!is.na(bz$index),]; bz$const <- 1L
load("/home/claude/repo/data/baranov_data.rda")
br <- haven::zap_labels(mother_data); br <- br[!is.na(br$dep_sample) & br$dep_sample==1,]
br <- br[!is.na(br$treat) & !is.na(br$grandmother) & !is.na(br$motherfinancial) &
         !is.na(br$grandmother_baseline),]
for (v in c("treat","grandmother","motherfinancial","grandmother_baseline"))
  br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial)); br$const <- 1L

cat(sprintf("%-10s %-22s %6s %8s %8s\n","data","spec","frac","LCB>0","med_LCB"))
for (cfg in list(list("Bursztyn",bz,"condition2","signed_up_number","applied_out_fl","const"),
                 list("Bursztyn",bz,"condition2","signed_up_number","applied_out_fl","college_deg"),
                 list("Baranov", br,"treat","grandmother","yb","const"),
                 list("Baranov", br,"treat","grandmother","yb","grandmother_baseline"))) {
  for (fr in fracs) {
    set.seed(3)
    v <- vapply(seq_len(400), function(i)
      lb_frac_affected_x(cfg[[2]],cfg[[3]],cfg[[4]],cfg[[5]], x=cfg[[6]],
                         at_group=0, method="split", split_frac=fr)$nu_lcb,
      numeric(1))
    cat(sprintf("%-10s %-22s %6.2f %8.2f %8.4f\n", cfg[[1]], cfg[[6]], fr,
                mean(v>0), median(v)))
  }
}
