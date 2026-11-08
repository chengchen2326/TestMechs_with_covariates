# Inference on the headline: raw Baranov, five-valued mediator, population cap.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_simulate_data_multiM_x.R","src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_theta_kk_min_defiers.R","src_eta_cond_joint.R",
            "src_lb_frac_affected_x_defiers.R")) source(file.path("../src",f))

load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
keep <- c("treat","relationship_husb","motherfinancial","grandmother_baseline","MIL","girl")
br <- br[stats::complete.cases(br[,keep]),]
for (v in keep) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
br$const <- 1L
cat(sprintf("n = %d\n\n", nrow(br)))

run <- function(xv, k, dbar, R=300) {
  set.seed(7)
  M <- t(vapply(seq_len(R), function(i) {
    r <- tryCatch(lb_frac_affected_x_defiers(br,"treat","mc","yb", x=xv, at_group=k,
                    max_defiers_share=dbar), error=function(e) NULL)
    if (is.null(r)) c(NA,NA,NA,NA) else c(r$eta, r$eta_lcb, r$dbar_used, r$n_cells_zeroed)
  }, numeric(4)))
  M <- M[!is.na(M[,1]),,drop=FALSE]
  if (!nrow(M)) return(NULL)
  c(reps=nrow(M), eta=median(M[,1]), lcb=median(M[,2]),
    share_pos=mean(M[,2]>0), dbar=median(M[,3]), zeroed=mean(M[,4]))
}

cat(sprintf("%-22s %2s %6s %5s %9s %9s %10s\n","covariate","k","dbar","reps","med eta","med LCB","LCB>0 share"))
for (xv in c("const","grandmother_baseline","MIL","girl")) {
  for (k in c(3,4)) for (dbar in c(0.03, 0.05)) {
    z <- run(xv, k, dbar)
    if (is.null(z)) { cat(sprintf("%-22s %2d %6.3f  infeasible\n", xv, k, dbar)); next }
    cat(sprintf("%-22s %2d %6.3f %5d %9.5f %9.5f %10.2f\n",
        xv, k, z["dbar"], z["reps"], z["eta"], z["lcb"], z["share_pos"]))
  }
}
