# =====================================================================
# kr_compare.R -- the comparison the deck was missing: our LCB criterion
# against KR's own test of the sharp null.
#
# Run as:  TM_SRC=../src Rscript kr_compare.R
# Output:  results/kr_compare.csv  and a printed table.
#
# Two calibrations, both from design_A.R:
#   sign      pooled bound exactly 0, conditional bound 0.06. If conditioning
#             is what buys power, KR's POOLED test should have none here either
#             -- it has the whole sample and aggregates over every k, so if it
#             still cannot see the violation, the violation is genuinely
#             invisible without covariates.
#   bursztyn  conditional and pooled bounds coincide. Conditioning buys no
#             identification, so this is where we should lose, and the size of
#             the loss is the price of splitting plus the price of using a
#             magnitude bound as an existence test.
#
# t = 0 is the sharp null exactly, so the t = 0 rows are size for every arm.
#
# Any output of this script from before 2026-11-10 is void: kr_reject did not
# then restore the RNG stream that KR's bootstrap pins, so every rep after the
# first redrew the same dataset. See the RNG note at the top of kr_arm.R.
# =====================================================================

suppressMessages({library(purrr); library(Rglpk); library(dplyr); library(magrittr)})
SRC <- Sys.getenv("TM_SRC", "../src")
source("violation_scale.R")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_simulate_data_binaryM.R", "src_simulate_data_binaryM_x.R",
            "src_lb_frac_affected_x.R")) source(file.path(SRC, f))
source("kr_arm.R")
source("design_A.R")

n_grid  <- as.numeric(strsplit(Sys.getenv("TM_N", "284,1000,5000"), ",")[[1]])
n_reps  <- as.integer(Sys.getenv("TM_REPS", "300"))
tgrid   <- as.numeric(strsplit(Sys.getenv("TM_T", "0,1"), ",")[[1]])
cals    <- strsplit(Sys.getenv("TM_CAL", "sign,bursztyn"), ",")[[1]]
outfile <- Sys.getenv("TM_OUT", "results/kr_compare.csv")
alpha   <- 0.05
at_group <- 0

# A full grid is about 20 minutes per (calibration, t) block, so a wall-clock
# limit can land mid-run. Each row is written the moment it is finished and
# every number in the file is machine-written: nothing is ever transcribed by
# hand from a log, and a kill costs one block rather than the whole study.
dir.create("results", showWarnings = FALSE)
emit <- function(row) {
  write.table(row, outfile, sep = ",", row.names = FALSE,
              col.names = !file.exists(outfile), append = file.exists(outfile),
              qmethod = "double")
}

set.seed(20261110)
rows <- list()
for (cal in cals) {
  for (tt in tgrid) {
    if (cal == "sign") { cells <- calib_sign; p_x <- c(0.5, 0.5)
    } else { bc <- bursztyn_cells(); cells <- bc$cells; p_x <- bc$p_x }
    cells <- scale_violation_binary(cells, tt)
    truth <- population_bound_binaryM_x(cells, p_x, at_group = at_group)
    cat(sprintf("\n[%s, t = %.1f]  E[eta(X)] = %.4f   eta_pooled = %.4f\n",
                cal, tt, truth$eta_conditional, truth$eta_pooled))
    for (n in n_grid) {
      ours <- pooled <- krrej <- numeric(0); krfail <- 0L
      for (r in seq_len(n_reps)) {
        df <- simulate_data_binaryM_x(cells, p_x = p_x, n = n, p_d = 0.5)
        df$const <- 1L
        f1 <- tryCatch(lb_frac_affected_x(df, "d", "m", "y", x = "x",
                 at_group = at_group, method = "split", alpha = alpha),
                 error = function(e) NULL)
        f0 <- tryCatch(lb_frac_affected_x(df, "d", "m", "y", x = "const",
                 at_group = at_group, method = "pooled_split", alpha = alpha),
                 error = function(e) NULL)
        if (!is.null(f1) && !is.na(f1$eta_lcb)) ours   <- c(ours,   f1$eta_lcb > 0)
        if (!is.null(f0) && !is.na(f0$eta_lcb)) pooled <- c(pooled, f0$eta_lcb > 0)
        kk <- kr_reject(df, "d", "m", "y", method = "CS", alpha = alpha)
        if (is.na(kk$reject)) krfail <- krfail + 1L else krrej <- c(krrej, kk$reject)
      }
      row <- data.frame(
        calibration = cal, t = tt, n = n, reps = n_reps,
        eta_cond = truth$eta_conditional, eta_pooled = truth$eta_pooled,
        ours = mean(ours), pooled = mean(pooled), kr_cs = mean(krrej),
        kr_failures = krfail)
      emit(row)
      rows[[length(rows) + 1L]] <- row
      cat(sprintf("   n = %5d   ours %.3f   pooled %.3f   kr_cs %.3f   (%d kr failures)\n",
                  n, mean(ours), mean(pooled), mean(krrej), krfail))
    }
  }
}

out <- do.call(rbind, rows)
cat("\n")
print(out[, c("calibration", "t", "n", "eta_cond", "eta_pooled",
              "ours", "pooled", "kr_cs")], row.names = FALSE, digits = 3)
cat("\nt = 0 rows are size; t = 1 rows are power.\n")
