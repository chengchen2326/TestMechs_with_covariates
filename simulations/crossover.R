# When does conditioning beat pooling?
#
# "Choice of test" in KR is a clean recommendation because their four tests all
# target the same object. Ours do not: the conditional and pooled procedures
# target different estimands, so a recommendation has to say WHEN, not WHICH.
#
# The economics of the choice are simple. Conditioning buys an identification
# gain, g = E[eta_k(X)] - eta_k^pooled, and costs extra sampling noise. Power
# favours conditioning once the gain outgrows the noise, and since the standard
# error falls like 1/sqrt(n), the crossover should be governed by g * sqrt(n).
# If that is right, the recommendation can be stated as a required sample size
# for a given gain -- the same kind of operational rule as KR's "at least 15
# independent observations per cell".
#
# Three heterogeneity levels from Design C, swept over the violation scale t
# and the sample size: homog (no gain by construction), split (the gain the
# Baranov covariate actually delivers) and signflip (a large constructed gain).
# Only the two arms that matter are run, so the grid is affordable.

suppressMessages({library(haven); library(purrr); library(Rglpk)})
SRC <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_discretize_y.R", "src_lb_frac_affected_x.R",
            "src_simulate_data_multiM_x.R")) source(file.path(SRC, f))
source("calib_C.R")
source("violation_scale.R")
source("design_C.R")          # build_config, signflip_cells

run_crossover <- function(configs = c("homog", "split", "signflip"),
                          t_grid = c(0.25, 0.5, 1),
                          n_grid = c(500, 1000, 2000, 4000, 8000),
                          k = 3, n_reps = 250, seed = 20261109) {
  set.seed(seed)
  rows <- list()
  for (cfg in configs) for (tt in t_grid) {
    B <- build_config(cfg, "maximal", 5, 5, t = tt)
    pb <- population_bound_multiM_x(B$cells, B$p_x, at_group = k)
    gap <- pb$eta_conditional - pb$eta_pooled
    cat(sprintf("\n[%s t=%.2f] cond %.5f  pooled %.5f  gap %.5f\n",
                cfg, tt, pb$eta_conditional, pb$eta_pooled, gap))
    for (n in n_grid) {
      po <- pp <- numeric(0)
      for (r in seq_len(n_reps)) {
        df <- simulate_data_multiM_x(B$cells, B$p_x, n = n, p_d = 0.5)
        df$const <- 1L
        a <- try(lb_frac_affected_x(df, "d","m","y", x = "x", at_group = k,
                                    method = "split"), silent = TRUE)
        b <- try(lb_frac_affected_x(df, "d","m","y", x = "const", at_group = k,
                                    method = "pooled_split"), silent = TRUE)
        if (!inherits(a, "try-error")) po <- c(po, a$eta_lcb)
        if (!inherits(b, "try-error")) pp <- c(pp, b$eta_lcb)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        config = cfg, t = tt, n = n, k = k,
        cond_bound = pb$eta_conditional, pooled_bound = pb$eta_pooled,
        gap = gap, gap_rootn = gap * sqrt(n),
        power_cond = mean(po > 0), power_pool = mean(pp > 0),
        diff = mean(po > 0) - mean(pp > 0),
        lcb_cond = mean(po), lcb_pool = mean(pp))
      cat(sprintf("   n=%5d  cond %.3f  pooled %.3f  diff %+.3f\n",
                  n, mean(po > 0), mean(pp > 0), mean(po > 0) - mean(pp > 0)))
    }
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  out <- run_crossover(k = as.integer(Sys.getenv("TM_K", "3")),
                       n_reps = as.integer(Sys.getenv("TM_REPS", "250")))
  dir.create("results", showWarnings = FALSE)
  write.csv(out, sprintf("results/crossover_k%s.csv", Sys.getenv("TM_K", "3")),
            row.names = FALSE)
  cat("\n")
  print(out[, c("config","t","n","gap","gap_rootn","power_cond","power_pool","diff")],
        row.names = FALSE, digits = 3)
}
