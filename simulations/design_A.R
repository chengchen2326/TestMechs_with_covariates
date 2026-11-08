# =====================================================================
# Design A: binary Y, binary M, one binary covariate.
#
# Two calibrations:
#   "sign"     the Section 4 example -- the binding outcome value differs
#              across cells, so the pooled bound is exactly 0 while the
#              conditional bound is 0.06.  Any power here is attributable to
#              conditioning.
#   "bursztyn" the empirical Bursztyn distribution split by college_deg, where
#              the conditional and pooled bounds nearly coincide.  This is the
#              honest calibration and it is expected to show little gain.
#
# Four arms:
#   ours      split-sample, covariates          (lb_frac_affected_x, "split")
#   pooled    split-sample, no covariates       ("pooled_split")
#   plugin_x  full-sample sup, covariates       ("plugin")
#   plugin_0  full-sample sup, no covariates    -- i.e. current practice
# =====================================================================

suppressMessages({library(purrr); library(Rglpk)})

SRC <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_simulate_data_binaryM.R", "src_simulate_data_binaryM_x.R",
            "src_lb_frac_affected_x.R")) source(file.path(SRC, f))

# ------------------------------------------------------------ calibrations
calib_sign <- list(
  list(p_y_m0d0 = c(0.65, 0.35), p_y_m0d1 = c(0.50, 0.50),
       p_y_m1d0 = c(0.50, 0.50), p_y_m1d1 = c(0.50, 0.50),
       p_m_0 = 0.6, p_m_1 = 0.6),
  list(p_y_m0d0 = c(0.50, 0.50), p_y_m0d1 = c(0.65, 0.35),
       p_y_m1d0 = c(0.50, 0.50), p_y_m1d1 = c(0.50, 0.50),
       p_m_0 = 0.6, p_m_1 = 0.6))

# Bursztyn cell counts by college_deg, (y, m) within arm, from the data.
bursztyn_cells <- function() {
  load("/home/claude/repo/data/burstzyn_data.rda")
  bz <- get(ls(pattern = "burstzyn|bursztyn")[1])
  bz <- bz[!is.na(bz$index) & !is.na(bz$college_deg), ]
  out <- list(); px <- numeric(0)
  for (v in sort(unique(bz$college_deg))) {
    s <- bz[bz$college_deg == v, ]
    f <- function(dv, mv) {
      z <- s[s$condition2 == dv & s$signed_up_number == mv, "applied_out_fl"]
      z <- z[!is.na(z)]
      p <- c(mean(z == 0), mean(z == 1))
      if (any(is.na(p))) c(0.5, 0.5) else p
    }
    out[[length(out) + 1]] <- list(
      p_y_m0d0 = f(0, 0), p_y_m0d1 = f(1, 0),
      p_y_m1d0 = f(0, 1), p_y_m1d1 = f(1, 1),
      p_m_0 = mean(s$signed_up_number[s$condition2 == 0] == 1),
      p_m_1 = mean(s$signed_up_number[s$condition2 == 1] == 1))
    px <- c(px, nrow(s))
  }
  list(cells = out, p_x = px / sum(px))
}

# ------------------------------------------------------------------- run
run_design_A <- function(calibration = c("sign", "bursztyn"),
                         n_grid = c(284, 500, 1000, 2000, 5000),
                         n_reps = 500, at_group = 0, alpha = 0.05,
                         seed = 20261108) {

  calibration <- match.arg(calibration)
  if (calibration == "sign") {
    cells <- calib_sign; p_x <- c(0.5, 0.5)
  } else {
    bc <- bursztyn_cells(); cells <- bc$cells; p_x <- bc$p_x
  }

  truth <- population_bound_binaryM_x(cells, p_x, at_group = at_group)
  cat(sprintf("\ncalibration = %s\n", calibration))
  cat(sprintf("  eta by cell      : %s\n",
              paste(sprintf("%.4f", truth$eta_by_cell), collapse = ", ")))
  cat(sprintf("  E[eta(X)] = %.4f   eta_pooled = %.4f   (gap %.4f)\n",
              truth$eta_conditional, truth$eta_pooled,
              truth$eta_conditional - truth$eta_pooled))
  cat(sprintf("  nu  cond  = %.4f   nu_pooled  = %.4f\n\n",
              truth$nu_conditional, truth$nu_pooled))

  set.seed(seed)
  rows <- list()
  for (n in n_grid) {
    acc <- list(ours = NULL, pooled = NULL, plugin_x = NULL, plugin_0 = NULL)
    for (r in seq_len(n_reps)) {
      df <- simulate_data_binaryM_x(cells, p_x = p_x, n = n, p_d = 0.5)
      df$const <- 1L
      fit <- function(xv, mth)
        tryCatch(lb_frac_affected_x(df, "d", "m", "y", x = xv,
                                    at_group = at_group, method = mth,
                                    alpha = alpha),
                 error = function(e) NULL)
      z <- list(ours     = fit("x",     "split"),
                pooled   = fit("const", "pooled_split"),
                plugin_x = fit("x",     "plugin"),
                plugin_0 = fit("const", "plugin"))
      for (nm in names(z)) if (!is.null(z[[nm]]))
        acc[[nm]] <- rbind(acc[[nm]], data.frame(
          eta = z[[nm]]$eta, lcb = z[[nm]]$eta_lcb, nu = z[[nm]]$nu,
          nu_lcb = z[[nm]]$nu_lcb,
          deg = z[[nm]]$n_cells_degenerate, mass = z[[nm]]$mass_degenerate))
    }
    for (nm in names(acc)) {
      A <- acc[[nm]]; if (is.null(A)) next
      target <- if (nm %in% c("pooled", "plugin_0")) truth$eta_pooled
                else truth$eta_conditional
      rows[[length(rows) + 1]] <- data.frame(
        calibration = calibration, n = n, arm = nm, reps = nrow(A),
        mean_eta   = mean(A$eta),
        bias       = mean(A$eta) - target,
        mean_lcb   = if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb, na.rm = TRUE),
        power      = if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb > 0, na.rm = TRUE),
        # Coverage is checked against the arm's OWN population bound. Using the
        # conditional bound for the pooled arm would make its job artificially
        # easy, since eta_pooled <= E[eta(X)].
        coverage   = if (all(is.na(A$lcb))) NA_real_
                     else mean(A$lcb <= target + 1e-12, na.rm = TRUE),
        mean_nu_lcb = if (all(is.na(A$nu_lcb))) NA_real_ else mean(A$nu_lcb, na.rm = TRUE),
        deg_cells  = mean(A$deg), deg_mass = mean(A$mass))
    }
    cat(sprintf("  n = %5d done\n", n))
  }
  res <- do.call(rbind, rows)
  attr(res, "truth") <- truth
  res
}

if (sys.nframe() == 0L) {
  n_reps <- as.integer(Sys.getenv("TM_REPS", "500"))
  out <- NULL
  for (cal in c("sign", "bursztyn")) {
    r <- run_design_A(calibration = cal, n_reps = n_reps)
    out <- rbind(out, r)
  }
  dir.create("results", showWarnings = FALSE)
  write.csv(out, "results/design_A.csv", row.names = FALSE)
  cat("\n")
  print(out[, c("calibration", "n", "arm", "mean_eta", "bias",
                "mean_lcb", "power", "coverage", "deg_mass")],
        row.names = FALSE, digits = 4)
}
