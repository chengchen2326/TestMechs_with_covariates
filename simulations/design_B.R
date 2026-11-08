# =====================================================================
# Design B: continuous covariates, where binning is not available.
#
# Binary Y, binary M, X ~ U[0,1]^p.  No compliers (the M marginals agree across
# arms), so monotonicity is trivial and eta_0(x) has a closed form:
#
#   P(M = 1 | D, X) = pm            (same in both arms)
#   q_d(x) = P(Y = 1 | D = d, M = 0, X = x) = 0.5 +- delta(x)/2
#
#   eta_0(x)     = (1 - pm) * |delta(x)|
#   eta_0^pooled = (1 - pm) * |E[delta(X)]|
#
# Two signals, both with E[delta] = 0 and E|delta| = amp * 2/pi, so the truth is
# the same for both and only the difficulty of recovering it changes:
#
#   "x1"    delta = amp * sin(2*pi*x_1)            -- one coordinate
#   "x1x2"  delta = amp * sin(2*pi*(x_1 + x_2))    -- needs two, jointly
#
# (x_1 + x_2 mod 1 is uniform, so E[delta] = 0 and E|delta| = amp*2/pi in both.)
# Coordinates beyond those are pure nuisance.
#
# Arms:
#   cts_glm        ours, logistic first stage in X          (misspecified for sin)
#   cts_poly       ours, logistic first stage in a cubic basis
#   binned_oracle  the discrete route, binning ONLY the signal coordinate(s).
#                  This is an oracle: it is told which covariates matter. It is
#                  included because it does well and the paper should say so.
#   binned_all     the discrete route as actually feasible: bin every covariate
#                  and interact.  Without knowing which covariate matters this
#                  is what KR's suggestion requires, and the cells collapse.
#   pooled         no covariates
# =====================================================================

suppressMessages({library(purrr); library(Rglpk)})
SRC <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_lb_frac_affected_x.R", "src_lb_frac_affected_x_cts.R"))
  source(file.path(SRC, f))

PM  <- 0.6
AMP <- 0.30
J_ORACLE <- 5    # bins per covariate for the oracle arm
J_ALL    <- 3    # bins per covariate when every covariate is binned

truth_B <- function(amp = AMP, pm = PM)
  list(eta_cond = (1 - pm) * amp * 2 / pi, eta_pooled = 0)

dgp_B <- function(n, p, signal = c("x1", "x1x2"), amp = AMP, pm = PM) {
  signal <- match.arg(signal)
  X <- matrix(stats::runif(n * p), n, p)
  d <- stats::rbinom(n, 1, 0.5)
  m <- stats::rbinom(n, 1, pm)
  delta <- if (signal == "x1") amp * sin(2 * pi * X[, 1])
           else                amp * sin(2 * pi * (X[, 1] + X[, 2]))
  q <- 0.5 + ifelse(d == 1, 1, -1) * delta / 2
  y <- integer(n)
  y[m == 0] <- stats::rbinom(sum(m == 0), 1, q[m == 0])
  y[m == 1] <- stats::rbinom(sum(m == 1), 1, 0.5)   # no violation in M = 1
  df <- data.frame(d = d, m = m, y = y, X)
  names(df)[-(1:3)] <- paste0("x", seq_len(p))
  df$const <- 1L
  df
}

cut_cols <- function(df, cols, J) {
  br <- seq(0, 1, length.out = J + 1)
  out <- lapply(cols, function(cc) cut(df[[cc]], breaks = br, include.lowest = TRUE))
  names(out) <- paste0("b_", cols)
  cbind(df, as.data.frame(out))
}

run_design_B <- function(signals = c("x1", "x1x2"),
                         n_grid = c(500, 2000, 8000), p_grid = c(2, 5),
                         n_reps = 400, seed = 20261108) {
  tr <- truth_B()
  cat(sprintf("truth: E[eta_0(X)] = %.4f, eta_0^pooled = 0 exactly\n", tr$eta_cond))
  cat(sprintf("oracle bins per covariate = %d; all-covariate bins = %d\n\n",
              J_ORACLE, J_ALL))
  set.seed(seed)
  rows <- list()
  for (sg in signals) for (n in n_grid) for (p in p_grid) {
    sig_cols <- if (sg == "x1") "x1" else c("x1", "x2")
    if (p < length(sig_cols)) next
    acc <- list()
    for (r in seq_len(n_reps)) {
      df <- dgp_B(n, p, signal = sg)
      xs <- paste0("x", seq_len(p))
      dfp <- df
      for (cc in xs) { dfp[[paste0(cc, "_2")]] <- df[[cc]]^2
                       dfp[[paste0(cc, "_3")]] <- df[[cc]]^3 }
      xp <- c(xs, paste0(xs, "_2"), paste0(xs, "_3"))

      z <- list()
      z$cts_glm  <- try(lb_frac_affected_x_cts(df,  "d","m","y", x = xs,
                          at_group = 0, pi_true = 0.5), silent = TRUE)
      z$cts_poly <- try(lb_frac_affected_x_cts(dfp, "d","m","y", x = xp,
                          at_group = 0, pi_true = 0.5), silent = TRUE)
      z$pooled   <- try(lb_frac_affected_x(df, "d","m","y", x = "const",
                          at_group = 0, method = "pooled_split"), silent = TRUE)
      do <- cut_cols(df, sig_cols, J_ORACLE)
      z$binned_oracle <- try(lb_frac_affected_x(do, "d","m","y",
                          x = paste0("b_", sig_cols), at_group = 0,
                          method = "split"), silent = TRUE)
      da <- cut_cols(df, xs, J_ALL)
      z$binned_all <- try(lb_frac_affected_x(da, "d","m","y",
                          x = paste0("b_", xs), at_group = 0,
                          method = "split"), silent = TRUE)

      for (nm in names(z)) {
        if (inherits(z[[nm]], "try-error")) next
        acc[[nm]] <- rbind(acc[[nm]], data.frame(
          eta = z[[nm]]$eta, lcb = z[[nm]]$eta_lcb,
          cells = if (is.null(z[[nm]]$n_cells)) 1 else z[[nm]]$n_cells,
          deg = if (is.null(z[[nm]]$mass_degenerate)) 0 else z[[nm]]$mass_degenerate))
      }
    }
    for (nm in names(acc)) {
      A <- acc[[nm]]
      target <- if (nm == "pooled") tr$eta_pooled else tr$eta_cond
      rows[[length(rows) + 1]] <- data.frame(
        signal = sg, n = n, p = p, arm = nm, reps = nrow(A),
        mean_eta = mean(A$eta), mean_lcb = mean(A$lcb),
        power = mean(A$lcb > 0), coverage = mean(A$lcb <= target + 1e-12),
        cells = mean(A$cells), deg_mass = mean(A$deg))
    }
    cat(sprintf("  signal = %-5s n = %5d p = %2d done\n", sg, n, p))
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  out <- run_design_B(n_reps = as.integer(Sys.getenv("TM_REPS", "400")))
  dir.create("results", showWarnings = FALSE)
  write.csv(out, "results/design_B.csv", row.names = FALSE)
  cat("\n")
  print(out[, c("signal","n","p","arm","mean_lcb","power","coverage","cells","deg_mass")],
        row.names = FALSE, digits = 3)
}
