# =====================================================================
# Design B: continuous covariates, where binning is not available.
#
# Binary Y, binary M, X ~ U[0,1]^p.  No compliers (the M marginals agree
# across arms), so monotonicity is trivial and eta_0(x) has a closed form.
#
#   P(M = 1 | D, X) = pm            (same in both arms)
#   q_d(x) = P(Y = 1 | D = d, M = 0, X = x) = 0.5 +- delta(x)/2
#   delta(x) = amp * sin(2*pi*x_1)
#
# Then
#   eta_0(x)      = (1 - pm) * |delta(x)|
#   E[eta_0(X)]   = (1 - pm) * amp * 2/pi
#   eta_0^pooled  = (1 - pm) * |E[delta(X)]| = 0   exactly, since E[sin] = 0.
#
# So the pooled bound is zero by construction while the conditional bound is
# positive -- the sign-changing structure of Design A, now in continuous X,
# where KR's suggested route of applying the test cell by cell has nothing to
# condition on.  Covariates x_2..x_p are pure nuisance: the signal is in x_1
# only, and the first stage has to find it.
#
# Arms:
#   cts_glm   ours, logistic first stage in X          (misspecified for sin)
#   cts_poly  ours, logistic first stage in a cubic basis of x_1
#   binned    lb_frac_affected_x with x_1 cut into J bins (the discrete route)
#   pooled    no covariates
# =====================================================================

suppressMessages({library(purrr); library(Rglpk)})
SRC <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_lb_frac_affected_x.R", "src_lb_frac_affected_x_cts.R"))
  source(file.path(SRC, f))

PM  <- 0.6     # P(M = 1), equal across arms
AMP <- 0.30    # amplitude of delta

truth_B <- function(amp = AMP, pm = PM)
  list(eta_cond = (1 - pm) * amp * 2 / pi, eta_pooled = 0)

dgp_B <- function(n, p, amp = AMP, pm = PM) {
  X <- matrix(stats::runif(n * p), n, p)
  d <- stats::rbinom(n, 1, 0.5)
  m <- stats::rbinom(n, 1, pm)
  delta <- amp * sin(2 * pi * X[, 1])
  q <- 0.5 + ifelse(d == 1, 1, -1) * delta / 2
  y <- integer(n)
  y[m == 0] <- stats::rbinom(sum(m == 0), 1, q[m == 0])
  # Within M = 1 the arms agree, so that stratum carries no violation.
  y[m == 1] <- stats::rbinom(sum(m == 1), 1, 0.5)
  df <- data.frame(d = d, m = m, y = y, X)
  names(df)[-(1:3)] <- paste0("x", seq_len(p))
  df$const <- 1L
  df
}

run_design_B <- function(n_grid = c(500, 2000), p_grid = c(1, 5),
                         bins = c(5, 20), n_reps = 400, seed = 20261108) {
  tr <- truth_B()
  cat(sprintf("truth: E[eta_0(X)] = %.4f, eta_0^pooled = %.4f (exactly 0)\n\n",
              tr$eta_cond, tr$eta_pooled))
  set.seed(seed)
  rows <- list()
  for (n in n_grid) for (p in p_grid) {
    acc <- list()
    for (r in seq_len(n_reps)) {
      df <- dgp_B(n, p)
      xs  <- paste0("x", seq_len(p))
      dfp <- df
      # cubic basis of x1, plus the nuisance covariates untransformed
      dfp$x1_2 <- df$x1^2; dfp$x1_3 <- df$x1^3
      xp <- c(xs, "x1_2", "x1_3")

      z <- list()
      z$cts_glm  <- try(lb_frac_affected_x_cts(df,  "d","m","y", x = xs,
                          at_group = 0, pi_true = 0.5), silent = TRUE)
      z$cts_poly <- try(lb_frac_affected_x_cts(dfp, "d","m","y", x = xp,
                          at_group = 0, pi_true = 0.5), silent = TRUE)
      z$pooled   <- try(lb_frac_affected_x(df, "d","m","y", x = "const",
                          at_group = 0, method = "pooled_split"), silent = TRUE)
      for (J in bins) {
        df$bin <- cut(df$x1, breaks = seq(0, 1, length.out = J + 1),
                      include.lowest = TRUE)
        z[[paste0("binned", J)]] <- try(lb_frac_affected_x(df, "d","m","y",
                          x = "bin", at_group = 0, method = "split"),
                          silent = TRUE)
      }
      for (nm in names(z)) {
        if (inherits(z[[nm]], "try-error")) next
        acc[[nm]] <- rbind(acc[[nm]], data.frame(
          eta = z[[nm]]$eta, lcb = z[[nm]]$eta_lcb,
          deg = if (is.null(z[[nm]]$mass_degenerate)) 0 else z[[nm]]$mass_degenerate))
      }
    }
    for (nm in names(acc)) {
      A <- acc[[nm]]
      target <- if (nm == "pooled") tr$eta_pooled else tr$eta_cond
      rows[[length(rows) + 1]] <- data.frame(
        n = n, p = p, arm = nm, reps = nrow(A),
        mean_eta = mean(A$eta), bias = mean(A$eta) - target,
        mean_lcb = mean(A$lcb), power = mean(A$lcb > 0),
        coverage = mean(A$lcb <= target + 1e-12),
        deg_mass = mean(A$deg))
    }
    cat(sprintf("  n = %5d  p = %2d  done\n", n, p))
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  out <- run_design_B(n_reps = as.integer(Sys.getenv("TM_REPS", "400")))
  dir.create("results", showWarnings = FALSE)
  write.csv(out, "results/design_B.csv", row.names = FALSE)
  cat("\n")
  print(out[, c("n","p","arm","mean_eta","bias","mean_lcb","power","coverage","deg_mass")],
        row.names = FALSE, digits = 4)
}
