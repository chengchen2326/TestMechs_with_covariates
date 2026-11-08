#' @title Covariate-assisted lower bound with continuous or many covariates
#' @description Same estimand as \code{lb_frac_affected_x}, but stage 1 uses a
#'   regression rather than cell means, so the covariates may be continuous or
#'   numerous. This is the case in which binning the covariates -- the route
#'   Kwon & Roth suggest for discrete X -- is not available at all.
#'
#'   Stage 1 fits \eqn{P(Y = v, M = m_k \mid D = d, X = x)} for each outcome
#'   value \eqn{v} and each arm on the first fold, and sets
#'   \eqn{\hat A_k(x) = \{v : \hat p_1(v \mid x) > \hat p_0(v \mid x)\}}.
#'   Nothing about this fit has to be correct: any \eqn{A} is dual-feasible, so
#'   a bad first stage costs tightness and not validity. Stage 1 is a
#'   classification problem, not an inference problem.
#'
#'   Stage 2 forms \eqn{g_k = 1\{Y \in \hat A_k(X), M = m_k\} - 1\{M \ge m_k\}}
#'   on the held-out fold and estimates its average treatment effect by AIPW.
#'   With \eqn{\pi} known -- the experimental case -- the AIPW estimator is
#'   unbiased for \emph{any} outcome regression, so the regression only buys
#'   variance reduction and cannot break validity.
#'
#'   Only \eqn{\eta_k} (the population-scale bound) is returned, not the ratio
#'   \eqn{\nu_k}. The denominator \eqn{E[\tilde\theta_{kk}^{\min}(X)]} involves a
#'   truncation inside the expectation, which is not a smooth functional of a
#'   regression; how to normalise in the continuous case is still open.
#'
#' @param df A data frame.
#' @param d,m,y Names of treatment, mediator and outcome.
#' @param x Character vector of covariate names. May be continuous.
#' @param at_group Value of M whose always-takers are of interest.
#' @param pi_true If the propensity is known by design, its value (a scalar or a
#'   vector of length \code{nrow(df)}). Strongly preferred in an experiment. If
#'   \code{NULL}, it is estimated by logistic regression on fold 1 and the
#'   estimator is only doubly robust rather than unconditionally valid.
#' @param outcome_reg If \code{TRUE} (default) the AIPW augmentation term is
#'   fitted; if \code{FALSE} the estimator is plain Horvitz--Thompson.
#' @param alpha One-sided level. Default 0.05.
#' @param split_frac Share of the sample for stage 1. Default 0.3; see
#'   \code{lb_frac_affected_x}.
#' @param folds Optional integer vector of 1s and 2s.
#' @param seed Optional seed for the split.
#' @param plugin If \code{TRUE}, stage 1 and stage 2 both use the whole sample.
#'   Upward biased; provided for comparison only, and no confidence statement
#'   is returned.
#' @return A list with \code{eta}, \code{eta_se}, \code{eta_lcb}, \code{n},
#'   \code{n_stage2}, \code{pi_known}, \code{method}.
#' @export
lb_frac_affected_x_cts <- function(df, d, m, y, x,
                                   at_group = NULL,
                                   pi_true = NULL,
                                   outcome_reg = TRUE,
                                   alpha = 0.05,
                                   split_frac = 0.3,
                                   folds = NULL,
                                   seed = NULL,
                                   plugin = FALSE) {

  if (is.null(at_group) || length(at_group) != 1L)
    stop("at_group must be a single value of m")
  if (is.null(x) || length(x) == 0L) stop("x is required")

  dd <- df[stats::complete.cases(df[, c(d, m, y, x), drop = FALSE]), , drop = FALSE]
  n <- nrow(dd)
  if (n == 0L) stop("no complete cases")

  D <- as.numeric(dd[[d]]); M <- dd[[m]]; Y <- dd[[y]]
  if (!all(D %in% c(0, 1))) stop("d must be binary 0/1")
  X <- dd[, x, drop = FALSE]
  k <- at_group
  yv <- sort(unique(Y))

  if (plugin) {
    f1 <- rep(TRUE, n); f2 <- rep(TRUE, n)
  } else {
    if (!is.null(folds)) {
      fv <- folds[stats::complete.cases(df[, c(d, m, y, x), drop = FALSE])]
    } else {
      if (!is.null(seed)) set.seed(seed)
      fv <- rep(2L, n); fv[sample.int(n, floor(split_frac * n))] <- 1L
    }
    f1 <- fv == 1L; f2 <- fv == 2L
  }

  fml <- stats::as.formula(paste("resp ~", paste(x, collapse = " + ")))

  # Logistic fit of P(resp = 1 | X) on a subset; falls back to the subset mean
  # when the fit is degenerate (no variation, separation, too few points). A
  # fallback costs tightness only.
  fit_prob <- function(resp, use) {
    train <- data.frame(resp = resp[use], X[use, , drop = FALSE])
    if (nrow(train) < 5L || length(unique(train$resp)) < 2L)
      return(rep(if (nrow(train) > 0) mean(train$resp) else 0, n))
    fo <- tryCatch(suppressWarnings(
            stats::glm(fml, data = train, family = stats::binomial())),
          error = function(e) NULL)
    if (is.null(fo)) return(rep(mean(train$resp), n))
    p <- tryCatch(as.numeric(stats::predict(fo, newdata = X, type = "response")),
                  error = function(e) rep(mean(train$resp), n))
    p[!is.finite(p)] <- mean(train$resp)
    p
  }

  # ------------------------------------------- stage 1: Ahat_k(x) by regression
  inA <- matrix(FALSE, nrow = n, ncol = length(yv))
  for (jj in seq_along(yv)) {
    resp <- as.numeric(Y == yv[jj] & M == k)
    p1 <- fit_prob(resp, f1 & D == 1)
    p0 <- fit_prob(resp, f1 & D == 0)
    inA[, jj] <- p1 > p0
  }
  inAi <- logical(n)
  for (jj in seq_along(yv)) inAi[Y == yv[jj]] <- inA[Y == yv[jj], jj]

  g <- as.numeric(inAi & M == k) - as.numeric(M >= k)

  # ------------------------------------------------------ stage 2: AIPW on g
  if (!is.null(pi_true)) {
    pi_hat <- if (length(pi_true) == 1L) rep(pi_true, n) else
              pi_true[stats::complete.cases(df[, c(d, m, y, x), drop = FALSE])]
    pi_known <- TRUE
  } else {
    pi_hat <- fit_prob(D, f1)
    pi_known <- FALSE
  }
  pi_hat <- pmin(pmax(pi_hat, 0.01), 0.99)

  if (outcome_reg) {
    # g is bounded in [-1, 0] for k at the bottom of the support and in
    # [-1, 1] generally, so a linear fit is adequate as an augmentation term.
    reg <- function(use) {
      tr <- data.frame(resp = g[use], X[use, , drop = FALSE])
      if (nrow(tr) < length(x) + 2L) return(rep(if (nrow(tr)) mean(tr$resp) else 0, n))
      fo <- tryCatch(stats::lm(fml, data = tr), error = function(e) NULL)
      if (is.null(fo)) return(rep(mean(tr$resp), n))
      p <- as.numeric(stats::predict(fo, newdata = X))
      p[!is.finite(p)] <- mean(tr$resp); p
    }
    mu1 <- reg(f1 & D == 1); mu0 <- reg(f1 & D == 0)
  } else {
    mu1 <- rep(0, n); mu0 <- rep(0, n)
  }

  idx <- which(f2); n2 <- length(idx)
  if (n2 == 0L) stop("stage-2 fold is empty")

  psi <- (D[idx] / pi_hat[idx]) * (g[idx] - mu1[idx]) -
         ((1 - D[idx]) / (1 - pi_hat[idx])) * (g[idx] - mu0[idx]) +
         mu1[idx] - mu0[idx]

  eta <- mean(psi)
  if (plugin) {
    return(list(eta = eta, eta_se = NA_real_, eta_lcb = NA_real_,
                n = n, n_stage2 = n2, pi_known = pi_known, method = "plugin_cts"))
  }
  eta_se  <- stats::sd(psi) / sqrt(n2)
  eta_lcb <- max(eta - stats::qnorm(1 - alpha) * eta_se, 0)

  list(eta = eta, eta_se = eta_se, eta_lcb = eta_lcb,
       n = n, n_stage2 = n2, pi_known = pi_known, method = "split_cts")
}
