#' @title Monotone projection of a pair of mediator marginals
#' @description Adjusts \code{pm1} so that \eqn{M \mid D = 1} stochastically
#'   dominates \eqn{M \mid D = 0}, which is KR's Assumption 2. Dominance is
#'   \eqn{F_1 \le F_0} pointwise, so the projection replaces \eqn{F_1} by
#'   \eqn{\min(F_1, F_0)}: still non-decreasing, still ends at one, and equal to
#'   \eqn{F_1} wherever dominance already held. Needed because Baranov's
#'   five-valued mediator violates dominance in sample (minimum defier share
#'   0.0154), so its marginals cannot be used as a monotone DGP unchanged.
#' @param pm0,pm1 Probability vectors over the ordered support of M.
#' @return A list with the projected \code{pm1} and the \code{shift}, the total
#'   mass moved.
#' @export
monotone_project <- function(pm0, pm1) {
  stopifnot(length(pm0) == length(pm1))
  c0 <- cumsum(pm0); c1 <- cumsum(pm1)
  c1n <- pmin(c1, c0)
  c1n[length(c1n)] <- 1
  pm1n <- diff(c(0, c1n))
  pm1n[pm1n < 0] <- 0
  pm1n <- pm1n / sum(pm1n)
  list(pm1 = pm1n, shift = sum(abs(pm1n - pm1)) / 2)
}


#' @title Comonotone coupling of two mediator marginals
#' @description The quantile coupling, built by the north-west corner rule.
#'   Under dominance it is supported on \eqn{l \le k}, so monotonicity holds by
#'   construction: comonotonicity means \eqn{M(d) = F_d^{-1}(U)}, and
#'   \eqn{F_1 \le F_0} gives \eqn{F_1^{-1} \ge F_0^{-1}}.
#' @param pm0,pm1 Probability vectors, with \code{pm1} dominating \code{pm0}.
#' @return A \code{K x K} matrix \code{theta} with \code{theta[l, k]} the mass of
#'   type \eqn{(m_l, m_k)}.
#' @export
comonotone_theta <- function(pm0, pm1) {
  K <- length(pm0)
  th <- matrix(0, K, K)
  i <- 1L; j <- 1L
  r0 <- pm0[1]; r1 <- pm1[1]
  guard <- 0L
  while (i <= K && j <= K && guard < 10L * K) {
    guard <- guard + 1L
    t <- min(r0, r1)
    th[i, j] <- th[i, j] + t
    r0 <- r0 - t; r1 <- r1 - t
    if (r0 <= 1e-12) { i <- i + 1L; if (i <= K) r0 <- pm0[i] }
    if (r1 <= 1e-12) { j <- j + 1L; if (j <= K) r1 <- pm1[j] }
  }
  th / sum(th)
}


#' @title Simulate data with a multivalued mediator and discrete covariates
#' @description Draws a covariate cell, then a compliance type
#'   \eqn{(M(0), M(1)) = (m_l, m_k)} from that cell's \code{theta}, then the
#'   treatment, then the outcome from the law attached to that type and arm.
#'   Because the type is drawn explicitly, the true \eqn{\theta_{kk}\nu_k} is
#'   known -- which the binary-mediator generator could not deliver, since it
#'   never pins the joint law of \eqn{(Y(0,m), Y(1,m))} for always-takers.
#'
#'   Monotonicity holds by construction: \code{theta} must be upper triangular.
#'
#' @param cells A list with one element per covariate cell, each a list with
#'   \code{theta} (\code{K x K}, upper triangular, sums to one), \code{J} (a
#'   list of \code{K} matrices of size \code{Q x Q}; \code{J[[k]]} is the joint
#'   law of \eqn{(Y(0,m_k), Y(1,m_k))} for always-takers at \eqn{m_k}), and
#'   \code{g0}, \code{g1} (each a \code{K x K} list-matrix of probability
#'   vectors over Y, used for complier types \eqn{l < k}: \code{g0[[l, k]]} is
#'   the law of \eqn{Y(0, m_l)} and \code{g1[[l, k]]} that of \eqn{Y(1, m_k)}).
#' @param p_x Cell probabilities; uniform by default.
#' @param n Sample size.
#' @param p_d P(D = 1).
#' @param seed Optional seed.
#' @param yvalues Support of Y; defaults to \code{1:Q}.
#' @return A data frame with \code{d}, \code{m}, \code{y}, \code{x}, and the
#'   latent \code{m0}, \code{m1} (the compliance type) for diagnostics.
#' @export
simulate_data_multiM_x <- function(cells, p_x = NULL, n, p_d = 0.5,
                                   seed = NULL, yvalues = NULL) {
  J0 <- length(cells)
  if (is.null(p_x)) p_x <- rep(1 / J0, J0)
  K <- nrow(cells[[1]]$theta)
  Q <- nrow(cells[[1]]$J[[1]])
  if (is.null(yvalues)) yvalues <- seq_len(Q)
  for (cj in cells)
    if (any(cj$theta[lower.tri(cj$theta)] > 1e-12))
      stop("theta must be upper triangular (monotonicity)")
  if (!is.null(seed)) set.seed(seed)

  xi <- sample.int(J0, n, replace = TRUE, prob = p_x)
  d  <- stats::rbinom(n, 1, p_d)
  m0 <- integer(n); m1 <- integer(n); y <- integer(n)

  for (jx in seq_len(J0)) {
    sel <- which(xi == jx)
    if (!length(sel)) next
    cj <- cells[[jx]]
    pairs <- which(cj$theta > 0, arr.ind = TRUE)
    pr <- cj$theta[pairs]
    pick <- sample.int(nrow(pairs), length(sel), replace = TRUE, prob = pr)
    m0[sel] <- pairs[pick, 1]; m1[sel] <- pairs[pick, 2]
    for (ii in seq_along(sel)) {
      i <- sel[ii]; l <- m0[i]; kk <- m1[i]
      if (l == kk) {
        # Draw the PAIR, reveal the arm's coordinate. This is what makes the
        # true nu_k known.
        cellp <- as.vector(cj$J[[l]])
        g <- sample.int(length(cellp), 1, prob = cellp)
        y0 <- ((g - 1L) %% Q) + 1L
        y1 <- ((g - 1L) %/% Q) + 1L
        y[i] <- if (d[i] == 1) y1 else y0
      } else {
        p <- if (d[i] == 1) cj$g1[[l, kk]] else cj$g0[[l, kk]]
        y[i] <- sample.int(Q, 1, prob = p)
      }
    }
  }
  data.frame(d = d, m = ifelse(d == 1, m1, m0), y = yvalues[y],
             x = xi, m0 = m0, m1 = m1)
}


#' @title Population bound and true affected mass for a multivalued design
#' @description Computes, exactly and without sampling, the observable
#'   \eqn{p_d(y, m \mid x)}, the conditional and pooled bounds \eqn{\eta_k}, and
#'   the \emph{true} \eqn{\theta_{kk}\nu_k}. Having the truth means coverage can
#'   be checked against the quantity of interest and not only against the bound,
#'   and it lets the slack of the bound be measured.
#' @param cells,p_x As in \code{simulate_data_multiM_x}.
#' @param at_group Index of the mediator value (1-based).
#' @export
population_bound_multiM_x <- function(cells, p_x = NULL, at_group) {
  J0 <- length(cells)
  if (is.null(p_x)) p_x <- rep(1 / J0, J0)
  K <- nrow(cells[[1]]$theta); Q <- nrow(cells[[1]]$J[[1]])
  k <- at_group

  obs <- function(cj) {
    p0 <- matrix(0, Q, K); p1 <- matrix(0, Q, K)
    for (mm in seq_len(K)) {
      p0[, mm] <- p0[, mm] + cj$theta[mm, mm] * rowSums(cj$J[[mm]])
      p1[, mm] <- p1[, mm] + cj$theta[mm, mm] * colSums(cj$J[[mm]])
      if (mm < K) for (kk in (mm + 1L):K)
        if (cj$theta[mm, kk] > 0) p0[, mm] <- p0[, mm] + cj$theta[mm, kk] * cj$g0[[mm, kk]]
      if (mm > 1L) for (l in 1L:(mm - 1L))
        if (cj$theta[l, mm] > 0) p1[, mm] <- p1[, mm] + cj$theta[l, mm] * cj$g1[[l, mm]]
    }
    list(p0 = p0, p1 = p1)
  }

  eta_of <- function(p0, p1) {
    sup_delta <- sum(pmax(p1[, k] - p0[, k], 0))
    above1 <- sum(p1[, k:K, drop = FALSE]); above0 <- sum(p0[, k:K, drop = FALSE])
    tkk <- max(sum(p1[, k]) - (above1 - above0), 0)
    list(eta = max(sup_delta - (sum(p1[, k]) - tkk), 0), theta = tkk)
  }

  O  <- lapply(cells, obs)
  pc <- lapply(O, function(o) eta_of(o$p0, o$p1))
  eta_cell   <- vapply(pc, function(z) z$eta,   numeric(1))
  theta_cell <- vapply(pc, function(z) z$theta, numeric(1))
  truth_cell <- vapply(cells, function(cj)
                   cj$theta[k, k] * (1 - sum(diag(cj$J[[k]]))), numeric(1))

  p0p <- Reduce(`+`, Map(function(w, o) w * o$p0, p_x, O))
  p1p <- Reduce(`+`, Map(function(w, o) w * o$p1, p_x, O))
  pooled <- eta_of(p0p, p1p)

  eta_cond <- sum(p_x * eta_cell)
  list(eta_by_cell = eta_cell, theta_by_cell = theta_cell,
       true_by_cell = truth_cell,
       eta_conditional = eta_cond, eta_pooled = pooled$eta,
       true_product = sum(p_x * truth_cell),
       slack_conditional = sum(p_x * truth_cell) - eta_cond,
       slack_pooled = sum(p_x * truth_cell) - pooled$eta)
}
