#' @title Minimum always-taker share allowing a bounded share of defiers
#' @description Solves
#'   \deqn{\min \theta_{kk} \quad \text{s.t.} \quad \sum_j \theta_{lj} = a_l,\
#'         \sum_l \theta_{lj} = b_j,\ \sum_{l > j} \theta_{lj} \le \bar d,\
#'         \theta \ge 0,}
#'   where \eqn{\theta_{lj} = P(M(0) = m_l, M(1) = m_j)}, \eqn{a} and \eqn{b}
#'   are the mediator marginals under control and treatment, and the third
#'   constraint caps the defier mass (types with \eqn{M(1) < M(0)}, the lower
#'   triangle). At \eqn{\bar d = 0} this is KR's monotone identified set and the
#'   value is their eq. (8); for \eqn{\bar d > 0} it is the relaxation their
#'   \code{max_defiers_share} argument implements for the pooled case.
#'
#'   KR's closed form does not survive the relaxation, so the value comes from
#'   an LP. The LP is of size \eqn{K^2} and does not involve the outcome at all,
#'   which is the separation that makes this affordable: its cost is independent
#'   of the size of the support of Y.
#'
#'   The dual is returned as well, and that is the point. Weak duality gives
#'   \eqn{a'u + b'v - \bar d w \le \tilde\theta_{kk}^{\min}} for every
#'   dual-feasible \eqn{(u, v, w)}, and \eqn{\tilde\theta_{kk}^{\min}} enters the
#'   bound with a positive sign, so substituting any dual-feasible triple can
#'   only lower the bound -- conservative, never invalid. Since the dual
#'   objective is linear in the marginals, a triple fixed on one fold turns the
#'   second stage back into a linear functional of the data, exactly as fixing
#'   the set \eqn{A} does.
#'
#' @param a,b Mediator marginals \eqn{P(M = m_l \mid D = 0)} and
#'   \eqn{P(M = m_j \mid D = 1)}, of equal length K.
#' @param k Index of the always-taker group, in 1..K.
#' @param dbar Maximum defier share. Default 0 (strict monotonicity).
#' @param dual If TRUE, also solve the dual and return \eqn{(u, v, w)}.
#' @return A list with \code{value}, \code{status}, and when \code{dual} is TRUE
#'   \code{u}, \code{v}, \code{w} and \code{dual_value}.
#' @export
theta_kk_min_defiers <- function(a, b, k, dbar = 0, dual = TRUE) {
  K <- length(a)
  stopifnot(length(b) == K, k >= 1, k <= K, dbar >= 0)

  idx <- function(l, j) (j - 1L) * K + l          # column-major, theta[l, j]
  nv  <- K * K
  cvec <- numeric(nv); cvec[idx(k, k)] <- 1

  # row sums = a, column sums = b, defier mass <= dbar
  A <- matrix(0, 2 * K + 1L, nv)
  for (l in 1:K) for (j in 1:K) A[l, idx(l, j)] <- 1
  for (j in 1:K) for (l in 1:K) A[K + j, idx(l, j)] <- 1
  for (l in 1:K) for (j in 1:K) if (l > j) A[2 * K + 1L, idx(l, j)] <- 1
  rhs <- c(a, b, dbar)
  dir <- c(rep("==", 2 * K), "<=")

  pr <- Rglpk::Rglpk_solve_LP(cvec, A, dir, rhs, max = FALSE)
  out <- list(value = pr$optimum, status = pr$status)
  if (!dual) return(out)

  # Dual:  max a'u + b'v - dbar*w
  #        s.t. u_l + v_j - w*1{l > j} <= c_lj,   w >= 0,  u, v free.
  nd <- 2 * K + 1L
  dobj <- c(a, b, -dbar)
  Ad <- matrix(0, nv, nd)
  rd <- numeric(nv)
  r <- 0L
  for (l in 1:K) for (j in 1:K) {
    r <- r + 1L
    Ad[r, l] <- 1; Ad[r, K + j] <- 1
    if (l > j) Ad[r, nd] <- -1
    rd[r] <- if (l == k && j == k) 1 else 0
  }
  bnds <- list(lower = list(ind = 1:(2 * K), val = rep(-Inf, 2 * K)),
               upper = list(ind = 1:(2 * K), val = rep(Inf, 2 * K)))
  du <- Rglpk::Rglpk_solve_LP(dobj, Ad, rep("<=", nv), rd, bounds = bnds,
                              max = TRUE)
  out$u <- du$solution[1:K]
  out$v <- du$solution[(K + 1):(2 * K)]
  out$w <- du$solution[nd]
  out$dual_value <- du$optimum
  out$dual_status <- du$status
  out
}


#' @title KR's closed form for the minimum always-taker share
#' @description \eqn{\max\{b_k - (\sum_{j \ge k} b_j - \sum_{l \ge k} a_l), 0\}},
#'   valid only under strict monotonicity. Kept for checking the LP at
#'   \eqn{\bar d = 0}.
#' @param a,b,k As in \code{theta_kk_min_defiers}.
#' @export
theta_kk_min_closed_form <- function(a, b, k) {
  K <- length(a)
  max(b[k] - (sum(b[k:K]) - sum(a[k:K])), 0)
}


#' @title Smallest defier share compatible with a pair of mediator marginals
#' @description The largest shortfall in the survivor function,
#'   \eqn{\max_m \{P(M \ge m \mid D=0) - P(M \ge m \mid D=1)\}_+}. Zero iff
#'   monotonicity is compatible with the marginals.
#' @param a,b As above.
#' @export
min_defier_share <- function(a, b) {
  K <- length(a)
  s0 <- rev(cumsum(rev(a))); s1 <- rev(cumsum(rev(b)))
  max(c(0, s0 - s1))
}
