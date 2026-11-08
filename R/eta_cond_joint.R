#' @title Conditional bound with a population-level defier cap
#' @description Minimises the conditional bound jointly over all covariate
#'   cells, subject to one defier cap on the population rather than one cap per
#'   cell:
#'   \deqn{\min_{\{\theta^{(x)}\}} \sum_x w(x)\Big(S(x) - b^{(x)}_k +
#'         \theta^{(x)}_{kk}\Big)_+}
#'   subject to each cell's mediator marginals and
#'   \eqn{\sum_x w(x) \sum_{l > j} \theta^{(x)}_{lj} \le \bar d}.
#'
#'   The positive part is handled with one auxiliary variable per cell, so this
#'   is a linear program with \eqn{J K^2 + J} variables and \eqn{2KJ + J + 1}
#'   constraints. The outcome enters only through the constants \eqn{S(x)}, so
#'   the program's size is independent of the support of Y -- the separation
#'   that makes the whole approach affordable survives the relaxation.
#'
#'   \strong{Why this imposition and not one cap per cell.} Dominance,
#'   \eqn{E[\eta_k(X)] \ge \eta_k^{\mathrm{pooled}}}, holds here by
#'   construction. Given any jointly feasible \eqn{\{\theta^{(x)}\}}, the
#'   aggregate \eqn{\bar\theta = \sum_x w(x)\theta^{(x)}} has the pooled
#'   marginals and defier mass at most \eqn{\bar d}, so it is feasible for the
#'   pooled program; convexity of \eqn{(\cdot)_+} together with
#'   \eqn{S_{\mathrm{pool}} \le \sum_x w(x) S(x)} makes the pooled objective at
#'   \eqn{\bar\theta} no larger than the conditional objective at
#'   \eqn{\{\theta^{(x)}\}}. Minimising the right-hand side gives the claim.
#'
#'   With a separate cap per cell the argument fails, because a cell whose own
#'   marginals need a larger defier share than the pooled marginals do faces a
#'   strictly larger identified set than the pooled problem. On Baranov that is
#'   not hypothetical: splitting by \code{girl} forces one cell to 0.074 against
#'   a pooled 0.0154, and the measured gain at \eqn{k=5} is \eqn{-0.032}.
#'
#' @param S Vector of \eqn{S(x) = \sup_A \Delta_k(A \mid x)}, one per cell.
#' @param A Matrix \eqn{J \times K} of control mediator marginals, by row.
#' @param B Matrix \eqn{J \times K} of treated mediator marginals, by row.
#' @param w Cell weights, summing to one.
#' @param k Index of the always-taker group.
#' @param dbar Population defier cap.
#' @return A list with \code{value}, \code{status}, \code{z} (the per-cell
#'   contributions) and \code{defier_by_cell}.
#' @export
eta_cond_joint <- function(S, A, B, w, k, dbar = 0) {
  J <- length(S); K <- ncol(A)
  stopifnot(nrow(A) == J, nrow(B) == J, length(w) == J, abs(sum(w) - 1) < 1e-8)

  nth <- J * K * K
  nv  <- nth + J                                  # theta's, then z's
  th  <- function(x, l, j) (x - 1L) * K * K + (j - 1L) * K + l
  zi  <- function(x) nth + x

  obj <- numeric(nv); obj[nth + seq_len(J)] <- w

  rows <- list(); dirs <- character(0); rhs <- numeric(0)
  add <- function(v, d, r) { rows[[length(rows) + 1L]] <<- v
                             dirs <<- c(dirs, d); rhs <<- c(rhs, r) }

  for (x in 1:J) {
    for (l in 1:K) { v <- numeric(nv)
      for (j in 1:K) v[th(x, l, j)] <- 1
      add(v, "==", A[x, l]) }
    for (j in 1:K) { v <- numeric(nv)
      for (l in 1:K) v[th(x, l, j)] <- 1
      add(v, "==", B[x, j]) }
    # z_x - theta^(x)_kk >= S(x) - b_k^(x)
    v <- numeric(nv); v[zi(x)] <- 1; v[th(x, k, k)] <- -1
    add(v, ">=", S[x] - B[x, k])
  }
  # population defier cap
  v <- numeric(nv)
  for (x in 1:J) for (l in 1:K) for (j in 1:K) if (l > j) v[th(x, l, j)] <- w[x]
  add(v, "<=", dbar)

  M <- do.call(rbind, rows)
  r <- Rglpk::Rglpk_solve_LP(obj, M, dirs, rhs, max = FALSE)
  if (r$status != 0) return(list(value = NA_real_, status = r$status))

  sol <- r$solution
  dfc <- vapply(1:J, function(x) {
    s <- 0
    for (l in 1:K) for (j in 1:K) if (l > j) s <- s + sol[th(x, l, j)]
    s }, numeric(1))
  list(value = r$optimum, status = r$status,
       z = sol[nth + seq_len(J)], defier_by_cell = dfc)
}
