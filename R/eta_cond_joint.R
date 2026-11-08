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


#' @title Dual of the population-capped conditional bound
#' @description Solves
#'   \deqn{\max \sum_{x,l} a^{(x)}_l u^{(x)}_l + \sum_{x,j} b^{(x)}_j v^{(x)}_j
#'         + \sum_x p_x\big(S(x) - b^{(x)}_k\big) - \bar d q}
#'   subject to
#'   \eqn{u^{(x)}_l + v^{(x)}_j - p_x 1\{(l,j)=(k,k)\} - q w(x) 1\{l>j\} \le 0},
#'   \eqn{0 \le p_x \le w(x)} and \eqn{q \ge 0}.
#'
#'   Weak duality makes every feasible point a valid lower bound on the
#'   conditional bound, hence on the truth. The objective is linear in the cell
#'   marginals and in \eqn{S(x)}, so a point fixed on one fold turns the second
#'   stage into the sample mean of a bounded transformed outcome.
#'
#'   The constraint \eqn{p_x \le w(x)} shows the dual scales with the cell
#'   weights, so the returned variables are also given normalised by \eqn{w(x)}:
#'   in that scaling \eqn{\tilde p_x \in [0,1]} and the stage-2 estimand is an
#'   ordinary inverse-probability-weighted mean.
#' @param S,A,B,w,k,dbar As in \code{eta_cond_joint}.
#' @return A list with \code{value}, \code{status}, the raw \code{u}, \code{v}
#'   (both \eqn{J \times K}), \code{p} (length J), \code{q}, and the normalised
#'   \code{u_t}, \code{v_t}, \code{p_t}.
#' @export
eta_cond_joint_dual <- function(S, A, B, w, k, dbar = 0, dual_box = 2) {
  J <- length(S); K <- ncol(A)
  ui <- function(x, l) (x - 1L) * K + l
  vi <- function(x, j) J * K + (x - 1L) * K + j
  pi_ <- function(x) 2L * J * K + x
  qi <- 2L * J * K + J + 1L
  nv <- qi

  obj <- numeric(nv)
  for (x in 1:J) {
    for (l in 1:K) obj[ui(x, l)] <- A[x, l]
    for (j in 1:K) obj[vi(x, j)] <- B[x, j]
    obj[pi_(x)] <- S[x] - B[x, k]
  }
  obj[qi] <- -dbar

  rows <- list(); rhs <- numeric(0)
  for (x in 1:J) for (l in 1:K) for (j in 1:K) {
    v <- numeric(nv)
    v[ui(x, l)] <- 1; v[vi(x, j)] <- v[vi(x, j)] + 1
    if (l == k && j == k) v[pi_(x)] <- -1
    if (l > j) v[qi] <- -w[x]
    rows[[length(rows) + 1L]] <- v; rhs <- c(rhs, 0)
  }
  for (x in 1:J) {                                   # p_x <= w(x)
    v <- numeric(nv); v[pi_(x)] <- 1
    rows[[length(rows) + 1L]] <- v; rhs <- c(rhs, w[x])
  }
  M <- do.call(rbind, rows)

  # u and v are free in the textbook dual, and that is a practical disaster:
  # u_l = +100 with v_j = -100 is feasible, so an LP vertex can carry enormous
  # components, and the second-stage estimate -- which is linear in the
  # marginals with these as coefficients -- then has exploding variance and sits
  # systematically below zero. Boxing them costs a little tightness and buys
  # conditioning. The box is not arbitrary: at dbar = 0 the point corresponding
  # to KR's closed form has coefficients 1{l >= k} on a and
  # 1{j = k} - 1{j >= k} on b, all within [-1, 1], so a box of +-2 (in the
  # normalised scaling, hence +-2*w(x) raw) still contains it.
  free <- 1:(2L * J * K)
  wrep <- rep(rep(w, each = 1L), times = 1L)
  wcol <- c(rep(w, each = K), rep(w, each = K))      # w(x) for each u, v entry
  bnds <- list(lower = list(ind = free, val = -dual_box * wcol),
               upper = list(ind = free, val =  dual_box * wcol))
  r <- Rglpk::Rglpk_solve_LP(obj, M, rep("<=", nrow(M)), rhs, bounds = bnds,
                             max = TRUE)
  if (r$status != 0) return(list(value = NA_real_, status = r$status))

  s <- r$solution
  U <- matrix(0, J, K); V <- matrix(0, J, K); P <- numeric(J)
  for (x in 1:J) {
    for (l in 1:K) U[x, l] <- s[ui(x, l)]
    for (j in 1:K) V[x, j] <- s[vi(x, j)]
    P[x] <- s[pi_(x)]
  }
  list(value = r$optimum, status = r$status, u = U, v = V, p = P, q = s[qi],
       u_t = U / w, v_t = V / w, p_t = P / w)
}
