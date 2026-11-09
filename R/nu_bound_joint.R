#' @title Lower bound on the affected SHARE, on KR's scale
#' @description Lower bound on
#'   \deqn{\bar\nu_k = \frac{\sum_x w(x)\,\theta^{(x)}_{kk}\nu^{(x)}_k}
#'                          {\sum_x w(x)\,\theta^{(x)}_{kk}},}
#'   the always-taker-weighted share of \eqn{k}-always-takers affected by the
#'   treatment. This is KR's \eqn{\bar\nu} aggregated over covariate cells.
#'
#'   \strong{Why this is not a division.} The numerator and the denominator are
#'   both bounded below over the identified set, and \eqn{a \le c} with
#'   \eqn{b \le d} says nothing about \eqn{a/b} versus \eqn{c/d}. Dividing a
#'   lower-bounded \eqn{\eta_k} by a lower-bounded
#'   \eqn{\tilde\theta^{\min}_{kk}} is therefore \emph{anti-conservative}: on
#'   Baranov's five-valued mediator it returns 0.0704 at \eqn{k=4} where the
#'   correct bound is 0.0495, and 0.2807 against 0.2670 at \eqn{k=5}. The two
#'   agree only when the denominator is pinned down, which happens at the
#'   bottom of the support of M -- monotonicity forces
#'   \eqn{\theta_{00} = P(M = m_0 \mid D = 1)} exactly -- and hence in every
#'   binary-mediator calculation at \eqn{k = 0}.
#'
#'   KR solve the pooled case as a linear-fractional program via the
#'   Charnes--Cooper transformation. The same transformation applies here with
#'   the cells carried along: write \eqn{\tilde\nu^{(x)}_k = \theta^{(x)}_{kk}
#'   \nu^{(x)}_k}, scale every variable by \eqn{z = 1/\sum_x w(x)
#'   \theta^{(x)}_{kk}}, and the program becomes linear in
#'   \eqn{(z\theta, z\tilde\nu, z)}.
#'
#' @param S Vector of \eqn{S(x) = \sup_A \Delta_k(A \mid x)}, one per cell.
#' @param A,B Matrices \eqn{J \times K} of control and treated mediator
#'   marginals, by row.
#' @param w Cell weights summing to one.
#' @param k Index of the always-taker group.
#' @param dbar Population defier cap.
#' @return A list with \code{value} (the lower bound on \eqn{\bar\nu_k}),
#'   \code{status}, \code{theta_kk} (the implied weighted always-taker share)
#'   and \code{numerator}.
#' @export
nu_bound_joint <- function(S, A, B, w, k, dbar = 0) {
  J <- length(S); K <- ncol(A)
  stopifnot(nrow(A) == J, nrow(B) == J, length(w) == J)

  # variables: Theta^x_{lj} (J*K*K), N^x (J), z (1)
  nth <- J * K * K
  iz  <- nth + J + 1L
  nv  <- iz
  th  <- function(x, l, j) (x - 1L) * K * K + (j - 1L) * K + l
  iN  <- function(x) nth + x

  obj <- numeric(nv)
  obj[nth + seq_len(J)] <- w                      # min sum_x w(x) N^x

  rows <- list(); dirs <- character(0); rhs <- numeric(0)
  add <- function(v, d, r) { rows[[length(rows) + 1L]] <<- v
                             dirs <<- c(dirs, d); rhs <<- c(rhs, r) }

  # normalisation: sum_x w(x) Theta^x_kk = 1
  v <- numeric(nv); for (x in 1:J) v[th(x, k, k)] <- w[x]
  add(v, "==", 1)

  for (x in 1:J) {
    for (l in 1:K) {                              # row sums = z * a^x_l
      v <- numeric(nv)
      for (j in 1:K) v[th(x, l, j)] <- 1
      v[iz] <- -A[x, l]
      add(v, "==", 0)
    }
    for (j in 1:K) {                              # col sums = z * b^x_j
      v <- numeric(nv)
      for (l in 1:K) v[th(x, l, j)] <- 1
      v[iz] <- -B[x, j]
      add(v, "==", 0)
    }
    # N^x >= z*S(x) - (z*b^x_k - Theta^x_kk)   [KR eq. (5), scaled]
    v <- numeric(nv)
    v[iN(x)] <- 1; v[th(x, k, k)] <- -1; v[iz] <- B[x, k] - S[x]
    add(v, ">=", 0)
    # N^x <= Theta^x_kk, since nu <= 1
    v <- numeric(nv); v[iN(x)] <- 1; v[th(x, k, k)] <- -1
    add(v, "<=", 0)
  }
  # population defier cap, scaled
  v <- numeric(nv)
  for (x in 1:J) for (l in 1:K) for (j in 1:K) if (l > j) v[th(x, l, j)] <- w[x]
  v[iz] <- -dbar
  add(v, "<=", 0)

  M <- do.call(rbind, rows)
  r <- Rglpk::Rglpk_solve_LP(obj, M, dirs, rhs, max = FALSE)
  if (r$status != 0) return(list(value = NA_real_, status = r$status))
  z <- r$solution[iz]
  list(value = r$optimum, status = r$status,
       theta_kk = if (z > 0) 1 / z else NA_real_,
       numerator = if (z > 0) r$optimum / z else NA_real_)
}
