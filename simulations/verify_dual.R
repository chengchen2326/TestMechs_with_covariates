# =====================================================================
# verify_dual.R -- is Proposition 1 the Kantorovich dual in disguise?
#
# Write Z = (Y, M). The sharp lower bound on theta_kk * nu_k is the value of
# the transport problem
#
#     min_gamma E_gamma[c_k],   c_k(z0, z1) = 1{m0 = m1 = m_k, y0 != y1},
#
# over couplings gamma of (Z(0), Z(1)) with the two identified marginals
# P_0, P_1 and support restricted to m1 >= m0 (monotonicity, A2). Its
# Kantorovich dual is
#
#     max  int nu_0 dP_0 + int nu_1 dP_1
#     s.t. nu_0(z0) + nu_1(z1) <= c_k(z0, z1)  for all (z0, z1) with m1 >= m0.
#
# CLAIM. For every A, the pair
#
#     nu_1 = g_k := 1{y in A, m = m_k} - 1{m >= m_k},   nu_0 = -g_k
#
# is dual-FEASIBLE, and its dual objective equals psi_k(A). So Proposition 1(i)
# -- psi_k(A) <= theta_kk nu_k for every A -- is exactly weak duality, and
# "choosing A" is "choosing a point in a one-parameter subfamily of the dual
# feasible set". The set formulation and the transport formulation are the same
# object; the slides' set language is the dual language with the names dropped.
#
# Four things are checked on random monotone-feasible laws, over ALL 2^Q subsets:
#   (1) dual feasibility of (-g, g), by enumeration over every admissible pair;
#   (2) the dual objective equals psi_k(A) -- the quantity estimated as an ATE;
#   (3) weak duality, psi_k(A) <= primal LP value;
#   (4) no loss at A* = {y : p_1(y, m_k) > p_0(y, m_k)} on the interior branch.
# Plus a falsification check: drop A2 and dual feasibility must break, which
# confirms that monotonicity is what makes the dual constraint hold rather than
# being decoration.
#
# Expected output (seed below): zero feasibility violation, identity exact,
# weak duality to 1e-16, and roughly 87% of draws failing without A2.
# =====================================================================

suppressMessages(library(Rglpk))
set.seed(20261110)

states <- function(K, Q) expand.grid(y = seq_len(Q), m = seq_len(K))

draw_monotone <- function(K, Q) {
  # Draw a coupling supported on m1 >= m0 and return its marginals, so the
  # primal is feasible by construction rather than by luck.
  Z <- states(K, Q); nz <- nrow(Z)
  G <- matrix(0, nz, nz)
  for (i in seq_len(nz)) for (j in seq_len(nz))
    if (Z$m[j] >= Z$m[i]) G[i, j] <- rgamma(1, shape = 0.6)
  G <- G / sum(G)
  list(Z = Z, p0 = rowSums(G), p1 = colSums(G))
}

cost_mat <- function(Z, k) {
  nz <- nrow(Z)
  outer(seq_len(nz), seq_len(nz), Vectorize(function(i, j)
    as.numeric(Z$m[i] == k && Z$m[j] == k && Z$y[i] != Z$y[j])))
}

g_vec <- function(Z, A, k) {
  as.numeric(Z$y %in% A & Z$m == k) - as.numeric(Z$m >= k)
}

primal_value <- function(Z, p0, p1, k, monotone = TRUE) {
  nz <- nrow(Z); nv <- nz * nz
  idx <- function(i, j) (i - 1L) * nz + j          # row-major: gamma[i, j]
  obj <- as.vector(t(cost_mat(Z, k)))
  rows <- list(); rhs <- numeric(0)
  for (i in seq_len(nz)) {                          # row marginals = p0
    v <- numeric(nv); v[idx(i, seq_len(nz))] <- 1
    rows[[length(rows) + 1L]] <- v; rhs <- c(rhs, p0[i])
  }
  for (j in seq_len(nz)) {                          # column marginals = p1
    v <- numeric(nv); v[idx(seq_len(nz), j)] <- 1
    rows[[length(rows) + 1L]] <- v; rhs <- c(rhs, p1[j])
  }
  ub <- rep(1, nv)
  if (monotone) for (i in seq_len(nz)) for (j in seq_len(nz))
    if (Z$m[j] < Z$m[i]) ub[idx(i, j)] <- 0         # no mass where m1 < m0
  r <- Rglpk_solve_LP(obj, do.call(rbind, rows), rep("==", length(rhs)), rhs,
                      bounds = list(upper = list(ind = seq_len(nv), val = ub)),
                      max = FALSE)
  if (r$status != 0) NA_real_ else r$optimum
}

feasibility_slack <- function(Z, A, k, monotone = TRUE) {
  # max over admissible pairs of  nu_0(z0) + nu_1(z1) - c(z0, z1); must be <= 0
  g <- g_vec(Z, A, k); C <- cost_mat(Z, k); nz <- nrow(Z)
  worst <- -Inf
  for (i in seq_len(nz)) for (j in seq_len(nz)) {
    if (monotone && Z$m[j] < Z$m[i]) next
    worst <- max(worst, (-g[i]) + g[j] - C[i, j])
  }
  worst
}

all_subsets <- function(Q) {
  out <- list(integer(0))
  for (r in seq_len(Q)) out <- c(out, combn(Q, r, simplify = FALSE))
  out
}

# --------------------------------------------------------------- main sweep
grid <- list(c(2, 3), c(3, 3), c(3, 4), c(4, 3), c(5, 2))
worst_slack <- -Inf; worst_gap <- 0; worst_excess <- -Inf; gapA <- numeric(0)

for (KQ in grid) {
  K <- KQ[1]; Q <- KQ[2]
  for (rep in 1:40) {
    d <- draw_monotone(K, Q); Z <- d$Z
    for (k in seq_len(K)) {
      prim <- primal_value(Z, d$p0, d$p1, k)
      if (!is.finite(prim)) next
      for (A in all_subsets(Q)) {
        g <- g_vec(Z, A, k)
        psi <- sum(g * d$p1) - sum(g * d$p0)          # E[g|D=1] - E[g|D=0]
        dual_obj <- sum(-g * d$p0) + sum(g * d$p1)    # int nu0 dP0 + int nu1 dP1
        worst_slack  <- max(worst_slack,  feasibility_slack(Z, A, k))
        worst_gap    <- max(worst_gap,    abs(dual_obj - psi))
        worst_excess <- max(worst_excess, psi - prim)
      }
      # A* = where the treated partial density exceeds the control one
      sl <- Z$m == k
      Astar <- Z$y[sl][d$p1[sl] > d$p0[sl]]
      gA <- g_vec(Z, Astar, k)
      gapA <- c(gapA, prim - (sum(gA * d$p1) - sum(gA * d$p0)))
    }
  }
  cat(sprintf("K = %d, Q = %d done\n", K, Q))
}

cat("\n")
cat(sprintf("dual feasibility   : max violation over admissible pairs = %.3e  (<= 0 required)\n", worst_slack))
cat(sprintf("objective identity : max |dual objective - psi_k(A)|     = %.3e\n", worst_gap))
cat(sprintf("weak duality       : max (psi_k(A) - primal value)       = %.3e  (<= 0 required)\n", worst_excess))
cat(sprintf("gap at A*          : min %.3e  median %.3e  share exactly 0: %.3f\n",
            min(gapA), median(gapA), mean(gapA < 1e-9)))

# ------------------------------------------- falsification: A2 is load-bearing
bad <- 0; tot <- 0
for (rep in 1:30) {
  d <- draw_monotone(3, 3); Z <- d$Z
  for (k in 1:3) for (A in list(integer(0), 1L, 2L, c(1L, 2L), c(1L, 2L, 3L))) {
    tot <- tot + 1
    if (feasibility_slack(Z, A, k, monotone = FALSE) > 1e-12) bad <- bad + 1
  }
}
cat(sprintf("\nwithout monotonicity: dual feasibility violated in %d of %d checks (%.0f%%)\n",
            bad, tot, 100 * bad / tot))
