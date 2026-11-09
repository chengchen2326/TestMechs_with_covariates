# C2 (Jon's question): with a fixed basis b(y), is stage 1 an LP in beta?
#
# 1. sup over SETS A equals sup over FUNCTIONS a: Y -> [0,1].  The objective is
#    linear in a and the extreme points of [0,1]^Y are the indicators, so the
#    relaxation is exact, not an approximation.
# 2. With a partition-of-unity basis the constraint 0 <= beta'b(y) <= 1 reduces
#    to beta in [0,1]^J, a box, so the LP has the closed form
#    beta_j = 1{coefficient_j > 0}.
# 3. Indicator bins are such a basis, and the closed form is then exactly KR's
#    discretisation. So "choice of bins" is a special case of "choice of basis".
suppressMessages(library(Rglpk))
set.seed(5)
cat(sprintf("%6s %14s %14s %14s %12s\n","Q","sup over sets","LP over [0,1]","basis closed form","max|diff|"))
worst <- 0
for (Q in c(2,5,10,25)) {
  for (rep in 1:200) {
    d <- rnorm(Q) * 0.05                       # d(y) = p1(y,m_k) - p0(y,m_k)
    # (1) brute force over all subsets, for small Q; else the known closed form
    if (Q <= 10) {
      best <- 0
      for (mask in 0:(2^Q - 1)) {
        sel <- as.logical(bitwAnd(mask, 2^(0:(Q-1))))
        best <- max(best, sum(d[sel]))
      }
    } else best <- sum(pmax(d,0))
    # (2) LP over a in [0,1]^Q
    lp <- Rglpk_solve_LP(d, diag(Q), rep("<=", Q), rep(1, Q), max = TRUE)
    # (3) partition-of-unity basis: contiguous bins of width 1 -> beta in [0,1]
    cf <- sum(d[d > 0])
    worst <- max(worst, abs(best - lp$optimum), abs(best - cf))
  }
  cat(sprintf("%6d %14.8f %14.8f %14.8f %12.2e\n", Q, best, lp$optimum, cf, worst))
}
cat(sprintf("\nmax discrepancy across all %d draws: %.2e\n", 4*200, worst))

cat("\n--- with a COARSER basis (bins of width > 1) the closed form still holds,\n")
cat("    and the loss relative to the finest basis is exactly what KR's\n")
cat("    discretisation gives up:\n")
set.seed(9); Q <- 24
cat(sprintf("%8s %14s %14s\n","bin width","sup attained","loss vs finest"))
d <- rnorm(Q) * 0.05
fine <- sum(pmax(d,0))
for (wdt in c(1,2,3,4,6,12,24)) {
  g <- rep(seq_len(Q/wdt), each = wdt)
  agg <- tapply(d, g, sum)
  cat(sprintf("%8d %14.6f %14.6f\n", wdt, sum(pmax(agg,0)), fine - sum(pmax(agg,0))))
}
