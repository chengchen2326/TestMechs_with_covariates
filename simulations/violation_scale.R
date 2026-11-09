# Scaling the violation of the sharp null, so every panel has a size row.
#
# KR index their designs by t in {0, 0.5, 1}, where t = 0 imposes the sharp
# null of full mediation. Our designs were each calibrated at one fixed
# violation, which left no clean size row: size could only be read off the
# pooled arm of the sign-changing calibration, where the pooled bound happens
# to be zero.
#
# The construction is the same in all three designs: interpolate each treated
# within-stratum outcome law toward the corresponding control law,
#
#     L1^t(. | m, x) = (1 - t) L0(. | m, x) + t L1(. | m, x).
#
# At t = 0 the two arms agree within every stratum, so with the maximal
# coupling nu_k = 0 and the sharp null holds exactly, not approximately. The
# mediator marginals are untouched, so monotonicity and the complier structure
# are unaffected and only the violation is scaled.
#
# That the bound is then exactly zero is worth checking rather than assuming.
# For the binary designs, eta_k = (sup_A Delta_k - C)_+ with
# C = sum_{j >= k} (b_j - a_j); at t = 0,
# sup_A Delta_k = (b_k - a_k)_+, and monotonicity gives
# sum_{j > k} (b_j - a_j) >= 0, hence C >= (b_k - a_k)_+ and eta_k = 0.

#' Interpolate the treated laws of a binary-mediator cell toward the control
#' laws. Used by Design A and anything else built on simulate_data_binaryM.
scale_violation_binary <- function(cells, t) {
  lapply(cells, function(cj) {
    cj$p_y_m0d1 <- (1 - t) * cj$p_y_m0d0 + t * cj$p_y_m0d1
    cj$p_y_m1d1 <- (1 - t) * cj$p_y_m1d0 + t * cj$p_y_m1d1
    cj
  })
}

#' Interpolate the treated laws of a multivalued-mediator cell toward the
#' control laws, rebuilding the always-taker couplings and the complier laws so
#' that the implied p_d stay consistent.
#' @param coupling "maximal" or "independent", as in calib_C.R.
scale_violation_multi <- function(cells, t, coupling = c("maximal", "independent")) {
  coupling <- match.arg(coupling)
  cpl <- if (coupling == "maximal") maximal_coupling else independent_coupling
  lapply(cells, function(cj) {
    K <- nrow(cj$theta)
    for (mm in 1:K) {
      L0 <- rowSums(cj$J[[mm]]); L0 <- L0 / sum(L0)
      L1 <- colSums(cj$J[[mm]]); L1 <- L1 / sum(L1)
      L1t <- (1 - t) * L0 + t * L1
      L1t <- pmax(L1t, 0); L1t <- L1t / sum(L1t)
      cj$J[[mm]] <- cpl(L0, L1t)
      for (l in 1:K) cj$g1[[l, mm]] <- L1t
    }
    cj
  })
}
