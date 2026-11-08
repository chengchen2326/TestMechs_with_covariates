# Baranov calibration for Design C: multivalued mediator.
#
# The DGP is built so that its implied p_d(y, m) equals Baranov's empirically,
# while the compliance type is drawn explicitly so the true theta_kk * nu_k is
# known.  Construction:
#
#   theta        = comonotone coupling of the (monotone-projected) M | D
#                  marginals, so monotonicity holds by construction.
#   Y laws       every type with M(0) = m_l shares the observed law of
#                  Y | M = m_l, D = 0, and every type with M(1) = m_k shares
#                  the observed law of Y | M = m_k, D = 1.  Summing over types
#                  then reproduces the observed p_d(., m) exactly.
#   J[[k]]       a coupling of those two laws for the always-takers at m_k.
#                  "maximal" puts min(p0, p1) on the diagonal, which minimises
#                  P(Y(1) != Y(0)) and makes nu_k equal the total variation
#                  distance -- the case where KR's bound is attained.
#                  "independent" makes nu_k larger, so the bound is slack.
#
# Baranov's five-valued mediator violates monotonicity in sample (minimum
# defier share 0.0154), so the projection is not cosmetic; the shift it makes
# is reported.

maximal_coupling <- function(p0, p1) {
  Q <- length(p0)
  d <- pmin(p0, p1)
  Jm <- diag(d, nrow = Q)
  r0 <- p0 - d; r1 <- p1 - d
  s <- sum(r0)
  if (s > 1e-12) Jm <- Jm + outer(r0 / s, r1 / s) * s
  Jm / sum(Jm)
}

independent_coupling <- function(p0, p1) outer(p0, p1)

#' Build Design C cells from Baranov.
#' @param K_target keep the top K_target mediator categories (merging the rest
#'   into the lowest), so K can be varied.
#' @param Q number of Y bins.
#' @param xvar covariate defining the cells, or NULL for a single cell.
#' @param coupling "maximal" (bound attained) or "independent" (bound slack).
baranov_cells <- function(K_target = 5, Q = 5, xvar = NULL,
                          coupling = c("maximal", "independent"),
                          mvar = "relationship_husb", yvar = "motherfinancial",
                          verbose = TRUE) {
  coupling <- match.arg(coupling)
  cpl <- if (coupling == "maximal") maximal_coupling else independent_coupling

  load("/home/claude/repo/data/baranov_data.rda")
  br <- haven::zap_labels(get("mother_data"))
  br <- br[!is.na(br$dep_sample) & br$dep_sample == 1, ]
  need <- c("treat", mvar, yvar, xvar)
  br <- br[stats::complete.cases(br[, need, drop = FALSE]), ]
  for (v in need) br[[v]] <- as.numeric(unclass(br[[v]]))

  # Mediator: collapse to K_target ordered categories.
  mv <- sort(unique(br[[mvar]]))
  if (length(mv) > K_target) {
    keep <- utils::tail(mv, K_target - 1L)
    br$mc <- ifelse(br[[mvar]] %in% keep, match(br[[mvar]], keep) + 1L, 1L)
  } else {
    br$mc <- match(br[[mvar]], mv)
  }
  K <- length(unique(br$mc))

  # Outcome: Q quantile bins.
  qs <- stats::quantile(br[[yvar]], probs = seq(0, 1, length.out = Q + 1), na.rm = TRUE)
  qs <- unique(qs); Qe <- length(qs) - 1L
  br$yc <- as.integer(cut(br[[yvar]], breaks = qs, include.lowest = TRUE))

  groups <- if (is.null(xvar)) list(seq_len(nrow(br))) else
            split(seq_len(nrow(br)), br[[xvar]])

  cells <- list(); p_x <- numeric(0)
  for (gi in seq_along(groups)) {
    s <- br[groups[[gi]], ]
    pm0 <- as.numeric(prop.table(table(factor(s$mc[s$treat == 0], levels = 1:K))))
    pm1 <- as.numeric(prop.table(table(factor(s$mc[s$treat == 1], levels = 1:K))))
    pr <- monotone_project(pm0, pm1)
    if (verbose) cat(sprintf("  cell %d: n = %4d, monotone projection shift = %.4f\n",
                             gi, nrow(s), pr$shift))
    th <- comonotone_theta(pm0, pr$pm1)

    ylaw <- function(dv, mm) {
      z <- s$yc[s$treat == dv & s$mc == mm]; z <- z[!is.na(z)]
      if (!length(z)) return(rep(1 / Qe, Qe))
      p <- as.numeric(prop.table(table(factor(z, levels = 1:Qe))))
      p[p < 0] <- 0; if (sum(p) == 0) rep(1 / Qe, Qe) else p / sum(p)
    }
    L0 <- lapply(1:K, function(mm) ylaw(0, mm))
    L1 <- lapply(1:K, function(mm) ylaw(1, mm))

    Jl <- lapply(1:K, function(mm) cpl(L0[[mm]], L1[[mm]]))
    g0 <- g1 <- matrix(list(), K, K)
    for (l in 1:K) for (kk in 1:K) {
      g0[[l, kk]] <- L0[[l]]; g1[[l, kk]] <- L1[[kk]]
    }
    cells[[length(cells) + 1]] <- list(theta = th, J = Jl, g0 = g0, g1 = g1)
    p_x <- c(p_x, nrow(s))
  }
  list(cells = cells, p_x = p_x / sum(p_x), K = K, Q = Qe)
}
