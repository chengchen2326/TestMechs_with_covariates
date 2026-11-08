#' @title Simulate data with a binary mediator and discrete covariates
#' @description Covariate-stratified version of \code{simulate_data_binaryM}.
#'   Draws a cell label X from \code{p_x}, then draws (D, M, Y) from the
#'   cell-specific partial densities. Each cell is generated exactly as in
#'   \code{simulate_data_binaryM}, so a single-cell call reproduces it.
#' @param cells A list with one element per covariate cell. Each element is a
#'   list with components \code{p_y_m1d1}, \code{p_y_m1d0}, \code{p_y_m0d1},
#'   \code{p_y_m0d0} (vectors over the support of Y), \code{p_m_1},
#'   \code{p_m_0} (scalars) and optionally \code{p_d} (scalar, default
#'   \code{p_d} below).
#' @param p_x Cell probabilities. Must be the same length as \code{cells}.
#'   Defaults to uniform.
#' @param n Total sample size.
#' @param p_d Default P(D = 1) within a cell, used for any cell that does not
#'   supply its own.
#' @param seed Optional seed.
#' @param yvalues Support of Y. Defaults to \code{1:length(p_y_m1d1)} of the
#'   first cell.
#' @return A data frame with columns \code{d}, \code{m}, \code{y}, \code{x}.
#' @export
simulate_data_binaryM_x <- function(cells,
                                    p_x = NULL,
                                    n,
                                    p_d = 0.5,
                                    seed = NULL,
                                    yvalues = NULL) {

  J <- length(cells)
  if (J < 1L) stop("cells must have at least one element")
  if (is.null(p_x)) p_x <- rep(1 / J, J)
  if (length(p_x) != J) stop("p_x must have the same length as cells")
  if (abs(sum(p_x) - 1) > 1e-8) stop("p_x must sum to one")
  if (is.null(yvalues)) yvalues <- seq_along(cells[[1]]$p_y_m1d1)

  required <- c("p_y_m1d1", "p_y_m1d0", "p_y_m0d1", "p_y_m0d0",
                "p_m_1", "p_m_0")
  for (j in seq_len(J)) {
    missing_fields <- setdiff(required, names(cells[[j]]))
    if (length(missing_fields) > 0)
      stop("cell ", j, " is missing: ", paste(missing_fields, collapse = ", "))
    if (cells[[j]]$p_m_1 < cells[[j]]$p_m_0)
      warning("cell ", j, " violates monotonicity in the M marginals: ",
              "p_m_1 = ", cells[[j]]$p_m_1, " < p_m_0 = ", cells[[j]]$p_m_0)
  }

  if (!is.null(seed)) set.seed(seed)

  # Multinomial cell assignment, then one simulate_data_binaryM draw per cell.
  n_by_cell <- as.vector(stats::rmultinom(1, size = n, prob = p_x))

  out <- vector("list", J)
  for (j in seq_len(J)) {
    nj <- n_by_cell[j]
    if (nj == 0L) next
    cj <- cells[[j]]
    pdj <- if (is.null(cj$p_d)) p_d else cj$p_d
    dfj <- simulate_data_binaryM(p_y_m1d1 = cj$p_y_m1d1,
                                 p_y_m1d0 = cj$p_y_m1d0,
                                 p_y_m0d1 = cj$p_y_m0d1,
                                 p_y_m0d0 = cj$p_y_m0d0,
                                 p_m_1    = cj$p_m_1,
                                 p_m_0    = cj$p_m_0,
                                 p_d      = pdj,
                                 n        = nj,
                                 seed     = NULL,
                                 yvalues  = yvalues)
    dfj$x <- j
    out[[j]] <- dfj
  }

  df <- do.call(rbind, out)
  rownames(df) <- NULL
  df[sample.int(nrow(df)), , drop = FALSE]
}


#' @title Population value of the bound in a simulated design
#' @description Computes eta_k and nu_k exactly from the partial densities used
#'   by \code{simulate_data_binaryM_x}, both cell by cell and pooled. Used to
#'   get the truth a Monte Carlo is measured against, with no sampling.
#' @param cells,p_x As in \code{simulate_data_binaryM_x}.
#' @param at_group Value of M (0 or 1) whose always-takers are of interest.
#' @return A list with \code{eta_by_cell}, \code{theta_by_cell},
#'   \code{eta_conditional} (= E[eta_k(X)]), \code{eta_pooled},
#'   \code{nu_conditional}, \code{nu_pooled}.
#' @export
population_bound_binaryM_x <- function(cells, p_x = NULL, at_group = 0) {

  J <- length(cells)
  if (is.null(p_x)) p_x <- rep(1 / J, J)
  k <- at_group
  if (!k %in% c(0, 1)) stop("at_group must be 0 or 1 for a binary mediator")

  # p_d[y, m] = P(Y = y, M = m | D = d, X = j) for each cell.
  joint <- function(cj, d) {
    pm1 <- if (d == 1) cj$p_m_1 else cj$p_m_0
    y_m1 <- if (d == 1) cj$p_y_m1d1 else cj$p_y_m1d0
    y_m0 <- if (d == 1) cj$p_y_m0d1 else cj$p_y_m0d0
    cbind(m0 = (1 - pm1) * y_m0 / sum(y_m0),
          m1 = pm1 * y_m1 / sum(y_m1))
  }

  # eta_k and theta_kk_min from a (p0, p1) pair, KR eq. (8) and Prop. 3.1.
  eta_of <- function(p0, p1) {
    col <- k + 1L
    sup_delta <- sum(pmax(p1[, col] - p0[, col], 0))
    above1 <- sum(p1[, col:2, drop = FALSE])
    above0 <- sum(p0[, col:2, drop = FALSE])
    tkk <- max(sum(p1[, col]) - (above1 - above0), 0)
    list(eta = max(sup_delta - (sum(p1[, col]) - tkk), 0), theta = tkk)
  }

  per <- lapply(cells, function(cj) eta_of(joint(cj, 0), joint(cj, 1)))
  eta_by_cell   <- vapply(per, function(z) z$eta,   numeric(1))
  theta_by_cell <- vapply(per, function(z) z$theta, numeric(1))

  # Pooled: mix the cell-level joints with the SAME weights p_x. Using arm-
  # specific weights P(X = x | D = d) breaks the mixture argument.
  p0_pool <- Reduce(`+`, Map(function(w, cj) w * joint(cj, 0), p_x, cells))
  p1_pool <- Reduce(`+`, Map(function(w, cj) w * joint(cj, 1), p_x, cells))
  pooled <- eta_of(p0_pool, p1_pool)

  eta_cond <- sum(p_x * eta_by_cell)
  theta_cond <- sum(p_x * theta_by_cell)

  list(eta_by_cell     = eta_by_cell,
       theta_by_cell   = theta_by_cell,
       eta_conditional = eta_cond,
       eta_pooled      = pooled$eta,
       nu_conditional  = if (theta_cond > 0) eta_cond / theta_cond else NA_real_,
       nu_pooled       = if (pooled$theta > 0) pooled$eta / pooled$theta else NA_real_)
}
