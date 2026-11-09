#' @title Covariate-assisted lower bound on the fraction of affected always-takers
#' @description Lower bound on the share of always-takers at \code{at_group} who
#'   are affected by the treatment despite no change in M, using discrete
#'   covariates.
#'
#'   The estimand is built from the identity that, for any fixed (possibly
#'   covariate-dependent) set \eqn{A(x)},
#'   \deqn{\psi_k(x;A) = \Delta_k(A(x)|x) - [P(M \ge m_k|D=1,x) - P(M \ge m_k|D=0,x)]}
#'   is a valid lower bound on \eqn{\eta_k(x)}, and equals Kwon & Roth's bound
#'   when \eqn{A = A^*_k(x)} and \eqn{\tilde\theta^{\min}_{kk}(x) > 0}. Writing
#'   \deqn{g_k(Y,M,X;A) = 1\{Y \in A(X), M = m_k\} - 1\{M \ge m_k\},}
#'   \eqn{E[\psi_k(X;A)]} is the average treatment effect of \eqn{g_k}. So once
#'   \eqn{A} is fixed, estimation is an ATE on a bounded transformed outcome and
#'   nothing more: no optimisation, no non-smooth functional, no bootstrap.
#'
#'   \code{method = "split"} picks \eqn{\hat A_k} on one fold and evaluates the
#'   ATE on the other, so the supremum is never estimated on the data used for
#'   inference. \code{method = "plugin"} does both on the full sample, which is
#'   the usual practice and is upward biased. \code{method = "pooled_split"} is
#'   \code{"split"} ignoring the covariates, which supplies a confidence
#'   statement for the pooled bound.
#'
#' @param df A data frame.
#' @param d Name of the treatment variable.
#' @param m Name of the mediator variable.
#' @param y Name of the outcome variable.
#' @param x Character vector of covariate names defining the cells. Ignored
#'   when \code{method = "pooled_split"}.
#' @param at_group Value of M whose always-takers are of interest. Required.
#' @param method One of \code{"split"}, \code{"plugin"}, \code{"pooled_split"}.
#' @param alpha One-sided level for the lower confidence bound. Default 0.05.
#' @param split_frac Share of the sample used for stage 1. Default 0.3. Stage 1
#'   only has to pick a set, and any set is dual-feasible, so a poor choice
#'   costs tightness rather than validity; stage 2 is where the central limit
#'   theorem is applied. That argues for giving stage 1 less than half, but not
#'   much less, since a very noisy \eqn{\hat A_k} shrinks the estimand itself.
#'   A sweep over 0.1--0.7 (see \code{simulations/split_sweep.R}) puts the
#'   optimum at 0.2--0.3: on Baranov, moving from 0.5 to 0.3 raises the share
#'   of splits with a positive 95\% bound from 0.85 to 0.97 and the median
#'   bound from 0.078 to 0.093.
#' @param folds Optional integer vector of 1s and 2s, length \code{nrow(df)},
#'   overriding the random split.
#' @param seed Optional seed for the split.
#' @return A list with \code{eta}, \code{eta_se}, \code{eta_lcb},
#'   \code{nu}, \code{nu_se}, \code{nu_lcb}, \code{theta}, \code{n},
#'   \code{n_stage2}, \code{n_cells}, \code{n_cells_degenerate},
#'   \code{mass_degenerate}, \code{n_cells_nonmonotone},
#'   \code{mass_nonmonotone}, \code{max_cell_defier_share},
#'   \code{defier_share_overall}, \code{method}.
#'
#'   Monotonicity is KR's Assumption 2 and is imposed, not tested. The returned
#'   defier shares are the smallest total defier mass that would reconcile the
#'   observed \eqn{M \mid D} marginals with \eqn{M(1) \ge M(0)}, overall and
#'   across cells -- the same descriptive quantity \code{lb_frac_affected}
#'   reports through \code{max_defiers_share}. They are reported rather than
#'   acted on: a sample shortfall need not be a population violation, and when
#'   the population sits on the boundary (no compliers) one appears in about
#'   half of all cells by chance, so zeroing those cells would be valid but
#'   would destroy power precisely in the designs of interest.
#' @export
lb_frac_affected_x <- function(df, d, m, y, x = NULL,
                               at_group = NULL,
                               method = c("split", "plugin", "pooled_split"),
                               alpha = 0.05,
                               split_frac = 0.3,
                               folds = NULL,
                               seed = NULL) {

  method <- match.arg(method)
  if (is.null(at_group) || length(at_group) != 1L)
    stop("at_group must be a single value of m")
  if (method != "pooled_split" && (is.null(x) || length(x) == 0L))
    stop("x is required unless method = 'pooled_split'")

  keep_cols <- c(d, m, y, if (method == "pooled_split") NULL else x)
  dd <- df[stats::complete.cases(df[, keep_cols, drop = FALSE]), , drop = FALSE]
  n <- nrow(dd)
  if (n == 0L) stop("no complete cases")

  D <- as.numeric(dd[[d]])
  M <- dd[[m]]
  Y <- dd[[y]]
  if (!all(D %in% c(0, 1))) stop("d must be binary 0/1")

  cell <- if (method == "pooled_split") {
    factor(rep("pooled", n))
  } else {
    interaction(dd[, x, drop = FALSE], drop = TRUE, sep = "|")
  }

  k  <- at_group
  yv <- sort(unique(Y))
  mvals_sorted <- sort(unique(M))

  # ---------------------------------------------------------------- folds
  if (method == "plugin") {
    f1 <- rep(TRUE, n)   # stage 1 and stage 2 both on everything
    f2 <- rep(TRUE, n)
  } else {
    if (!is.null(folds)) {
      if (length(folds) != nrow(df)) stop("folds must have length nrow(df)")
      fv <- folds[stats::complete.cases(df[, keep_cols, drop = FALSE])]
    } else {
      if (!is.null(seed)) set.seed(seed)
      fv <- rep(2L, n)
      fv[sample.int(n, floor(split_frac * n))] <- 1L
    }
    f1 <- fv == 1L
    f2 <- fv == 2L
  }

  # ------------------------------------------- stage 1: Ahat_k(x) on fold 1
  # Ahat_k(x) = { y : phat_1(y, m_k | x) > phat_0(y, m_k | x) }.
  # Any set is dual-feasible, so a cell with no stage-1 data simply gets the
  # empty set: that costs tightness, never validity.
  Ahat <- stats::setNames(vector("list", nlevels(cell)), levels(cell))
  for (cx in levels(cell)) {
    i1 <- f1 & cell == cx
    t1 <- i1 & D == 1
    t0 <- i1 & D == 0
    if (!any(t1) || !any(t0)) { Ahat[[cx]] <- logical(length(yv)); next }
    p1 <- vapply(yv, function(v) mean(Y[t1] == v & M[t1] == k), numeric(1))
    p0 <- vapply(yv, function(v) mean(Y[t0] == v & M[t0] == k), numeric(1))
    Ahat[[cx]] <- p1 > p0
  }

  # --------------------------------- stage 2: the transformed outcome g_k
  inA <- logical(n)
  for (cx in levels(cell)) {
    sel <- cell == cx
    inA[sel] <- Y[sel] %in% yv[Ahat[[cx]]]
  }
  g <- as.numeric(inA & M == k) - as.numeric(M >= k)

  # Pieces of theta_kk_min(x): a on the treated arm, b on the control arm,
  # so that theta_kk_min(x) = max( E[a|D=1,x] + E[b|D=0,x], 0 ).
  a <- as.numeric(M == k) - as.numeric(M >= k)
  b <- as.numeric(M >= k)

  # Monotonicity is not optional here. KR's eq. (8) closed form for
  # theta_kk_min is derived under M(1) >= M(0), i.e. stochastic dominance of
  # M | D = 1 over M | D = 0. When the sample violates it, the expression is no
  # longer the minimum over Theta_I and psi_k stops being a valid lower bound:
  # on the Baranov five-valued mediator it returns -0.222 at k = 2.
  # A violating cell therefore contributes zero, which is valid because
  # theta_kk(x) nu_k(x) >= 0, and is counted so the cost is visible.
  # Smallest total defier mass that would make the M | D marginals compatible
  # with M(1) >= M(0): the largest shortfall in the survivor function. Zero
  # means monotonicity holds in this subsample. Reported rather than acted on,
  # the way lb_frac_affected reports max_defiers_share, because a sample
  # shortfall need not be a population violation -- when the population sits on
  # the boundary (the no-complier case) one appears about half the time.
  min_defier <- function(j) {
    if (length(j) == 0L) return(NA_real_)
    j1 <- j[D[j] == 1]; j0 <- j[D[j] == 0]
    if (length(j1) == 0L || length(j0) == 0L) return(NA_real_)
    lv <- sort(unique(M[j]))
    s1 <- vapply(lv, function(v) mean(M[j1] >= v), numeric(1))
    s0 <- vapply(lv, function(v) mean(M[j0] >= v), numeric(1))
    max(c(0, s0 - s1))
  }

  idx <- which(f2)
  n2  <- length(idx)
  if (n2 == 0L) stop("stage-2 fold is empty")

  cl2 <- droplevels(cell[idx])
  lv  <- levels(cl2)

  tau <- stats::setNames(numeric(length(lv)), lv)   # cell contribution to eta
  thc <- stats::setNames(numeric(length(lv)), lv)   # cell theta_kk_min
  w   <- stats::setNames(numeric(length(lv)), lv)   # cell weight
  degenerate <- stats::setNames(logical(length(lv)), lv)
  defier     <- stats::setNames(rep(NA_real_, length(lv)), lv)
  mu1 <- mu0 <- pih <- a1 <- b0 <- tau            # per-cell nuisances

  for (cx in lv) {
    j  <- idx[cl2 == cx]
    j1 <- j[D[j] == 1]
    j0 <- j[D[j] == 0]
    w[cx] <- length(j) / n2
    if (length(j1) == 0L || length(j0) == 0L) {
      # No contrast available in this cell. Contributing zero is valid, since
      # theta_kk(x) nu_k(x) >= 0; it is conservative, and it is exactly the
      # price of thin cells that motivates not binning at all.
      degenerate[cx] <- TRUE
      tau[cx] <- 0; thc[cx] <- 0
      next
    }
    # A cell can fail stochastic dominance in sample while satisfying it in the
    # population -- and does so about half the time when the population sits on
    # the boundary, which is the no-complier case. Zeroing such a cell would be
    # valid but would destroy power exactly where the design is most
    # interesting, so the violation is recorded and reported, not acted on.
    defier[cx] <- min_defier(j)
    pih[cx] <- length(j1) / length(j)
    mu1[cx] <- mean(g[j1]); mu0[cx] <- mean(g[j0])
    tau[cx] <- mu1[cx] - mu0[cx]
    a1[cx]  <- mean(a[j1]); b0[cx]  <- mean(b[j0])
    thc[cx] <- max(a1[cx] + b0[cx], 0)
  }

  eta   <- sum(w * tau)
  theta <- sum(w * thc)

  defier_overall <- min_defier(idx)

  # ------------------------------------------------- influence functions
  # Saturated AIPW: within a cell the treated and control residuals each sum to
  # zero, and the cell terms are centred at the estimate, so mean(IF) = 0.
  ci <- as.integer(cl2)                         # cell index per stage-2 unit
  Di <- D[idx]; gi <- g[idx]; ai <- a[idx]; bi <- b[idx]
  degi <- degenerate[ci]
  pi_i  <- pih[ci];  mu1i <- mu1[ci];  mu0i <- mu0[ci];  taui <- tau[ci]
  a1i   <- a1[ci];   b0i  <- b0[ci];   thci <- thc[ci]
  pi_i[degi] <- 0.5                             # keep the unused branch finite

  if_eta <- ifelse(degi, 0 - eta,
                   (Di / pi_i) * (gi - mu1i) -
                   ((1 - Di) / (1 - pi_i)) * (gi - mu0i) +
                   taui - eta)
  if_th  <- ifelse(degi | thci <= 0, 0 - theta,
                   (Di / pi_i) * (ai - a1i) +
                   ((1 - Di) / (1 - pi_i)) * (bi - b0i) +
                   thci - theta)

  z <- stats::qnorm(1 - alpha)
  eta_se  <- stats::sd(if_eta) / sqrt(n2)
  eta_lcb <- max(eta - z * eta_se, 0)

  # The ratio divides by an estimated theta that can be very small -- for the
  # second category of Baranov's five-valued mediator it is 0.040 -- and a
  # slightly negative eta then becomes a large negative nu (-0.222 there). The
  # numerator is truncated first, since nu >= 0, and the ratio is withheld when
  # the denominator is too small to support it.
  # The ratio eta / theta is a valid lower bound on nu_k only when the
  # denominator is pinned down by the identified set. Monotonicity pins it at
  # BOTH ends of the support of M -- treated units at the bottom category were
  # there under control too, so theta_00 = P(M = m_0 | D = 1), and controls at
  # the top stay there, so theta_KK = P(M = m_K | D = 0) -- and nowhere else.
  # In the interior both numerator and denominator are only bounded below, and
  # a <= c with b <= d says nothing about a/b versus c/d: the division is
  # ANTI-CONSERVATIVE. On Baranov's five-valued mediator it gives 0.0704 at
  # k = 4 where the correct bound is 0.0495. Use nu_bound_joint(), which solves
  # KR's linear-fractional program by the Charnes-Cooper transformation.
  k_pinned <- (kk_index <- match(at_group, mvals_sorted)) %in% c(1L, length(mvals_sorted))
  theta_floor <- 5 / n2
  if (k_pinned && theta > theta_floor) {
    nu <- max(eta, 0) / theta
    if_nu  <- (if_eta - nu * if_th) / theta        # delta method
    nu_se  <- stats::sd(if_nu) / sqrt(n2)
    nu_lcb <- max(eta - z * eta_se, 0) / theta
  } else {
    nu <- nu_se <- NA_real_
    nu_lcb <- NA_real_
  }
  if (!k_pinned)
    attr(nu, "reason") <- paste("theta_kk is not pinned down at an interior",
      "mediator value; use nu_bound_joint() for the share scale")

  # The plug-in point estimate is the one practitioners compute: the supremum
  # is taken on the same data used to report the number, so it is biased up.
  # No confidence statement is attached, because a normal one would not be
  # valid here -- that is the point of the comparison.
  if (method == "plugin") {
    eta_se <- NA_real_; eta_lcb <- NA_real_
    nu_se  <- NA_real_; nu_lcb  <- NA_real_
  }

  list(eta = eta, eta_se = eta_se, eta_lcb = eta_lcb,
       nu = nu, nu_se = nu_se, nu_lcb = nu_lcb,
       theta = theta,
       n = n, n_stage2 = n2,
       n_cells = length(lv),
       n_cells_degenerate = sum(degenerate),
       mass_degenerate = sum(w[degenerate]),
       n_cells_nonmonotone = sum(defier > 1e-9, na.rm = TRUE),
       mass_nonmonotone = sum(w[which(defier > 1e-9)]),
       max_cell_defier_share = if (all(is.na(defier))) NA_real_ else
                                 max(defier, na.rm = TRUE),
       defier_share_overall = defier_overall,
       method = method)
}
