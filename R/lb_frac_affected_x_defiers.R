#' @title Inference on the conditional bound with a population defier cap
#' @description Split-sample lower confidence bound on
#'   \eqn{E[\theta_{kk}(X)\nu_k(X)]} when up to \eqn{\bar d} of the population
#'   may be defiers, with the cap imposed once on the population rather than
#'   cell by cell.
#'
#'   Stage 1, on the first fold, estimates each cell's mediator marginals and
#'   \eqn{\hat A_k(x)}, then solves the dual of the joint program. Stage 2
#'   evaluates the resulting linear functional on the held-out fold.
#'
#'   Validity rests on weak duality twice over and on nothing else. Any
#'   dual-feasible point gives a value at or below the conditional bound, which
#'   is itself at or below the truth, so a badly chosen first stage costs
#'   tightness and not coverage. In the normalised scaling the dual constraints
#'   \eqn{\tilde u^{(x)}_l + \tilde v^{(x)}_j - \tilde p_x 1\{(l,j)=(k,k)\} -
#'   q 1\{l > j\} \le 0}, \eqn{\tilde p_x \in [0,1]}, \eqn{q \ge 0} do not
#'   involve the cell weights at all, so a point solved with fold-1 weights
#'   stays feasible for the population program -- the estimate is valid even
#'   though the weights are estimated.
#'
#'   A cell with no treated or no control unit in the second fold has its dual
#'   variables set to zero, which is feasible and contributes nothing. That is
#'   conservative rather than undefined.
#'
#' @param df A data frame.
#' @param d,m,y Names of treatment, mediator and outcome.
#' @param x Character vector of covariate names defining the cells.
#' @param at_group Value of m whose always-takers are of interest.
#' @param max_defiers_share The population cap \eqn{\bar d}. Default 0.
#' @param alpha One-sided level. Default 0.05.
#' @param split_frac Share of the sample for stage 1. Default 0.3.
#' @param folds,seed Optional fold assignment or seed.
#' @return A list with \code{eta}, \code{eta_se}, \code{eta_lcb},
#'   \code{dbar_used}, \code{dbar_min_feasible}, \code{n}, \code{n_stage2},
#'   \code{n_cells}, \code{n_cells_zeroed}, \code{stage1_bound}.
#' @export
lb_frac_affected_x_defiers <- function(df, d, m, y, x,
                                       at_group = NULL,
                                       max_defiers_share = 0,
                                       alpha = 0.05,
                                       split_frac = 0.3,
                                       folds = NULL,
                                       seed = NULL) {

  if (is.null(at_group) || length(at_group) != 1L)
    stop("at_group must be a single value of m")
  cols <- c(d, m, y, x)
  dd <- df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]
  n <- nrow(dd)
  D <- as.numeric(dd[[d]]); M <- dd[[m]]; Y <- dd[[y]]
  if (!all(D %in% c(0, 1))) stop("d must be binary 0/1")

  mvals <- sort(unique(M)); K <- length(mvals)
  kk <- match(at_group, mvals)
  if (is.na(kk)) stop("at_group is not a value taken by m")
  yv <- sort(unique(Y))
  cell <- interaction(dd[, x, drop = FALSE], drop = TRUE, sep = "|")
  lv <- levels(cell); J <- length(lv)

  if (!is.null(folds)) {
    fv <- folds[stats::complete.cases(df[, cols, drop = FALSE])]
  } else {
    if (!is.null(seed)) set.seed(seed)
    fv <- rep(2L, n); fv[sample.int(n, floor(split_frac * n))] <- 1L
  }
  f1 <- fv == 1L; f2 <- fv == 2L
  if (!any(f1) || !any(f2)) stop("one of the folds is empty")

  # ------------------------------------------------- stage 1, on fold 1
  Am <- Bm <- matrix(0, J, K)
  Sv <- numeric(J); wv <- numeric(J)
  Ahat <- stats::setNames(vector("list", J), lv)
  ok1 <- logical(J)
  for (jx in seq_len(J)) {
    cx <- lv[jx]
    i1 <- which(f1 & cell == cx)
    wv[jx] <- length(i1)
    t1 <- i1[D[i1] == 1]; t0 <- i1[D[i1] == 0]
    Ahat[[cx]] <- logical(length(yv))
    if (!length(t1) || !length(t0)) next
    ok1[jx] <- TRUE
    Am[jx, ] <- vapply(mvals, function(v) mean(M[t0] == v), numeric(1))
    Bm[jx, ] <- vapply(mvals, function(v) mean(M[t1] == v), numeric(1))
    p1 <- vapply(yv, function(v) mean(Y[t1] == v & M[t1] == at_group), numeric(1))
    p0 <- vapply(yv, function(v) mean(Y[t0] == v & M[t0] == at_group), numeric(1))
    Ahat[[cx]] <- p1 > p0
    Sv[jx] <- sum(pmax(p1 - p0, 0))
  }
  wv <- if (sum(wv) > 0) wv / sum(wv) else rep(1 / J, J)
  # Cells with no stage-1 contrast get uniform marginals, which keeps the
  # program feasible; their dual variables are zeroed below in any case.
  for (jx in seq_len(J)) if (!ok1[jx]) {
    Am[jx, ] <- rep(1 / K, K); Bm[jx, ] <- rep(1 / K, K); Sv[jx] <- 0
  }

  # dbar is a sensitivity parameter, not something to estimate. The minimum
  # feasible defier share is a maximum of survivor-function shortfalls -- a
  # supremum of noisy estimates -- so it is biased UPWARD, the mirror image of
  # the winner's curse in sup_A Delta. Raising dbar to the fold-1 estimate of it
  # inflates the identified set and destroys the bound: on Baranov a fold of 30%
  # puts that estimate at 0.045-0.052 against a full-sample 0.020, and at 0.05
  # the bound is zero. So dbar is taken as given, exactly as KR take
  # max_defiers_share, and the full-sample minimum is reported as a diagnostic.
  #
  # When the fold-1 marginals are themselves incompatible with dbar, they are
  # pulled toward monotonicity until the program is feasible. This cannot harm
  # validity: in the normalised scaling the dual constraints involve only
  # p_tilde and q, never the marginals, so the marginals enter the objective
  # alone and a projection changes which feasible point is chosen, not whether
  # it is feasible.
  dmin <- sum(wv * vapply(seq_len(J), function(jx)
            min_defier_share(Am[jx, ], Bm[jx, ]), numeric(1)))
  dbar <- max_defiers_share

  Bp <- Bm
  if (dmin > dbar) {
    Bmono <- Bm
    for (jx in seq_len(J))
      Bmono[jx, ] <- monotone_project(Am[jx, ], Bm[jx, ])$pm1
    lam_lo <- 0; lam_hi <- 1
    for (it in 1:40) {
      lam <- (lam_lo + lam_hi) / 2
      Btry <- (1 - lam) * Bm + lam * Bmono
      dtry <- sum(wv * vapply(seq_len(J), function(jx)
                min_defier_share(Am[jx, ], Btry[jx, ]), numeric(1)))
      if (dtry > dbar) lam_lo <- lam else lam_hi <- lam
    }
    Bp <- (1 - lam_hi) * Bm + lam_hi * Bmono
  }

  du <- eta_cond_joint_dual(Sv, Am, Bp, wv, kk, dbar = dbar)
  if (is.na(du$value)) stop("the stage-1 dual program is infeasible")
  Ut <- du$u_t; Vt <- du$v_t; Pt <- du$p_t; q <- du$q

  # ------------------------------------------------- stage 2, on fold 2
  idx <- which(f2); n2 <- length(idx)
  cl2 <- cell[idx]
  inA <- logical(n2)
  for (jj in seq_along(yv)) {
    sel <- Y[idx] == yv[jj]
    inA[sel] <- vapply(as.character(cl2[sel]), function(cx) Ahat[[cx]][jj], logical(1))
  }
  hitA <- as.numeric(inA & M[idx] == at_group)
  isk  <- as.numeric(M[idx] == at_group)
  mpos <- match(M[idx], mvals)
  jpos <- match(as.character(cl2), lv)

  zeroed <- logical(J)
  m1 <- m0 <- pih <- numeric(J); w2 <- numeric(J)
  for (jx in seq_len(J)) {
    s <- which(jpos == jx); w2[jx] <- length(s)
    s1 <- s[D[idx][s] == 1]; s0 <- s[D[idx][s] == 0]
    if (!length(s1) || !length(s0) || !ok1[jx]) { zeroed[jx] <- TRUE; next }
    pih[jx] <- length(s1) / length(s)
    H1 <- Vt[jx, mpos[s1]] + Pt[jx] * (hitA[s1] - isk[s1])
    H0 <- Ut[jx, mpos[s0]] - Pt[jx] * hitA[s0]
    m1[jx] <- mean(H1); m0[jx] <- mean(H0)
  }
  w2 <- if (sum(w2) > 0) w2 / sum(w2) else rep(1 / J, J)
  contrib <- ifelse(zeroed, 0, m1 + m0)
  eta <- sum(w2 * contrib) - dbar * q

  # Saturated influence function: within a cell each arm's residuals sum to
  # zero, and the cell terms are centred, so the mean is zero by construction.
  ifv <- numeric(n2)
  for (jx in seq_len(J)) {
    s <- which(jpos == jx)
    if (zeroed[jx]) { ifv[s] <- 0 - sum(w2 * contrib); next }
    Di <- D[idx][s]
    H1 <- Vt[jx, mpos[s]] + Pt[jx] * (hitA[s] - isk[s])
    H0 <- Ut[jx, mpos[s]] - Pt[jx] * hitA[s]
    ifv[s] <- (Di / pih[jx]) * (H1 - m1[jx]) +
              ((1 - Di) / (1 - pih[jx])) * (H0 - m0[jx]) +
              contrib[jx] - sum(w2 * contrib)
  }
  eta_se <- stats::sd(ifv) / sqrt(n2)
  eta_lcb <- max(eta - stats::qnorm(1 - alpha) * eta_se, 0)

  list(eta = eta, eta_se = eta_se, eta_lcb = eta_lcb,
       dbar_used = dbar, dbar_min_feasible_fold1 = dmin,
       projection_used = dmin > dbar,
       n = n, n_stage2 = n2, n_cells = J, n_cells_zeroed = sum(zeroed),
       stage1_bound = du$value)
}
