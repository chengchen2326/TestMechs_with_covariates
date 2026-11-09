#' @title Cross-fitted conditional bound
#' @description Runs the split-sample procedure twice with the folds exchanged
#'   and averages, so every observation contributes to the estimate rather than
#'   half of them being spent on choosing \eqn{\hat A_k}.
#'
#'   \strong{Validity of the point estimate.} Fold 1 selects \eqn{\hat A_1} and
#'   fold 2 evaluates it; fold 2 selects \eqn{\hat A_2} and fold 1 evaluates it.
#'   Each half-estimate is conditionally unbiased for \eqn{\eta^{\hat A_f}_k},
#'   and any fixed \eqn{A} gives a valid lower bound, so the average estimates
#'   \eqn{(\eta^{\hat A_1}_k + \eta^{\hat A_2}_k)/2 \le \eta_k}. Averaging two
#'   valid lower bounds is a valid lower bound; nothing needs to converge.
#'
#'   \strong{Variance.} The two halves are dependent, because the nuisance in
#'   one is a function of the data in the other. Standard cross-fitting ignores
#'   this and pools the influence functions, which is justified when the
#'   nuisance settles down -- and \eqn{\hat A_k} need not. Cauchy--Schwarz gives
#'   an alternative that never relies on that:
#'   \deqn{\mathrm{Var}\Big(\tfrac{\hat\eta_a + \hat\eta_b}{2}\Big)
#'         \le \Big(\tfrac{\mathrm{sd}(\hat\eta_a) + \mathrm{sd}(\hat\eta_b)}{2}\Big)^2 .}
#'   With equal folds that conservative standard error equals the one a single
#'   split would carry, so cross-fitting is free: the same interval width around
#'   a better-centred and less dispersed estimate. The pooled-influence-function
#'   standard error is also returned, as \code{eta_se_pooled}; whether it is
#'   trustworthy here is an empirical question and
#'   \code{simulations/crossfit_mc.R} answers it.
#'
#' @param df,d,m,y,x,at_group,alpha As in \code{lb_frac_affected_x}.
#' @param split_frac Share of the sample in the first fold. Default 0.5, which
#'   is the variance-minimising choice when both folds are used.
#' @param seed Optional seed for the split.
#' @return A list with \code{eta}, \code{eta_se} (conservative),
#'   \code{eta_se_pooled} (optimistic), \code{eta_lcb},
#'   \code{eta_lcb_pooled}, and the two half-fits in \code{fits}.
#' @export
lb_frac_affected_x_cf <- function(df, d, m, y, x,
                                  at_group = NULL,
                                  alpha = 0.05,
                                  split_frac = 0.5,
                                  seed = NULL) {
  n <- nrow(df)
  if (!is.null(seed)) set.seed(seed)
  fv <- rep(2L, n)
  fv[sample.int(n, floor(split_frac * n))] <- 1L

  fa <- lb_frac_affected_x(df, d, m, y, x, at_group = at_group,
                           method = "split", alpha = alpha, folds = fv)
  fb <- lb_frac_affected_x(df, d, m, y, x, at_group = at_group,
                           method = "split", alpha = alpha,
                           folds = ifelse(fv == 1L, 2L, 1L))

  eta <- (fa$eta + fb$eta) / 2
  se_cons <- (fa$eta_se + fb$eta_se) / 2            # Cauchy-Schwarz
  # Pooled influence function: each half's residuals, rescaled so the two
  # halves together have the weight of one sample.
  se_pool <- sqrt((fa$eta_se^2 + fb$eta_se^2)) / 2  # assumes zero covariance

  z <- stats::qnorm(1 - alpha)
  list(eta = eta,
       eta_se = se_cons, eta_lcb = max(eta - z * se_cons, 0),
       eta_se_pooled = se_pool, eta_lcb_pooled = max(eta - z * se_pool, 0),
       n = n, fits = list(a = fa, b = fb))
}
