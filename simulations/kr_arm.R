# =====================================================================
# kr_arm.R -- KR's own test as a Monte Carlo arm.
#
# WHY THIS EXISTS. Every other arm in these simulations reports a lower
# CONFIDENCE bound and calls the sharp null rejected when that bound clears
# zero. KR do something different: they test the sharp null directly, as a
# system of moment inequalities, and report a p-value. The two are not
# comparable as estimators, but they ARE comparable as tests of the same null,
# and an audience will ask. On the real data KR's test wins
# (simulations/kr_test.R: Bursztyn p = 0.0428, Baranov p = 0.00146, against our
# LCB clearing zero in 42% and 83% of splits), for two structural reasons:
#
#   1. we spend half the sample choosing Ahat, and
#   2. "is the lower bound on ONE k detectably positive" is a blunter question
#      than "is the inequality system across ALL k and all A infeasible" --
#      KR aggregate evidence we discard.
#
# What this arm settles is the question those two reasons leave open: is our
# power gain under heterogeneity a real gain, or an artefact of comparing
# against a pooled bound rather than against a proper test? KR's test is
# pooled, so if conditioning is doing the work it should beat KR's test in the
# sign-changing design and lose to it in the homogeneous one.
#
# SCOPE. The binary-M branch of KR's Cox-Shi implementation has no nuisance
# parameters and solves its QP with quadprog. The multivalued branch profiles
# out theta and needs osqp, which is not installable in this environment
# (not on CRAN mirrors reachable here, no Debian package). Design C therefore
# falls back to ARP, which takes the same A matrices through a different
# critical value; where neither is available the arm returns NA rather than
# silently substituting a different solver into KR's own procedure.
# =====================================================================

SRC_KR <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_discretize_y.R", "src_cox_shi_nonuisance.R",
            "src_test_sharp_null_coxandshi_binary_m.R",
            "src_nonparametric_bootstrap.R", "src_test_sharp_null.R")) {
  p <- file.path(SRC_KR, f)
  if (file.exists(p)) source(p)
}

# ---------------------------------------------------------------------
# RNG HAZARD. KR's variance step calls compute_bootstrap_draws_clustered,
# which runs set.seed(b) for b = 1, ..., B (src_nonparametric_bootstrap.R).
# The global stream is therefore left in a FIXED, data-independent state by
# every call. In a Monte Carlo loop that is fatal: the next rep's
# simulate_data_* call starts from that fixed state and redraws the SAME
# dataset, so a 300-rep study collapses to one draw repeated 300 times. It
# does not error and it does not look obviously wrong -- it shows up as
# rejection rates pinned at 0.000 and 1.000 and as non-monotone size
# (observed 2026-11-10: pooled power 0.993 at n = 284 then 0.000 at n = 1000,
# KR size 0.330 at one n and 0.000 either side).
#
# So kr_reject saves the stream and restores it on exit. The guard lives here
# rather than in each driver so that no caller has to remember. KR's own
# bootstrap is unaffected: it seeds itself.
# ---------------------------------------------------------------------
with_rng_preserved <- function(expr) {
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1L)
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  force(expr)
}

#' Reject the sharp null by KR's test?
#'
#' The global RNG state is restored before returning, so the caller's stream
#' advances only through its own draws. See the RNG note above.
#'
#' @param df,d,m,y as elsewhere.
#' @param method "CS" (Cox-Shi, KR's headline) or "ARP".
#' @param alpha level.
#' @param num_Ybins passed through; NULL leaves Y alone.
#' @return list(reject, pval, status). status is "ok", or the error, so a
#'   failure is visible in the output instead of being counted as a
#'   non-rejection -- which would bias the arm downward exactly where the
#'   design is hardest.
kr_reject <- function(df, d, m, y, method = "CS", alpha = 0.05,
                      num_Ybins = NULL, cluster = NULL) {
  binary_M <- length(unique(stats::na.omit(df[[m]]))) == 2L
  out <- with_rng_preserved(tryCatch({
    if (binary_M) {
      r <- test_sharp_null_coxandshi_binary_m(df = df, d = d, m = m, y = y,
                                              num_Ybins = num_Ybins,
                                              cluster = cluster,
                                              alpha = alpha)
    } else {
      r <- test_sharp_null(df = df, d = d, m = m, y = y, method = method,
                           num_Ybins = num_Ybins, cluster = cluster,
                           alpha = alpha, analytic_variance = TRUE)
    }
    r
  }, error = function(e) structure(list(), class = "kr_err",
                                   msg = conditionMessage(e))))
  if (inherits(out, "kr_err"))
    return(list(reject = NA, pval = NA_real_, status = attr(out, "msg")))
  p <- suppressWarnings(as.numeric(unlist(out)[["pval"]]))
  if (is.null(p) || is.na(p)) {
    nm <- names(unlist(out))
    cand <- nm[grepl("p.?val|pvalue", nm, ignore.case = TRUE)]
    if (length(cand)) p <- as.numeric(unlist(out)[[cand[1]]])
  }
  list(reject = as.integer(p <= alpha), pval = p, status = "ok")
}
