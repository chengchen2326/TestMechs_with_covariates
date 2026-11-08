# Does the conditioning gain survive the defier relaxation on RAW Baranov?
# Identification level only: empirical distributions treated as populations, no
# sampling, no splitting. If the gain does not survive here, no amount of
# inference machinery will recover it.
#
# Two ways to impose the relaxation, the N2 choice:
#   cell   each cell gets its own dbar, at least its own minimum feasible share
#   pop    one dbar for the pooled marginals; cells inherit it, raised to their
#          own minimum feasible share where necessary (otherwise infeasible)
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_theta_kk_min_defiers.R")) source(file.path("../src",f))

load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
keep <- c("treat","relationship_husb","motherfinancial","grandmother_baseline","MIL","girl")
br <- br[stats::complete.cases(br[,keep]),]
for (v in keep) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
K <- length(mv); yv <- sort(unique(br$yb))
cat(sprintf("Baranov n = %d, K = %d, Q = %d\n\n", nrow(br), K, length(yv)))

marg <- function(s) list(a = as.numeric(prop.table(table(factor(s$mc[s$treat==0], levels=1:K)))),
                         b = as.numeric(prop.table(table(factor(s$mc[s$treat==1], levels=1:K)))))
supd <- function(s, k) {
  p1 <- sapply(yv, function(yy) mean(s$yb==yy & s$mc==k & s$treat==1))/mean(s$treat==1)
  p0 <- sapply(yv, function(yy) mean(s$yb==yy & s$mc==k & s$treat==0))/mean(s$treat==0)
  sum(pmax(p1-p0,0))
}
eta_of <- function(s, k, dbar) {
  m <- marg(s); d <- max(dbar, min_defier_share(m$a, m$b))
  r <- theta_kk_min_defiers(m$a, m$b, k, dbar = d, dual = FALSE)
  if (r$status != 0) return(c(eta = NA_real_, theta = NA_real_, dbar = d))
  c(eta = max(supd(s,k) - (m$b[k] - r$value), 0), theta = r$value, dbar = d)
}

for (xv in c("grandmother_baseline","MIL","girl")) {
  gs <- split(seq_len(nrow(br)), br[[xv]]); w <- vapply(gs, length, 1L); w <- w/sum(w)
  cat(sprintf("--- cells by %s (sizes %s) ---\n", xv, paste(vapply(gs,length,1L), collapse="/")))
  cat(sprintf("  min feasible dbar: pooled %.4f | cells %s\n",
      min_defier_share(marg(br)$a, marg(br)$b),
      paste(sprintf("%.4f", vapply(gs, function(g) { m <- marg(br[g,]); min_defier_share(m$a,m$b) }, 0)), collapse=" ")))
  cat(sprintf("%4s %8s %11s %11s %9s %s\n","k","dbar","E[eta(X)]","eta_pooled","gain",""))
  for (k in 3:K) for (dbar in c(0, 0.0154, 0.05, 0.10)) {
    pc <- vapply(gs, function(g) eta_of(br[g,], k, dbar), numeric(3))
    pl <- eta_of(br, k, dbar)
    ec <- sum(w * pc["eta",])
    cat(sprintf("%4d %8.4f %11.6f %11.6f %9.6f  (cell dbar used: %s)\n",
        k, dbar, ec, pl["eta"], ec - pl["eta"],
        paste(sprintf("%.3f", pc["dbar",]), collapse=" ")))
  }
  cat("\n")
}
