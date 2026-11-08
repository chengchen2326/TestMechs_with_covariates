suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_theta_kk_min_defiers.R")) source(file.path("../src",f))

cat("=== CHECK 1: dbar = 0 must reproduce KR's closed form, and strong duality ===\n")
set.seed(4); worst_cf <- 0; worst_dg <- 0; ntest <- 0
for (K in 2:6) for (rep in 1:300) {
  a <- rgamma(K,1); a <- a/sum(a); b <- rgamma(K,1); b <- b/sum(b)
  # make monotone so dbar = 0 is feasible
  s0 <- rev(cumsum(rev(a))); s1 <- rev(cumsum(rev(b)))
  if (any(s1 < s0 - 1e-12)) next
  for (k in 1:K) {
    r <- theta_kk_min_defiers(a,b,k,dbar=0)
    if (r$status != 0) next
    cf <- theta_kk_min_closed_form(a,b,k)
    worst_cf <- max(worst_cf, abs(r$value - cf))
    worst_dg <- max(worst_dg, abs(r$value - r$dual_value))
    ntest <- ntest + 1
  }
}
cat(sprintf("  %d tests | max |LP - eq.(8)| = %.2e | max duality gap = %.2e\n\n",
            ntest, worst_cf, worst_dg))

cat("=== CHECK 2: dbar > 0 must match lb_frac_affected's max_defiers_share ===\n")
load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
br <- br[stats::complete.cases(br[,c("treat","relationship_husb","motherfinancial")]),]
for (v in c("treat","relationship_husb","motherfinancial")) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
a <- as.numeric(prop.table(table(br$mc[br$treat==0])))
b <- as.numeric(prop.table(table(br$mc[br$treat==1])))
cat(sprintf("  minimum feasible defier share = %.6f\n", min_defier_share(a,b)))
cat(sprintf("%6s %4s %14s %16s %12s %12s\n","dbar","k","theta_kk LP","via lb_frac_aff","nu_LP","nu_pkg"))
for (dbar in c(0.0154, 0.05, 0.10, 0.25)) {
  for (k in 3:5) {
    r <- theta_kk_min_defiers(a,b,k,dbar=dbar)
    sup_delta <- { p1 <- sapply(sort(unique(br$yb)), function(yy) mean(br$yb==yy & br$mc==k & br$treat==1))/mean(br$treat==1)
                   p0 <- sapply(sort(unique(br$yb)), function(yy) mean(br$yb==yy & br$mc==k & br$treat==0))/mean(br$treat==0)
                   sum(pmax(p1-p0,0)) }
    eta <- max(sup_delta - (b[k] - r$value), 0)
    nu_lp <- if (r$value > 1e-12) eta/r$value else NA_real_
    pkg <- tryCatch(as.numeric(lb_frac_affected(br, d="treat", m="mc", y="yb", at_group=k,
                      max_defiers_share=dbar)), error=function(e) NA_real_, warning=function(w) NA_real_)
    cat(sprintf("%6.4f %4d %14.6f %16s %12.6f %12.6f\n", dbar, k, r$value, "-", nu_lp, pkg))
  }
}
