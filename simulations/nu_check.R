# Does the shortcut nu = eta / theta_kk_min agree with KR's linear-fractional
# program? Run on monotone-projected Baranov marginals so both are well defined.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_simulate_data_multiM_x.R","src_theta_kk_min_defiers.R",
            "src_lb_frac_affected_x.R")) source(file.path("../src",f))
load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
br <- br[stats::complete.cases(br[,c("treat","relationship_husb","motherfinancial")]),]
for (v in c("treat","relationship_husb","motherfinancial")) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
K <- length(mv); br$const <- 1L

# Reweight the treated arm so the M marginals are monotone, then rebuild a
# dataset with those marginals (weights keep lb_frac_affected applicable).
a <- as.numeric(prop.table(table(br$mc[br$treat==0])))
b <- as.numeric(prop.table(table(br$mc[br$treat==1])))
bp <- monotone_project(a,b)$pm1
cat(sprintf("projection shift = %.4f\n", monotone_project(a,b)$shift))
br$w <- ifelse(br$treat==1, bp[br$mc]/b[br$mc], 1)

cat(sprintf("\n%4s %14s %14s %12s %10s\n","k","eta/theta (mine)","lb_frac_affected","|diff|","theta_min"))
for (k in 1:K) {
  r <- lb_frac_affected_x(br,"treat","mc","yb", x="const", at_group=k, method="plugin")
  pkg <- tryCatch(as.numeric(lb_frac_affected(br, d="treat", m="mc", y="yb", w="w", at_group=k)),
                  error=function(e) NA_real_, warning=function(w) NA_real_)
  cat(sprintf("%4d %14s %14s %12s %10.4f\n", k,
      ifelse(is.na(r$nu),"NA",sprintf("%.6f",r$nu)),
      ifelse(is.na(pkg),"NA",sprintf("%.6f",pkg)),
      ifelse(is.na(r$nu)||is.na(pkg),"-",sprintf("%.2e",abs(r$nu-pkg))), r$theta))
}

cat("\n=== Charnes-Cooper program vs lb_frac_affected (pooled, J = 1) ===\n")
source("../src/src_nu_bound_joint.R")
yv <- sort(unique(br$yb))
am <- matrix(as.numeric(prop.table(table(br$mc[br$treat==0]))), 1)
bmv <- tapply(br$w[br$treat==1], br$mc[br$treat==1], sum); bmv <- bmv/sum(bmv)
bm <- matrix(as.numeric(bmv), 1)
cat(sprintf("%4s %16s %16s %12s\n","k","Charnes-Cooper","lb_frac_affected","|diff|"))
for (k in 1:K) {
  t1 <- br$treat==1; t0 <- br$treat==0
  p1 <- vapply(yv, function(v) sum(br$w[t1 & br$yb==v & br$mc==k])/sum(br$w[t1]), 0)
  p0 <- vapply(yv, function(v) mean(br$yb[t0]==v & br$mc[t0]==k), 0)
  Sx <- sum(pmax(p1-p0,0))
  cc <- nu_bound_joint(Sx, am, bm, 1, k, dbar=0)
  pkg <- tryCatch(as.numeric(lb_frac_affected(br, d="treat", m="mc", y="yb", w="w", at_group=k)),
                  error=function(e) NA_real_, warning=function(w) NA_real_)
  cat(sprintf("%4d %16s %16s %12s\n", k,
      ifelse(is.na(cc$value),"NA",sprintf("%.8f",cc$value)),
      ifelse(is.na(pkg),"NA",sprintf("%.8f",pkg)),
      ifelse(is.na(cc$value)||is.na(pkg),"-",sprintf("%.2e",abs(cc$value-pkg)))))
}
