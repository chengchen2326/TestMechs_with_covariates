# Cell-level versus population-level defier cap. Identification level only.
#
# NOTE, and this trap has now bitten twice: the pooled comparison must be built
# from the SAME covariate distribution as the conditional side, i.e.
# sum_x w(x) a^(x) and sum_x w(x) b^(x), NOT the raw pooled tabulation. The raw
# tabulation weights p_d(. | x) by P(X = x | D = d), which differs across arms
# in finite samples even under randomisation, and the comparison then appears to
# violate dominance.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_theta_kk_min_defiers.R","src_eta_cond_joint.R")) source(file.path("../src",f))

cell_parts <- function(s, k, K, yv) {
  a <- as.numeric(prop.table(table(factor(s$mc[s$treat==0], levels=1:K))))
  b <- as.numeric(prop.table(table(factor(s$mc[s$treat==1], levels=1:K))))
  p1 <- sapply(yv, function(yy) mean(s$yb==yy & s$mc==k & s$treat==1))/mean(s$treat==1)
  p0 <- sapply(yv, function(yy) mean(s$yb==yy & s$mc==k & s$treat==0))/mean(s$treat==0)
  list(a=a, b=b, p1=p1, p0=p0, S=sum(pmax(p1-p0,0)))
}
eta_from <- function(a, b, S, k, dbar, raise=TRUE) {
  d <- if (raise) max(dbar, min_defier_share(a,b)) else dbar
  r <- theta_kk_min_defiers(a,b,k,dbar=d,dual=FALSE)
  if (r$status!=0) return(c(eta=NA_real_, dbar=d))
  c(eta=max(S - (b[k]-r$value), 0), dbar=d)
}
# Pooled built by MIXING the cells with the common weights w.
pooled_parts <- function(P, w, k) {
  a <- Reduce(`+`, Map(function(wi,p) wi*p$a, w, P))
  b <- Reduce(`+`, Map(function(wi,p) wi*p$b, w, P))
  d <- Reduce(`+`, Map(function(wi,p) wi*(p$p1-p$p0), w, P))
  list(a=a, b=b, S=sum(pmax(d,0)))
}

cat("=== CHECK: synthetic monotone cells, dbar = 0: joint == cell-by-cell ===\n")
set.seed(21); worst <- 0; nt <- 0
for (rep in 1:200) {
  K <- 4; J <- 3
  A <- B <- matrix(0,J,K); Sv <- numeric(J)
  ok <- TRUE
  for (x in 1:J) {
    a <- rgamma(K,1); a <- a/sum(a)
    b <- rgamma(K,1); b <- b/sum(b)
    if (min_defier_share(a,b) > 1e-12) { ok <- FALSE; break }
    A[x,] <- a; B[x,] <- b; Sv[x] <- runif(1,0,0.1)
  }
  if (!ok) next
  w <- rep(1/J,J)
  for (k in 1:K) {
    j <- eta_cond_joint(Sv,A,B,w,k,dbar=0)
    c0 <- sum(w*vapply(1:J, function(x) eta_from(A[x,],B[x,],Sv[x],k,0,raise=FALSE)["eta"], 0))
    if (!is.na(j$value)) { worst <- max(worst, abs(j$value-c0)); nt <- nt+1 }
  }
}
cat(sprintf("  %d tests | max |joint - cell-by-cell| = %.2e\n\n", nt, worst))

load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
keep <- c("treat","relationship_husb","motherfinancial","grandmother_baseline","MIL","girl")
br <- br[stats::complete.cases(br[,keep]),]
for (v in keep) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
K <- length(mv); yv <- sort(unique(br$yb))

cat("=== Baranov: cell-level vs population-level cap (pooled mixed correctly) ===\n")
for (xv in c("grandmother_baseline","MIL","girl")) {
  gs <- split(seq_len(nrow(br)), br[[xv]]); w <- as.numeric(vapply(gs,length,1L)); w <- w/sum(w)
  P3 <- lapply(gs, function(g) cell_parts(br[g,],3,K,yv))
  cat(sprintf("\n--- %s | cell min feasible dbar: %s | mixed-pooled min: %.4f ---\n", xv,
      paste(sprintf("%.4f", vapply(P3, function(p) min_defier_share(p$a,p$b), 0)), collapse=" "),
      min_defier_share(pooled_parts(P3,w,3)$a, pooled_parts(P3,w,3)$b)))
  cat(sprintf("%3s %7s %11s %11s %11s %10s %10s\n","k","dbar","cond(cell)","cond(pop)","pooled","gain cell","gain pop"))
  for (k in 3:K) for (dbar in c(0.02, 0.03, 0.05, 0.08)) {
    P <- lapply(gs, function(g) cell_parts(br[g,],k,K,yv))
    Sv <- vapply(P, `[[`, numeric(1), "S")
    Am <- do.call(rbind, lapply(P,`[[`,"a")); Bm <- do.call(rbind, lapply(P,`[[`,"b"))
    pp <- pooled_parts(P, w, k)
    jp <- eta_cond_joint(Sv, Am, Bm, w, k, dbar=dbar)
    cc <- sum(w*vapply(seq_along(P), function(x)
            eta_from(P[[x]]$a,P[[x]]$b,P[[x]]$S,k,dbar)["eta"], 0))
    pl <- eta_from(pp$a, pp$b, pp$S, k, dbar)["eta"]
    f1 <- if (!is.na(jp$value) && jp$value < pl-1e-9) "  POP-FAIL" else ""
    f2 <- if (!is.na(cc) && cc < pl-1e-9) "  CELL-FAIL" else ""
    cat(sprintf("%3d %7.4f %11.6f %11.6f %11.6f %+10.6f %+10.6f%s%s\n",
      k, dbar, cc, ifelse(is.na(jp$value),NA,jp$value), pl, cc-pl,
      ifelse(is.na(jp$value),NA,jp$value-pl), f1, f2))
  }
}
