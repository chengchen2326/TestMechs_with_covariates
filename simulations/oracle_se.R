# Oracle stage 1: fix the dual point at the FULL-SAMPLE optimum, then ask what
# standard error stage 2 carries. This separates the sampling error of the
# estimator from the noise in choosing the dual point, and so answers whether
# the identification result is detectable at this sample size at all.
suppressMessages({library(haven);library(Rglpk)})
for (f in c("src_simulate_data_multiM_x.R","src_theta_kk_min_defiers.R",
            "src_eta_cond_joint.R")) source(file.path("../src",f))
load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
keep <- c("treat","relationship_husb","motherfinancial","grandmother_baseline","MIL")
br <- br[stats::complete.cases(br[,keep]),]
for (v in keep) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
K <- length(mv); yv <- sort(unique(br$yb)); br$const <- 1L

oracle <- function(xv, kk, dbar, box = 2) {
  gs <- split(seq_len(nrow(br)), br[[xv]]); J <- length(gs)
  w <- as.numeric(vapply(gs,length,1L)); w <- w/sum(w)
  A <- B <- matrix(0,J,K); S <- numeric(J); Ah <- vector("list",J)
  for (x in 1:J) { s <- br[gs[[x]],]; t1 <- s$treat==1; t0 <- s$treat==0
    A[x,] <- vapply(1:K, function(v) mean(s$mc[t0]==v), 0)
    B[x,] <- vapply(1:K, function(v) mean(s$mc[t1]==v), 0)
    p1 <- vapply(yv, function(v) mean(s$yb[t1]==v & s$mc[t1]==kk), 0)
    p0 <- vapply(yv, function(v) mean(s$yb[t0]==v & s$mc[t0]==kk), 0)
    Ah[[x]] <- p1 > p0; S[x] <- sum(pmax(p1-p0,0)) }
  dmin <- sum(w*vapply(1:J, function(x) min_defier_share(A[x,],B[x,]), 0))
  if (dmin > dbar) return(NULL)
  du <- eta_cond_joint_dual(S,A,B,w,kk,dbar,dual_box=box)
  if (is.na(du$value)) return(NULL)
  Ut <- du$u_t; Vt <- du$v_t; Pt <- du$p_t
  # influence function of the SAME functional on the full sample
  n <- nrow(br); ifv <- numeric(n); m1 <- m0 <- pih <- numeric(J); cont <- numeric(J)
  for (x in 1:J) { g <- gs[[x]]; s <- br[g,]
    hitA <- as.numeric(s$yb %in% yv[Ah[[x]]] & s$mc==kk); isk <- as.numeric(s$mc==kk)
    H1 <- Vt[x, s$mc] + Pt[x]*(hitA-isk); H0 <- Ut[x, s$mc] - Pt[x]*hitA
    t1 <- s$treat==1
    pih[x] <- mean(t1); m1[x] <- mean(H1[t1]); m0[x] <- mean(H0[!t1])
    cont[x] <- m1[x]+m0[x] }
  est <- sum(w*cont) - dbar*du$q
  for (x in 1:J) { g <- gs[[x]]; s <- br[g,]
    hitA <- as.numeric(s$yb %in% yv[Ah[[x]]] & s$mc==kk); isk <- as.numeric(s$mc==kk)
    H1 <- Vt[x, s$mc] + Pt[x]*(hitA-isk); H0 <- Ut[x, s$mc] - Pt[x]*hitA
    Di <- s$treat
    ifv[g] <- (Di/pih[x])*(H1-m1[x]) - 0 + ((1-Di)/(1-pih[x]))*(H0-m0[x]) + cont[x] - sum(w*cont) }
  se <- sd(ifv)/sqrt(n)
  c(bound=du$value, est=est, se=se, lcb=max(est-qnorm(.95)*se,0), dmin=dmin,
    ratio=du$value/se)
}
cat(sprintf("%-22s %2s %6s %6s %9s %9s %9s %7s\n","covariate","k","dbar","box","bound","se","LCB","b/se"))
for (xv in c("const","grandmother_baseline")) for (kk in c(3,4,5)) {
  dbar <- 0.021
  for (box in c(0.1,0.25,0.5,1,2)) {
    z <- oracle(xv,kk,dbar,box)
    if (is.null(z)) { cat(sprintf("%-22s %2d %6.3f %6.2f  infeasible\n", xv, kk, dbar, box)); next }
    cat(sprintf("%-22s %2d %6.3f %6.2f %9.5f %9.5f %9.5f %7.2f\n",
        xv, kk, dbar, box, z["bound"], z["se"], z["lcb"], z["ratio"]))
  }
  cat("\n")
}
