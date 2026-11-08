suppressMessages({library(haven);library(Rglpk)})
for (f in c("src_simulate_data_multiM_x.R","src_theta_kk_min_defiers.R","src_eta_cond_joint.R")) source(file.path("../src",f))
load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
keep <- c("treat","relationship_husb","motherfinancial","grandmother_baseline")
br <- br[stats::complete.cases(br[,keep]),]
for (v in keep) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
mv <- sort(unique(br$relationship_husb)); br$mc <- match(br$relationship_husb, mv)
K <- length(mv); yv <- sort(unique(br$yb)); kk <- 4; dbar <- 0.05

gs <- split(seq_len(nrow(br)), br$grandmother_baseline); J <- length(gs)
w <- as.numeric(vapply(gs,length,1L)); w <- w/sum(w)
A <- B <- matrix(0,J,K); S <- numeric(J); Ah <- vector("list",J)
for (x in 1:J) { s <- br[gs[[x]],]
  t1 <- s$treat==1; t0 <- s$treat==0
  A[x,] <- vapply(1:K, function(v) mean(s$mc[t0]==v), 0)
  B[x,] <- vapply(1:K, function(v) mean(s$mc[t1]==v), 0)
  p1 <- vapply(yv, function(v) mean(s$yb[t1]==v & s$mc[t1]==kk), 0)
  p0 <- vapply(yv, function(v) mean(s$yb[t0]==v & s$mc[t0]==kk), 0)
  Ah[[x]] <- p1 > p0; S[x] <- sum(pmax(p1-p0,0)) }

pr <- eta_cond_joint(S,A,B,w,kk,dbar)
du <- eta_cond_joint_dual(S,A,B,w,kk,dbar,dual_box=2)
cat(sprintf("primal %.8f | dual %.8f | gap %.2e\n", pr$value, du$value, abs(pr$value-du$value)))

# Reconstruct the dual objective from the raw duals, as a sanity check
obj <- sum(vapply(1:J, function(x) sum(A[x,]*du$u[x,]) + sum(B[x,]*du$v[x,]) +
                                   du$p[x]*(S[x]-B[x,kk]), 0)) - dbar*du$q
cat(sprintf("objective rebuilt from raw duals: %.8f\n", obj))

# Now the normalised / transformed-outcome form, evaluated on the SAME data
Ut <- du$u_t; Vt <- du$v_t; Pt <- du$p_t
tot <- 0
for (x in 1:J) { s <- br[gs[[x]],]
  hitA <- as.numeric(s$yb %in% yv[Ah[[x]]] & s$mc==kk)
  isk  <- as.numeric(s$mc==kk)
  t1 <- s$treat==1; t0 <- s$treat==0
  H1 <- Vt[x, s$mc[t1]] + Pt[x]*(hitA[t1]-isk[t1])
  H0 <- Ut[x, s$mc[t0]] - Pt[x]*hitA[t0]
  tot <- tot + w[x]*(mean(H1) + mean(H0))
  cat(sprintf("  cell %d: E[H1|D=1]=%.6f  E[H0|D=0]=%.6f  sum=%.6f\n", x, mean(H1), mean(H0), mean(H1)+mean(H0)))
}
cat(sprintf("transformed-outcome form: %.8f  (should equal dual + dbar*q = %.8f)\n",
            tot - dbar*du$q, du$value))
cat(sprintf("  dbar*q = %.8f,  q = %.6f,  p_t = %s\n", dbar*du$q, du$q, paste(sprintf("%.4f",Pt),collapse=" ")))
cat(sprintf("  range u_t = [%.3f, %.3f], v_t = [%.3f, %.3f]\n", min(Ut),max(Ut),min(Vt),max(Vt)))
