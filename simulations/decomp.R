# Numerical check of the two-channel decomposition of the conditioning gain.
#
#   S(x) = sum_y ( p1(y,m_k|x) - p0(y,m_k|x) )_+
#   C(x) = P(M >= m_k | D=1, x) - P(M >= m_k | D=0, x)
#   eta_k(x) = ( S(x) - C(x) )_+
#
#   gain = E[eta_k(X)] - eta_k^pooled
#        = [ E[(S-C)_+] - (E S - E C)_+ ]        = G_B, outer truncation
#        + [ (E S - E C)_+ - (S_pool - E C)_+ ]  = G_A, outcome convexity
#
# with S_pool = sum_y ( E[p1 - p0] )_+ <= E[S].  Both terms are non-negative.
# Corollary: if m_k is the smallest value in the support of M then C == 0, so
# G_B == 0 and the whole gain must come from G_A.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R",
            "src_simulate_data_multiM_x.R")) source(file.path("../src",f))
source("calib_C.R")

obs_of <- function(cj, K, Q) {
  p0 <- matrix(0,Q,K); p1 <- matrix(0,Q,K)
  for (mm in 1:K) {
    p0[,mm] <- p0[,mm] + cj$theta[mm,mm]*rowSums(cj$J[[mm]])
    p1[,mm] <- p1[,mm] + cj$theta[mm,mm]*colSums(cj$J[[mm]])
    if (mm<K) for (kk in (mm+1):K) if (cj$theta[mm,kk]>0) p0[,mm] <- p0[,mm]+cj$theta[mm,kk]*cj$g0[[mm,kk]]
    if (mm>1) for (l in 1:(mm-1))  if (cj$theta[l,mm]>0) p1[,mm] <- p1[,mm]+cj$theta[l,mm]*cj$g1[[l,mm]]
  }
  list(p0=p0,p1=p1)
}

decompose <- function(P0, P1, w, k) {
  # P0, P1: lists of Q x K matrices, one per cell;  w: cell weights
  K <- ncol(P0[[1]])
  S <- vapply(seq_along(P0), function(j) sum(pmax(P1[[j]][,k]-P0[[j]][,k],0)), numeric(1))
  C <- vapply(seq_along(P0), function(j)
         sum(P1[[j]][,k:K,drop=FALSE]) - sum(P0[[j]][,k:K,drop=FALSE]), numeric(1))
  ES <- sum(w*S); EC <- sum(w*C)
  d_pool <- Reduce(`+`, Map(function(wi,a,b) wi*(b[,k]-a[,k]), w, P0, P1))
  S_pool <- sum(pmax(d_pool,0))
  eta_cond <- sum(w*pmax(S-C,0))
  eta_pool <- max(S_pool-EC,0)
  G_B <- eta_cond - max(ES-EC,0)
  G_A <- max(ES-EC,0) - eta_pool
  list(eta_cond=eta_cond, eta_pool=eta_pool, gain=eta_cond-eta_pool,
       G_A=G_A, G_B=G_B, S=S, C=C, ES=ES, EC=EC, S_pool=S_pool)
}

cat("=== Design C, split configuration (cells = grandmother_baseline) ===\n")
bc <- baranov_cells(K_target=5, Q=5, xvar="grandmother_baseline", coupling="maximal", verbose=FALSE)
O <- lapply(bc$cells, obs_of, K=bc$K, Q=bc$Q)
P0 <- lapply(O, `[[`, "p0"); P1 <- lapply(O, `[[`, "p1")
cat(sprintf("%3s %9s %9s %9s %9s %9s %10s %8s %8s\n","k","eta_cond","eta_pool","gain","G_A","G_B","resid","E[S]","E[C]"))
for (k in 1:bc$K) {
  z <- decompose(P0, P1, bc$p_x, k)
  pb <- population_bound_multiM_x(bc$cells, bc$p_x, at_group=k)
  resid <- (z$G_A + z$G_B) - z$gain
  chk <- abs(z$eta_cond - pb$eta_conditional) + abs(z$eta_pool - pb$eta_pooled)
  cat(sprintf("%3d %9.5f %9.5f %9.5f %9.5f %9.5f %10.1e %8.4f %8.4f%s\n",
      k, z$eta_cond, z$eta_pool, z$gain, z$G_A, z$G_B, resid, z$ES, z$EC,
      if (chk > 1e-12) sprintf("  [mismatch %.1e]", chk) else ""))
}

cat("\n  per-cell S(x) and C(x) at k = 3 and k = 5:\n")
for (k in c(3,5)) { z <- decompose(P0,P1,bc$p_x,k)
  cat(sprintf("    k=%d:  S = %s   C = %s   S-C = %s\n", k,
      paste(sprintf("%.4f",z$S),collapse=" "), paste(sprintf("%.4f",z$C),collapse=" "),
      paste(sprintf("%+.4f",z$S-z$C),collapse=" "))) }

cat("\n=== Corollary: k at the bottom of the support => C == 0 and G_B == 0 ===\n")
emp_cells <- function(df, dv, mv, yv, xv) {
  mvals <- sort(unique(df[[mv]])); yvals <- sort(unique(df[[yv]]))
  gs <- split(seq_len(nrow(df)), df[[xv]]); P0 <- P1 <- list(); w <- numeric(0)
  for (g in gs) { s <- df[g,]
    f <- function(d) { M <- matrix(0,length(yvals),length(mvals))
      ss <- s[s[[dv]]==d,]
      for (a in seq_along(yvals)) for (b in seq_along(mvals))
        M[a,b] <- mean(ss[[yv]]==yvals[a] & ss[[mv]]==mvals[b]); M }
    P0[[length(P0)+1]] <- f(0); P1[[length(P1)+1]] <- f(1); w <- c(w, nrow(s)) }
  list(P0=P0, P1=P1, w=w/sum(w))
}
load("/home/claude/repo/data/burstzyn_data.rda"); bz <- get(ls(pattern="burstzyn|bursztyn")[1])
bz <- bz[!is.na(bz$index) & !is.na(bz$college_deg),]
e <- emp_cells(bz,"condition2","signed_up_number","applied_out_fl","college_deg")
for (k in 1:2) { z <- decompose(e$P0,e$P1,e$w,k)
  cat(sprintf("  Bursztyn, k=%d (m=%d): gain=%.6f  G_A=%.6f  G_B=%.6f  E[C]=%.6f\n",
      k, k-1, z$gain, z$G_A, z$G_B, z$EC)) }
load("/home/claude/repo/data/baranov_data.rda"); br <- haven::zap_labels(mother_data)
br <- br[!is.na(br$dep_sample)&br$dep_sample==1,]
br <- br[stats::complete.cases(br[,c("treat","grandmother","motherfinancial","grandmother_baseline")]),]
for (v in c("treat","grandmother","motherfinancial","grandmother_baseline")) br[[v]] <- as.numeric(unclass(br[[v]]))
br$yb <- as.integer(br$motherfinancial > median(br$motherfinancial))
e2 <- emp_cells(br,"treat","grandmother","yb","grandmother_baseline")
for (k in 1:2) { z <- decompose(e2$P0,e2$P1,e2$w,k)
  cat(sprintf("  Baranov,  k=%d (m=%d): gain=%.6f  G_A=%.6f  G_B=%.6f  E[C]=%.6f\n",
      k, k-1, z$gain, z$G_A, z$G_B, z$EC)) }
