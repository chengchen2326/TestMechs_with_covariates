# Exact check: feed lb_frac_affected the population distribution as weights,
# so the comparison is machine precision rather than Monte Carlo error.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_simulate_data_multiM_x.R")) source(file.path("../src",f))
source("calib_C.R")

exact_df <- function(cells, p_x, K, Q, p_d = 0.5) {
  obs <- function(cj) {
    p0 <- matrix(0,Q,K); p1 <- matrix(0,Q,K)
    for (mm in 1:K) {
      p0[,mm] <- p0[,mm] + cj$theta[mm,mm]*rowSums(cj$J[[mm]])
      p1[,mm] <- p1[,mm] + cj$theta[mm,mm]*colSums(cj$J[[mm]])
      if (mm<K) for (kk in (mm+1):K) if (cj$theta[mm,kk]>0) p0[,mm] <- p0[,mm]+cj$theta[mm,kk]*cj$g0[[mm,kk]]
      if (mm>1) for (l in 1:(mm-1))  if (cj$theta[l,mm]>0) p1[,mm] <- p1[,mm]+cj$theta[l,mm]*cj$g1[[l,mm]]
    }
    list(p0=p0,p1=p1)
  }
  O <- lapply(cells, obs)
  p0 <- Reduce(`+`, Map(function(w,o) w*o$p0, p_x, O))
  p1 <- Reduce(`+`, Map(function(w,o) w*o$p1, p_x, O))
  g <- expand.grid(y=1:Q, m=1:K, d=0:1)
  g$w <- ifelse(g$d==1, p1[cbind(g$y,g$m)]*p_d, p0[cbind(g$y,g$m)]*(1-p_d))
  g[g$w > 1e-14, ]
}

for (cpl in c("maximal","independent")) {
  bc <- baranov_cells(K_target=5, Q=5, coupling=cpl, verbose=FALSE)
  ed <- exact_df(bc$cells, bc$p_x, bc$K, bc$Q)
  cat(sprintf("\n--- %s coupling: exact weighted comparison (support points: %d) ---\n", cpl, nrow(ed)))
  cat(sprintf("%4s %16s %16s %12s\n","k","population nu","lb_frac_aff nu","|diff|"))
  for (k in 1:bc$K) {
    pb <- population_bound_multiM_x(bc$cells, bc$p_x, at_group=k)
    th <- sum(bc$p_x*pb$theta_by_cell)
    nu_pop <- if (th > 1e-12) pb$eta_pooled/th else NA_real_
    a <- tryCatch(as.numeric(lb_frac_affected(ed, d="d", m="m", y="y", w="w", at_group=k)),
                  error=function(e) NA_real_, warning=function(w) NA_real_)
    cat(sprintf("%4d %16.10f %16.10f %12.2e\n", k, nu_pop, a, abs(a-nu_pop)))
  }
}
