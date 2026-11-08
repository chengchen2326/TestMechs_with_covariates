suppressMessages({library(haven);library(purrr);library(Rglpk)})
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_lb_frac_affected_x.R","src_simulate_data_multiM_x.R")) source(file.path("../src",f))
source("calib_C.R")

cat("=== Baranov calibration, K=5, Q=5, maximal coupling ===\n")
bc <- baranov_cells(K_target=5, Q=5, coupling="maximal")
cat(sprintf("  K = %d, Q = %d, cells = %d\n", bc$K, bc$Q, length(bc$cells)))
cat(sprintf("  theta diagonal (always-taker shares): %s\n",
    paste(sprintf("%.3f", diag(bc$cells[[1]]$theta)), collapse=" ")))

cat("\n=== CHECK 1: bound <= truth (validity), and slack ===\n")
cat(sprintf("%4s %12s %12s %12s %10s\n","k","eta_cond","eta_pooled","true_prod","slack"))
for (k in 1:bc$K) {
  pb <- population_bound_multiM_x(bc$cells, bc$p_x, at_group=k)
  cat(sprintf("%4d %12.6f %12.6f %12.6f %10.6f%s\n", k, pb$eta_conditional,
      pb$eta_pooled, pb$true_product, pb$slack_conditional,
      if (pb$slack_conditional < -1e-9) "  <-- INVALID" else ""))
}

cat("\n=== CHECK 2: population bound vs lb_frac_affected on a large sample ===\n")
set.seed(99)
df <- simulate_data_multiM_x(bc$cells, bc$p_x, n=400000, p_d=0.5)
df$const <- 1L
cat(sprintf("  simulated P(M=m|D=0) = %s\n", paste(sprintf("%.4f", prop.table(table(df$m[df$d==0]))), collapse=" ")))
cat(sprintf("  simulated P(M=m|D=1) = %s\n", paste(sprintf("%.4f", prop.table(table(df$m[df$d==1]))), collapse=" ")))
cat(sprintf("%4s %14s %14s %12s %12s\n","k","population eta","lb_frac_aff nu","eta/theta","|diff|"))
for (k in 1:bc$K) {
  pb <- population_bound_multiM_x(bc$cells, bc$p_x, at_group=k)
  nu_pop <- if (sum(bc$p_x*pb$theta_by_cell) > 0) pb$eta_pooled/sum(bc$p_x*pb$theta_by_cell) else NA
  a <- tryCatch(as.numeric(lb_frac_affected(df, d="d", m="m", y="y", at_group=k)),
                error=function(e) NA_real_, warning=function(w) NA_real_)
  cat(sprintf("%4d %14.6f %14.6f %12.6f %12.2e\n", k, pb$eta_pooled, a, nu_pop, abs(a-nu_pop)))
}

cat("\n=== CHECK 3: independent coupling should leave the bound slack ===\n")
bi <- baranov_cells(K_target=5, Q=5, coupling="independent", verbose=FALSE)
cat(sprintf("%4s %12s %12s %10s\n","k","eta_cond","true_prod","slack"))
for (k in 1:bi$K) {
  pb <- population_bound_multiM_x(bi$cells, bi$p_x, at_group=k)
  cat(sprintf("%4d %12.6f %12.6f %10.6f\n", k, pb$eta_conditional, pb$true_product, pb$slack_conditional))
}
