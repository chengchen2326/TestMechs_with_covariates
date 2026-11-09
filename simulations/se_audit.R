# Is the reported standard error right in the no-complier boundary case?
# Direct test: compare the empirical dispersion of eta_hat across replications
# with the mean of the standard error the estimator reports. If the reported SE
# sits below the empirical SD, the influence function is missing a term.
suppressMessages({library(haven);library(purrr);library(Rglpk)})
SRC <- "../src"
for (f in c("src_remove_missing_from_df.R","src_lb_frac_affected.R","src_discretize_y.R",
            "src_lb_frac_affected_x.R","src_simulate_data_multiM_x.R")) source(file.path(SRC,f))
source("calib_C.R"); source("violation_scale.R"); source("design_C.R")

REPS <- as.integer(Sys.getenv("TM_REPS","2000"))
cat(sprintf("%-9s %2s %6s %-7s %10s %10s %8s %8s\n",
            "config","k","n","arm","emp sd","mean SE","ratio","size"))
for (cfg in c("signflip","homog")) {
  B <- build_config(cfg, "maximal", 5, 5, t = 0)
  for (k in c(3,5)) for (n in c(1000, 4000, 16000)) {
    set.seed(4242)
    M <- matrix(NA_real_, REPS, 4)
    for (r in seq_len(REPS)) {
      df <- simulate_data_multiM_x(B$cells, B$p_x, n = n, p_d = 0.5)
      df$const <- 1L
      a <- try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=k,method="split"), silent=TRUE)
      b <- try(lb_frac_affected_x(df,"d","m","y",x="const",at_group=k,method="pooled_split"), silent=TRUE)
      if (!inherits(a,"try-error")) M[r,1:2] <- c(a$eta, a$eta_se)
      if (!inherits(b,"try-error")) M[r,3:4] <- c(b$eta, b$eta_se)
    }
    for (j in c(1,3)) {
      nm <- if (j==1) "ours" else "pooled"
      es <- sd(M[,j], na.rm=TRUE); ms <- mean(M[,j+1], na.rm=TRUE)
      lcb <- pmax(M[,j] - qnorm(.95)*M[,j+1], 0)
      cat(sprintf("%-9s %2d %6d %-7s %10.6f %10.6f %8.3f %8.4f%s\n",
          cfg, k, n, nm, es, ms, ms/es, mean(lcb > 0, na.rm=TRUE),
          if (ms/es < 0.95) "  <-- SE too small" else ""))
    }
  }
}
