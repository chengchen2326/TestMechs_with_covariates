# =====================================================================
# Design C: multivalued mediator, calibrated to Baranov.
#
# The upgrade over Design A is that the compliance type is drawn explicitly, so
# the TRUE theta_kk * nu_k is known. Coverage can therefore be checked against
# the quantity of interest and not only against the bound, and the slack of the
# bound can be measured.
#
# Configurations:
#   homog     one Baranov-calibrated cell replicated, so the cells are
#             identical. Conditioning should gain nothing; this measures what it
#             costs.
#   split     Baranov-calibrated cells by a real covariate: whatever
#             heterogeneity the data happen to have.
#   signflip  two cells whose Y | D = 1 laws are perturbed in opposite
#             directions, so the per-cell bound is positive while the pooled
#             bound is exactly zero. The multivalued analogue of Design A's
#             sign-changing calibration.
#
# Couplings:
#   maximal      always-takers' potential outcomes are coupled to maximise
#                agreement, so nu_k is the total variation distance and KR's
#                bound is attained where no compliers enter m_k.
#   independent  nu_k is larger, so the bound is slack. Reported because a
#                procedure can only be judged against the truth, and the truth
#                is usually not at the bound.
# =====================================================================

suppressMessages({library(haven); library(purrr); library(Rglpk)})
SRC <- Sys.getenv("TM_SRC", "../src")
for (f in c("src_remove_missing_from_df.R", "src_lb_frac_affected.R",
            "src_discretize_y.R", "src_lb_frac_affected_x.R",
            "src_simulate_data_multiM_x.R")) source(file.path(SRC, f))
source("calib_C.R")

# ---------------------------------------------------- sign-flipping cells
# Shift mass between the extreme Y bins of the treated law, in opposite
# directions in the two cells. The M marginals and theta are untouched, so
# monotonicity is unaffected; averaging the two cells restores the control law
# exactly, so the pooled bound is zero by construction.
signflip_cells <- function(base, k, delta = 0.08, coupling = "maximal") {
  cpl <- if (coupling == "maximal") maximal_coupling else independent_coupling
  K <- nrow(base$theta); Q <- nrow(base$J[[1]])
  L0 <- lapply(1:K, function(mm) rowSums(base$J[[mm]]) / sum(base$J[[mm]]))
  # Move mass between the two HEAVIEST bins of Y | M = m_k, D = 0. The
  # perturbation is applied with both signs, so each of the two bins must be
  # able to give up delta as well as receive it: feasibility is bounded by the
  # smaller of the two masses. Picking the extreme or the lightest bin would
  # therefore cap delta at whatever little mass that bin carries -- 0.05 in
  # this calibration -- and leave the design with almost no signal.
  ord <- order(L0[[k]], decreasing = TRUE)
  a <- ord[1]; b <- ord[2]
  v <- numeric(Q); v[a] <- 1; v[b] <- -1          # zero-sum perturbation
  feas <- min(L0[[k]][a], L0[[k]][b], 1 - L0[[k]][a], 1 - L0[[k]][b])
  delta_eff <- min(delta, 0.98 * feas)
  if (delta_eff < delta)
    message(sprintf("signflip: delta capped from %.3f to %.3f by the support of Y | M = m_%d, D = 0",
                    delta, delta_eff, k))
  out <- list()
  for (sgn in c(1, -1)) {
    cj <- base
    L1k <- L0[[k]] + sgn * delta_eff * v
    stopifnot(all(L1k >= -1e-12), abs(sum(L1k) - 1) < 1e-10)
    L1k <- pmax(L1k, 0); L1k <- L1k / sum(L1k)
    cj$J[[k]] <- cpl(L0[[k]], L1k)
    for (l in 1:K) cj$g1[[l, k]] <- L1k
    out[[length(out) + 1]] <- cj
  }
  out
}

build_config <- function(config, coupling, K_target = 5, Q = 5,
                         xvar = "grandmother_baseline", delta = 0.20) {
  if (config == "split") {
    bc <- baranov_cells(K_target, Q, xvar = xvar, coupling = coupling, verbose = FALSE)
    return(bc)
  }
  bc <- baranov_cells(K_target, Q, xvar = NULL, coupling = coupling, verbose = FALSE)
  base <- bc$cells[[1]]
  if (config == "homog")
    return(list(cells = list(base, base), p_x = c(0.5, 0.5), K = bc$K, Q = bc$Q))
  if (config == "signflip") {
    # Remove compliers, so that the M | D marginals coincide and the pooled
    # bound is exactly zero rather than merely small: with theta diagonal
    # everyone is an always-taker, the complier-mass term vanishes, and
    # averaging the two cells' treated outcome laws restores the control law.
    base$theta <- diag(rowSums(base$theta))
    cl <- signflip_cells(base, k = K_target, delta = delta, coupling = coupling)
    return(list(cells = cl, p_x = c(0.5, 0.5), K = bc$K, Q = bc$Q))
  }
  stop("unknown config")
}

run_design_C <- function(configs = c("homog", "split", "signflip"),
                         couplings = c("maximal", "independent"),
                         focal_k = NULL, n_grid = c(1000, 4000),
                         n_reps = 300, K_target = 5, Q = 5, seed = 20261108) {
  set.seed(seed)
  rows <- list()
  for (cfg in configs) for (cpl in couplings) {
    B <- build_config(cfg, cpl, K_target, Q)
    ks <- if (is.null(focal_k)) B$K else focal_k
    for (k in ks) {
      pb <- population_bound_multiM_x(B$cells, B$p_x, at_group = k)
      cat(sprintf("\n[%s / %s] k = %d: E[eta(X)] = %.5f  eta_pooled = %.5f  truth = %.5f  slack = %.5f\n",
                  cfg, cpl, k, pb$eta_conditional, pb$eta_pooled,
                  pb$true_product, pb$slack_conditional))
      for (n in n_grid) {
        acc <- list()
        for (r in seq_len(n_reps)) {
          df <- simulate_data_multiM_x(B$cells, B$p_x, n = n, p_d = 0.5)
          df$const <- 1L
          z <- list(
            ours     = try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=k,method="split"), silent=TRUE),
            pooled   = try(lb_frac_affected_x(df,"d","m","y",x="const",at_group=k,method="pooled_split"), silent=TRUE),
            plugin_x = try(lb_frac_affected_x(df,"d","m","y",x="x",at_group=k,method="plugin"), silent=TRUE),
            plugin_0 = try(lb_frac_affected_x(df,"d","m","y",x="const",at_group=k,method="plugin"), silent=TRUE))
          for (nm in names(z)) { if (inherits(z[[nm]],"try-error")) next
            acc[[nm]] <- rbind(acc[[nm]], data.frame(eta=z[[nm]]$eta, lcb=z[[nm]]$eta_lcb,
                                 defier=z[[nm]]$max_cell_defier_share)) }
        }
        for (nm in names(acc)) {
          A <- acc[[nm]]
          bound <- if (nm %in% c("pooled","plugin_0")) pb$eta_pooled else pb$eta_conditional
          rows[[length(rows)+1]] <- data.frame(
            config=cfg, coupling=cpl, k=k, n=n, arm=nm, reps=nrow(A),
            bound=bound, truth=pb$true_product,
            mean_eta=mean(A$eta), bias_vs_bound=mean(A$eta)-bound,
            mean_lcb=if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb, na.rm=TRUE),
            power=if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb>0, na.rm=TRUE),
            cover_bound=if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb<=bound+1e-12, na.rm=TRUE),
            cover_truth=if (all(is.na(A$lcb))) NA_real_ else mean(A$lcb<=pb$true_product+1e-12, na.rm=TRUE),
            max_defier=mean(A$defier, na.rm=TRUE))
        }
        cat(sprintf("    n = %5d done\n", n))
      }
    }
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  out <- run_design_C(n_reps = as.integer(Sys.getenv("TM_REPS","300")),
                      focal_k = as.integer(strsplit(Sys.getenv("TM_K","5"), ",")[[1]]))
  dir.create("results", showWarnings = FALSE)
  write.csv(out, "results/design_C.csv", row.names = FALSE)
  cat("\n")
  print(out[, c("config","coupling","k","n","arm","bound","truth","mean_lcb",
                "power","cover_bound","cover_truth")],
        row.names = FALSE, digits = 3)
}
