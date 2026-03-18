#' β < d/2 Smoothness Regime Simulations
#'
#' Tests theoretical β > d/2 condition for tree-based DML validity.
#' Explores performance across smoothness regimes: β ∈ {3, 2, 1} with d = 4.
#'
#' **Purpose:** Empirically characterize what happens when the smoothness
#' condition is satisfied (β=3 > d/2), at the boundary (β=2 = d/2), and
#' violated (β=1 < d/2). Constitution §9: exploratory stress testing with
#' honest reporting regardless of outcome.
#'
#' **Grid:**
#' - DGPs: beta_high (β=3), beta_boundary (β=2), beta_low (β=1)
#' - Methods: tree-DML, forest-DML, linear-DML (skip rashomon for speed)
#' - Sample sizes: n ∈ {400, 800, 1600}
#' - Replications: 500 per configuration
#' - Total: 3 DGPs × 3 methods × 3 sizes × 500 reps = 13,500 runs
#'
#' **Output:**
#' - results/beta_study_YYYY-MM-DD/simulation_results.rds
#' - results/beta_study_YYYY-MM-DD/summary_stats.csv
#'
#' **Runtime:** ~10-12 hours with 4-core parallelization
#'
#' **Output suppression:** All tree JSON and verbose output suppressed to prevent
#' massive log files (would be 500GB+ without suppression). Progress printed every
#' 100 simulations.
#'
#' **Three-way fidelity:**
#' - Paper Section 4 or Appendix C (supplementary simulations)
#' - DGPs constructed with piecewise polynomials (exact β control)
#' - Results demonstrate empirical β > d/2 boundary effects

# Load packages and functions
suppressMessages({
  library(parallel)
  library(dplyr)
  devtools::load_all("../../../optimaltrees")
})

# Load simulation helpers (prevents log bloat)
source("../simulation_helpers.R")

# Source doubletree functions SILENTLY
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/dml_att_repeated.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_beta_continuous.R",
  "methods/method_forest_dml.R",
  "methods/method_linear_dml.R"
), safe_source))

# Helper: Count leaves in a tree
count_tree_leaves <- function(tree) {
  if (is.null(tree) || !is.list(tree)) return(0)

  # If it's a leaf node (has name == "class" and no "feature")
  if (!is.null(tree$name) && tree$name == "class" && is.null(tree$feature)) {
    return(1)
  }

  # If it's an internal node, recurse on children
  left <- if (!is.null(tree$false)) count_tree_leaves(tree$false) else 0
  right <- if (!is.null(tree$true)) count_tree_leaves(tree$true) else 0

  return(left + right)
}

# Configuration
N_VALUES <- c(400, 800, 1600)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 500  # Production quality; use 50-100 for testing
SEED_OFFSET <- 20000  # Different offset from primary sims

# DGP functions (beta regime study)
DGPS <- list(
  beta_high = generate_dgp_beta_high,       # β=3 > d/2 (control)
  beta_boundary = generate_dgp_beta_boundary,  # β=2 = d/2 (boundary)
  beta_low = generate_dgp_beta_low          # β=1 < d/2 (FAILS condition)
)

# Methods (skip rashomon for speed)
METHODS <- c("tree", "forest", "linear")

# Parallelization - reduced to 2 cores to limit memory usage
N_CORES <- 2
cat(sprintf("Using %d cores for parallelization (memory-limited mode)\n\n", N_CORES))

# Create output directory
output_dir <- sprintf("results/beta_study_%s", Sys.Date())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create simulation grid
sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(strrep("=", 70), "\n")
cat("β < d/2 Smoothness Regime Simulation Study\n")
cat(strrep("=", 70), "\n\n")

cat("Simulation grid:\n")
cat(sprintf("  DGPs: %s\n", paste(names(DGPS), collapse = ", ")))
cat(sprintf("  Methods: %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("  Sample sizes: %s\n", paste(N_VALUES, collapse = ", ")))
cat(sprintf("  Replications: %d per configuration\n", N_REPS))
cat(sprintf("  Total simulations: %d\n\n", nrow(sim_grid)))

cat(sprintf("Estimated runtime: %.1f hours (with %d cores)\n\n",
            nrow(sim_grid) * 3 / 3600 / N_CORES, N_CORES))

# Simulation function for single replication
run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  # Progress tracking: console only, infrequent (prevents log bloat)
  progress_msg(sim_id, nrow(grid), every = 100, label = "Simulations")

  # Memory monitoring every 10 simulations
  if (sim_id %% 10 == 0) {
    monitor_memory(threshold_mb = 12000, verbose = FALSE)
  }

  # Generate data
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  # Fit model based on method - ALL wrapped in suppress_all()
  result <- tryCatch({

    if (row$method == "tree") {
      # Tree-DML (fold-specific regularization)
      # NEW: Pass continuous X with adaptive discretization
      # Suppress ALL output to prevent massive logs
      fit <- suppress_all({
        estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = FALSE,
          verbose = FALSE,
          discretize_method = "quantiles",   # Quantile-based binning
          discretize_bins = "adaptive"       # b_n = max(2, ceil(log(n)/3))
        )
      })

      # Extract tree complexity (average n_leaves across folds)
      # fit$nuisance_fits is a list of K elements, each with e_model and m0_model
      n_leaves_e <- if (!is.null(fit$nuisance_fits) && length(fit$nuisance_fits) > 0) {
        leaves_vec <- sapply(fit$nuisance_fits, function(fold_fit) {
          if (!is.null(fold_fit$e_model) && !is.null(fold_fit$e_model$model) &&
              !is.null(fold_fit$e_model$model$tree_json)) {
            count_tree_leaves(fold_fit$e_model$model$tree_json)
          } else {
            NA
          }
        })
        mean(leaves_vec, na.rm = TRUE)
      } else {
        NA
      }

      n_leaves_m0 <- if (!is.null(fit$nuisance_fits) && length(fit$nuisance_fits) > 0) {
        leaves_vec <- sapply(fit$nuisance_fits, function(fold_fit) {
          if (!is.null(fold_fit$m0_model) && !is.null(fold_fit$m0_model$model) &&
              !is.null(fold_fit$m0_model$model$tree_json)) {
            count_tree_leaves(fold_fit$m0_model$model$tree_json)
          } else {
            NA
          }
        })
        mean(leaves_vec, na.rm = TRUE)
      } else {
        NA
      }

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        n_leaves_e = n_leaves_e,       # NEW: propensity tree complexity
        n_leaves_m0 = n_leaves_m0      # NEW: outcome tree complexity
      )

    } else if (row$method == "forest") {
      # Forest-DML (ranger)
      fit <- suppress_all({
        att_forest(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          num.trees = 500,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        n_leaves_e = NA,    # Not applicable for forest
        n_leaves_m0 = NA
      )

    } else if (row$method == "linear") {
      # Linear-DML (logistic regression)
      fit <- suppress_all({
        att_linear(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          interactions = FALSE,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        n_leaves_e = NA,    # Not applicable for linear
        n_leaves_m0 = NA
      )

    } else {
      stop("Unknown method: ", row$method)
    }

  }, error = function(e) {
    list(
      theta = NA,
      sigma = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      n_leaves_e = NA,
      n_leaves_m0 = NA
    )
  })

  # Return results with metadata
  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = d$true_att,
    beta = d$diagnostics$beta,
    rate_regime = d$diagnostics$rate_regime,
    theoretical_sn = d$diagnostics$theoretical_sn,  # NEW: for s_n verification
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    n_leaves_e = result$n_leaves_e,      # NEW: fitted tree complexity
    n_leaves_m0 = result$n_leaves_m0,    # NEW: fitted tree complexity
    stringsAsFactors = FALSE
  )
}

# Run simulations in parallel
cat("Starting simulations...\n")
start_time <- Sys.time()

results_list <- mclapply(
  sim_grid$sim_id,
  function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET),
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Combine results
results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\nSimulations complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

# Save full results (atomic write)
safe_save(results, file.path(output_dir, "simulation_results.rds"))

# Compute summary statistics by β regime
summary_stats <- results %>%
  filter(converged) %>%
  group_by(dgp, beta, rate_regime, method, n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att, na.rm = TRUE),
    rmse = sqrt(mean((theta - true_att)^2, na.rm = TRUE)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    median_ci_width = median(ci_upper - ci_lower, na.rm = TRUE),
    .groups = "drop"
  )

# Save summary statistics
write.csv(summary_stats,
          file.path(output_dir, "summary_stats.csv"),
          row.names = FALSE)
cat(sprintf("Summary stats saved to: %s\n", file.path(output_dir, "summary_stats.csv")))

# Print summary
cat("\n")
cat(strrep("=", 70), "\n")
cat("Summary Statistics by β Regime\n")
cat("=" %R% 70, "\n\n")
print(summary_stats, n = Inf)

# Quality checks: Focus on β regime differences
cat("\n")
cat(strrep("=", 70), "\n")
cat("β Regime Analysis\n")
cat("=" %R% 70, "\n\n")

# Check if beta_low shows degradation vs beta_high
beta_comparison <- summary_stats %>%
  filter(method == "tree", n == 800) %>%
  select(dgp, beta, rate_regime, coverage, rmse, mean_ci_width)

cat("Tree-DML performance at n=800 by β regime:\n\n")
print(beta_comparison)

# Statistical test: Is coverage in beta_low significantly different from 95%?
beta_low_results <- results %>%
  filter(converged, dgp == "beta_low", method == "tree", n == 800)

if (nrow(beta_low_results) > 0) {
  coverage_beta_low <- mean(beta_low_results$ci_lower <= beta_low_results$true_att &
                            beta_low_results$ci_upper >= beta_low_results$true_att)
  n_beta_low <- nrow(beta_low_results)

  # Z-test for proportion
  z_stat <- (coverage_beta_low - 0.95) / sqrt(0.95 * 0.05 / n_beta_low)
  p_value <- 2 * pnorm(-abs(z_stat))

  cat(sprintf("\n\nStatistical test: β=1 < d/2 coverage vs 95%%:\n"))
  cat(sprintf("  Observed coverage: %.3f\n", coverage_beta_low))
  cat(sprintf("  Z-statistic: %.2f\n", z_stat))
  cat(sprintf("  P-value: %.4f\n", p_value))
  if (p_value < 0.05) {
    cat(sprintf("  ✓ Significantly different from 95%% (p < 0.05)\n"))
  } else {
    cat(sprintf("  ✗ Not significantly different from 95%% (p ≥ 0.05)\n"))
  }
}

# Compare beta_high vs beta_low
beta_high_results <- results %>%
  filter(converged, dgp == "beta_high", method == "tree", n == 800)

if (nrow(beta_high_results) > 0 && nrow(beta_low_results) > 0) {
  coverage_beta_high <- mean(beta_high_results$ci_lower <= beta_high_results$true_att &
                             beta_high_results$ci_upper >= beta_high_results$true_att)

  cat(sprintf("\n\nComparison: β=3 (high) vs β=1 (low) at n=800:\n"))
  cat(sprintf("  Coverage β=3: %.3f\n", coverage_beta_high))
  cat(sprintf("  Coverage β=1: %.3f\n", coverage_beta_low))
  cat(sprintf("  Difference: %.3f (%.1f%% relative)\n",
              coverage_beta_high - coverage_beta_low,
              100 * (coverage_beta_high - coverage_beta_low) / coverage_beta_high))

  # Two-proportion z-test
  n_high <- nrow(beta_high_results)
  n_low <- nrow(beta_low_results)
  successes_high <- sum(beta_high_results$ci_lower <= beta_high_results$true_att &
                        beta_high_results$ci_upper >= beta_high_results$true_att)
  successes_low <- sum(beta_low_results$ci_lower <= beta_low_results$true_att &
                       beta_low_results$ci_upper >= beta_low_results$true_att)

  pooled_p <- (successes_high + successes_low) / (n_high + n_low)
  se_diff <- sqrt(pooled_p * (1 - pooled_p) * (1/n_high + 1/n_low))
  z_diff <- (coverage_beta_high - coverage_beta_low) / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))

  cat(sprintf("  Z-statistic: %.2f\n", z_diff))
  cat(sprintf("  P-value: %.4f\n", p_diff))
  if (p_diff < 0.05) {
    cat(sprintf("  ✓ Significantly different (p < 0.05)\n"))
  } else {
    cat(sprintf("  ✗ Not significantly different (p ≥ 0.05)\n"))
  }
}

# Flag any concerning patterns
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Quality Checks\n")
cat("=" %R% 70, "\n\n")

# Check convergence by regime
conv_by_regime <- results %>%
  group_by(dgp, rate_regime) %>%
  summarize(
    convergence_rate = mean(converged),
    .groups = "drop"
  )

cat("Convergence rates by β regime:\n")
print(conv_by_regime)

if (any(conv_by_regime$convergence_rate < 0.98)) {
  cat("\n⚠️  Some regimes have convergence < 98%\n")
} else {
  cat("\n✓ All regimes achieve ≥ 98% convergence\n")
}

# Overall summary
cat("\n")
cat(strrep("=", 70), "\n")
cat("Simulation Complete\n")
cat("=" %R% 70, "\n\n")

cat("Next steps:\n")
cat("  1. Run analyze_beta_study.R to generate figures and tables\n")
cat("  2. Review results for manuscript Section 4 or Appendix C\n")
cat("  3. Update session notes with key findings\n\n")

cat(sprintf("Output directory: %s\n", output_dir))

# Final check: warn if large files were created
check_large_files(output_dir, min_mb = 10)
