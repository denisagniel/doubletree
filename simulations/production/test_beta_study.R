#' Test Run: β < d/2 Smoothness Regime Simulations
#'
#' Runs a small test (10 reps per configuration) to verify:
#' 1. Continuous features pass through to optimaltrees correctly
#' 2. Adaptive discretization works
#' 3. Tree complexity (n_leaves) is captured in results
#' 4. Fitted s_n ≈ theoretical s_n (within reasonable range)
#'
#' Runtime: ~2-3 minutes

# Load packages and functions
library(parallel)
library(dplyr)

devtools::load_all("../../../optimaltrees")

source("../../R/estimate_att.R")
source("../../R/dml_att_repeated.R")
source("../../R/nuisance_trees.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/utils.R")

source("dgps/dgps_beta_continuous.R")  # NEW: continuous features
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

# Helper: Count leaves in a tree
count_tree_leaves <- function(tree) {
  if (is.null(tree) || !is.list(tree)) return(0)

  # If it's a leaf node
  if (!is.null(tree$name) && tree$name == "class" && is.null(tree$feature)) {
    return(1)
  }

  # If it's an internal node, recurse
  left <- if (!is.null(tree$false)) count_tree_leaves(tree$false) else 0
  right <- if (!is.null(tree$true)) count_tree_leaves(tree$true) else 0

  return(left + right)
}

# Test configuration
N_VALUES <- c(800)  # Single sample size for quick test
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 50  # 50 reps for better coverage estimates (was 10)
SEED_OFFSET <- 20000

# DGP functions
DGPS <- list(
  beta_high = generate_dgp_beta_high,
  beta_boundary = generate_dgp_beta_boundary,
  beta_low = generate_dgp_beta_low
)

# Test tree method only for speed
METHODS <- c("tree")

# Single core for debugging
N_CORES <- 1

cat(strrep("=", 70), "\n")
cat("TEST RUN: β < d/2 Smoothness Regime Simulation\n")
cat(strrep("=", 70), "\n\n")

cat("Test configuration:\n")
cat(sprintf("  DGPs: %s\n", paste(names(DGPS), collapse = ", ")))
cat(sprintf("  Methods: %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("  Sample sizes: %s\n", paste(N_VALUES, collapse = ", ")))
cat(sprintf("  Replications: %d per configuration\n", N_REPS))
cat(sprintf("  Total simulations: %d\n\n", length(DGPS) * length(METHODS) * length(N_VALUES) * N_REPS))

# Create simulation grid
sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

# Simulation function
run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  result <- tryCatch({
    # Tree-DML with continuous features + adaptive discretization
    capture.output({
      fit <- estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        K = k_folds,
        regularization = log(row$n) / row$n,
        cv_regularization = FALSE,
        use_rashomon = FALSE,
        verbose = FALSE,
        discretize_method = "quantiles",
        discretize_bins = "adaptive"
      )
    }, file = tempfile())

    # Extract tree complexity (average across folds)
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
      n_leaves_e = n_leaves_e,
      n_leaves_m0 = n_leaves_m0
    )

  }, error = function(e) {
    cat(sprintf("ERROR in sim %d: %s\n", sim_id, e$message))
    list(
      theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, n_leaves_e = NA, n_leaves_m0 = NA
    )
  })

  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = d$true_att,
    beta = d$diagnostics$beta,
    rate_regime = d$diagnostics$rate_regime,
    theoretical_sn = d$diagnostics$theoretical_sn,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    n_leaves_e = result$n_leaves_e,
    n_leaves_m0 = result$n_leaves_m0,
    stringsAsFactors = FALSE
  )
}

# Run test simulations
cat("Starting test simulations...\n")
start_time <- Sys.time()

results_list <- lapply(sim_grid$sim_id, function(id) {
  cat(sprintf("  Sim %d/%d\n", id, nrow(sim_grid)))
  run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
})

results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat(sprintf("\nTest complete in %.1f minutes\n\n", elapsed))

# Analyze results
cat(strrep("=", 70), "\n")
cat("Test Results\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("Convergence rate: %.1f%% (%d/%d)\n",
            100 * mean(results$converged),
            sum(results$converged),
            nrow(results)))

if (all(results$converged)) {
  cat("✓ All simulations converged\n\n")
} else {
  cat("⚠️  Some simulations failed\n\n")
}

# Tree complexity verification
cat("Tree Complexity (s_n) Verification:\n")
cat(strrep("-", 70), "\n")

sn_summary <- results %>%
  filter(converged) %>%
  group_by(dgp, beta, theoretical_sn) %>%
  summarize(
    mean_n_leaves_e = mean(n_leaves_e, na.rm = TRUE),
    mean_n_leaves_m0 = mean(n_leaves_m0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ratio_e = mean_n_leaves_e / theoretical_sn,
    ratio_m0 = mean_n_leaves_m0 / theoretical_sn
  )

print(sn_summary, n = Inf)

cat("\nInterpretation:\n")
cat("  • Ratio ≈ 1.0: Fitted trees match theoretical prediction ✓\n")
cat("  • Ratio < 0.5: Over-regularized (trees too small) ⚠️\n")
cat("  • Ratio > 2.0: Under-regularized (trees too large) ⚠️\n\n")

# Check if any ratios are concerning
if (any(sn_summary$ratio_e < 0.5 | sn_summary$ratio_e > 2.0, na.rm = TRUE)) {
  cat("⚠️  WARNING: Some fitted trees deviate significantly from theory\n")
  cat("   Consider adjusting regularization parameter\n\n")
} else {
  cat("✓ Fitted tree complexity aligns with theoretical predictions\n\n")
}

# Coverage
cat("Coverage (n=800):\n")
cat(strrep("-", 70), "\n")

coverage_summary <- results %>%
  filter(converged) %>%
  group_by(dgp, beta) %>%
  summarize(
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att),
    mean_ci_width = mean(ci_upper - ci_lower),
    .groups = "drop"
  )

print(coverage_summary, n = Inf)

cat("\n")
if (all(coverage_summary$coverage >= 0.80)) {
  cat("✓ Coverage looks reasonable for a small test\n")
} else {
  cat("⚠️  Some coverage < 80% (but this is only 10 reps)\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("Test Summary\n")
cat(strrep("=", 70), "\n\n")

cat("Checks:\n")
cat(sprintf("  [%s] All simulations converged\n",
            ifelse(all(results$converged), "✓", "✗")))
cat(sprintf("  [%s] Tree complexity captured (n_leaves not NA)\n",
            ifelse(all(!is.na(results$n_leaves_e[results$converged])), "✓", "✗")))
cat(sprintf("  [%s] Fitted s_n aligns with theory (ratio ∈ [0.5, 2.0])\n",
            ifelse(all(sn_summary$ratio_e >= 0.5 & sn_summary$ratio_e <= 2.0, na.rm = TRUE), "✓", "✗")))

if (all(results$converged) &&
    all(!is.na(results$n_leaves_e[results$converged])) &&
    all(sn_summary$ratio_e >= 0.5 & sn_summary$ratio_e <= 2.0, na.rm = TRUE)) {
  cat("\n✓ Test PASSED - Ready for full simulation\n")
  cat("\nNext step: Run full simulation with 500 reps\n")
  cat("  Rscript run_beta_study.R\n")
} else {
  cat("\n✗ Test FAILED - Review implementation before full run\n")
}

cat("\n")
