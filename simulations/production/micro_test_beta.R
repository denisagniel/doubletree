#' Micro Test: β Smoothness Regime Study
#'
#' Quick verification (3 reps, ~2 minutes) to confirm:
#' 1. Predict bug fix working (predictions vary, not constant 0.5)
#' 2. No memory issues
#' 3. All three β regimes run successfully
#'
#' Run this BEFORE committing to full 500-rep simulation.

library(dplyr)

# Load packages
devtools::load_all("../../../optimaltrees")

# Source dmltree functions
source("../../R/estimate_att.R")
source("../../R/nuisance_trees.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/utils.R")

source("dgps/dgps_beta_continuous.R")

# Helper: Count leaves in a tree
count_tree_leaves <- function(tree) {
  if (is.null(tree) || !is.list(tree)) return(0)

  if (!is.null(tree$name) && tree$name == "class" && is.null(tree$feature)) {
    return(1)
  }

  left <- if (!is.null(tree$false)) count_tree_leaves(tree$false) else 0
  right <- if (!is.null(tree$true)) count_tree_leaves(tree$true) else 0

  return(left + right)
}

# Micro configuration
N_REPS <- 3
N <- 400  # Small sample
TAU <- 0.10
K_FOLDS <- 2  # Reduced for speed
SEED_OFFSET <- 90000

# DGPs
DGPS <- list(
  beta_high = generate_dgp_beta_high,       # β=3 > d/2
  beta_boundary = generate_dgp_beta_boundary,  # β=2 = d/2
  beta_low = generate_dgp_beta_low          # β=1 < d/2
)

cat(strrep("=", 70), "\n")
cat("MICRO TEST: β Smoothness Regime Study\n")
cat(strrep("=", 70), "\n\n")

cat("Configuration:\n")
cat(sprintf("  Reps per DGP: %d\n", N_REPS))
cat(sprintf("  Sample size: %d\n", N))
cat(sprintf("  Folds: %d\n", K_FOLDS))
cat(sprintf("  Total runs: %d (3 DGPs × %d reps)\n\n", 3 * N_REPS, N_REPS))

cat("Expected runtime: ~2 minutes\n\n")

# Track memory before
mem_before <- gc()
cat("Memory before:\n")
print(mem_before)
cat("\n")

# Run micro test
results_list <- list()
sim_id <- 0

for (dgp_name in names(DGPS)) {
  cat(sprintf("Testing %s (β=%s)...\n",
              dgp_name,
              switch(dgp_name,
                     beta_high = "3",
                     beta_boundary = "2",
                     beta_low = "1")))

  dgp_func <- DGPS[[dgp_name]]

  for (rep in 1:N_REPS) {
    sim_id <- sim_id + 1
    seed <- SEED_OFFSET + sim_id

    # Generate data
    d <- dgp_func(n = N, tau = TAU, seed = seed)

    # Fit tree-DML
    result <- tryCatch({
      capture.output({
        fit <- estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = K_FOLDS,
          regularization = log(N) / N,
          cv_regularization = FALSE,
          use_rashomon = FALSE,
          verbose = FALSE,
          discretize_method = "quantiles",
          discretize_bins = "adaptive"
        )
      }, file = tempfile())

      # Extract tree complexity
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

      # Check predictions variance (predict bug check)
      pred_sd_e <- NA
      pred_sd_m0 <- NA

      if (!is.null(fit$nuisance_fits) && length(fit$nuisance_fits) > 0) {
        # Get predictions from first fold as diagnostic
        fold1 <- fit$nuisance_fits[[1]]
        if (!is.null(fold1$e_pred)) {
          pred_sd_e <- sd(fold1$e_pred)
        }
        if (!is.null(fold1$m0_pred)) {
          pred_sd_m0 <- sd(fold1$m0_pred)
        }
      }

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        n_leaves_e = n_leaves_e,
        n_leaves_m0 = n_leaves_m0,
        pred_sd_e = pred_sd_e,      # NEW: for predict bug check
        pred_sd_m0 = pred_sd_m0
      )

    }, error = function(e) {
      list(
        theta = NA,
        sigma = NA,
        ci_lower = NA,
        ci_upper = NA,
        converged = FALSE,
        n_leaves_e = NA,
        n_leaves_m0 = NA,
        pred_sd_e = NA,
        pred_sd_m0 = NA,
        error_msg = as.character(e)
      )
    })

    # Store result
    results_list[[sim_id]] <- data.frame(
      sim_id = sim_id,
      dgp = dgp_name,
      n = N,
      rep = rep,
      true_att = d$true_att,
      beta = d$diagnostics$beta,
      theta = result$theta,
      sigma = result$sigma,
      ci_lower = result$ci_lower,
      ci_upper = result$ci_upper,
      converged = result$converged,
      n_leaves_e = result$n_leaves_e,
      n_leaves_m0 = result$n_leaves_m0,
      pred_sd_e = result$pred_sd_e,
      pred_sd_m0 = result$pred_sd_m0,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  Rep %d: theta=%.3f, sigma=%.3f, converged=%s\n",
                rep, result$theta, result$sigma, result$converged))
  }

  cat("\n")
}

# Combine results
results <- do.call(rbind, results_list)

# Track memory after
mem_after <- gc()
cat("Memory after:\n")
print(mem_after)
cat("\n")

cat(strrep("=", 70), "\n")
cat("MICRO TEST RESULTS\n")
cat(strrep("=", 70), "\n\n")

# Check 1: Convergence
cat("Check 1: Convergence\n")
convergence_rate <- mean(results$converged)
cat(sprintf("  Convergence: %.0f%% (%d/%d)\n",
            100 * convergence_rate,
            sum(results$converged),
            nrow(results)))

if (convergence_rate == 1.0) {
  cat("  ✓ All simulations converged\n\n")
} else {
  cat("  ✗ Some simulations failed\n\n")
}

# Check 2: Predict bug (CRITICAL)
cat("Check 2: Predict Bug Fix (predictions should vary, not constant 0.5)\n")

pred_check <- results %>%
  filter(converged) %>%
  summarize(
    mean_pred_sd_e = mean(pred_sd_e, na.rm = TRUE),
    mean_pred_sd_m0 = mean(pred_sd_m0, na.rm = TRUE),
    min_pred_sd_e = min(pred_sd_e, na.rm = TRUE),
    min_pred_sd_m0 = min(pred_sd_m0, na.rm = TRUE)
  )

cat(sprintf("  Propensity predictions SD: mean=%.4f, min=%.4f\n",
            pred_check$mean_pred_sd_e, pred_check$min_pred_sd_e))
cat(sprintf("  Outcome predictions SD:    mean=%.4f, min=%.4f\n",
            pred_check$mean_pred_sd_m0, pred_check$min_pred_sd_m0))

# If SD < 1e-6, predictions are constant (bug not fixed)
if (pred_check$min_pred_sd_e < 1e-6 || pred_check$min_pred_sd_m0 < 1e-6) {
  cat("  ✗ PREDICT BUG DETECTED: Predictions are constant!\n")
  cat("  ✗ DO NOT RUN FULL SIMULATION - FIX PREDICT BUG FIRST\n\n")
} else {
  cat("  ✓ Predictions vary correctly (bug fixed)\n\n")
}

# Check 3: Tree complexity
cat("Check 3: Tree Complexity\n")

tree_stats <- results %>%
  filter(converged) %>%
  group_by(dgp, beta) %>%
  summarize(
    mean_leaves_e = mean(n_leaves_e, na.rm = TRUE),
    mean_leaves_m0 = mean(n_leaves_m0, na.rm = TRUE),
    .groups = "drop"
  )

print(tree_stats)

if (any(tree_stats$mean_leaves_e < 2 | tree_stats$mean_leaves_m0 < 2)) {
  cat("\n  ⚠️  Some trees have <2 leaves (constant predictions)\n\n")
} else {
  cat("\n  ✓ All trees have reasonable complexity\n\n")
}

# Check 4: Coverage
cat("Check 4: Coverage (very rough with n=3)\n")

coverage_stats <- results %>%
  filter(converged) %>%
  group_by(dgp, beta) %>%
  summarize(
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

print(coverage_stats)
cat("\n  Note: With only 3 reps, coverage estimates are unreliable\n\n")

# Overall verdict
cat(strrep("=", 70), "\n")
cat("VERDICT\n")
cat(strrep("=", 70), "\n\n")

all_converged <- convergence_rate == 1.0
predict_ok <- pred_check$min_pred_sd_e >= 1e-6 && pred_check$min_pred_sd_m0 >= 1e-6
trees_ok <- all(tree_stats$mean_leaves_e >= 2 & tree_stats$mean_leaves_m0 >= 2)

if (all_converged && predict_ok && trees_ok) {
  cat("✓ MICRO TEST PASSED\n\n")
  cat("Ready to run full simulation:\n")
  cat("  1. Close Chrome to free memory (saves ~4GB)\n")
  cat("  2. Pause OneDrive to prevent sync conflicts\n")
  cat("  3. Run: Rscript run_beta_study.R\n")
  cat("  4. Expected runtime: ~10-12 hours\n\n")
} else {
  cat("✗ MICRO TEST FAILED\n\n")
  cat("Issues detected:\n")
  if (!all_converged) cat("  - Not all simulations converged\n")
  if (!predict_ok) cat("  - Predict bug still present (constant predictions)\n")
  if (!trees_ok) cat("  - Tree complexity too low\n")
  cat("\nFix issues before running full simulation.\n\n")
}

# Save micro test results for reference
saveRDS(results, "results/micro_test_beta_results.rds")
cat("Micro test results saved to: results/micro_test_beta_results.rds\n")
