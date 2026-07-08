#!/usr/bin/env Rscript

# Propensity Score Diagnostics
# Phase 1.1 of diagnostic plan
#
# Analyzes quality of tree-based propensity score estimation:
# - Bias, RMSE, max error
# - Calibration quality
# - Overlap and positivity
# - Tree complexity
#
# Created: 2026-05-27

# ============================================================================
# Setup
# ============================================================================

library(optimaltrees)
library(doubletree)
library(ggplot2)
library(gridExtra)

# Source utilities
source("diagnostics/utils/tree_diagnostics.R")
source("diagnostics/utils/eif_components.R")
source("diagnostics/utils/plotting.R")

# Source DGPs
source("code/dgps.R")

# ============================================================================
# Configuration
# ============================================================================

# Simulation parameters
n_reps <- 100
sample_sizes <- c(500, 1000, 2000)
dgps <- 1:4  # All four DGPs
dgp_names <- c("simple", "moderate", "complex", "continuous")

# Regularization
use_cv <- FALSE  # Start with fixed lambda
fixed_lambda_multiplier <- 1.0  # log(n)/n

# Random seed
set.seed(20260527)

# Output directory
output_dir <- "diagnostics/results/propensity"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Helper Functions
# ============================================================================

#' Fit propensity tree and diagnose
#'
#' @param X Covariates
#' @param A Treatment
#' @param e_true True propensity scores
#' @param regularization Lambda value
#' @param use_cv Whether to use CV for lambda selection
#' @return List with diagnostics
fit_and_diagnose_ps <- function(X, A, e_true, regularization, use_cv = FALSE) {

  # Fit propensity tree
  if (use_cv) {
    # Use CV to select lambda
    cv_result <- optimaltrees::cv_regularization(
      X = X,
      y = A,
      K = 5,
      lambda_grid = regularization * c(0.25, 0.5, 1, 2, 4),
      loss_function = "log_loss"
    )
    lambda_selected <- cv_result$best_lambda
    lambda_relative <- lambda_selected / regularization
    cv_losses <- cv_result$cv_loss
  } else {
    # Use fixed lambda
    lambda_selected <- regularization
    lambda_relative <- 1.0
    cv_losses <- NA
  }

  # Fit tree with selected lambda
  ps_tree <- optimaltrees::fit_tree(
    X = X,
    y = A,
    loss = "log_loss",
    regularization = lambda_selected
  )

  # Get predictions
  e_hat <- predict(ps_tree, X)

  # Compute metrics
  n <- nrow(X)

  # 1. Bias and RMSE
  bias_e <- mean(e_hat - e_true)
  rmse_e <- sqrt(mean((e_hat - e_true)^2))
  mae_e <- mean(abs(e_hat - e_true))
  max_abs_error_e <- max(abs(e_hat - e_true))

  # 2. Prediction metrics (log loss, calibration)
  pred_metrics <- compute_prediction_metrics(e_hat, A, type = "classification")

  # 3. Overlap diagnostics
  overlap <- analyze_overlap(e_hat, A)

  # 4. Tree complexity
  n_leaves_ps <- count_leaves(ps_tree)
  n_splits_ps <- n_leaves_ps - 1
  tree_depth_ps <- max_depth(ps_tree)

  # 5. True vs estimated comparison
  cor_e <- cor(e_hat, e_true)

  return(list(
    # Predictions
    e_hat = e_hat,
    e_true = e_true,
    A = A,

    # Error metrics
    bias_e = bias_e,
    rmse_e = rmse_e,
    mae_e = mae_e,
    max_error_e = max_abs_error_e,
    cor_e = cor_e,

    # Prediction quality
    log_loss = pred_metrics$log_loss,
    brier_score = pred_metrics$brier_score,
    calibration_slope = pred_metrics$calibration_slope,
    calibration_intercept = pred_metrics$calibration_intercept,

    # Overlap
    min_e_hat = overlap$min_e,
    max_e_hat = overlap$max_e,
    extreme_low = overlap$extreme_low,
    extreme_high = overlap$extreme_high,
    extreme_total = overlap$extreme_total,
    ess_ratio = overlap$ess_ratio,
    positivity_violations = overlap$positivity_violations,

    # Tree complexity
    n_leaves = n_leaves_ps,
    n_splits = n_splits_ps,
    tree_depth = tree_depth_ps,

    # Regularization
    lambda = lambda_selected,
    lambda_relative = lambda_relative,
    lambda_theory = regularization,

    # CV (if used)
    cv_used = use_cv,
    cv_losses = cv_losses,

    # Tree object
    tree = ps_tree
  ))
}

# ============================================================================
# Main Simulation Loop
# ============================================================================

cat("=================================================\n")
cat("Propensity Score Diagnostics\n")
cat("=================================================\n")
cat(sprintf("Replications: %d\n", n_reps))
cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("DGPs: %s\n", paste(dgp_names, collapse = ", ")))
cat(sprintf("CV for lambda: %s\n", ifelse(use_cv, "YES", "NO")))
cat(sprintf("Fixed lambda multiplier: %.2f\n", fixed_lambda_multiplier))
cat("=================================================\n\n")

# Storage
results <- list()
result_counter <- 1

# Loop over conditions
for (dgp_id in dgps) {
  dgp_name <- dgp_names[dgp_id]

  # Get DGP function
  dgp_fun <- switch(
    dgp_id,
    generate_dgp_simple,
    generate_dgp_moderate,
    generate_dgp_complex,
    generate_dgp_continuous
  )

  cat(sprintf("\n--- DGP %d: %s ---\n", dgp_id, dgp_name))

  for (n in sample_sizes) {
    cat(sprintf("\n  Sample size: n = %d\n", n))

    # Regularization parameter
    lambda_theory <- log(n) / n
    lambda_used <- fixed_lambda_multiplier * lambda_theory

    cat(sprintf("    Theory lambda: %.6f\n", lambda_theory))
    cat(sprintf("    Used lambda: %.6f\n", lambda_used))

    # Run replications
    for (rep in 1:n_reps) {
      if (rep %% 25 == 0) {
        cat(sprintf("    Rep %d / %d (%.0f%%)\n", rep, n_reps, 100 * rep / n_reps))
      }

      # Generate data
      data <- tryCatch({
        dgp_fun(n = n)
      }, error = function(e) {
        cat(sprintf("    ERROR in DGP generation (rep %d): %s\n", rep, conditionMessage(e)))
        return(NULL)
      })

      if (is.null(data)) {
        next
      }

      # Fit and diagnose
      diag <- tryCatch({
        fit_and_diagnose_ps(
          X = data$X,
          A = data$A,
          e_true = data$e_true,
          regularization = lambda_used,
          use_cv = use_cv
        )
      }, error = function(e) {
        cat(sprintf("    ERROR in fitting (rep %d): %s\n", rep, conditionMessage(e)))
        return(NULL)
      })

      if (is.null(diag)) {
        next
      }

      # Store results (exclude heavy objects like tree and predictions)
      results[[result_counter]] <- list(
        dgp = dgp_id,
        dgp_name = dgp_name,
        n = n,
        rep = rep,

        # Error metrics
        bias_e = diag$bias_e,
        rmse_e = diag$rmse_e,
        mae_e = diag$mae_e,
        max_error_e = diag$max_error_e,
        cor_e = diag$cor_e,

        # Prediction quality
        log_loss = diag$log_loss,
        brier_score = diag$brier_score,
        calibration_slope = diag$calibration_slope,
        calibration_intercept = diag$calibration_intercept,

        # Overlap
        min_e_hat = diag$min_e_hat,
        max_e_hat = diag$max_e_hat,
        extreme_low = diag$extreme_low,
        extreme_high = diag$extreme_high,
        extreme_total = diag$extreme_total,
        ess_ratio = diag$ess_ratio,
        positivity_violations = diag$positivity_violations,

        # Tree complexity
        n_leaves = diag$n_leaves,
        n_splits = diag$n_splits,
        tree_depth = diag$tree_depth,

        # Regularization
        lambda = diag$lambda,
        lambda_relative = diag$lambda_relative,
        lambda_theory = diag$lambda_theory
      )

      result_counter <- result_counter + 1
    }

    cat(sprintf("    Completed %d replications\n", n_reps))
  }
}

# ============================================================================
# Convert to Data Frame and Save
# ============================================================================

cat("\n=================================================\n")
cat("Converting results to data frame...\n")

results_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

cat(sprintf("Total results: %d\n", nrow(results_df)))

# Save raw results
output_file <- file.path(output_dir, "propensity_diagnostics.rds")
saveRDS(results_df, file = output_file)
cat(sprintf("Saved to: %s\n", output_file))

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=================================================\n")
cat("Summary Statistics\n")
cat("=================================================\n\n")

# Overall summary
cat("Overall Summary:\n")
cat(sprintf("  Mean RMSE: %.4f (SD: %.4f)\n",
            mean(results_df$rmse_e), sd(results_df$rmse_e)))
cat(sprintf("  Mean MAE: %.4f (SD: %.4f)\n",
            mean(results_df$mae_e), sd(results_df$mae_e)))
cat(sprintf("  Mean bias: %.4f (SD: %.4f)\n",
            mean(results_df$bias_e), sd(results_df$bias_e)))
cat(sprintf("  Mean calibration slope: %.3f (SD: %.3f)\n",
            mean(results_df$calibration_slope, na.rm = TRUE),
            sd(results_df$calibration_slope, na.rm = TRUE)))
cat(sprintf("  Mean correlation: %.3f (SD: %.3f)\n",
            mean(results_df$cor_e), sd(results_df$cor_e)))
cat(sprintf("  Mean tree size: %.1f leaves (SD: %.1f)\n",
            mean(results_df$n_leaves), sd(results_df$n_leaves)))

# By DGP
cat("\n\nBy DGP:\n")
for (dgp_id in dgps) {
  dgp_name <- dgp_names[dgp_id]
  subset <- results_df[results_df$dgp == dgp_id, ]

  cat(sprintf("\n%s (DGP %d):\n", dgp_name, dgp_id))
  cat(sprintf("  RMSE: %.4f (%.4f)\n", mean(subset$rmse_e), sd(subset$rmse_e)))
  cat(sprintf("  Calibration slope: %.3f (%.3f)\n",
              mean(subset$calibration_slope, na.rm = TRUE),
              sd(subset$calibration_slope, na.rm = TRUE)))
  cat(sprintf("  Tree size: %.1f (%.1f) leaves\n",
              mean(subset$n_leaves), sd(subset$n_leaves)))
  cat(sprintf("  Extreme weights: %.1f%%\n",
              100 * mean(subset$extreme_total)))
}

# By sample size (for complex DGP only)
cat("\n\nComplex DGP by Sample Size:\n")
subset_complex <- results_df[results_df$dgp == 3, ]
for (n in sample_sizes) {
  subset_n <- subset_complex[subset_complex$n == n, ]
  cat(sprintf("\nn = %d:\n", n))
  cat(sprintf("  RMSE: %.4f (%.4f)\n", mean(subset_n$rmse_e), sd(subset_n$rmse_e)))
  cat(sprintf("  Calibration slope: %.3f (%.3f)\n",
              mean(subset_n$calibration_slope, na.rm = TRUE),
              sd(subset_n$calibration_slope, na.rm = TRUE)))
  cat(sprintf("  Tree size: %.1f (%.1f) leaves\n",
              mean(subset_n$n_leaves), sd(subset_n$n_leaves)))
}

# ============================================================================
# Generate Plots
# ============================================================================

cat("\n=================================================\n")
cat("Generating diagnostic plots...\n")
cat("=================================================\n\n")

# Focus on complex DGP, n=1000 for detailed plots
cat("Generating example plots (complex DGP, n=1000, rep=1)...\n")

# Re-generate one example for plotting
data_example <- generate_dgp_complex(n = 1000)
set.seed(1)
diag_example <- fit_and_diagnose_ps(
  X = data_example$X,
  A = data_example$A,
  e_true = data_example$e_true,
  regularization = log(1000) / 1000 * fixed_lambda_multiplier,
  use_cv = use_cv
)

# Plot 1: Distribution
p1 <- plot_propensity_distribution(
  diag_example$e_hat,
  diag_example$e_true,
  diag_example$A
)
ggsave(
  file.path(output_dir, "ps_distribution_example.png"),
  plot = p1,
  width = 8,
  height = 6
)

# Plot 2: Calibration
p2 <- plot_calibration(diag_example$e_hat, diag_example$A)
ggsave(
  file.path(output_dir, "ps_calibration_example.png"),
  plot = p2,
  width = 6,
  height = 6
)

# Plot 3: Error vs true
p3 <- plot_prediction_error(
  diag_example$e_hat,
  diag_example$e_true,
  title = "Propensity Score: Prediction Error vs True Value"
)
ggsave(
  file.path(output_dir, "ps_error_example.png"),
  plot = p3,
  width = 8,
  height = 6
)

# Plot 4: RMSE by DGP
p4 <- ggplot(results_df, aes(x = factor(n), y = rmse_e, fill = dgp_name)) +
  geom_boxplot() +
  facet_wrap(~ dgp_name, scales = "free_y") +
  labs(
    title = "Propensity Score RMSE by DGP and Sample Size",
    x = "Sample Size",
    y = "RMSE"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "ps_rmse_by_dgp.png"),
  plot = p4,
  width = 10,
  height = 8
)

# Plot 5: Tree size by DGP
p5 <- ggplot(results_df, aes(x = factor(n), y = n_leaves, fill = dgp_name)) +
  geom_boxplot() +
  facet_wrap(~ dgp_name) +
  labs(
    title = "Tree Complexity by DGP and Sample Size",
    x = "Sample Size",
    y = "Number of Leaves"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "tree_size_by_dgp.png"),
  plot = p5,
  width = 10,
  height = 8
)

cat("Plots saved to:\n")
cat(sprintf("  %s\n", output_dir))

# ============================================================================
# Done
# ============================================================================

cat("\n=================================================\n")
cat("Propensity Diagnostics Complete\n")
cat("=================================================\n")
cat(sprintf("\nResults saved to: %s\n", output_file))
cat("Next steps:\n")
cat("  1. Review summary statistics above\n")
cat("  2. Examine plots in diagnostics/results/propensity/\n")
cat("  3. Run 02_outcome_diagnostics.R for outcome model analysis\n")
cat("\n")
