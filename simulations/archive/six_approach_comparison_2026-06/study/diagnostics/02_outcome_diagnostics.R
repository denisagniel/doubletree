#!/usr/bin/env Rscript

# Outcome Model Diagnostics
# Phase 1.2 of diagnostic plan
#
# Analyzes quality of tree-based outcome model estimation:
# - Bias, RMSE on control units
# - Extrapolation error to treated units
# - Tree complexity
# - Oracle performance (fit to true outcome function)
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
use_cv <- FALSE
fixed_lambda_multiplier <- 1.0

# Random seed
set.seed(20260527)

# Output directory
output_dir <- "diagnostics/results/outcome"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Helper Functions
# ============================================================================

#' Fit outcome tree and diagnose
#'
#' @param X Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param mu0_true True outcome function E[Y(0)|X]
#' @param regularization Lambda value
#' @param use_cv Whether to use CV for lambda selection
#' @return List with diagnostics
fit_and_diagnose_outcome <- function(X, A, Y, mu0_true, regularization, use_cv = FALSE) {

  # Fit outcome tree on control units only
  controls <- A == 0
  X_controls <- X[controls, , drop = FALSE]
  Y_controls <- Y[controls]

  n_control <- sum(controls)

  if (n_control < 20) {
    return(NULL)  # Too few controls to fit tree
  }

  # Regularization
  if (use_cv) {
    cv_result <- optimaltrees::cv_regularization(
      X = X_controls,
      y = Y_controls,
      K = 5,
      lambda_grid = regularization * c(0.25, 0.5, 1, 2, 4),
      loss_function = "log_loss"  # Binary outcome
    )
    lambda_selected <- cv_result$best_lambda
    lambda_relative <- lambda_selected / regularization
  } else {
    lambda_selected <- regularization
    lambda_relative <- 1.0
  }

  # Fit tree
  outcome_tree <- optimaltrees::fit_tree(
    X = X_controls,
    y = Y_controls,
    loss = "log_loss",  # Binary outcome
    regularization = lambda_selected
  )

  # Predictions
  mu0_hat_all <- predict(outcome_tree, X)
  mu0_hat_control <- mu0_hat_all[controls]
  mu0_hat_treated <- mu0_hat_all[!controls]

  # True values
  mu0_true_all <- mu0_true
  mu0_true_control <- mu0_true[controls]
  mu0_true_treated <- mu0_true[!controls]

  # ============================================================================
  # Metrics on control units (in-sample)
  # ============================================================================

  bias_control <- mean(mu0_hat_control - mu0_true_control)
  rmse_control <- sqrt(mean((mu0_hat_control - mu0_true_control)^2))
  mae_control <- mean(abs(mu0_hat_control - mu0_true_control))
  max_error_control <- max(abs(mu0_hat_control - mu0_true_control))
  cor_control <- cor(mu0_hat_control, mu0_true_control)

  # Prediction metrics (log loss, calibration)
  pred_metrics_control <- compute_prediction_metrics(
    mu0_hat_control,
    Y_controls,
    type = "classification"
  )

  # ============================================================================
  # Metrics on treated units (extrapolation)
  # ============================================================================

  n_treated <- sum(!controls)

  if (n_treated > 0) {
    bias_treated <- mean(mu0_hat_treated - mu0_true_treated)
    rmse_treated <- sqrt(mean((mu0_hat_treated - mu0_true_treated)^2))
    mae_treated <- mean(abs(mu0_hat_treated - mu0_true_treated))
    max_error_treated <- max(abs(mu0_hat_treated - mu0_true_treated))
    cor_treated <- cor(mu0_hat_treated, mu0_true_treated)

    # Extrapolation error (treated vs control)
    extrapolation_rmse_increase <- rmse_treated - rmse_control
    extrapolation_rmse_ratio <- rmse_treated / rmse_control
  } else {
    bias_treated <- NA_real_
    rmse_treated <- NA_real_
    mae_treated <- NA_real_
    max_error_treated <- NA_real_
    cor_treated <- NA_real_
    extrapolation_rmse_increase <- NA_real_
    extrapolation_rmse_ratio <- NA_real_
  }

  # ============================================================================
  # Tree complexity
  # ============================================================================

  n_leaves <- count_leaves(outcome_tree)
  n_splits <- n_leaves - 1
  tree_depth <- max_depth(outcome_tree)

  # ============================================================================
  # Oracle performance: fit tree to true function
  # ============================================================================

  oracle_result <- compute_oracle_performance(
    X = X_controls,
    y_true = mu0_true_control,
    loss = "squared_error",  # Fit to continuous true values
    regularization = lambda_selected
  )

  oracle_rmse <- oracle_result$metrics$rmse
  oracle_n_leaves <- oracle_result$metrics$n_leaves

  # Oracle vs data: Can tree represent function?
  expressiveness_gap <- rmse_control - oracle_rmse

  return(list(
    # Sample info
    n = nrow(X),
    n_control = n_control,
    n_treated = n_treated,

    # Control units (in-sample)
    bias_control = bias_control,
    rmse_control = rmse_control,
    mae_control = mae_control,
    max_error_control = max_error_control,
    cor_control = cor_control,
    log_loss_control = pred_metrics_control$log_loss,
    calibration_slope_control = pred_metrics_control$calibration_slope,

    # Treated units (extrapolation)
    bias_treated = bias_treated,
    rmse_treated = rmse_treated,
    mae_treated = mae_treated,
    max_error_treated = max_error_treated,
    cor_treated = cor_treated,
    extrapolation_rmse_increase = extrapolation_rmse_increase,
    extrapolation_rmse_ratio = extrapolation_rmse_ratio,

    # Tree complexity
    n_leaves = n_leaves,
    n_splits = n_splits,
    tree_depth = tree_depth,

    # Oracle
    oracle_rmse = oracle_rmse,
    oracle_n_leaves = oracle_n_leaves,
    expressiveness_gap = expressiveness_gap,

    # Regularization
    lambda = lambda_selected,
    lambda_relative = lambda_relative,

    # Predictions (for plotting)
    mu0_hat_all = mu0_hat_all,
    mu0_true_all = mu0_true_all,
    A = A,

    # Tree
    tree = outcome_tree
  ))
}

# ============================================================================
# Main Simulation Loop
# ============================================================================

cat("=================================================\n")
cat("Outcome Model Diagnostics\n")
cat("=================================================\n")
cat(sprintf("Replications: %d\n", n_reps))
cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("DGPs: %s\n", paste(dgp_names, collapse = ", ")))
cat(sprintf("CV for lambda: %s\n", ifelse(use_cv, "YES", "NO")))
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

    # Regularization
    lambda_theory <- log(n) / n
    lambda_used <- fixed_lambda_multiplier * lambda_theory

    # Run replications
    for (rep in 1:n_reps) {
      if (rep %% 25 == 0) {
        cat(sprintf("    Rep %d / %d (%.0f%%)\n", rep, n_reps, 100 * rep / n_reps))
      }

      # Generate data
      data <- tryCatch({
        dgp_fun(n = n)
      }, error = function(e) {
        cat(sprintf("    ERROR in DGP (rep %d): %s\n", rep, conditionMessage(e)))
        return(NULL)
      })

      if (is.null(data)) {
        next
      }

      # Fit and diagnose
      diag <- tryCatch({
        fit_and_diagnose_outcome(
          X = data$X,
          A = data$A,
          Y = data$Y,
          mu0_true = data$mu0_true,
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

      # Store results (exclude heavy objects)
      results[[result_counter]] <- list(
        dgp = dgp_id,
        dgp_name = dgp_name,
        n = n,
        rep = rep,

        # Sample
        n_control = diag$n_control,
        n_treated = diag$n_treated,

        # Control units
        bias_control = diag$bias_control,
        rmse_control = diag$rmse_control,
        mae_control = diag$mae_control,
        max_error_control = diag$max_error_control,
        cor_control = diag$cor_control,
        log_loss_control = diag$log_loss_control,
        calibration_slope_control = diag$calibration_slope_control,

        # Treated units
        bias_treated = diag$bias_treated,
        rmse_treated = diag$rmse_treated,
        mae_treated = diag$mae_treated,
        max_error_treated = diag$max_error_treated,
        cor_treated = diag$cor_treated,
        extrapolation_rmse_increase = diag$extrapolation_rmse_increase,
        extrapolation_rmse_ratio = diag$extrapolation_rmse_ratio,

        # Tree
        n_leaves = diag$n_leaves,
        n_splits = diag$n_splits,
        tree_depth = diag$tree_depth,

        # Oracle
        oracle_rmse = diag$oracle_rmse,
        oracle_n_leaves = diag$oracle_n_leaves,
        expressiveness_gap = diag$expressiveness_gap,

        # Regularization
        lambda = diag$lambda,
        lambda_relative = diag$lambda_relative
      )

      result_counter <- result_counter + 1
    }
  }
}

# ============================================================================
# Convert to Data Frame and Save
# ============================================================================

cat("\n=================================================\n")
cat("Converting results...\n")

results_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

cat(sprintf("Total results: %d\n", nrow(results_df)))

# Save
output_file <- file.path(output_dir, "outcome_diagnostics.rds")
saveRDS(results_df, file = output_file)
cat(sprintf("Saved to: %s\n", output_file))

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=================================================\n")
cat("Summary Statistics\n")
cat("=================================================\n\n")

# Overall
cat("Overall Summary:\n")
cat(sprintf("  Mean RMSE (control): %.4f (SD: %.4f)\n",
            mean(results_df$rmse_control), sd(results_df$rmse_control)))
cat(sprintf("  Mean RMSE (treated): %.4f (SD: %.4f)\n",
            mean(results_df$rmse_treated, na.rm = TRUE),
            sd(results_df$rmse_treated, na.rm = TRUE)))
cat(sprintf("  Mean extrapolation ratio: %.3f (SD: %.3f)\n",
            mean(results_df$extrapolation_rmse_ratio, na.rm = TRUE),
            sd(results_df$extrapolation_rmse_ratio, na.rm = TRUE)))
cat(sprintf("  Mean oracle RMSE: %.4f (SD: %.4f)\n",
            mean(results_df$oracle_rmse), sd(results_df$oracle_rmse)))
cat(sprintf("  Mean expressiveness gap: %.4f (SD: %.4f)\n",
            mean(results_df$expressiveness_gap), sd(results_df$expressiveness_gap)))
cat(sprintf("  Mean tree size: %.1f leaves (SD: %.1f)\n",
            mean(results_df$n_leaves), sd(results_df$n_leaves)))

# By DGP
cat("\n\nBy DGP:\n")
for (dgp_id in dgps) {
  dgp_name <- dgp_names[dgp_id]
  subset <- results_df[results_df$dgp == dgp_id, ]

  cat(sprintf("\n%s (DGP %d):\n", dgp_name, dgp_id))
  cat(sprintf("  RMSE (control): %.4f (%.4f)\n",
              mean(subset$rmse_control), sd(subset$rmse_control)))
  cat(sprintf("  RMSE (treated): %.4f (%.4f)\n",
              mean(subset$rmse_treated, na.rm = TRUE),
              sd(subset$rmse_treated, na.rm = TRUE)))
  cat(sprintf("  Extrapolation ratio: %.3f (%.3f)\n",
              mean(subset$extrapolation_rmse_ratio, na.rm = TRUE),
              sd(subset$extrapolation_rmse_ratio, na.rm = TRUE)))
  cat(sprintf("  Oracle RMSE: %.4f (%.4f)\n",
              mean(subset$oracle_rmse), sd(subset$oracle_rmse)))
  cat(sprintf("  Tree size: %.1f (%.1f) leaves\n",
              mean(subset$n_leaves), sd(subset$n_leaves)))
}

# Complex DGP by sample size
cat("\n\nComplex DGP by Sample Size:\n")
subset_complex <- results_df[results_df$dgp == 3, ]
for (n in sample_sizes) {
  subset_n <- subset_complex[subset_complex$n == n, ]
  cat(sprintf("\nn = %d:\n", n))
  cat(sprintf("  RMSE (control): %.4f (%.4f)\n",
              mean(subset_n$rmse_control), sd(subset_n$rmse_control)))
  cat(sprintf("  Oracle RMSE: %.4f (%.4f)\n",
              mean(subset_n$oracle_rmse), sd(subset_n$oracle_rmse)))
  cat(sprintf("  Tree size: %.1f (%.1f) leaves\n",
              mean(subset_n$n_leaves), sd(subset_n$n_leaves)))
}

# ============================================================================
# Key Finding: Is oracle RMSE large?
# ============================================================================

cat("\n=================================================\n")
cat("KEY DIAGNOSTIC: Oracle Performance\n")
cat("=================================================\n\n")

cat("Oracle trees are fit to the TRUE outcome function (no noise).\n")
cat("If oracle RMSE is large, trees cannot represent the function well.\n\n")

for (dgp_id in dgps) {
  dgp_name <- dgp_names[dgp_id]
  subset <- results_df[results_df$dgp == dgp_id, ]

  oracle_rmse_mean <- mean(subset$oracle_rmse)

  cat(sprintf("%s: Oracle RMSE = %.4f\n", dgp_name, oracle_rmse_mean))

  if (oracle_rmse_mean > 0.10) {
    cat("  ⚠ HIGH: Trees may struggle to represent this function\n")
  } else if (oracle_rmse_mean > 0.05) {
    cat("  ⚠ MODERATE: Some representation difficulty\n")
  } else {
    cat("  ✓ LOW: Trees can represent function well\n")
  }
}

cat("\n")

# ============================================================================
# Generate Plots
# ============================================================================

cat("\n=================================================\n")
cat("Generating plots...\n")
cat("=================================================\n\n")

# Example: complex DGP, n=1000
data_example <- generate_dgp_complex(n = 1000)
set.seed(1)
diag_example <- fit_and_diagnose_outcome(
  X = data_example$X,
  A = data_example$A,
  Y = data_example$Y,
  mu0_true = data_example$mu0_true,
  regularization = log(1000) / 1000,
  use_cv = use_cv
)

# Plot 1: Prediction error
p1 <- plot_prediction_error(
  diag_example$mu0_hat_all,
  diag_example$mu0_true_all,
  title = "Outcome Model: Prediction Error"
)
ggsave(
  file.path(output_dir, "outcome_error_example.png"),
  plot = p1,
  width = 8,
  height = 6
)

# Plot 2: By treatment group
df_example <- data.frame(
  mu0_true = diag_example$mu0_true_all,
  mu0_hat = diag_example$mu0_hat_all,
  A = factor(diag_example$A, levels = c(0, 1), labels = c("Control", "Treated"))
)

p2 <- ggplot(df_example, aes(x = mu0_true, y = mu0_hat, color = A)) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ A) +
  labs(
    title = "Outcome Model: True vs Predicted by Treatment Group",
    x = "True E[Y(0)|X]",
    y = "Predicted E[Y(0)|X]"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "outcome_by_treatment_example.png"),
  plot = p2,
  width = 10,
  height = 5
)

# Plot 3: RMSE comparison
p3 <- ggplot(results_df, aes(x = factor(n), y = rmse_control, fill = dgp_name)) +
  geom_boxplot() +
  facet_wrap(~ dgp_name, scales = "free_y") +
  labs(
    title = "Outcome Model RMSE (Control Units)",
    x = "Sample Size",
    y = "RMSE"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "outcome_rmse_by_dgp.png"),
  plot = p3,
  width = 10,
  height = 8
)

# Plot 4: Oracle vs data RMSE
results_df$rmse_type_data <- results_df$rmse_control
results_df$rmse_type_oracle <- results_df$oracle_rmse

df_long <- rbind(
  data.frame(
    dgp_name = results_df$dgp_name,
    n = results_df$n,
    rmse = results_df$rmse_type_data,
    type = "Data"
  ),
  data.frame(
    dgp_name = results_df$dgp_name,
    n = results_df$n,
    rmse = results_df$rmse_type_oracle,
    type = "Oracle"
  )
)

p4 <- ggplot(df_long, aes(x = factor(n), y = rmse, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ dgp_name, scales = "free_y") +
  labs(
    title = "Data RMSE vs Oracle RMSE",
    subtitle = "Oracle = fit to true function (no noise). Gap = estimation error.",
    x = "Sample Size",
    y = "RMSE",
    fill = ""
  ) +
  theme_minimal()

ggsave(
  file.path(output_dir, "oracle_vs_data_rmse.png"),
  plot = p4,
  width = 12,
  height = 8
)

cat("Plots saved.\n")

# ============================================================================
# Done
# ============================================================================

cat("\n=================================================\n")
cat("Outcome Model Diagnostics Complete\n")
cat("=================================================\n")
cat(sprintf("\nResults: %s\n", output_file))
cat("Next: Run 03_rashomon_diagnostics.R\n\n")
