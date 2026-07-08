#!/usr/bin/env Rscript

# EIF Component Decomposition
# Phase 5 of diagnostic plan
#
# Decomposes ATT estimation bias into components:
# - Component 1: Outcome model error on treated units
# - Component 2: Propensity-weighted residuals on control units
#
# This identifies whether propensity score or outcome model is the primary
# source of bias.
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

# Source DGPs and estimators
source("code/dgps.R")

# ============================================================================
# Configuration
# ============================================================================

# Simulation parameters
n_reps <- 100
sample_sizes <- c(500, 1000, 2000)
dgps <- 1:4
dgp_names <- c("simple", "moderate", "complex", "continuous")

# Fixed regularization
fixed_lambda_multiplier <- 1.0

# Random seed
set.seed(20260527)

# Output directory
output_dir <- "diagnostics/results/eif_decomposition"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Helper Function: Fit Full Pipeline and Decompose
# ============================================================================

#' Fit propensity and outcome trees, compute ATT, decompose EIF
#'
#' @param X Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param dgp DGP number (for computing true functions)
#' @param regularization Lambda value
#' @return List with full decomposition
fit_and_decompose <- function(X, A, Y, dgp, regularization) {

  n <- nrow(X)

  # Compute true nuisance functions
  e_true <- compute_true_propensity(X, dgp)
  mu0_true <- compute_true_outcome(X, dgp)

  # ------------------------------------------
  # Fit propensity score tree
  # ------------------------------------------
  ps_tree <- optimaltrees::fit_tree(
    X = X,
    y = A,
    loss = "log_loss",
    regularization = regularization
  )
  e_hat <- predict(ps_tree, X)

  # ------------------------------------------
  # Fit outcome model tree (on controls)
  # ------------------------------------------
  controls <- A == 0
  X_controls <- X[controls, , drop = FALSE]
  Y_controls <- Y[controls]

  outcome_tree <- optimaltrees::fit_tree(
    X = X_controls,
    y = Y_controls,
    loss = "log_loss",
    regularization = regularization
  )
  mu0_hat <- predict(outcome_tree, X)

  # ------------------------------------------
  # Compute ATT estimate
  # ------------------------------------------
  P_A1 <- mean(A)
  theta_hat <- mean(A / P_A1 * (Y - mu0_hat)) -
               mean((A / P_A1 - 1) * (e_hat / (1 - e_hat)) * (Y - mu0_hat))

  # ------------------------------------------
  # EIF decomposition
  # ------------------------------------------
  decomp <- decompose_eif_components(
    X = X,
    A = A,
    Y = Y,
    e_hat = e_hat,
    mu0_hat = mu0_hat,
    e_true = e_true,
    mu0_true = mu0_true,
    theta_true = 0.15
  )

  # ------------------------------------------
  # Tree complexity
  # ------------------------------------------
  n_leaves_ps <- count_leaves(ps_tree)
  n_leaves_outcome <- count_leaves(outcome_tree)

  return(list(
    # ATT
    theta_hat = theta_hat,
    theta_true = 0.15,
    theta_bias = theta_hat - 0.15,

    # EIF components
    comp1_true = decomp$comp1_true,
    comp1_est = decomp$comp1_est,
    comp1_bias = decomp$comp1_bias,

    comp2_true = decomp$comp2_true,
    comp2_est = decomp$comp2_est,
    comp2_bias = decomp$comp2_bias,

    total_bias_explained = decomp$total_bias_explained,

    # Nuisance function errors
    ps_rmse_control = decomp$ps_rmse_control,
    ps_bias_control = decomp$ps_bias_control,
    outcome_rmse_treated = decomp$outcome_rmse_treated,
    outcome_bias_treated = decomp$outcome_bias_treated,
    outcome_rmse_control = decomp$outcome_rmse_control,
    outcome_bias_control = decomp$outcome_bias_control,

    # Weight diagnostics
    mean_weight_est = decomp$mean_weight_est,
    max_weight_est = decomp$max_weight_est,
    extreme_weights = decomp$extreme_weights,
    weight_bias = decomp$weight_bias,

    # Tree complexity
    n_leaves_ps = n_leaves_ps,
    n_leaves_outcome = n_leaves_outcome,

    # Sample info
    n = decomp$n,
    n_treated = decomp$n_treated,
    n_control = decomp$n_control
  ))
}

# ============================================================================
# Main Simulation Loop
# ============================================================================

cat("=================================================\n")
cat("EIF Component Decomposition\n")
cat("=================================================\n")
cat(sprintf("Replications: %d\n", n_reps))
cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("DGPs: %s\n", paste(dgp_names, collapse = ", ")))
cat("=================================================\n\n")

# Storage
results <- list()
result_counter <- 1

# Loop
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
        cat(sprintf("    ERROR (rep %d): %s\n", rep, conditionMessage(e)))
        return(NULL)
      })

      if (is.null(data)) {
        next
      }

      # Fit and decompose
      decomp <- tryCatch({
        fit_and_decompose(
          X = data$X,
          A = data$A,
          Y = data$Y,
          dgp = dgp_id,
          regularization = lambda_used
        )
      }, error = function(e) {
        cat(sprintf("    ERROR in decomposition (rep %d): %s\n", rep, conditionMessage(e)))
        return(NULL)
      })

      if (is.null(decomp)) {
        next
      }

      # Store
      results[[result_counter]] <- list(
        dgp = dgp_id,
        dgp_name = dgp_name,
        n = n,
        rep = rep,

        # ATT
        theta_hat = decomp$theta_hat,
        theta_bias = decomp$theta_bias,

        # Components
        comp1_bias = decomp$comp1_bias,
        comp2_bias = decomp$comp2_bias,
        total_bias_explained = decomp$total_bias_explained,

        # Nuisance errors
        ps_rmse_control = decomp$ps_rmse_control,
        ps_bias_control = decomp$ps_bias_control,
        outcome_rmse_treated = decomp$outcome_rmse_treated,
        outcome_bias_treated = decomp$outcome_bias_treated,
        outcome_rmse_control = decomp$outcome_rmse_control,
        outcome_bias_control = decomp$outcome_bias_control,

        # Weights
        mean_weight_est = decomp$mean_weight_est,
        max_weight_est = decomp$max_weight_est,
        extreme_weights = decomp$extreme_weights,
        weight_bias = decomp$weight_bias,

        # Trees
        n_leaves_ps = decomp$n_leaves_ps,
        n_leaves_outcome = decomp$n_leaves_outcome,

        # Sample
        n_treated = decomp$n_treated,
        n_control = decomp$n_control
      )

      result_counter <- result_counter + 1
    }
  }
}

# ============================================================================
# Convert and Save
# ============================================================================

cat("\n=================================================\n")
cat("Converting results...\n")

results_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

cat(sprintf("Total: %d\n", nrow(results_df)))

output_file <- file.path(output_dir, "eif_decomposition.rds")
saveRDS(results_df, file = output_file)
cat(sprintf("Saved: %s\n", output_file))

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=================================================\n")
cat("BIAS DECOMPOSITION SUMMARY\n")
cat("=================================================\n\n")

cat("Component 1: Outcome model error on treated units\n")
cat("Component 2: Propensity-weighted residuals on control units\n\n")

# Overall
cat("Overall:\n")
cat(sprintf("  Total bias: %.4f (SD: %.4f)\n",
            mean(results_df$theta_bias), sd(results_df$theta_bias)))
cat(sprintf("  Component 1 bias: %.4f (SD: %.4f)\n",
            mean(results_df$comp1_bias), sd(results_df$comp1_bias)))
cat(sprintf("  Component 2 bias: %.4f (SD: %.4f)\n",
            mean(results_df$comp2_bias), sd(results_df$comp2_bias)))
cat(sprintf("  Bias explained: %.4f (SD: %.4f)\n",
            mean(results_df$total_bias_explained), sd(results_df$total_bias_explained)))

# Percent contribution
mean_comp1 <- mean(abs(results_df$comp1_bias))
mean_comp2 <- mean(abs(results_df$comp2_bias))
total_contrib <- mean_comp1 + mean_comp2

cat(sprintf("\nRelative contribution:\n"))
cat(sprintf("  Component 1: %.1f%%\n", 100 * mean_comp1 / total_contrib))
cat(sprintf("  Component 2: %.1f%%\n", 100 * mean_comp2 / total_contrib))

# By DGP
cat("\n\nBy DGP:\n")
for (dgp_id in dgps) {
  dgp_name <- dgp_names[dgp_id]
  subset <- results_df[results_df$dgp == dgp_id, ]

  cat(sprintf("\n%s:\n", dgp_name))
  cat(sprintf("  Total bias: %.4f (%.4f)\n",
              mean(subset$theta_bias), sd(subset$theta_bias)))
  cat(sprintf("  Component 1 bias: %.4f (%.4f)\n",
              mean(subset$comp1_bias), sd(subset$comp1_bias)))
  cat(sprintf("  Component 2 bias: %.4f (%.4f)\n",
              mean(subset$comp2_bias), sd(subset$comp2_bias)))

  # Which dominates?
  mean_c1 <- mean(abs(subset$comp1_bias))
  mean_c2 <- mean(abs(subset$comp2_bias))

  if (mean_c1 > 2 * mean_c2) {
    cat("  → PRIMARY ISSUE: Outcome model (Component 1)\n")
  } else if (mean_c2 > 2 * mean_c1) {
    cat("  → PRIMARY ISSUE: Propensity score (Component 2)\n")
  } else {
    cat("  → BOTH components contribute\n")
  }
}

# Complex DGP by sample size
cat("\n\nComplex DGP by Sample Size:\n")
subset_complex <- results_df[results_df$dgp == 3, ]

for (n in sample_sizes) {
  subset_n <- subset_complex[subset_complex$n == n, ]

  cat(sprintf("\nn = %d:\n", n))
  cat(sprintf("  Total bias: %.4f (%.4f)\n",
              mean(subset_n$theta_bias), sd(subset_n$theta_bias)))
  cat(sprintf("  Component 1: %.4f (%.4f)\n",
              mean(subset_n$comp1_bias), sd(subset_n$comp1_bias)))
  cat(sprintf("  Component 2: %.4f (%.4f)\n",
              mean(subset_n$comp2_bias), sd(subset_n$comp2_bias)))
}

# ============================================================================
# Nuisance Function Quality
# ============================================================================

cat("\n=================================================\n")
cat("NUISANCE FUNCTION QUALITY\n")
cat("=================================================\n\n")

# Focus on complex DGP
subset_complex <- results_df[results_df$dgp == 3, ]

cat("Complex DGP (where bias is largest):\n\n")

cat("Propensity Score (control units):\n")
cat(sprintf("  RMSE: %.4f (%.4f)\n",
            mean(subset_complex$ps_rmse_control),
            sd(subset_complex$ps_rmse_control)))
cat(sprintf("  Bias: %.4f (%.4f)\n",
            mean(subset_complex$ps_bias_control),
            sd(subset_complex$ps_bias_control)))

cat("\nOutcome Model (treated units - extrapolation):\n")
cat(sprintf("  RMSE: %.4f (%.4f)\n",
            mean(subset_complex$outcome_rmse_treated),
            sd(subset_complex$outcome_rmse_treated)))
cat(sprintf("  Bias: %.4f (%.4f)\n",
            mean(subset_complex$outcome_bias_treated),
            sd(subset_complex$outcome_bias_treated)))

cat("\nOutcome Model (control units - in-sample):\n")
cat(sprintf("  RMSE: %.4f (%.4f)\n",
            mean(subset_complex$outcome_rmse_control),
            sd(subset_complex$outcome_rmse_control)))

cat("\nWeights:\n")
cat(sprintf("  Mean weight: %.2f (%.2f)\n",
            mean(subset_complex$mean_weight_est),
            sd(subset_complex$mean_weight_est)))
cat(sprintf("  Max weight: %.2f (%.2f)\n",
            mean(subset_complex$max_weight_est),
            sd(subset_complex$max_weight_est)))
cat(sprintf("  Extreme weights (>10): %.1f%%\n",
            100 * mean(subset_complex$extreme_weights)))

# ============================================================================
# Key Finding
# ============================================================================

cat("\n=================================================\n")
cat("KEY FINDING\n")
cat("=================================================\n\n")

# Overall assessment
overall_comp1 <- mean(abs(results_df$comp1_bias))
overall_comp2 <- mean(abs(results_df$comp2_bias))

if (overall_comp1 > 1.5 * overall_comp2) {
  cat("PRIMARY SOURCE OF BIAS: Outcome model (Component 1)\n\n")
  cat("The outcome model is poorly estimating E[Y(0)|X] on treated units.\n")
  cat("This suggests:\n")
  cat("  - Trees fitted on controls don't extrapolate well to treated\n")
  cat("  - OR trees are too simple to capture outcome function\n")
  cat("  - OR regularization is too strong\n\n")
  cat("Next steps:\n")
  cat("  1. Check oracle outcome tree performance (02_outcome_diagnostics.R)\n")
  cat("  2. Test weaker regularization for outcome trees\n")
  cat("  3. Consider separate lambda for outcome vs propensity\n")

} else if (overall_comp2 > 1.5 * overall_comp1) {
  cat("PRIMARY SOURCE OF BIAS: Propensity score (Component 2)\n\n")
  cat("The propensity weights are creating bias in the residual term.\n")
  cat("This suggests:\n")
  cat("  - Propensity trees are mis-estimating e(X)\n")
  cat("  - Leading to biased weights on control units\n")
  cat("  - Trees may be too simple to capture propensity function\n\n")
  cat("Next steps:\n")
  cat("  1. Check propensity tree diagnostics (01_propensity_diagnostics.R)\n")
  cat("  2. Test weaker regularization for propensity trees\n")
  cat("  3. Check overlap quality\n")

} else {
  cat("BOTH components contribute roughly equally\n\n")
  cat("Both propensity and outcome models have estimation errors.\n")
  cat("This suggests:\n")
  cat("  - General issue: trees too simple for DGP complexity\n")
  cat("  - Regularization may be too strong for both\n")
  cat("  - OR sample size too small for required tree depth\n\n")
  cat("Next steps:\n")
  cat("  1. Check both propensity and outcome diagnostics\n")
  cat("  2. Test weaker regularization globally\n")
  cat("  3. Consider lambda grid calibration experiments\n")
}

# ============================================================================
# Generate Plots
# ============================================================================

cat("\n=================================================\n")
cat("Generating plots...\n")

# Plot 1: Bias decomposition by DGP
df_plot <- data.frame(
  dgp_name = rep(results_df$dgp_name, 2),
  component = rep(c("Component 1 (Outcome)", "Component 2 (PS)"), each = nrow(results_df)),
  bias = c(results_df$comp1_bias, results_df$comp2_bias)
)

p1 <- ggplot(df_plot, aes(x = component, y = bias, fill = component)) +
  geom_boxplot() +
  facet_wrap(~ dgp_name) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "EIF Bias Decomposition by DGP",
    subtitle = "Component 1 = Outcome model on treated, Component 2 = PS-weighted residuals on control",
    x = "",
    y = "Bias Contribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "bias_decomposition_by_dgp.png"),
  plot = p1,
  width = 12,
  height = 8
)

# Plot 2: Complex DGP - components over sample size
subset_complex <- results_df[results_df$dgp == 3, ]

df_complex <- data.frame(
  n = rep(subset_complex$n, 2),
  component = rep(c("Component 1", "Component 2"), each = nrow(subset_complex)),
  bias = c(subset_complex$comp1_bias, subset_complex$comp2_bias)
)

p2 <- ggplot(df_complex, aes(x = factor(n), y = bias, fill = component)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Complex DGP: Bias Components by Sample Size",
    x = "Sample Size",
    y = "Bias Contribution",
    fill = ""
  ) +
  theme_minimal()

ggsave(
  file.path(output_dir, "complex_bias_by_n.png"),
  plot = p2,
  width = 10,
  height = 6
)

# Plot 3: Scatterplot - comp1 vs comp2 bias
p3 <- ggplot(results_df, aes(x = comp1_bias, y = comp2_bias, color = dgp_name)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    title = "Component 1 vs Component 2 Bias",
    subtitle = "Points off diagonal indicate one component dominates",
    x = "Component 1 Bias (Outcome Model)",
    y = "Component 2 Bias (Propensity Score)",
    color = "DGP"
  ) +
  theme_minimal()

ggsave(
  file.path(output_dir, "comp1_vs_comp2_scatter.png"),
  plot = p3,
  width = 10,
  height = 8
)

cat("Plots saved.\n")

# ============================================================================
# Done
# ============================================================================

cat("\n=================================================\n")
cat("EIF Decomposition Complete\n")
cat("=================================================\n")
cat(sprintf("\nResults: %s\n", output_file))
cat("\nUse findings to guide calibration experiments (06_calibration_experiments.R)\n\n")
