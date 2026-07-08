#!/usr/bin/env Rscript

# Quick Test of Diagnostic Infrastructure
# Reduced settings for fast validation
# Created: 2026-05-27

cat("\n================================================================\n")
cat("  DIAGNOSTIC TEST RUN\n")
cat("  Quick validation with reduced settings\n")
cat("================================================================\n\n")

# ============================================================================
# Setup
# ============================================================================

library(optimaltrees)
library(doubletree)
library(ggplot2)

# Source utilities
source("diagnostics/utils/tree_diagnostics.R")
source("diagnostics/utils/eif_components.R")
source("diagnostics/utils/plotting.R")

# Source DGPs
source("code/dgps.R")

# ============================================================================
# REDUCED Configuration for Testing
# ============================================================================

cat("Test configuration:\n")
cat("  - 10 replications (instead of 100)\n")
cat("  - n = 1000 only (instead of 500, 1000, 2000)\n")
cat("  - Complex DGP only (dgp = 3)\n")
cat("  - Fixed lambda (no CV)\n\n")

n_reps <- 10
sample_sizes <- c(1000)
dgps <- 3  # Complex DGP only
dgp_names <- c("complex")

fixed_lambda_multiplier <- 1.0

set.seed(20260527)

# Output directory
output_dir <- "diagnostics/results/test"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Test 1: Propensity Score Diagnostics
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 1: Propensity Score Diagnostics\n")
cat("----------------------------------------------------------------\n\n")

test1_start <- Sys.time()

ps_results <- list()

for (rep in 1:n_reps) {
  cat(sprintf("  Rep %d / %d\r", rep, n_reps))

  # Generate data
  data <- generate_dgp_complex(n = 1000)

  # Fit propensity tree
  lambda <- log(1000) / 1000
  ps_tree <- optimaltrees::fit_tree(
    X = data$X,
    y = data$A,
    loss_function = "log_loss",
    regularization = lambda
  )

  # Get predictions
  e_hat <- predict(ps_tree, data$X)

  # Compute metrics
  bias_e <- mean(e_hat - data$e_true)
  rmse_e <- sqrt(mean((e_hat - data$e_true)^2))
  n_leaves <- count_leaves(ps_tree)

  ps_results[[rep]] <- list(
    rep = rep,
    bias_e = bias_e,
    rmse_e = rmse_e,
    n_leaves = n_leaves
  )
}

cat("\n\n")

ps_df <- do.call(rbind, lapply(ps_results, as.data.frame))

cat(sprintf("Results:\n"))
cat(sprintf("  Mean RMSE: %.4f (SD: %.4f)\n", mean(ps_df$rmse_e), sd(ps_df$rmse_e)))
cat(sprintf("  Mean bias: %.4f (SD: %.4f)\n", mean(ps_df$bias_e), sd(ps_df$bias_e)))
cat(sprintf("  Mean tree size: %.1f leaves (SD: %.1f)\n", mean(ps_df$n_leaves), sd(ps_df$n_leaves)))

test1_time <- difftime(Sys.time(), test1_start, units = "secs")
cat(sprintf("  Time: %.1f seconds\n", as.numeric(test1_time)))

if (mean(ps_df$rmse_e) > 0.15) {
  cat("  ⚠ HIGH RMSE: Propensity trees struggling\n")
} else {
  cat("  ✓ RMSE reasonable\n")
}

# ============================================================================
# Test 2: Outcome Model Diagnostics
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 2: Outcome Model Diagnostics\n")
cat("----------------------------------------------------------------\n\n")

test2_start <- Sys.time()

outcome_results <- list()

for (rep in 1:n_reps) {
  cat(sprintf("  Rep %d / %d\r", rep, n_reps))

  # Generate data
  data <- generate_dgp_complex(n = 1000)

  # Fit outcome tree (on controls)
  controls <- data$A == 0
  X_controls <- data$X[controls, , drop = FALSE]
  Y_controls <- data$Y[controls]

  lambda <- log(sum(controls)) / sum(controls)
  outcome_tree <- optimaltrees::fit_tree(
    X = X_controls,
    y = Y_controls,
    loss_function = "log_loss",
    regularization = lambda
  )

  # Predictions
  mu0_hat <- predict(outcome_tree, data$X)

  # Metrics on controls
  rmse_control <- sqrt(mean((mu0_hat[controls] - data$mu0_true[controls])^2))

  # Metrics on treated (extrapolation)
  treated <- data$A == 1
  rmse_treated <- sqrt(mean((mu0_hat[treated] - data$mu0_true[treated])^2))

  # Oracle: fit to true function
  oracle_tree <- optimaltrees::fit_tree(
    X = X_controls,
    y = data$mu0_true[controls],
    loss = "squared_error",
    regularization = lambda
  )
  oracle_pred <- predict(oracle_tree, data$X[controls, , drop = FALSE])
  oracle_rmse <- sqrt(mean((oracle_pred - data$mu0_true[controls])^2))

  n_leaves <- count_leaves(outcome_tree)

  outcome_results[[rep]] <- list(
    rep = rep,
    rmse_control = rmse_control,
    rmse_treated = rmse_treated,
    oracle_rmse = oracle_rmse,
    n_leaves = n_leaves
  )
}

cat("\n\n")

outcome_df <- do.call(rbind, lapply(outcome_results, as.data.frame))

cat(sprintf("Results:\n"))
cat(sprintf("  Mean RMSE (control): %.4f (SD: %.4f)\n",
            mean(outcome_df$rmse_control), sd(outcome_df$rmse_control)))
cat(sprintf("  Mean RMSE (treated): %.4f (SD: %.4f)\n",
            mean(outcome_df$rmse_treated), sd(outcome_df$rmse_treated)))
cat(sprintf("  Mean oracle RMSE: %.4f (SD: %.4f)\n",
            mean(outcome_df$oracle_rmse), sd(outcome_df$oracle_rmse)))
cat(sprintf("  Mean tree size: %.1f leaves (SD: %.1f)\n",
            mean(outcome_df$n_leaves), sd(outcome_df$n_leaves)))

test2_time <- difftime(Sys.time(), test2_start, units = "secs")
cat(sprintf("  Time: %.1f seconds\n", as.numeric(test2_time)))

oracle_rmse_mean <- mean(outcome_df$oracle_rmse)
if (oracle_rmse_mean > 0.10) {
  cat("  ⚠ HIGH ORACLE RMSE: Trees cannot represent function well\n")
} else if (oracle_rmse_mean > 0.05) {
  cat("  ⚠ MODERATE ORACLE RMSE: Some expressiveness issues\n")
} else {
  cat("  ✓ LOW ORACLE RMSE: Trees can represent function\n")
}

# ============================================================================
# Test 3: EIF Decomposition
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 3: EIF Bias Decomposition\n")
cat("----------------------------------------------------------------\n\n")

test3_start <- Sys.time()

eif_results <- list()

for (rep in 1:n_reps) {
  cat(sprintf("  Rep %d / %d\r", rep, n_reps))

  # Generate data
  data <- generate_dgp_complex(n = 1000)

  # Fit propensity tree
  lambda <- log(1000) / 1000
  ps_tree <- optimaltrees::fit_tree(
    X = data$X,
    y = data$A,
    loss_function = "log_loss",
    regularization = lambda
  )
  e_hat <- predict(ps_tree, data$X)

  # Fit outcome tree
  controls <- data$A == 0
  X_controls <- data$X[controls, , drop = FALSE]
  Y_controls <- data$Y[controls]

  outcome_tree <- optimaltrees::fit_tree(
    X = X_controls,
    y = Y_controls,
    loss_function = "log_loss",
    regularization = lambda
  )
  mu0_hat <- predict(outcome_tree, data$X)

  # Compute ATT
  P_A1 <- mean(data$A)
  theta_hat <- mean(data$A / P_A1 * (data$Y - mu0_hat)) -
               mean((data$A / P_A1 - 1) * (e_hat / (1 - e_hat)) * (data$Y - mu0_hat))

  # Decompose
  decomp <- decompose_eif_components(
    X = data$X,
    A = data$A,
    Y = data$Y,
    e_hat = e_hat,
    mu0_hat = mu0_hat,
    e_true = data$e_true,
    mu0_true = data$mu0_true,
    theta_true = 0.15
  )

  eif_results[[rep]] <- list(
    rep = rep,
    theta_hat = theta_hat,
    theta_bias = theta_hat - 0.15,
    comp1_bias = decomp$comp1_bias,
    comp2_bias = decomp$comp2_bias
  )
}

cat("\n\n")

eif_df <- do.call(rbind, lapply(eif_results, as.data.frame))

cat(sprintf("Results:\n"))
cat(sprintf("  Mean ATT bias: %.4f (SD: %.4f)\n",
            mean(eif_df$theta_bias), sd(eif_df$theta_bias)))
cat(sprintf("  Component 1 (Outcome): %.4f (SD: %.4f)\n",
            mean(eif_df$comp1_bias), sd(eif_df$comp1_bias)))
cat(sprintf("  Component 2 (Propensity): %.4f (SD: %.4f)\n",
            mean(eif_df$comp2_bias), sd(eif_df$comp2_bias)))

test3_time <- difftime(Sys.time(), test3_start, units = "secs")
cat(sprintf("  Time: %.1f seconds\n", as.numeric(test3_time)))

comp1_abs <- mean(abs(eif_df$comp1_bias))
comp2_abs <- mean(abs(eif_df$comp2_bias))

if (comp1_abs > 2 * comp2_abs) {
  cat("  → PRIMARY ISSUE: Outcome model (Component 1)\n")
} else if (comp2_abs > 2 * comp1_abs) {
  cat("  → PRIMARY ISSUE: Propensity score (Component 2)\n")
} else {
  cat("  → BOTH components contribute\n")
}

# ============================================================================
# Test 4: Plotting
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 4: Plotting Functions\n")
cat("----------------------------------------------------------------\n\n")

test4_start <- Sys.time()

# Generate one example for plotting
data_example <- generate_dgp_complex(n = 1000)

# Fit trees
ps_tree <- optimaltrees::fit_tree(
  X = data_example$X,
  y = data_example$A,
  loss_function = "log_loss",
  regularization = log(1000) / 1000
)
e_hat <- predict(ps_tree, data_example$X)

cat("  Creating plots...\n")

# Test plot 1: Propensity distribution
p1 <- plot_propensity_distribution(e_hat, data_example$e_true, data_example$A)
ggsave(
  file.path(output_dir, "test_ps_distribution.png"),
  plot = p1,
  width = 8,
  height = 6
)
cat("    ✓ Propensity distribution plot\n")

# Test plot 2: Calibration
p2 <- plot_calibration(e_hat, data_example$A)
ggsave(
  file.path(output_dir, "test_calibration.png"),
  plot = p2,
  width = 6,
  height = 6
)
cat("    ✓ Calibration plot\n")

# Test plot 3: Prediction error
p3 <- plot_prediction_error(e_hat, data_example$e_true, "Test: Propensity Error")
ggsave(
  file.path(output_dir, "test_ps_error.png"),
  plot = p3,
  width = 8,
  height = 6
)
cat("    ✓ Prediction error plot\n")

test4_time <- difftime(Sys.time(), test4_start, units = "secs")
cat(sprintf("  Time: %.1f seconds\n", as.numeric(test4_time)))

# ============================================================================
# Summary
# ============================================================================

cat("\n================================================================\n")
cat("  TEST SUMMARY\n")
cat("================================================================\n\n")

total_time <- test1_time + test2_time + test3_time + test4_time

cat("All tests completed successfully!\n\n")

cat("Timing:\n")
cat(sprintf("  Test 1 (Propensity): %.1f sec\n", as.numeric(test1_time)))
cat(sprintf("  Test 2 (Outcome): %.1f sec\n", as.numeric(test2_time)))
cat(sprintf("  Test 3 (EIF): %.1f sec\n", as.numeric(test3_time)))
cat(sprintf("  Test 4 (Plotting): %.1f sec\n", as.numeric(test4_time)))
cat(sprintf("  Total: %.1f sec (%.1f min)\n",
            as.numeric(total_time), as.numeric(total_time) / 60))

cat("\nExtrapolation to full run:\n")
full_reps <- 100
full_dgps <- 4
full_n_sizes <- 3
scaling_factor <- (full_reps / n_reps) * full_dgps * full_n_sizes

cat(sprintf("  Full run would take: ~%.0f minutes (%.1f hours)\n",
            as.numeric(total_time) / 60 * scaling_factor,
            as.numeric(total_time) / 3600 * scaling_factor))

cat("\nKey Findings (from test):\n")
cat(sprintf("  - Propensity RMSE: %.4f\n", mean(ps_df$rmse_e)))
cat(sprintf("  - Outcome oracle RMSE: %.4f\n", mean(outcome_df$oracle_rmse)))
cat(sprintf("  - ATT bias: %.4f\n", mean(eif_df$theta_bias)))

cat("\nOutput:\n")
cat(sprintf("  - Test plots saved to: %s\n", output_dir))

cat("\n✓ Infrastructure validated and ready for full run\n\n")
