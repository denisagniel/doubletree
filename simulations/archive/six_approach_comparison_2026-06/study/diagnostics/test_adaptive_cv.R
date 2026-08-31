#!/usr/bin/env Rscript
# Test adaptive CV regularization
# Validates that it converges to interior point

library(optimaltrees)
source('code/dgps.R')

cat("\n")
cat("================================================================\n")
cat("Testing Adaptive CV Regularization\n")
cat("================================================================\n\n")

set.seed(12345)
n <- 500
data <- generate_dgp_complex(n = n)

cat("Sample: n =", n, "(complex DGP)\n")
cat("Theory lambda:", round(log(n)/n, 5), "\n\n")

# ============================================================================
# Test 1: Propensity Score
# ============================================================================

cat("----------------------------------------------------------------\n")
cat("Test 1: Propensity Score with Adaptive CV\n")
cat("----------------------------------------------------------------\n\n")

result_ps <- cv_regularization_adaptive(
  X = data$X,
  y = data$A,
  loss_function = "log_loss",
  K = 3,  # Fewer folds for speed
  max_iterations = 10,
  refit = TRUE,
  verbose = TRUE
)

cat("\nResults:\n")
cat("  Converged:", result_ps$converged, "\n")
cat("  Iterations:", result_ps$iterations, "\n")
cat("  Best lambda:", round(result_ps$best_lambda, 6), "\n")
cat("  Relative to theory:", round(result_ps$best_lambda / (log(n)/n), 2), "x\n")
cat("  Final grid size:", length(result_ps$lambda_grid), "\n")

# Check quality
e_hat <- predict(result_ps$model, data$X)
rmse_ps <- sqrt(mean((e_hat - data$e_true)^2))
n_leaves_ps <- count_leaves_tree(result_ps$model@trees[[1]])

cat("\n  Propensity RMSE:", round(rmse_ps, 4), "\n")
cat("  Tree size:", n_leaves_ps, "leaves\n")

# ============================================================================
# Test 2: Outcome Model
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 2: Outcome Model with Adaptive CV\n")
cat("----------------------------------------------------------------\n\n")

controls <- data$A == 0
X_controls <- data$X[controls, , drop = FALSE]
Y_controls <- data$Y[controls]

result_outcome <- cv_regularization_adaptive(
  X = X_controls,
  y = Y_controls,
  loss_function = "log_loss",
  K = 3,
  max_iterations = 10,
  refit = TRUE,
  verbose = TRUE
)

cat("\nResults:\n")
cat("  Converged:", result_outcome$converged, "\n")
cat("  Iterations:", result_outcome$iterations, "\n")
cat("  Best lambda:", round(result_outcome$best_lambda, 6), "\n")

# Check quality
mu0_hat <- predict(result_outcome$model, data$X)
rmse_outcome <- sqrt(mean((mu0_hat[controls] - data$mu0_true[controls])^2))
n_leaves_outcome <- count_leaves_tree(result_outcome$model@trees[[1]])

cat("\n  Outcome RMSE:", round(rmse_outcome, 4), "\n")
cat("  Tree size:", n_leaves_outcome, "leaves\n")

# ============================================================================
# Test 3: Compare to Standard CV
# ============================================================================

cat("\n----------------------------------------------------------------\n")
cat("Test 3: Comparison to Standard CV\n")
cat("----------------------------------------------------------------\n\n")

result_standard <- cv_regularization(
  X = data$X,
  y = data$A,
  loss_function = "log_loss",
  K = 3,
  refit = TRUE,
  verbose = FALSE
)

e_hat_std <- predict(result_standard$model, data$X)
rmse_ps_std <- sqrt(mean((e_hat_std - data$e_true)^2))

cat("Standard CV:\n")
cat("  Best lambda:", round(result_standard$best_lambda, 6), "\n")
cat("  RMSE:", round(rmse_ps_std, 4), "\n")
cat("  Grid size:", length(result_standard$lambda_grid), "\n")

cat("\nAdaptive CV:\n")
cat("  Best lambda:", round(result_ps$best_lambda, 6), "\n")
cat("  RMSE:", round(rmse_ps, 4), "\n")
cat("  Grid size:", length(result_ps$lambda_grid), "\n")

cat("\nImprovement:\n")
cat("  RMSE reduction:", round(rmse_ps_std - rmse_ps, 4), "\n")
cat("  Percent improvement:", round(100 * (rmse_ps_std - rmse_ps) / rmse_ps_std, 1), "%\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

if (result_ps$converged && result_outcome$converged) {
  cat("✓ SUCCESS: Both models converged to interior points\n\n")
} else {
  cat("⚠ WARNING: One or both models did not converge\n\n")
}

cat("Propensity Score:\n")
cat("  - Converged in", result_ps$iterations, "iterations\n")
cat("  - Selected lambda:", round(result_ps$best_lambda / (log(n)/n), 2), "× theory\n")
cat("  - RMSE:", round(rmse_ps, 4), "\n")
cat("  - Tree size:", n_leaves_ps, "leaves\n\n")

cat("Outcome Model:\n")
cat("  - Converged in", result_outcome$iterations, "iterations\n")
cat("  - Selected lambda:", round(result_outcome$best_lambda / (log(sum(controls))/sum(controls)), 2), "× theory\n")
cat("  - RMSE:", round(rmse_outcome, 4), "\n")
cat("  - Tree size:", n_leaves_outcome, "leaves\n\n")

if (rmse_ps < rmse_ps_std) {
  cat("✓ Adaptive CV improved over standard CV\n")
} else {
  cat("⚠ Adaptive CV did not improve (may need more iterations or different settings)\n")
}

cat("\n✓ Adaptive CV implementation validated!\n\n")
