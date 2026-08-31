#!/usr/bin/env Rscript
# Quick check: What lambda does CV actually select?
# Run single replication with CV diagnostics
library(magrittr)

library(optimaltrees)
source('code/dgps.R')

cat("Checking CV lambda selection...\n\n")

# Test on complex DGP (where bias is worst)
set.seed(12345)
n <- 1000
data <- generate_dgp_complex(n = n)

cat("Sample size: n =", n, "\n")
cat("Theory lambda: log(n)/n =", log(n)/n, "\n")
cat("Default CV grid: (log n / n) * c(0.25, 0.5, 1, 2, 4)\n")
cat("              = [", paste(round((log(n)/n) * c(0.25, 0.5, 1, 2, 4), 5), collapse=", "), "]\n\n")

# ============================================================================
# Test 1: Propensity Score
# ============================================================================

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("Test 1: Propensity Score CV\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

cv_result_ps <- optimaltrees::cv_regularization(
  X = data$X,
  y = data$A,
  loss_function = "log_loss",
  K = 5,
  refit = TRUE,
  verbose = FALSE  # We'll extract info ourselves
)

cat("CV Results:\n")
cat("  Best lambda:", cv_result_ps$best_lambda, "\n")
cat("  Relative to theory:", round(cv_result_ps$best_lambda / (log(n)/n), 2), "x\n")
cat("  Lambda grid:", paste(round(cv_result_ps$lambda_grid, 5), collapse=", "), "\n")
cat("  CV losses:  ", paste(round(cv_result_ps$cv_loss, 4), collapse=", "), "\n")

# Which lambda was selected?
best_idx <- which.min(cv_result_ps$cv_loss)
cat("\n  Selected index:", best_idx, "out of", length(cv_result_ps$lambda_grid), "\n")

if (best_idx == 1) {
  cat("  → CV selected WEAKEST lambda (smallest regularization)\n")
  cat("  → Grid may not go weak enough!\n")
} else if (best_idx == length(cv_result_ps$lambda_grid)) {
  cat("  → CV selected STRONGEST lambda (largest regularization)\n")
  cat("  → Trees being over-regularized!\n")
} else {
  cat("  → CV selected interior point (reasonable)\n")
}

# Check propensity quality
e_hat <- predict(cv_result_ps$model, data$X)
rmse_ps <- sqrt(mean((e_hat - data$e_true)^2))
cat("\n  Propensity RMSE:", round(rmse_ps, 4), "\n")

# Count leaves
n_leaves_ps <- optimaltrees::count_leaves_tree(cv_result_ps$model@trees[[1]])
cat("  Tree size:", n_leaves_ps, "leaves\n")

# ============================================================================
# Test 2: Outcome Model
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("Test 2: Outcome Model CV\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

controls <- data$A == 0
X_controls <- data$X[controls, , drop = FALSE]
Y_controls <- data$Y[controls]

cv_result_outcome <- optimaltrees::cv_regularization(
  X = X_controls,
  y = Y_controls,
  loss_function = "log_loss",
  K = 5,
  refit = TRUE,
  verbose = FALSE
)

n_control <- sum(controls)
cat("Training size:", n_control, "(controls only)\n")
cat("Theory lambda for controls:", round(log(n_control)/n_control, 5), "\n\n")

cat("CV Results:\n")
cat("  Best lambda:", cv_result_outcome$best_lambda, "\n")
cat("  Relative to theory (full n):", round(cv_result_outcome$best_lambda / (log(n)/n), 2), "x\n")
cat("  Relative to theory (n_control):", round(cv_result_outcome$best_lambda / (log(n_control)/n_control), 2), "x\n")
cat("  Lambda grid:", paste(round(cv_result_outcome$lambda_grid, 5), collapse=", "), "\n")
cat("  CV losses:  ", paste(round(cv_result_outcome$cv_loss, 4), collapse=", "), "\n")

best_idx <- which.min(cv_result_outcome$cv_loss)
cat("\n  Selected index:", best_idx, "out of", length(cv_result_outcome$lambda_grid), "\n")

if (best_idx == 1) {
  cat("  → CV selected WEAKEST lambda\n")
  cat("  → Grid may not go weak enough!\n")
} else if (best_idx == length(cv_result_outcome$lambda_grid)) {
  cat("  → CV selected STRONGEST lambda\n")
  cat("  → Trees being over-regularized!\n")
} else {
  cat("  → CV selected interior point\n")
}

# Check outcome quality
mu0_hat <- predict(cv_result_outcome$model, data$X)
rmse_outcome_control <- sqrt(mean((mu0_hat[controls] - data$mu0_true[controls])^2))
rmse_outcome_treated <- sqrt(mean((mu0_hat[!controls] - data$mu0_true[!controls])^2))

cat("\n  Outcome RMSE (control):", round(rmse_outcome_control, 4), "\n")
cat("  Outcome RMSE (treated):", round(rmse_outcome_treated, 4), "\n")

n_leaves_outcome <- optimaltrees::count_leaves_tree(cv_result_outcome$model@trees[[1]])
cat("  Tree size:", n_leaves_outcome, "leaves\n")

# ============================================================================
# Test 3: Oracle comparison
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("Test 3: Oracle Tree (True Function)\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Fit oracle tree to TRUE outcome function
oracle_tree <- optimaltrees::fit_tree(
  X = X_controls,
  y = data$mu0_true[controls],
  loss_function = "squared_error",
  regularization = cv_result_outcome$best_lambda  # Use same lambda
)

oracle_pred <- predict(oracle_tree, X_controls)
oracle_rmse <- sqrt(mean((oracle_pred - data$mu0_true[controls])^2))

cat("Oracle tree (same lambda as data):\n")
cat("  RMSE:", round(oracle_rmse, 4), "\n")
cat("  Gap (data - oracle):", round(rmse_outcome_control - oracle_rmse, 4), "\n")

n_leaves_oracle <- optimaltrees::count_leaves_tree(oracle_tree@trees[[1]])
cat("  Tree size:", n_leaves_oracle, "leaves\n")

if (oracle_rmse < 0.01) {
  cat("\n  → Trees CAN represent function perfectly\n")
  cat("  → Problem is estimation noise / regularization\n")
} else if (oracle_rmse < 0.10) {
  cat("\n  → Trees can represent function reasonably well\n")
  cat("  → Some expressiveness limitation\n")
} else {
  cat("\n  → Trees CANNOT represent function well\n")
  cat("  → Fundamental expressiveness problem\n")
}

# ============================================================================
# Test 4: Try weaker lambda
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("Test 4: What if we use weaker lambda?\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Test with extended grid going to weaker lambda
extended_grid <- (log(n)/n) * c(0.1, 0.25, 0.5, 1, 2, 4)

cv_extended <- optimaltrees::cv_regularization(
  X = data$X,
  y = data$A,
  loss_function = "log_loss",
  lambda_grid = extended_grid,
  K = 5,
  refit = TRUE,
  verbose = FALSE
)

cat("Extended grid (added λ×0.1):\n")
cat("  Lambda grid:", paste(round(cv_extended$lambda_grid, 5), collapse=", "), "\n")
cat("  CV losses:  ", paste(round(cv_extended$cv_loss, 4), collapse=", "), "\n")
cat("\n  Best lambda:", cv_extended$best_lambda, "\n")
cat("  Relative to theory:", round(cv_extended$best_lambda / (log(n)/n), 2), "x\n")

best_idx_ext <- which.min(cv_extended$cv_loss)
cat("  Selected index:", best_idx_ext, "out of", length(cv_extended$lambda_grid), "\n")

if (best_idx_ext == 1) {
  cat("\n  → STILL selecting weakest lambda!\n")
  cat("  → Grid needs to go even weaker (try λ×0.05, λ×0.01)\n")
} else {
  cat("\n  → Interior point selected\n")

  # Compare quality
  e_hat_ext <- predict(cv_extended$model, data$X)
  rmse_ps_ext <- sqrt(mean((e_hat_ext - data$e_true)^2))

  cat("\n  Comparison:\n")
  cat("    Default grid RMSE:", round(rmse_ps, 4), "\n")
  cat("    Extended grid RMSE:", round(rmse_ps_ext, 4), "\n")
  cat("    Improvement:", round(rmse_ps - rmse_ps_ext, 4), "\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

cat("Key findings:\n\n")

cat("1. CV-selected lambda:\n")
cat("   - Propensity:", round(cv_result_ps$best_lambda / (log(n)/n), 2), "× theory value\n")
cat("   - Outcome:   ", round(cv_result_outcome$best_lambda / (log(n_control)/n_control), 2), "× theory value\n\n")

cat("2. Quality metrics:\n")
cat("   - Propensity RMSE:", round(rmse_ps, 4), "\n")
cat("   - Outcome RMSE (data):", round(rmse_outcome_control, 4), "\n")
cat("   - Outcome RMSE (oracle):", round(oracle_rmse, 4), "\n")
cat("   - Gap:", round(rmse_outcome_control - oracle_rmse, 4), "\n\n")

cat("3. Tree complexity:\n")
cat("   - Propensity:", n_leaves_ps, "leaves\n")
cat("   - Outcome:   ", n_leaves_outcome, "leaves\n")
cat("   - Oracle:    ", n_leaves_oracle, "leaves\n\n")

# Diagnosis
if (best_idx == 1 || best_idx_ext == 1) {
  cat("DIAGNOSIS: CV grid too strong!\n")
  cat("  - CV consistently selects weakest available lambda\n")
  cat("  - Held-out loss continues improving with weaker λ\n")
  cat("  - Need grid: (log n / n) * c(0.01, 0.05, 0.1, 0.25, 0.5)\n\n")

  cat("ACTION: Extend lambda grid to weaker values\n")
} else if (best_idx == length(cv_result_ps$lambda_grid)) {
  cat("DIAGNOSIS: CV selecting too strong lambda!\n")
  cat("  - CV prefers strongest regularization\n")
  cat("  - May be overfitting to held-out loss\n")
  cat("  - Or trees genuinely can't improve fit\n\n")

  cat("ACTION: Check if trees can represent functions (oracle test)\n")
} else {
  cat("DIAGNOSIS: CV selecting reasonable lambda\n")
  cat("  - Interior point selected\n")
  cat("  - Grid appears appropriate\n")
  cat("  - Problem may not be CV selection\n\n")

  if (oracle_rmse < 0.01 && rmse_outcome_control > 0.20) {
    cat("ACTION: Large gap between oracle and data suggests:\n")
    cat("  - Estimation noise from small samples\n")
    cat("  - Or need different loss function\n")
    cat("  - Or inference loss ≠ prediction loss\n")
  }
}

cat("\n")
