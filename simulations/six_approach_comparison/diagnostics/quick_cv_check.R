#!/usr/bin/env Rscript
# Ultra-simple CV check - just report what lambda is selected
library(optimaltrees)
source('code/dgps.R')

set.seed(12345)
data <- generate_dgp_complex(n = 500)  # Smaller n for speed

cat("Quick CV check (n=500, complex DGP)\n\n")

# Propensity
cat("Propensity Score:\n")
cv_ps <- optimaltrees::cv_regularization(
  X = data$X, y = data$A,
  loss_function = "log_loss",
  K = 3,  # Fewer folds
  refit = FALSE,  # Don't refit
  verbose = FALSE
)

cat("  Best lambda:", cv_ps$best_lambda, "\n")
cat("  Theory lambda:", log(500)/500, "\n")
cat("  Ratio:", round(cv_ps$best_lambda / (log(500)/500), 2), "x\n")
cat("  Grid:", paste(round(cv_ps$lambda_grid, 5), collapse=", "), "\n")
cat("  Losses:", paste(round(cv_ps$cv_loss, 4), collapse=", "), "\n")

best_idx <- which.min(cv_ps$cv_loss)
cat("  Selected index:", best_idx, "/", length(cv_ps$lambda_grid), "\n")

if (best_idx == 1) {
  cat("  → BOUNDARY: Selected weakest lambda!\n")
  cat("  → Grid doesn't go weak enough\n\n")
} else if (best_idx == length(cv_ps$lambda_grid)) {
  cat("  → BOUNDARY: Selected strongest lambda!\n")
  cat("  → Over-regularizing\n\n")
} else {
  cat("  → Interior point (OK)\n\n")
}

# Outcome
controls <- data$A == 0
cat("\nOutcome Model:\n")
cv_outcome <- optimaltrees::cv_regularization(
  X = data$X[controls,], y = data$Y[controls],
  loss_function = "log_loss",
  K = 3,
  refit = FALSE,
  verbose = FALSE
)

cat("  Best lambda:", cv_outcome$best_lambda, "\n")
cat("  Ratio:", round(cv_outcome$best_lambda / (log(500)/500), 2), "x\n")
cat("  Grid:", paste(round(cv_outcome$lambda_grid, 5), collapse=", "), "\n")
cat("  Losses:", paste(round(cv_outcome$cv_loss, 4), collapse=", "), "\n")

best_idx2 <- which.min(cv_outcome$cv_loss)
cat("  Selected index:", best_idx2, "/", length(cv_outcome$lambda_grid), "\n")

if (best_idx2 == 1) {
  cat("  → BOUNDARY: Selected weakest lambda!\n\n")
} else if (best_idx2 == length(cv_outcome$lambda_grid)) {
  cat("  → BOUNDARY: Selected strongest lambda!\n\n")
} else {
  cat("  → Interior point (OK)\n\n")
}

cat("FINDING: ")
if (best_idx == 1 || best_idx2 == 1) {
  cat("CV hits lower boundary - grid too strong!\n")
  cat("Need: (log n / n) * c(0.05, 0.1, 0.25, 0.5, 1)\n")
} else if (best_idx == length(cv_ps$lambda_grid) || best_idx2 == length(cv_outcome$lambda_grid)) {
  cat("CV hits upper boundary - trees can't improve!\n")
} else {
  cat("CV selecting interior points - grid is reasonable\n")
}
