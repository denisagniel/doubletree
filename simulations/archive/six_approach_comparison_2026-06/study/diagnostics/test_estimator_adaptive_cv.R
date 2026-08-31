#!/usr/bin/env Rscript
# Quick test that estimators work with adaptive CV
# Single replication to validate deployment

library(doubletree)
library(optimaltrees)
source('code/dgps.R')

cat("\n")
cat("================================================================\n")
cat("Testing Estimators with Adaptive CV\n")
cat("================================================================\n\n")

set.seed(42)
n <- 500
data <- generate_dgp_complex(n = n)

cat("Sample: n =", n, "(complex DGP)\n")
cat("True ATT:", 0.15, "\n\n")

# Test 1: estimate_att with cv_regularization=TRUE (should use adaptive now)
cat("----------------------------------------------------------------\n")
cat("Test: estimate_att with CV regularization\n")
cat("----------------------------------------------------------------\n\n")

result <- estimate_att(
  X = data$X,
  A = data$A,
  Y = data$Y,
  K = 3,  # Fewer folds for speed
  outcome_type = "binary",
  cv_regularization = TRUE,  # This should now use adaptive CV
  cv_K = 3,
  verbose = TRUE,
  seed = 42
)

cat("\n")
cat("Results:\n")
cat("  ATT estimate:", round(result$theta, 4), "\n")
cat("  Standard error:", round(result$sigma, 4), "\n")
cat("  95% CI: [", round(result$ci_95[1], 4), ",", round(result$ci_95[2], 4), "]\n")
cat("  True ATT:", 0.15, "\n")
cat("  Bias:", round(result$theta - 0.15, 4), "\n")
cat("  CI covers truth:", result$ci_95[1] <= 0.15 && 0.15 <= result$ci_95[2], "\n")

cat("\n")
cat("================================================================\n")
cat("✓ Estimator works with adaptive CV!\n")
cat("================================================================\n\n")
