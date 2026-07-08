#!/usr/bin/env Rscript

# Detailed test of Approach 3 (doubletree with CV)

library(optimaltrees)
library(doubletree)

source("code/dgps.R")
source("code/estimators.R")

cat("\n=== Testing Approach 3 (Doubletree with CV) ===\n\n")

# Try with larger n and more complex DGP
set.seed(123)
cat("Generating complex DGP with n=1000...\n")
data <- generate_dgp_complex(n = 1000)

cat(sprintf("  n = %d\n", nrow(data$X)))
cat(sprintf("  p = %d features\n", ncol(data$X)))
cat(sprintf("  true ATT = %.4f\n", data$true_att))
cat(sprintf("  treatment prop = %.3f\n\n", mean(data$A)))

cat("Running approach 3...\n")
start_time <- Sys.time()

result <- tryCatch({
  estimate_att_doubletree(X = data$X, A = data$A, Y = data$Y, K = 5)
}, error = function(e) {
  cat("\nERROR CAUGHT:\n")
  cat(as.character(e), "\n\n")
  list(theta = NA, se = NA, error = as.character(e))
})

elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

if (!is.null(result$error)) {
  cat(sprintf("\n✗ FAILED after %.1f sec\n", elapsed))
  cat(sprintf("Error: %s\n", result$error))
} else {
  bias <- result$theta - data$true_att
  z_score <- abs(bias / result$se)
  covered <- z_score <= 1.96

  cat(sprintf("\n✓ SUCCESS after %.1f sec\n", elapsed))
  cat(sprintf("  theta = %.4f (true = %.4f)\n", result$theta, data$true_att))
  cat(sprintf("  bias = %.4f\n", bias))
  cat(sprintf("  se = %.4f\n", result$se))
  cat(sprintf("  z-score = %.2f\n", z_score))
  cat(sprintf("  95%% CI covers truth: %s\n", covered))
}
