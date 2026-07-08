#!/usr/bin/env Rscript

# Test approaches 4 and 6 with CV implementation

library(optimaltrees)
devtools::load_all("../..", quiet = TRUE)  # Load doubletree (go up 2 levels from six_approach_comparison)

source("code/dgps.R")
source("code/estimators.R")

cat("\n=== Testing Approaches 4 & 6 with CV ===\n\n")

# Generate test data
set.seed(456)
cat("Generating simple DGP with n=500...\n")
data <- generate_dgp_simple(n = 500)

cat(sprintf("  n = %d\n", nrow(data$X)))
cat(sprintf("  p = %d features\n", ncol(data$X)))
cat(sprintf("  true ATT = %.4f\n", data$true_att))
cat(sprintf("  treatment prop = %.3f\n\n", mean(data$A)))

# Test Approach 4
cat("\n--- Approach 4: Doubletree Averaged (with CV) ---\n")
start_time4 <- Sys.time()
result4 <- tryCatch({
  estimate_att_doubletree_averaged(X = data$X, A = data$A, Y = data$Y, K = 5)
}, error = function(e) {
  cat("\nERROR CAUGHT:\n")
  cat(as.character(e), "\n\n")
  list(theta = NA, se = NA, error = as.character(e))
})
elapsed4 <- as.numeric(Sys.time() - start_time4, units = "secs")

if (!is.null(result4$error)) {
  cat(sprintf("✗ FAILED after %.1f sec\n", elapsed4))
  cat(sprintf("Error: %s\n", result4$error))
} else {
  bias4 <- result4$theta - data$true_att
  z_score4 <- abs(bias4 / result4$se)
  covered4 <- z_score4 <= 1.96

  cat(sprintf("✓ SUCCESS after %.1f sec\n", elapsed4))
  cat(sprintf("  theta = %.4f (true = %.4f)\n", result4$theta, data$true_att))
  cat(sprintf("  bias = %.4f\n", bias4))
  cat(sprintf("  se = %.4f\n", result4$se))
  cat(sprintf("  z-score = %.2f\n", z_score4))
  cat(sprintf("  95%% CI covers truth: %s\n", covered4))
}

# Test Approach 6
cat("\n--- Approach 6: M-Split Averaged (with CV) ---\n")
start_time6 <- Sys.time()
result6 <- tryCatch({
  estimate_att_msplit_averaged(X = data$X, A = data$A, Y = data$Y, M = 5, K = 3)
}, error = function(e) {
  cat("\nERROR CAUGHT:\n")
  cat(as.character(e), "\n\n")
  list(theta = NA, se = NA, error = as.character(e))
})
elapsed6 <- as.numeric(Sys.time() - start_time6, units = "secs")

if (!is.null(result6$error)) {
  cat(sprintf("✗ FAILED after %.1f sec\n", elapsed6))
  cat(sprintf("Error: %s\n", result6$error))
} else {
  bias6 <- result6$theta - data$true_att
  z_score6 <- abs(bias6 / result6$se)
  covered6 <- z_score6 <= 1.96

  cat(sprintf("✓ SUCCESS after %.1f sec\n", elapsed6))
  cat(sprintf("  theta = %.4f (true = %.4f)\n", result6$theta, data$true_att))
  cat(sprintf("  bias = %.4f\n", bias6))
  cat(sprintf("  se = %.4f\n", result6$se))
  cat(sprintf("  z-score = %.2f\n", z_score6))
  cat(sprintf("  95%% CI covers truth: %s\n", covered6))
}

# Summary
cat("\n\n=== SUMMARY ===\n")
success_count <- sum(!is.na(c(result4$theta, result6$theta)))
cat(sprintf("Successful: %d/2\n", success_count))

if (success_count == 2) {
  cat("\n✓ Both approaches 4 and 6 are working with CV!\n")
  cat("  Ready to include in full simulation.\n")
} else {
  cat("\n✗ Some approaches still failing.\n")
  cat("  Need further debugging before deployment.\n")
}
