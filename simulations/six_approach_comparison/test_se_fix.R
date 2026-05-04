#!/usr/bin/env Rscript
#
# Test SE Fix for Approach 3 (doubletree)
#
# Verifies that fixing att_se() resolves the 23x SE inflation
#
# Created: 2026-05-04

suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

source("code/dgps.R")
source("code/estimators.R")

cat("=== Testing SE Fix ===\n\n")

set.seed(123)
data <- generate_dgp_simple(n = 500)
cat(sprintf("Data: n=%d, true ATT=%.4f\n\n", nrow(data$X), data$true_att))

# Test all approaches
cat("Approach 1 (full_sample): ")
result1 <- estimate_att_fullsample(data$X, data$A, data$Y)
cat(sprintf("theta=%.4f, se=%.4f\n", result1$theta, result1$se))

cat("Approach 2 (crossfit): ")
result2 <- estimate_att_crossfit(data$X, data$A, data$Y)
cat(sprintf("theta=%.4f, se=%.4f\n", result2$theta, result2$se))

cat("Approach 3 (doubletree): ")
result3 <- estimate_att_doubletree(data$X, data$A, data$Y)
cat(sprintf("theta=%.4f, se=%.4f\n", result3$theta, result3$se))

cat("Approach 5 (msplit): ")
result5 <- estimate_att_msplit(data$X, data$A, data$Y, M = 10)
cat(sprintf("theta=%.4f, se=%.4f\n", result5$theta, result5$se))

# Check ratios
ratio_3_to_1 <- result3$se / result1$se
ratio_5_to_1 <- result5$se / result1$se

cat("\n=== SE Ratios ===\n")
cat(sprintf("Approach 3 / Approach 1: %.2fx\n", ratio_3_to_1))
cat(sprintf("Approach 5 / Approach 1: %.2fx\n", ratio_5_to_1))

cat("\n=== Assessment ===\n")
if (ratio_3_to_1 < 2) {
  cat("✓ Approach 3 SE is now reasonable (<2x approach 1)\n")
} else if (ratio_3_to_1 < 5) {
  cat("⚠ Approach 3 SE is moderately high (2-5x approach 1)\n")
} else {
  cat("✗ Approach 3 SE is still suspiciously high (>5x approach 1)\n")
}

if (ratio_5_to_1 < 2) {
  cat("✓ Approach 5 SE is now reasonable (<2x approach 1)\n")
} else if (ratio_5_to_1 < 5) {
  cat("⚠ Approach 5 SE is moderately high (2-5x approach 1)\n")
} else {
  cat("✗ Approach 5 SE is still suspiciously high (>5x approach 1)\n")
}

cat("\n")
