#!/usr/bin/env Rscript
#
# Test: Does Sample Size Help Rashomon Intersection?
#
# Hypothesis: Larger n → more stable structures → higher overlap
#
# Created: 2026-05-04

suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

source("code/dgps.R")
source("code/estimators.R")

cat("======================================\n")
cat("Rashomon Intersection vs Sample Size\n")
cat("======================================\n\n")

test_sizes <- c(500, 1000, 2000, 4000)
results <- list()

for (n in test_sizes) {
  cat(sprintf("\n--- Testing n = %d ---\n", n))

  set.seed(12345)  # Same seed for fair comparison
  data <- generate_dgp_simple(n = n)

  cat("Running doubletree_singlefit...\n")
  result <- tryCatch({
    estimate_att_doubletree_singlefit(
      data$X, data$A, data$Y,
      K = 5,
      regularization = 0.1
    )
  }, error = function(e) {
    list(theta = NA, error = e$message)
  })

  success <- !is.null(result$error) == FALSE && !is.na(result$theta)

  cat(sprintf("  Success: %s\n", success))
  if (success) {
    cat(sprintf("  theta = %.4f (true = %.4f)\n", result$theta, data$true_att))
    cat(sprintf("  se = %.4f\n", result$se))
    cat(sprintf("  e_structure: %s\n", !is.null(result$structures$e)))
    cat(sprintf("  m0_structure: %s\n", !is.null(result$structures$m0)))
  } else {
    cat(sprintf("  Error: %s\n", result$error))
  }

  results[[as.character(n)]] <- list(
    n = n,
    success = success,
    theta = if (success) result$theta else NA,
    has_e_struct = success && !is.null(result$structures$e),
    has_m0_struct = success && !is.null(result$structures$m0)
  )
}

# Summary
cat("\n\n======================================\n")
cat("Summary\n")
cat("======================================\n\n")

cat("Sample Size | Success | e_struct | m0_struct\n")
cat("----------- | ------- | -------- | ---------\n")
for (n in test_sizes) {
  r <- results[[as.character(n)]]
  cat(sprintf("%11d | %7s | %8s | %9s\n",
              r$n,
              ifelse(r$success, "YES", "NO"),
              ifelse(r$has_e_struct, "YES", "NO"),
              ifelse(r$has_m0_struct, "YES", "NO")))
}

# Conclusion
cat("\n\nConclusion:\n")
success_count <- sum(sapply(results, function(r) r$success))
if (success_count == 0) {
  cat("✗ No sample size succeeded in finding Rashomon intersection\n")
  cat("  → Suggests structural issue, not just sample size\n")
} else if (success_count < length(test_sizes)) {
  cat(sprintf("⚠ Only %d/%d sample sizes succeeded\n", success_count, length(test_sizes)))
  cat("  → Sample size helps but doesn't guarantee success\n")
} else {
  cat("✓ All sample sizes succeeded\n")
  cat("  → Sample size resolves the intersection problem\n")
}

cat("\n")
