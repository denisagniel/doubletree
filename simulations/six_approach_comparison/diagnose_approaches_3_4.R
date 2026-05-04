#!/usr/bin/env Rscript
#
# Diagnostic Script for Approaches 3 and 4
#
# Investigates:
# 1. Why approach 4 (doubletree_singlefit) fails with Rashomon intersection
# 2. Why approach 3 (doubletree) has very high SE (~1.0 vs ~0.04)
#
# Created: 2026-05-04

suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

source("code/dgps.R")
source("code/estimators.R")

cat("===========================================\n")
cat("Diagnostic Analysis: Approaches 3 and 4\n")
cat("===========================================\n\n")

# Generate test data
set.seed(123)
data <- generate_dgp_simple(n = 500)
cat(sprintf("Data: n=%d, true ATT=%.4f\n\n", nrow(data$X), data$true_att))

# ============================================================================
# Issue 1: Approach 4 - Rashomon Intersection Failure
# ============================================================================

cat("--- Issue 1: Approach 4 Rashomon Intersection ---\n\n")

# Run doubletree with detailed diagnostics
cat("Running doubletree::estimate_att with use_rashomon=TRUE...\n")
result_rashomon <- doubletree::estimate_att(
  X = data$X,
  A = data$A,
  Y = data$Y,
  K = 5,
  regularization = 0.1,
  outcome_type = "binary",
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 0.05,
  verbose = TRUE  # Turn on verbosity
)

cat("\n")
cat("Result structure:\n")
cat(sprintf("  theta: %s\n", result_rashomon$theta))
cat(sprintf("  sigma: %s\n", result_rashomon$sigma))
cat(sprintf("  nuisance_fits available: %s\n", !is.null(result_rashomon$nuisance_fits)))

if (!is.null(result_rashomon$nuisance_fits)) {
  cat(sprintf("  cf_e available: %s\n", !is.null(result_rashomon$nuisance_fits$cf_e)))
  cat(sprintf("  cf_m0 available: %s\n", !is.null(result_rashomon$nuisance_fits$cf_m0)))

  if (!is.null(result_rashomon$nuisance_fits$cf_e)) {
    e_struct <- result_rashomon$nuisance_fits$cf_e$structure
    cat(sprintf("  e structure: %s\n", if (is.null(e_struct)) "NULL" else "Available"))
  }

  if (!is.null(result_rashomon$nuisance_fits$cf_m0)) {
    m0_struct <- result_rashomon$nuisance_fits$cf_m0$structure
    cat(sprintf("  m0 structure: %s\n", if (is.null(m0_struct)) "NULL" else "Available"))
  }
}

# Test with different rashomon_bound_multiplier values
cat("\n\nTesting different Rashomon bound multipliers:\n")
multipliers <- c(0.01, 0.05, 0.1, 0.2)

for (mult in multipliers) {
  cat(sprintf("\nmultiplier = %.2f: ", mult))

  result_test <- tryCatch({
    doubletree::estimate_att(
      X = data$X,
      A = data$A,
      Y = data$Y,
      K = 5,
      regularization = 0.1,
      outcome_type = "binary",
      use_rashomon = TRUE,
      rashomon_bound_multiplier = mult,
      verbose = FALSE
    )
  }, error = function(e) list(theta = NA, error = e$message))

  if (is.na(result_test$theta)) {
    cat("FAIL\n")
  } else {
    e_ok <- !is.null(result_test$nuisance_fits$cf_e$structure)
    m0_ok <- !is.null(result_test$nuisance_fits$cf_m0$structure)
    cat(sprintf("theta=%.4f, e_struct=%s, m0_struct=%s\n",
                result_test$theta, e_ok, m0_ok))
  }
}

# ============================================================================
# Issue 2: Approach 3 - High Standard Errors
# ============================================================================

cat("\n\n--- Issue 2: Approach 3 High Standard Errors ---\n\n")

# Compare SE calculation methods
cat("Running approach 1 (full_sample)...\n")
result1 <- estimate_att_fullsample(data$X, data$A, data$Y)
cat(sprintf("  theta=%.4f, se=%.4f\n", result1$theta, result1$se))

cat("\nRunning approach 2 (crossfit)...\n")
result2 <- estimate_att_crossfit(data$X, data$A, data$Y)
cat(sprintf("  theta=%.4f, se=%.4f\n", result2$theta, result2$se))

cat("\nRunning approach 3 (doubletree)...\n")
result3 <- estimate_att_doubletree(data$X, data$A, data$Y)
cat(sprintf("  theta=%.4f, se=%.4f\n", result3$theta, result3$se))

# Examine what doubletree returns
cat("\n\nInvestigating doubletree SE calculation:\n")
cat(sprintf("  doubletree returns: sigma = %.4f\n", result_rashomon$sigma))

# Check if this is the SE or something else
if (!is.null(result_rashomon$nuisance_fits)) {
  cat("\nChecking nuisance_fits structure:\n")
  cat(sprintf("  Names: %s\n", paste(names(result_rashomon$nuisance_fits), collapse=", ")))
}

# Check all return values from doubletree
cat("\nAll return values from doubletree::estimate_att:\n")
cat(sprintf("  Names: %s\n", paste(names(result_rashomon), collapse=", ")))

# Try to extract EIF scores if available
if (!is.null(result_rashomon$eif_scores)) {
  cat("\nEIF scores available!\n")
  cat(sprintf("  Length: %d\n", length(result_rashomon$eif_scores)))
  cat(sprintf("  Mean: %.4f (should be theta)\n", mean(result_rashomon$eif_scores)))
  cat(sprintf("  SD: %.4f\n", sd(result_rashomon$eif_scores)))

  # Compute SE manually
  se_manual <- sqrt(mean(result_rashomon$eif_scores^2) / length(result_rashomon$eif_scores))
  cat(sprintf("  Manual SE: %.4f\n", se_manual))
}

# ============================================================================
# Issue 3: Compare Predictions
# ============================================================================

cat("\n\n--- Issue 3: Compare Predictions Across Approaches ---\n\n")

cat("Prediction ranges:\n")
cat(sprintf("  Approach 1 - e_hat: [%.4f, %.4f]\n",
            min(result1$e_hat), max(result1$e_hat)))
cat(sprintf("  Approach 1 - m0_hat: [%.4f, %.4f]\n",
            min(result1$m0_hat), max(result1$m0_hat)))

cat(sprintf("\n  Approach 2 - e_hat: [%.4f, %.4f]\n",
            min(result2$e_hat), max(result2$e_hat)))
cat(sprintf("  Approach 2 - m0_hat: [%.4f, %.4f]\n",
            min(result2$m0_hat), max(result2$m0_hat)))

if (!all(is.na(result3$e_hat))) {
  cat(sprintf("\n  Approach 3 - e_hat: [%.4f, %.4f]\n",
              min(result3$e_hat, na.rm=TRUE), max(result3$e_hat, na.rm=TRUE)))
  cat(sprintf("  Approach 3 - m0_hat: [%.4f, %.4f]\n",
              min(result3$m0_hat, na.rm=TRUE), max(result3$m0_hat, na.rm=TRUE)))
}

# ============================================================================
# Summary
# ============================================================================

cat("\n\n===========================================\n")
cat("Summary\n")
cat("===========================================\n\n")

cat("Issue 1 (Approach 4 Rashomon Failure):\n")
if (is.na(result_rashomon$theta)) {
  cat("  ✗ Rashomon intersection failed even with multiplier=0.05\n")
  cat("  → May need larger multiplier or different DGP\n")
} else {
  e_has_struct <- !is.null(result_rashomon$nuisance_fits$cf_e$structure)
  m0_has_struct <- !is.null(result_rashomon$nuisance_fits$cf_m0$structure)

  if (e_has_struct && m0_has_struct) {
    cat("  ✓ Both structures found successfully\n")
  } else {
    cat("  ✗ Missing structures:\n")
    if (!e_has_struct) cat("    - Propensity structure missing\n")
    if (!m0_has_struct) cat("    - Outcome structure missing\n")
  }
}

cat("\nIssue 2 (High SEs):\n")
cat(sprintf("  Approach 1 SE: %.4f\n", result1$se))
cat(sprintf("  Approach 2 SE: %.4f\n", result2$se))
cat(sprintf("  Approach 3 SE: %.4f (%.1fx higher)\n",
            result3$se, result3$se / result1$se))

if (result3$se > 10 * result1$se) {
  cat("  → Suspiciously high (>10x) - likely bug in doubletree package\n")
} else if (result3$se > 2 * result1$se) {
  cat("  → Moderately high (2-10x) - investigate variance estimation\n")
} else {
  cat("  → Reasonable range - may be honest uncertainty\n")
}

cat("\n")
