# Quick verification test for approaches 4 and 6
# Created: 2026-05-26

# Load development version
devtools::load_all(".", quiet = TRUE)

# Generate simple test data
set.seed(123)
n <- 300
X <- data.frame(
  x1 = rbinom(n, 1, 0.5),
  x2 = rbinom(n, 1, 0.5),
  x3 = rbinom(n, 1, 0.5)
)

# Treatment depends on X
A <- rbinom(n, 1, plogis(X$x1 + 0.5 * X$x2))

# Outcome depends on treatment and X
Y <- rbinom(n, 1, plogis(A * 0.3 + X$x1 + 0.5 * X$x3))

cat("Test data generated: n =", n, ", n_treated =", sum(A), "\n\n")

# ============================================================================
# Test Approach 4: Doubletree Averaged
# ============================================================================

cat("=== Testing Approach 4: Doubletree Averaged ===\n")

approach4_result <- tryCatch({
  estimate_att_doubletree_averaged(
    X = X,
    A = A,
    Y = Y,
    K = 3,  # Small K for quick test
    outcome_type = "binary",
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR in approach 4:", e$message, "\n")
  NULL
})

if (!is.null(approach4_result)) {
  cat("\nApproach 4 Results:\n")
  cat(sprintf("  ATT: %.4f\n", approach4_result$theta))
  cat(sprintf("  SE: %.4f\n", approach4_result$sigma))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
              approach4_result$ci_95[1], approach4_result$ci_95[2]))
  cat(sprintf("  K trees averaged: %d\n", approach4_result$n_trees_averaged))
  cat(sprintf("  Converged: %s\n", approach4_result$converged))
  cat("  Approach 4 PASSED\n")
} else {
  cat("  Approach 4 FAILED\n")
}

cat("\n")

# ============================================================================
# Test Approach 6: M-Split Averaged
# ============================================================================

cat("=== Testing Approach 6: M-Split Averaged ===\n")

approach6_result <- tryCatch({
  estimate_att_msplit_averaged(
    X = X,
    A = A,
    Y = Y,
    M = 3,  # Small M for quick test
    K = 3,  # Small K for quick test
    outcome_type = "binary",
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR in approach 6:", e$message, "\n")
  NULL
})

if (!is.null(approach6_result)) {
  cat("\nApproach 6 Results:\n")
  cat(sprintf("  ATT: %.4f\n", approach6_result$theta))
  cat(sprintf("  SE: %.4f\n", approach6_result$sigma))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
              approach6_result$ci_95[1], approach6_result$ci_95[2]))
  cat(sprintf("  M×K trees averaged: %d\n", approach6_result$n_trees_averaged))
  cat("  Approach 6 PASSED\n")
} else {
  cat("  Approach 6 FAILED\n")
}

cat("\n")

# ============================================================================
# Summary
# ============================================================================

cat("=== Summary ===\n")
if (!is.null(approach4_result) && !is.null(approach6_result)) {
  cat("✓ Both approaches 4 and 6 completed successfully\n")
  cat("\nComparison:\n")
  cat(sprintf("  Approach 4 ATT: %.4f (K=%d trees)\n",
              approach4_result$theta, approach4_result$K))
  cat(sprintf("  Approach 6 ATT: %.4f (M×K=%d trees)\n",
              approach6_result$theta, approach6_result$n_trees_averaged))
  cat(sprintf("  Difference: %.4f\n",
              abs(approach4_result$theta - approach6_result$theta)))
} else {
  cat("✗ One or both approaches failed\n")
  if (is.null(approach4_result)) cat("  - Approach 4 failed\n")
  if (is.null(approach6_result)) cat("  - Approach 6 failed\n")
}
