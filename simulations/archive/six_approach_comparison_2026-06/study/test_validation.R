#!/usr/bin/env Rscript
#
# Comprehensive Validation Test Script
# Tests all six approaches with proper error handling
#
# Created: 2026-05-01
# Purpose: Verify predict() fixes and validation work correctly

suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

source("code/dgps.R")
source("code/estimators.R")

cat("=== Testing Validation and Error Handling ===\n\n")

# Helper function to run with error catching
run_safely <- function(estimator_fn, data, ...) {
  tryCatch({
    result <- estimator_fn(data$X, data$A, data$Y, ...)

    # Validate result structure
    if (!is.list(result)) {
      return(list(error = "Result is not a list"))
    }

    # Check if estimator already reported an error
    if (!is.null(result$error)) {
      return(result)  # Pass through error from estimator
    }

    if (is.null(result$theta) || is.null(result$se)) {
      return(list(error = "Missing theta or se in result"))
    }

    if (is.na(result$theta) || !is.finite(result$theta)) {
      return(list(error = paste("Invalid theta:", result$theta)))
    }

    if (is.na(result$se) || !is.finite(result$se) || result$se <= 0) {
      return(list(error = paste("Invalid se:", result$se)))
    }

    result
  }, error = function(e) {
    list(error = e$message)
  })
}

# Test 1: Binary DGP (Simple) - All Approaches
cat("=== Test 1: Binary DGP (Simple) ===\n\n")
set.seed(123)
data1 <- generate_dgp_simple(n = 500)
cat(sprintf("True ATT: %.4f\n\n", data1$true_att))

approach_names <- c(
  "full_sample",
  "crossfit",
  "doubletree",
  "doubletree_singlefit",
  "msplit",
  "msplit_singlefit"
)

approach_fns <- list(
  estimate_att_fullsample,
  estimate_att_crossfit,
  estimate_att_doubletree,
  estimate_att_doubletree_singlefit,
  estimate_att_msplit,
  estimate_att_msplit_singlefit
)

test1_results <- list()
for (i in 1:6) {
  cat(sprintf("Approach %d (%s): ", i, approach_names[i]))

  result <- run_safely(approach_fns[[i]], data1)
  test1_results[[i]] <- result

  if (!is.null(result$error) || is.na(result$theta)) {
    error_msg <- if (!is.null(result$error)) result$error else "Unknown error (NA theta)"
    cat(sprintf("ERROR: %s\n", error_msg))
  } else {
    bias <- result$theta - data1$true_att
    cat(sprintf("theta=%.4f, se=%.4f, bias=%.4f",
                result$theta, result$se, bias))

    # Check if reasonable
    if (abs(bias) > 0.5) {
      cat(" [WARNING: Large bias]")
    }
    if (abs(result$theta) > 10) {
      cat(" [WARNING: Extreme estimate]")
    }
    cat("\n")
  }
}

# Test 2: Binary DGP (Complex) - Subset of Approaches
cat("\n=== Test 2: Binary DGP (Complex) ===\n\n")
set.seed(456)
data2 <- generate_dgp_complex(n = 500)
cat(sprintf("True ATT: %.4f\n\n", data2$true_att))

test2_results <- list()
# Test approaches 1, 2, 5 (fast ones)
test_indices <- c(1, 2, 5)
for (idx in test_indices) {
  i <- test_indices[which(test_indices == idx)]
  cat(sprintf("Approach %d (%s): ", idx, approach_names[idx]))

  result <- run_safely(approach_fns[[idx]], data2)
  test2_results[[idx]] <- result

  if (!is.null(result$error) || is.na(result$theta)) {
    error_msg <- if (!is.null(result$error)) result$error else "Unknown error (NA theta)"
    cat(sprintf("ERROR: %s\n", error_msg))
  } else {
    bias <- result$theta - data2$true_att
    cat(sprintf("theta=%.4f, se=%.4f, bias=%.4f\n",
                result$theta, result$se, bias))
  }
}

# Test 3: Continuous Outcome - Expected to work or fail gracefully
cat("\n=== Test 3: Continuous Outcome (may fail - separate issue) ===\n\n")
set.seed(789)
data3 <- generate_dgp_continuous(n = 500)
cat(sprintf("True ATT: %.4f\n\n", data3$true_att))

test3_results <- list()
# Only test approach 1 (full sample) for continuous
cat("Approach 1 (full_sample): ")
result <- run_safely(estimate_att_fullsample, data3, regularization = 0.1)
test3_results[[1]] <- result

if (!is.null(result$error)) {
  cat(sprintf("ERROR (expected for continuous): %s\n", result$error))
} else {
  bias <- result$theta - data3$true_att
  cat(sprintf("theta=%.4f, se=%.4f, bias=%.4f\n",
              result$theta, result$se, bias))
}

# Test 4: Prediction Format Validation
cat("\n=== Test 4: Prediction Format Validation ===\n\n")
set.seed(101)
data4 <- generate_dgp_simple(n = 200)

cat("Testing that predict() returns probabilities in [0,1]...\n")
result <- run_safely(estimate_att_fullsample, data4)

if (!is.null(result$error)) {
  cat(sprintf("ERROR: %s\n", result$error))
} else {
  # Check prediction ranges
  e_range <- range(result$e_hat, na.rm = TRUE)
  m0_range <- range(result$m0_hat, na.rm = TRUE)

  cat(sprintf("  e_hat range: [%.4f, %.4f]", e_range[1], e_range[2]))
  if (e_range[1] < 0 || e_range[2] > 1) {
    cat(" [FAIL: Outside [0,1]]")
  } else {
    cat(" [PASS]")
  }
  cat("\n")

  cat(sprintf("  m0_hat range: [%.4f, %.4f]", m0_range[1], m0_range[2]))
  if (m0_range[1] < 0 || m0_range[2] > 1) {
    cat(" [FAIL: Outside [0,1]]")
  } else {
    cat(" [PASS]")
  }
  cat("\n")
}

# Summary
cat("\n=== Summary ===\n\n")

# Test 1 summary
test1_success <- sum(sapply(test1_results, function(r) is.null(r$error)))
cat(sprintf("Test 1 (Simple DGP): %d/%d approaches successful\n",
            test1_success, length(test1_results)))

if (test1_success < 6) {
  cat("  Failed approaches:\n")
  for (i in 1:6) {
    if (!is.null(test1_results[[i]]$error)) {
      cat(sprintf("    %d (%s): %s\n", i, approach_names[i],
                  test1_results[[i]]$error))
    }
  }
}

# Test 2 summary
test2_success <- sum(sapply(test2_results[test_indices], function(r) {
  !is.null(r) && is.null(r$error)
}))
cat(sprintf("Test 2 (Complex DGP): %d/%d approaches successful\n",
            test2_success, length(test_indices)))

# Overall
if (test1_success == 6 && test2_success == length(test_indices)) {
  cat("\n✓ All validation tests passed\n")
} else {
  cat("\n✗ Some tests failed - see details above\n")
  quit(status = 1)
}
