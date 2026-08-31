#!/usr/bin/env Rscript
# Test approach 4 with auto_tune_intersecting enabled
# Before resubmitting to cluster

library(optimaltrees)
library(doubletree)

# Source simulation code
source("code/dgps.R")
source("code/estimators.R")

cat("==============================================\n")
cat("Testing Approach 4 with Auto-Tuning\n")
cat("==============================================\n\n")

# Test with DGP 3 (complex), n=500 - the one that was failing
set.seed(12345)
n <- 500

cat("Generating data: DGP 3 (complex), n=500\n")
data <- generate_dgp_complex(n)

cat("Running approach 4 (doubletree_averaged) with auto_tune_intersecting=TRUE\n\n")

# Time it
start_time <- Sys.time()

result <- tryCatch({
  estimate_att_doubletree_averaged(
    X = data$X,
    A = data$A,
    Y = data$Y,
    K = 5,
    regularization = 0.1  # unused but passed for compatibility
  )
}, error = function(e) {
  list(error = e$message)
})

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\n==============================================\n")
cat("Results\n")
cat("==============================================\n")

if (!is.null(result$error)) {
  cat("✗ FAILED:\n")
  cat("  Error:", result$error, "\n")
} else {
  cat("✓ SUCCESS:\n")
  cat(sprintf("  ATT estimate: %.4f\n", result$theta))
  cat(sprintf("  Standard error: %.4f\n", result$se))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
              result$theta - 1.96*result$se,
              result$theta + 1.96*result$se))
  cat(sprintf("  True ATT: %.4f\n", data$true_att))
  cat(sprintf("  Bias: %.4f\n", result$theta - data$true_att))
  cat(sprintf("  Time: %.2f seconds\n", elapsed))

  if (!is.null(result$averaged_trees)) {
    cat("\n  Averaged trees created successfully\n")
  }
}

cat("\n==============================================\n")
cat("Test complete\n")
cat("==============================================\n")

# Return invisibly for interactive use
invisible(result)
