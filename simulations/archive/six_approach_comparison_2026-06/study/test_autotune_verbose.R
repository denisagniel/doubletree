#!/usr/bin/env Rscript
# Test auto-tuning with VERBOSE output to see what epsilon values are tried

library(optimaltrees)
library(doubletree)

# Source simulation code
source("code/dgps.R")

cat("==============================================\n")
cat("Testing Auto-Tuning with Verbose Output\n")
cat("==============================================\n\n")

# Test with DGP 3 (complex), n=500
set.seed(12345)
n <- 500

cat("Generating data: DGP 3 (complex), n=500\n\n")
data <- generate_dgp_complex(n)

cat("Running with auto_tune_intersecting=TRUE, verbose=TRUE\n")
cat("This will show all epsilon values attempted\n\n")

# Call package function directly (not wrapper) to control verbosity
start_time <- Sys.time()

result <- tryCatch({
  doubletree::estimate_att_doubletree_averaged(
    X = data$X,
    A = data$A,
    Y = data$Y,
    K = 5,
    outcome_type = "binary",
    cv_regularization = TRUE,
    auto_tune_intersecting = TRUE,  # Enable auto-tuning
    verbose = TRUE  # VERBOSE to see what happens
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
  cat(sprintf("  Standard error: %.4f\n", result$sigma))
  cat(sprintf("  True ATT: %.4f\n", data$true_att))
  cat(sprintf("  Time: %.2f seconds\n", elapsed))

  # Check if epsilon_n is in result
  if (!is.null(result$epsilon_n)) {
    cat(sprintf("  Final epsilon_n used: %.4f\n", result$epsilon_n))
  }
}

cat("\n==============================================\n")
