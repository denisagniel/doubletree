#!/usr/bin/env Rscript

#' Test Continuous Outcome Support in Baseline Methods
#'
#' Verifies that forest and linear methods now work with DGP5 (continuous Y)

suppressMessages({
  library(doubletree)
})

# Source baseline methods
source("methods/method_forest.R")
source("methods/method_linear.R")

# Source DGPs
source("dgps/dgps_continuous.R")

cat("Testing Continuous Outcome Support\n")
cat(strrep("=", 60), "\n\n", sep = "")

# Generate DGP5 data (continuous outcome)
set.seed(123)
n <- 400
data <- generate_dgp_continuous_continuous(n, tau = 0.10)

cat("Generated DGP5 data:\n")
cat(sprintf("  n = %d\n", n))
cat(sprintf("  Y range: [%.2f, %.2f]\n", min(data$Y), max(data$Y)))
cat(sprintf("  Y is continuous: %s\n", !all(data$Y %in% c(0, 1))))
cat(sprintf("  True ATT: %.3f\n\n", data$tau))

# Test forest method
cat("Testing method_forest...\n")
cat(strrep("-", 60), "\n", sep = "")

fit_forest <- tryCatch({
  att_forest(data$X, data$A, data$Y, K = 5, seed = 456, verbose = TRUE)
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  NULL
})

if (!is.null(fit_forest)) {
  cat("\n✓ Forest method succeeded!\n")
  cat(sprintf("  Estimate: %.4f (True: %.3f)\n", fit_forest$theta, data$tau))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", fit_forest$ci[1], fit_forest$ci[2]))
  cat(sprintf("  Covers truth: %s\n",
              ifelse(data$tau >= fit_forest$ci[1] && data$tau <= fit_forest$ci[2],
                     "Yes", "No")))
} else {
  cat("\n✗ Forest method failed!\n")
}

cat("\n")

# Test linear method
cat("Testing method_linear...\n")
cat(strrep("-", 60), "\n", sep = "")

fit_linear <- tryCatch({
  att_linear(data$X, data$A, data$Y, K = 5, seed = 456, verbose = TRUE)
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  NULL
})

if (!is.null(fit_linear)) {
  cat("\n✓ Linear method succeeded!\n")
  cat(sprintf("  Estimate: %.4f (True: %.3f)\n", fit_linear$theta, data$tau))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", fit_linear$ci[1], fit_linear$ci[2]))
  cat(sprintf("  Covers truth: %s\n",
              ifelse(data$tau >= fit_linear$ci[1] && data$tau <= fit_linear$ci[2],
                     "Yes", "No")))
} else {
  cat("\n✗ Linear method failed!\n")
}

# Summary
cat("\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Summary\n")
cat(strrep("=", 60), "\n", sep = "")

if (!is.null(fit_forest) && !is.null(fit_linear)) {
  cat("✓ Both methods now support continuous outcomes!\n")
  cat("✓ Ready to resubmit DGP5 jobs to O2\n")
} else {
  cat("✗ Fix incomplete - please debug\n")
}
