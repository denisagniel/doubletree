#!/usr/bin/env Rscript

#' Quick test that Phase 2 DGPs work with batch script infrastructure

suppressMessages({
  library(optimaltrees)
  library(doubletree)
})

# Source all DGPs
source("dgps/dgps_smooth.R")
source("dgps/dgps_continuous.R")
source("dgps/dgps_phase2.R")

cat("Testing Phase 2 DGP Integration\n")
cat(strrep("=", 60), "\n\n", sep = "")

# Test each DGP
dgps <- list(
  list(name = "dgp7", func = generate_dgp7, outcome = "binary"),
  list(name = "dgp8", func = generate_dgp8, outcome = "continuous"),
  list(name = "dgp9", func = generate_dgp9, outcome = "binary")
)

n <- 400
tau <- 0.10
seed <- 123

all_passed <- TRUE

for (dgp_info in dgps) {
  cat(sprintf("Testing %s (expected outcome: %s)...\n",
              dgp_info$name, dgp_info$outcome))

  # Generate data
  d <- tryCatch({
    dgp_info$func(n = n, tau = tau, seed = seed)
  }, error = function(e) {
    cat(sprintf("  ERROR generating data: %s\n", e$message))
    NULL
  })

  if (is.null(d)) {
    all_passed <- FALSE
    next
  }

  # Check outcome type matches expectation
  is_binary <- all(d$Y %in% c(0, 1))
  actual_outcome <- if (is_binary) "binary" else "continuous"

  if (actual_outcome != dgp_info$outcome) {
    cat(sprintf("  ERROR: Expected %s outcome, got %s\n",
                dgp_info$outcome, actual_outcome))
    all_passed <- FALSE
  } else {
    cat(sprintf("  ✓ Outcome type correct (%s)\n", actual_outcome))
  }

  # Check data structure
  if (!all(c("X", "A", "Y", "tau", "true_att") %in% names(d))) {
    cat("  ERROR: Missing required fields\n")
    all_passed <- FALSE
  } else {
    cat("  ✓ Data structure correct\n")
  }

  # Check sample size
  if (nrow(d$X) != n || length(d$A) != n || length(d$Y) != n) {
    cat("  ERROR: Sample size mismatch\n")
    all_passed <- FALSE
  } else {
    cat("  ✓ Sample size correct\n")
  }

  cat("\n")
}

cat(strrep("=", 60), "\n", sep = "")
if (all_passed) {
  cat("✓ All Phase 2 DGPs pass integration tests\n")
  cat("✓ Ready for O2 batch processing\n")
  quit(status = 0)
} else {
  cat("✗ Some tests failed - fix before deploying\n")
  quit(status = 1)
}
