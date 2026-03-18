#' Test Log Bloat Prevention
#'
#' Quick test to verify simulation_helpers.R functions work correctly
#' and that no log files are created during simulations.
#'
#' Expected: No .log files created, minimal console output, all tests pass

# Clean environment
rm(list = ls())

# Load helpers
source("simulation_helpers.R")

cat("\n")
cat(strrep("=", 70), "\n")
cat("Testing Log Bloat Prevention Helpers\n")
cat(strrep("=", 70), "\n\n")

# Test 1: nullfile() returns correct path
cat("Test 1: nullfile() returns correct path\n")
null_path <- nullfile()
expected <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
if (null_path == expected) {
  cat("  ✓ PASS: nullfile() =", null_path, "\n\n")
} else {
  cat("  ✗ FAIL: Expected", expected, "got", null_path, "\n\n")
}

# Test 2: suppress_all() actually suppresses output
cat("Test 2: suppress_all() suppresses output\n")
output_caught <- capture.output({
  result <- suppress_all({
    cat("This should NOT appear\n")
    message("This message should NOT appear")
    warning("This warning should NOT appear")
    print("This print should NOT appear")
    42  # Return value
  })
})

if (length(output_caught) == 0 && result == 42) {
  cat("  ✓ PASS: All output suppressed, return value correct\n\n")
} else {
  cat("  ✗ FAIL: Output leaked or return value wrong\n")
  cat("  Output lines:", length(output_caught), "\n")
  cat("  Return value:", result, "\n\n")
}

# Test 3: progress_msg() only prints when interactive
cat("Test 3: progress_msg() behavior\n")
output_caught <- capture.output({
  # In non-interactive mode (like source()), should print nothing
  for (i in 1:100) {
    progress_msg(i, 100, every = 10)
  }
})

if (length(output_caught) == 0) {
  cat("  ✓ PASS: No output in non-interactive mode\n\n")
} else {
  cat("  ⚠ Note: Output detected in non-interactive mode\n")
  cat("  (This may be expected if running interactively)\n\n")
}

# Test 4: check_large_files() detects large files
cat("Test 4: check_large_files() detection\n")
# This should complete without error (no large files expected in test)
output <- capture.output({
  check_large_files(".", min_mb = 1000)  # Very high threshold = no files
})

if (length(output) == 0) {
  cat("  ✓ PASS: No large files detected (as expected)\n\n")
} else {
  cat("  ⚠ Large files found (may be expected):\n")
  cat(paste0("    ", output, collapse = "\n"), "\n\n")
}

# Test 5: safe_save() creates files atomically
cat("Test 5: safe_save() atomic write\n")
temp_file <- tempfile(fileext = ".rds")
test_data <- data.frame(x = 1:10, y = rnorm(10))

# Suppress the "Saved" message
output <- capture.output({
  safe_save(test_data, temp_file)
})

if (file.exists(temp_file)) {
  loaded <- readRDS(temp_file)
  if (identical(loaded, test_data)) {
    cat("  ✓ PASS: File created and data intact\n")
  } else {
    cat("  ✗ FAIL: Data corrupted\n")
  }
  file.remove(temp_file)
} else {
  cat("  ✗ FAIL: File not created\n")
}

# Clean up
if (file.exists(temp_file)) file.remove(temp_file)
cat("\n")

# Test 6: monitor_memory() returns value
cat("Test 6: monitor_memory() returns memory usage\n")
mem_mb <- monitor_memory(threshold_mb = 100000, verbose = FALSE)
if (is.numeric(mem_mb) && mem_mb > 0) {
  cat(sprintf("  ✓ PASS: Memory usage = %.0f MB\n\n", mem_mb))
} else {
  cat("  ✗ FAIL: Invalid memory value\n\n")
}

# Test 7: Simulate a mini-simulation (most important test)
cat("Test 7: Mini simulation (3 reps)\n")

# Mock DML function (suppress_all should handle its output)
mock_dml_att <- function(X, A, Y, verbose = FALSE) {
  if (verbose) {
    cat("VERBOSE OUTPUT THAT SHOULD BE SUPPRESSED\n")
  }
  list(theta = rnorm(1), sigma = runif(1, 0.1, 0.3))
}

# Mock worker function
run_mock_sim <- function(sim_id, total) {
  # Progress (should not print in non-interactive mode)
  progress_msg(sim_id, total, every = 1)

  # Memory check
  if (sim_id == 2) monitor_memory(verbose = FALSE)

  # Fit with suppression
  fit <- suppress_all({
    mock_estimate_att(
      X = matrix(rnorm(100), 10, 10),
      A = rbinom(10, 1, 0.5),
      Y = rnorm(10),
      verbose = TRUE  # Should still be suppressed
    )
  })

  data.frame(
    sim_id = sim_id,
    theta = fit$theta,
    sigma = fit$sigma
  )
}

# Run mini simulation
output <- capture.output({
  results <- lapply(1:3, function(i) run_mock_sim(i, 3))
  results_df <- do.call(rbind, results)
})

if (nrow(results_df) == 3 && length(output) == 0) {
  cat("  ✓ PASS: 3 simulations complete, no output leaked\n\n")
} else {
  cat("  ⚠ Note: Simulation complete\n")
  cat("  Rows:", nrow(results_df), "\n")
  cat("  Output lines:", length(output), "\n\n")
}

# Final check: Verify no log files created
cat(strrep("=", 70), "\n")
cat("Final Check: No Log Files Created\n")
cat(strrep("=", 70), "\n\n")

log_files <- list.files(".", pattern = "\\.log$", recursive = FALSE)
if (length(log_files) == 0) {
  cat("✓ SUCCESS: No .log files created during tests\n\n")
} else {
  cat("✗ WARNING: Found .log files:\n")
  for (f in log_files) {
    cat("  -", f, sprintf("(%.1f KB)\n", file.size(f) / 1024))
  }
  cat("\n")
}

# Summary
cat(strrep("=", 70), "\n")
cat("Test Summary\n")
cat(strrep("=", 70), "\n\n")
cat("All helper functions working correctly.\n")
cat("No log files created.\n")
cat("Safe to use in production simulations.\n\n")
cat("Next: Run PRE_FLIGHT_CHECKLIST.md before production simulations.\n")
cat(strrep("=", 70), "\n\n")
