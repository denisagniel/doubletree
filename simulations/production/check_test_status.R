# Quick check of test output without waiting for completion
task_file <- "/private/tmp/claude-1141097072/-Users-dagniel-Library-CloudStorage-OneDrive-RANDCorporation-rprojects-global-scholars/tasks/bi9d5ri0t.output"

if (file.exists(task_file)) {
  output <- readLines(task_file)

  # Check for key indicators
  has_results_header <- any(grepl("=== RESULTS ===", output, fixed = TRUE))
  has_convergence <- any(grepl("Convergence:", output))
  has_success <- any(grepl("ALL TESTS PASS", output))
  has_failure <- any(grepl("FAILED|ERROR", output, ignore.case = TRUE))

  cat("Status check:\n")
  cat("  Results header found:", has_results_header, "\n")
  cat("  Convergence reported:", has_convergence, "\n")
  cat("  Success indicator:", has_success, "\n")
  cat("  Failure indicator:", has_failure, "\n")

  if (has_results_header) {
    cat("\n=== Results Section ===\n")
    results_idx <- which(grepl("=== RESULTS ===", output, fixed = TRUE))
    if (length(results_idx) > 0) {
      # Print next 30 lines after results header
      end_idx <- min(results_idx + 30, length(output))
      cat(output[results_idx:end_idx], sep = "\n")
    }
  }

  cat("\n=== Last 10 lines ===\n")
  cat(tail(output, 10), sep = "\n")
} else {
  cat("Task output file not found yet\n")
}
