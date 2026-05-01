#!/usr/bin/env Rscript

# Combine results from all 120 jobs into single dataset

library(dplyr)
library(readr)

cat("==============================================\n")
cat("Combining Six-Approach Comparison Results\n")
cat("==============================================\n\n")

# Find all result files
fast_files <- list.files("results/raw", pattern = "^fast_approach_", full.names = TRUE)
medium_files <- list.files("results/raw", pattern = "^medium_approach_", full.names = TRUE)
msplit_files <- list.files("results/raw", pattern = "^msplit_approach_", full.names = TRUE)

all_files <- c(fast_files, medium_files, msplit_files)

cat(sprintf("Found %d result files:\n", length(all_files)))
cat(sprintf("  Fast approaches:   %d / 36\n", length(fast_files)))
cat(sprintf("  Medium approaches: %d / 24\n", length(medium_files)))
cat(sprintf("  M-split approach:  %d / 60\n", length(msplit_files)))
cat("\n")

if (length(all_files) != 120) {
  warning(sprintf("Expected 120 files, found %d. Some jobs may have failed.\n", length(all_files)))
}

if (length(all_files) == 0) {
  stop("No result files found. Check that jobs completed successfully.")
}

# Read and combine
cat("Reading and combining results...\n")
all_results <- lapply(all_files, function(f) {
  tryCatch({
    readRDS(f)
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s\n", f, e$message))
    NULL
  })
})

# Remove NULLs (failed reads)
all_results <- Filter(Negate(is.null), all_results)

combined <- bind_rows(all_results)

cat(sprintf("\nCombined %d replications\n", nrow(combined)))
cat(sprintf("  Approaches: %s\n", paste(sort(unique(combined$approach)), collapse = ", ")))
cat(sprintf("  DGPs: %s\n", paste(sort(unique(combined$dgp)), collapse = ", ")))
cat(sprintf("  Sample sizes: %s\n", paste(sort(unique(combined$n)), collapse = ", ")))
cat("\n")

# Check for errors
n_errors <- sum(is.na(combined$theta_hat))
if (n_errors > 0) {
  cat(sprintf("⚠️  Warning: %d replications have NA estimates\n", n_errors))
  error_summary <- combined %>%
    filter(is.na(theta_hat)) %>%
    group_by(approach_name, dgp_name, n) %>%
    summarise(count = n(), .groups = "drop")
  print(error_summary)
  cat("\n")
}

# Save combined results
cat("Saving combined results...\n")
saveRDS(combined, "results/combined/all_results.rds")
cat("  ✓ results/combined/all_results.rds\n")

# Generate summary tables
cat("\nGenerating summary tables...\n")

# 1. Inference summary
inference_summary <- combined %>%
  group_by(approach, approach_name, dgp, dgp_name, n) %>%
  summarise(
    reps = n(),
    mean_theta = mean(theta_hat, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    sd_theta = sd(theta_hat, na.rm = TRUE),
    rmse = sqrt(mean(bias^2, na.rm = TRUE)),
    mean_se = mean(se, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    coverage_lower = mean(coverage_lower, na.rm = TRUE),
    coverage_upper = mean(coverage_upper, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    median_ci_width = median(ci_width, na.rm = TRUE),
    n_errors = sum(is.na(theta_hat)),
    .groups = "drop"
  )

write_csv(inference_summary, "results/combined/summary_inference.csv")
cat("  ✓ results/combined/summary_inference.csv\n")

# 2. Timing summary
timing_summary <- combined %>%
  group_by(approach, approach_name) %>%
  summarise(
    total_reps = n(),
    total_time_hours = sum(elapsed_time, na.rm = TRUE) / 3600,
    mean_time_sec = mean(elapsed_time, na.rm = TRUE),
    median_time_sec = median(elapsed_time, na.rm = TRUE),
    sd_time_sec = sd(elapsed_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(approach)

write_csv(timing_summary, "results/combined/summary_timing.csv")
cat("  ✓ results/combined/summary_timing.csv\n")

# 3. Quick summary by approach
approach_summary <- inference_summary %>%
  group_by(approach, approach_name) %>%
  summarise(
    total_reps = sum(reps),
    mean_bias = mean(bias),
    mean_rmse = mean(rmse),
    mean_coverage = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(approach)

cat("\n==============================================\n")
cat("Summary by Approach\n")
cat("==============================================\n\n")
print(approach_summary, n = 6)

cat("\n==============================================\n")
cat("✓ Results combined successfully!\n")
cat("==============================================\n\n")

cat("Output files created:\n")
cat("  - results/combined/all_results.rds\n")
cat("  - results/combined/summary_inference.csv\n")
cat("  - results/combined/summary_timing.csv\n\n")

cat("Next step: Analyze results\n")
cat("  Rscript code/analyze_results.R\n\n")
