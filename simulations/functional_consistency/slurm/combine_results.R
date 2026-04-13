#!/usr/bin/env Rscript
# Combine results from all simulation replications

cat("========================================\n")
cat("Combining Functional Consistency Results\n")
cat("========================================\n\n")

# Find result files
result_dir <- "results"
rds_files <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

cat("Found", length(rds_files), "result files\n\n")

if (length(rds_files) == 0) {
  stop("No result files found in ", result_dir)
}

# Load and combine
cat("Loading results...\n")
results_list <- lapply(rds_files, function(f) {
  tryCatch(
    readRDS(f),
    error = function(e) {
      warning("Failed to load: ", f, " - ", e$message)
      NULL
    }
  )
})

# Remove NULLs (failed loads)
results_list <- results_list[!sapply(results_list, is.null)]

cat("Successfully loaded", length(results_list), "replications\n\n")

# Combine into single data frame
results <- do.call(rbind, results_list)

# Check completeness
cat("Checking completeness...\n")
expected_configs <- 5 * 3 * 3 * 3  # n × dgp × method × K = 135
expected_reps_per_config <- 500
expected_total <- expected_configs * expected_reps_per_config

config_summary <- aggregate(
  seed ~ n + dgp + method + K,
  data = results,
  FUN = length
)
names(config_summary)[5] <- "n_reps"

cat("\nReplications per configuration:\n")
print(config_summary)

cat("\nTotal replications:", nrow(results), "/", expected_total, "\n")

incomplete <- config_summary[config_summary$n_reps < expected_reps_per_config, ]
if (nrow(incomplete) > 0) {
  cat("\nWARNING: Incomplete configurations:\n")
  print(incomplete)
}

# Save combined results
output_file <- "results/combined_fc_simulations.rds"
saveRDS(results, file = output_file)
cat("\nCombined results saved to:", output_file, "\n")

# Also save as CSV
csv_file <- "results/combined_fc_simulations.csv"
write.csv(results, file = csv_file, row.names = FALSE)
cat("CSV version saved to:", csv_file, "\n")

# Summary statistics
cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n\n")

# Overall coverage by method
cat("Coverage by method:\n")
coverage_by_method <- aggregate(coverage ~ method, data = results, FUN = mean)
print(coverage_by_method)
cat("\n")

# Coverage by method and n
cat("Coverage by method and n:\n")
coverage_by_n <- aggregate(coverage ~ method + n, data = results, FUN = mean)
print(coverage_by_n)
cat("\n")

# Functional consistency by method
cat("Functional consistency (max_diff) by method:\n")
fc_summary <- aggregate(
  cbind(max_diff_e, max_diff_m0) ~ method,
  data = results,
  FUN = function(x) c(mean = mean(x), max = max(x))
)
print(fc_summary)
cat("\n")

# Bias and standardized bias by method and n
cat("Bias and standardized bias by method and n:\n")
bias_summary <- aggregate(
  cbind(bias, standardized_bias) ~ method + n,
  data = results,
  FUN = mean
)
print(bias_summary)
cat("\n")

cat("========================================\n")
cat("Done!\n")
cat("========================================\n")
