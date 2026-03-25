#!/usr/bin/env Rscript

#' Combine individual replication results from O2 array jobs
#'
#' Reads all .rds files from scratch directory and combines into single dataset
#' Saves combined results to final output directory

library(dplyr)

# Output directory on scratch (where individual replications are saved)
SCRATCH_DIR <- file.path("/n/scratch/users",
                         substr(Sys.getenv("USER"), 1, 1),
                         Sys.getenv("USER"),
                         "global-scholars/results/o2_primary")

# Final output directory (can be on home or data directory)
OUTPUT_DIR <- "results/o2_primary_combined"

cat("\n========================================\n")
cat("Combining DML-ATT Simulation Results\n")
cat("========================================\n\n")

cat("Reading from:", SCRATCH_DIR, "\n")
cat("Writing to:", OUTPUT_DIR, "\n\n")

# Check if scratch directory exists
if (!dir.exists(SCRATCH_DIR)) {
  stop("Scratch directory does not exist: ", SCRATCH_DIR)
}

# Find all .rds files
rds_files <- list.files(SCRATCH_DIR, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No .rds files found in ", SCRATCH_DIR)
}

cat("Found", length(rds_files), "replication files\n")
cat("Expected: 18,000 replications\n")
cat("Completion:", round(100 * length(rds_files) / 18000, 1), "%\n\n")

# Read all files and combine
cat("Reading files...\n")
results_list <- lapply(seq_along(rds_files), function(i) {
  if (i %% 1000 == 0) {
    cat("  Read", i, "/", length(rds_files), "\n")
  }
  tryCatch({
    readRDS(rds_files[i])
  }, error = function(e) {
    warning("Failed to read ", rds_files[i], ": ", e$message)
    NULL
  })
})

# Remove failed reads
results_list <- Filter(Negate(is.null), results_list)

cat("\nSuccessfully read", length(results_list), "files\n")

# Combine into single data frame
cat("Combining...\n")
results <- do.call(rbind, results_list)

cat("Combined dataset:\n")
cat("  Rows:", nrow(results), "\n")
cat("  Columns:", ncol(results), "\n\n")

# Check for duplicates
duplicates <- results %>%
  group_by(dgp, n, method, replication) %>%
  filter(n() > 1) %>%
  nrow()

if (duplicates > 0) {
  warning("Found ", duplicates, " duplicate replications")
  cat("Removing duplicates...\n")
  results <- results %>%
    distinct(dgp, n, method, replication, .keep_all = TRUE)
  cat("After deduplication:", nrow(results), "rows\n\n")
}

# Convergence summary
cat("Convergence summary:\n")
conv_summary <- results %>%
  group_by(dgp, method) %>%
  summarize(
    total = n(),
    converged = sum(converged),
    pct = 100 * mean(converged),
    .groups = "drop"
  ) %>%
  arrange(dgp, method)

print(as.data.frame(conv_summary))
cat("\n")

# Overall convergence
overall_conv <- mean(results$converged)
cat(sprintf("Overall convergence: %.1f%%\n\n", 100 * overall_conv))

# Performance summary (converged only)
if (sum(results$converged) > 0) {
  cat("Performance summary (converged replications):\n")
  perf_summary <- results %>%
    filter(converged) %>%
    group_by(dgp, method, n) %>%
    summarize(
      reps = n(),
      bias = mean(theta - true_att),
      rmse = sqrt(mean((theta - true_att)^2)),
      coverage = 100 * mean(ci_lower <= true_att & ci_upper >= true_att),
      mean_ci_width = mean(ci_upper - ci_lower),
      .groups = "drop"
    )

  print(as.data.frame(perf_summary), row.names = FALSE)
  cat("\n")
}

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Save combined results
output_file <- file.path(OUTPUT_DIR, sprintf("combined_results_%s.rds", Sys.Date()))
saveRDS(results, output_file)
cat("Saved combined results to:", output_file, "\n")

# Save CSV for inspection
csv_file <- file.path(OUTPUT_DIR, sprintf("combined_results_%s.csv", Sys.Date()))
write.csv(results, csv_file, row.names = FALSE)
cat("Saved CSV to:", csv_file, "\n")

# Save summaries
summary_file <- file.path(OUTPUT_DIR, sprintf("summary_stats_%s.csv", Sys.Date()))
if (exists("perf_summary")) {
  write.csv(perf_summary, summary_file, row.names = FALSE)
  cat("Saved summary statistics to:", summary_file, "\n")
}

# Save metadata
metadata <- list(
  timestamp = Sys.time(),
  n_replications = nrow(results),
  n_converged = sum(results$converged),
  convergence_rate = mean(results$converged),
  n_files_read = length(rds_files),
  scratch_dir = SCRATCH_DIR,
  output_dir = OUTPUT_DIR
)
metadata_file <- file.path(OUTPUT_DIR, sprintf("metadata_%s.rds", Sys.Date()))
saveRDS(metadata, metadata_file)
cat("Saved metadata to:", metadata_file, "\n")

cat("\n========================================\n")
cat("Combination complete!\n")
cat("========================================\n\n")
