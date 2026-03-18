# combine_batch_results.R
# Combine results from all batches into single summary table

results_dir <- "simulations/results_extended"

message("=== Combining Batch Results ===\n")

# Find all result files
result_files <- list.files(results_dir, pattern = "^result_.*\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in ", results_dir)
}

message("Found ", length(result_files), " result files")

# Expected: 64 configs (4 DGPs × 4 n × 4 epsilon)
expected_configs <- 4 * 4 * 4
if (length(result_files) < expected_configs) {
  warning("Expected ", expected_configs, " result files, found only ", length(result_files))
}

# Load and combine all results
message("\nLoading results...")
summary_rows <- list()

for (i in seq_along(result_files)) {
  if (i %% 10 == 0) message("  Loaded ", i, "/", length(result_files))

  result <- readRDS(result_files[i])

  # Extract config and metrics
  cfg <- result$config

  # One row per method
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "fold_specific", as.data.frame(result$fold_specific)
  )
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "rashomon", as.data.frame(result$rashomon)
  )
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "oracle", as.data.frame(result$oracle)
  )
}

message("Loaded ", length(result_files), " configs → ", length(summary_rows), " rows")

# Combine into single table
message("\nCombining into summary table...")
if (requireNamespace("dplyr", quietly = TRUE)) {
  summary_table <- dplyr::bind_rows(summary_rows)
} else {
  # Fallback: handle different column counts
  all_cols <- unique(unlist(lapply(summary_rows, names)))
  summary_rows_filled <- lapply(summary_rows, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }
    df[, all_cols]
  })
  summary_table <- do.call(rbind, summary_rows_filled)
}

message("Summary table: ", nrow(summary_table), " rows × ", ncol(summary_table), " columns")

# Save combined results
saveRDS(summary_table, file.path(results_dir, "simulation_summary.rds"))
write.csv(summary_table, file.path(results_dir, "simulation_summary.csv"), row.names = FALSE)

message("\n=== Results Combined ===")
message("Saved to:")
message("  ", file.path(results_dir, "simulation_summary.rds"))
message("  ", file.path(results_dir, "simulation_summary.csv"))

# Quick summary
message("\n=== Quick Summary ===")
cat("\nConfigurations by DGP:\n")
print(table(summary_table$dgp[summary_table$method == "rashomon"]))

cat("\nConfigurations by sample size:\n")
print(table(summary_table$n[summary_table$method == "rashomon"]))

cat("\nMean coverage by method:\n")
print(aggregate(coverage_95 ~ method, data = summary_table, FUN = mean, na.rm = TRUE))

message("\nRun analysis with:")
message("  Rscript simulations/analyze_results.R")
