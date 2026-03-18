# Combine results from individual DGP batches into single dataset

library(dplyr)

cat("Combining batch results...\n\n")

# Find batch result directories
batch_dirs <- list.dirs("results", recursive = FALSE, full.names = TRUE)
batch_dirs <- batch_dirs[grepl("dgp[123]_batch", batch_dirs)]

if (length(batch_dirs) == 0) {
  stop("No batch results found in results/ directory")
}

cat("Found", length(batch_dirs), "batch directories:\n")
for (dir in batch_dirs) {
  cat("  -", dir, "\n")
}
cat("\n")

# Load and combine results
all_results <- list()
for (dir in batch_dirs) {
  dgp_name <- sub(".*/(dgp[123])_batch.*", "\\1", dir)
  rds_file <- file.path(dir, paste0(dgp_name, "_results.rds"))

  if (file.exists(rds_file)) {
    cat("Loading", rds_file, "...")
    results <- readRDS(rds_file)
    all_results[[dgp_name]] <- results
    cat(" ", nrow(results), "rows\n")
  } else {
    warning("Missing: ", rds_file)
  }
}

if (length(all_results) == 0) {
  stop("No results loaded")
}

# Combine into single dataframe
combined <- dplyr::bind_rows(all_results)

cat("\nCombined results:\n")
cat("  Total rows:", nrow(combined), "\n")
cat("  DGPs:", paste(unique(combined$dgp), collapse = ", "), "\n")
cat("  Methods:", paste(unique(combined$method), collapse = ", "), "\n")
cat("  Sample sizes:", paste(sort(unique(combined$n)), collapse = ", "), "\n")
cat("  Convergence rate:", sprintf("%.1f%%", 100 * mean(combined$converged, na.rm = TRUE)), "\n")

# Save combined results
output_dir <- sprintf("results/primary_%s", Sys.Date())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(combined, file.path(output_dir, "simulation_results.rds"))
cat("\nSaved to:", file.path(output_dir, "simulation_results.rds"), "\n")

# Summary statistics
summary_stats <- combined %>%
  filter(converged) %>%
  group_by(dgp, method, n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att, na.rm = TRUE),
    rmse = sqrt(mean((theta - true_att)^2, na.rm = TRUE)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(output_dir, "summary_stats.csv"), row.names = FALSE)
cat("Saved summary to:", file.path(output_dir, "summary_stats.csv"), "\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("Summary Statistics\n")
cat(strrep("=", 70), "\n\n")
print(summary_stats, n = Inf)

cat("\n✓ Batch results combined successfully\n")
