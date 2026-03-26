#!/usr/bin/env Rscript

#' Analyze Current DML-ATT Simulation Results
#'
#' Combines all completed .rds files and generates summary statistics
#' for current progress.

suppressMessages({
  library(dplyr)
  library(tidyr)
})

# Configuration
results_dir <- "/n/scratch/users/d/dma12/global-scholars/results/o2_primary"
output_file <- "current_results_summary.rds"

cat("DML-ATT Simulation Results Analysis\n")
cat("====================================\n\n")

# Find all result files
cat("Loading result files...\n")
result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
n_files <- length(result_files)

cat(sprintf("Found %d result files\n", n_files))

if (n_files == 0) {
  stop("No result files found in ", results_dir)
}

# Load and combine all results
cat("Combining results...\n")
all_results <- lapply(result_files, readRDS)
combined <- bind_rows(all_results)

cat(sprintf("Total replications: %d\n\n", nrow(combined)))

# Summary by configuration
cat("Results by Configuration:\n")
cat("-------------------------\n")
config_summary <- combined %>%
  group_by(dgp, n, method) %>%
  summarise(
    n_reps = n(),
    n_converged = sum(converged),
    convergence_rate = mean(converged),
    mean_theta = mean(theta, na.rm = TRUE),
    sd_theta = sd(theta, na.rm = TRUE),
    mean_sigma = mean(sigma, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dgp, n, method)

print(as.data.frame(config_summary), row.names = FALSE)

cat("\n\nDetailed Statistics (converged replications only):\n")
cat("--------------------------------------------------\n")

detailed_summary <- combined %>%
  filter(converged) %>%
  mutate(
    bias = theta - true_att,
    covered = (true_att >= ci_lower) & (true_att <= ci_upper),
    ci_width = ci_upper - ci_lower
  ) %>%
  group_by(dgp, n, method) %>%
  summarise(
    n_converged = n(),
    mean_theta = mean(theta),
    bias = mean(bias),
    rmse = sqrt(mean(bias^2)),
    coverage = mean(covered),
    mean_ci_width = mean(ci_width),
    .groups = "drop"
  ) %>%
  arrange(dgp, n, method)

print(as.data.frame(detailed_summary), row.names = FALSE)

# Save combined results
cat("\n\nSaving combined results to:", output_file, "\n")
saveRDS(list(
  combined = combined,
  config_summary = config_summary,
  detailed_summary = detailed_summary,
  n_files = n_files,
  timestamp = Sys.time()
), output_file)

# Key findings
cat("\n\nKey Findings:\n")
cat("-------------\n")

# Complete configurations (1000 reps)
complete_configs <- config_summary %>%
  filter(n_reps == 1000)
cat(sprintf("Complete configurations (1000 reps): %d/36\n", nrow(complete_configs)))

# Convergence issues
low_convergence <- config_summary %>%
  filter(convergence_rate < 0.95)
if (nrow(low_convergence) > 0) {
  cat("\nConfigurations with convergence < 95%:\n")
  print(as.data.frame(low_convergence %>% select(dgp, n, method, n_reps, convergence_rate)),
        row.names = FALSE)
}

# Coverage statistics (converged only)
coverage_stats <- detailed_summary %>%
  group_by(method) %>%
  summarise(
    mean_coverage = mean(coverage),
    min_coverage = min(coverage),
    max_coverage = max(coverage),
    .groups = "drop"
  )

cat("\n\nCoverage by Method (converged replications):\n")
print(as.data.frame(coverage_stats), row.names = FALSE)

# RMSE by method
rmse_stats <- detailed_summary %>%
  group_by(method) %>%
  summarise(
    mean_rmse = mean(rmse),
    .groups = "drop"
  ) %>%
  arrange(mean_rmse)

cat("\nRMSE by Method (converged replications):\n")
print(as.data.frame(rmse_stats), row.names = FALSE)

cat("\n\nAnalysis complete!\n")
cat("Results saved to:", output_file, "\n")
