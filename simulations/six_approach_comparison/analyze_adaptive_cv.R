#!/usr/bin/env Rscript
# Analyze Adaptive CV Validation Results
# Combines batch results and compares to baseline

library(dplyr)
library(tidyr)

# Configuration
true_att <- 0.15
results_dir <- "/n/scratch/users/d/dagniel/global-scholars/adaptive_cv_validation"

# If running locally, can specify different directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  results_dir <- args[1]
}

cat("\n")
cat("================================================================\n")
cat("Adaptive CV Validation - Results Analysis\n")
cat("================================================================\n\n")

cat("Results directory:", results_dir, "\n\n")

# Find all result files
result_files <- list.files(results_dir, pattern = "adaptive_cv_n.*\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in: ", results_dir)
}

cat("Found", length(result_files), "batch result files\n\n")

# Load and combine all results
all_results <- do.call(rbind, lapply(result_files, function(f) {
  batch <- readRDS(f)
  # Extract n from filename
  n <- as.integer(sub(".*_n([0-9]+)_.*", "\\1", basename(f)))
  batch$n <- n
  batch
}))

cat("Total replications loaded:", nrow(all_results), "\n\n")

# Summarize by sample size
summary_by_n <- all_results %>%
  group_by(n) %>%
  summarise(
    n_reps = n(),
    mean_theta = mean(theta),
    se_theta = sd(theta),
    mean_bias = mean(bias),
    rmse = sqrt(mean(bias^2)),
    coverage = mean(covers) * 100,
    mean_se = mean(sigma),
    se_ratio = mean(sigma) / sd(theta),
    mean_time = mean(time_sec),
    .groups = 'drop'
  )

cat("================================================================\n")
cat("RESULTS BY SAMPLE SIZE\n")
cat("================================================================\n\n")

for (i in 1:nrow(summary_by_n)) {
  row <- summary_by_n[i, ]

  cat("n =", row$n, "\n")
  cat("  Replications:", row$n_reps, "\n\n")

  cat("  Point estimation:\n")
  cat("    Mean estimate:", round(row$mean_theta, 4), "\n")
  cat("    True ATT:", true_att, "\n")
  cat("    Bias:", round(row$mean_bias, 4), "\n")
  cat("    RMSE:", round(row$rmse, 4), "\n")
  cat("    Std dev:", round(row$se_theta, 4), "\n\n")

  cat("  Inference:\n")
  cat("    Coverage:", round(row$coverage, 1), "%\n")
  cat("    Mean SE:", round(row$mean_se, 4), "\n")
  cat("    SE/SD ratio:", round(row$se_ratio, 3), "\n\n")

  cat("  Computation:\n")
  cat("    Mean time:", round(row$mean_time, 1), "sec\n\n")

  cat("  ----------------------------------------------------------------\n\n")
}

cat("================================================================\n")
cat("COMPARISON TO BASELINE (Standard CV)\n")
cat("================================================================\n\n")

# Baseline values from original results
baseline <- data.frame(
  n = c(1000, 2000),
  bias = c(-0.020, -0.020),
  coverage = c(85, 88)
)

comparison <- summary_by_n %>%
  select(n, mean_bias, coverage) %>%
  left_join(baseline, by = "n", suffix = c("_adaptive", "_baseline"))

cat("Sample Size | Metric   | Baseline | Adaptive | Change\n")
cat("-----------:|----------|----------|----------|--------\n")

for (i in 1:nrow(comparison)) {
  row <- comparison[i, ]

  bias_change <- row$mean_bias_adaptive - row$bias_baseline
  cov_change <- row$coverage_adaptive - row$coverage_baseline

  cat(sprintf("n = %4d   | Bias     | %7.4f  | %7.4f  | %+.4f\n",
              row$n, row$bias_baseline, row$mean_bias_adaptive, bias_change))
  cat(sprintf("           | Coverage | %7.1f%% | %7.1f%% | %+.1f%%\n",
              row$coverage_baseline, row$coverage_adaptive, cov_change))
  cat("\n")
}

cat("================================================================\n")
cat("ASSESSMENT\n")
cat("================================================================\n\n")

# Check if improvements achieved
for (i in 1:nrow(comparison)) {
  row <- comparison[i, ]

  bias_improved <- abs(row$mean_bias_adaptive) < abs(row$bias_baseline)
  coverage_improved <- row$coverage_adaptive > row$coverage_baseline
  coverage_target <- row$coverage_adaptive >= 94 && row$coverage_adaptive <= 96

  cat("n =", row$n, ":\n")

  if (bias_improved) {
    cat("  ✓ Bias improved: |", round(row$mean_bias_adaptive, 4), "| < |",
        round(row$bias_baseline, 4), "|\n")
  } else {
    cat("  ✗ Bias not improved: |", round(row$mean_bias_adaptive, 4), "| >= |",
        round(row$bias_baseline, 4), "|\n")
  }

  if (coverage_target) {
    cat("  ✓ Coverage in target range: ", round(row$coverage_adaptive, 1), "% (94-96%)\n")
  } else if (coverage_improved) {
    cat("  ⚠ Coverage improved but not in target: ", round(row$coverage_adaptive, 1),
        "% (target 94-96%)\n")
  } else {
    cat("  ✗ Coverage not improved: ", round(row$coverage_adaptive, 1),
        "% vs ", round(row$coverage_baseline, 1), "%\n")
  }

  cat("\n")
}

cat("================================================================\n\n")

# Save combined results
combined_file <- file.path(results_dir, "adaptive_cv_combined_results.rds")
saveRDS(list(
  all_results = all_results,
  summary = summary_by_n,
  comparison = comparison
), combined_file)

cat("Combined results saved to:", combined_file, "\n\n")
