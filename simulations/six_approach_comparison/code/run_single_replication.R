#!/usr/bin/env Rscript

# Main simulation worker for six-approach comparison
# Handles all 6 approaches with common interface
# Called by SLURM job arrays

# Parse command line arguments
library(optparse)

option_list <- list(
  make_option("--approach", type="integer", help="Approach number (1-6)"),
  make_option("--dgp", type="integer", help="DGP number (1-4)"),
  make_option("--n", type="integer", help="Sample size"),
  make_option("--reps", type="integer", default=500, help="Number of replications"),
  make_option("--rep_start", type="integer", default=1, help="Starting replication"),
  make_option("--rep_end", type="integer", default=NULL, help="Ending replication"),
  make_option("--output", type="character", help="Output file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$approach) || is.null(opt$dgp) || is.null(opt$n) || is.null(opt$output)) {
  stop("Must specify --approach, --dgp, --n, and --output")
}

# Set rep range
if (is.null(opt$rep_end)) {
  opt$rep_end <- opt$rep_start + opt$reps - 1
}

# Load required packages
suppressPackageStartupMessages({
  library(optimaltrees)
  library(doubletree)
})

# Source helper functions
code_dir <- "code"
source(file.path(code_dir, "dgps.R"))
source(file.path(code_dir, "estimators.R"))
source(file.path(code_dir, "metrics.R"))

# Map approach number to function
approach_map <- list(
  `1` = estimate_att_fullsample,
  `2` = estimate_att_crossfit,
  `3` = estimate_att_doubletree,
  `4` = estimate_att_doubletree_singlefit,
  `5` = estimate_att_msplit,
  `6` = estimate_att_msplit_singlefit
)

approach_names <- c(
  "full_sample", "crossfit_separate", "doubletree",
  "doubletree_singlefit", "msplit", "msplit_singlefit"
)

# Map DGP number to function (now 4 DGPs)
dgp_map <- list(
  `1` = generate_dgp_simple,
  `2` = generate_dgp_moderate,
  `3` = generate_dgp_complex,
  `4` = generate_dgp_continuous
)

dgp_names <- c("simple", "moderate", "complex", "continuous")

# Get functions
estimator <- approach_map[[as.character(opt$approach)]]
dgp_fun <- dgp_map[[as.character(opt$dgp)]]
approach_name <- approach_names[opt$approach]
dgp_name <- dgp_names[opt$dgp]

# Print job info
cat(sprintf("==========================================\n"))
cat(sprintf("Six-Approach Comparison - Simulation Job\n"))
cat(sprintf("==========================================\n"))
cat(sprintf("Approach: %d (%s)\n", opt$approach, approach_name))
cat(sprintf("DGP: %d (%s)\n", opt$dgp, dgp_name))
cat(sprintf("Sample size: %d\n", opt$n))
cat(sprintf("Replications: %d to %d (%d total)\n",
            opt$rep_start, opt$rep_end, opt$rep_end - opt$rep_start + 1))
cat(sprintf("Output: %s\n", opt$output))
cat(sprintf("==========================================\n\n"))

# Run replications
results <- list()
n_reps <- opt$rep_end - opt$rep_start + 1

for (i in 1:n_reps) {
  rep <- opt$rep_start + i - 1

  # Progress reporting
  if (i %% 25 == 0 || i == n_reps) {
    cat(sprintf("[%s] Rep %d / %d (%.1f%% complete)\n",
                Sys.time(), i, n_reps, 100 * i / n_reps))
  }

  # Set seed for reproducibility
  # Seed formula: base + approach*1e6 + dgp*1e5 + n + rep
  seed <- 1000000 + opt$approach * 1000000 + opt$dgp * 100000 + opt$n + rep
  set.seed(seed)

  # Generate data
  data <- tryCatch({
    dgp_fun(n = opt$n)
  }, error = function(e) {
    list(error = paste("DGP error:", as.character(e)))
  })

  if (!is.null(data$error)) {
    results[[i]] <- list(
      rep = rep,
      approach = opt$approach,
      approach_name = approach_name,
      dgp = opt$dgp,
      dgp_name = dgp_name,
      n = opt$n,
      theta_hat = NA_real_,
      se = NA_real_,
      theta_true = NA_real_,
      elapsed_time = NA_real_,
      error = data$error
    )
    next
  }

  # Time the estimation
  start_time <- Sys.time()

  # Run estimator
  result <- tryCatch({
    estimator(X = data$X, A = data$A, Y = data$Y)
  }, error = function(e) {
    list(theta = NA_real_, se = NA_real_, error = as.character(e))
  })

  elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")

  # Compute coverage metrics
  if (!is.na(result$theta) && !is.na(result$se)) {
    coverage_metrics <- compute_coverage(
      theta_hat = result$theta,
      se = result$se,
      theta_true = data$true_att
    )
  } else {
    coverage_metrics <- c(
      coverage = NA, coverage_lower = NA, coverage_upper = NA,
      ci_width = NA, ci_lower = NA, ci_upper = NA
    )
  }

  # Store results
  results[[i]] <- list(
    rep = rep,
    seed = seed,
    approach = opt$approach,
    approach_name = approach_name,
    dgp = opt$dgp,
    dgp_name = dgp_name,
    n = opt$n,
    theta_hat = result$theta,
    se = result$se,
    theta_true = data$true_att,
    bias = result$theta - data$true_att,
    coverage = coverage_metrics["coverage"],
    coverage_lower = coverage_metrics["coverage_lower"],
    coverage_upper = coverage_metrics["coverage_upper"],
    ci_width = coverage_metrics["ci_width"],
    ci_lower = coverage_metrics["ci_lower"],
    ci_upper = coverage_metrics["ci_upper"],
    elapsed_time = elapsed_time,
    error = if(is.null(result$error)) NA_character_ else result$error
  )
}

# Convert to data frame
cat(sprintf("\nCombining %d results...\n", length(results)))
results_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

# Add summary stats
cat(sprintf("\nResults summary:\n"))
cat(sprintf("  Successful: %d / %d\n", sum(!is.na(results_df$theta_hat)), nrow(results_df)))
cat(sprintf("  Mean bias: %.4f\n", mean(results_df$bias, na.rm = TRUE)))
cat(sprintf("  Mean RMSE: %.4f\n", sqrt(mean(results_df$bias^2, na.rm = TRUE))))
cat(sprintf("  Coverage: %.3f\n", mean(results_df$coverage, na.rm = TRUE)))
cat(sprintf("  Mean time: %.2f sec\n", mean(results_df$elapsed_time, na.rm = TRUE)))

# Save results
cat(sprintf("\nSaving results to: %s\n", opt$output))
saveRDS(results_df, file = opt$output)

cat("\n✓ Job complete!\n")
