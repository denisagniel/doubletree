#!/usr/bin/env Rscript
# Adaptive CV Validation - Worker Script for Cluster
# Runs batch of replications for one configuration

library(doubletree)
library(optimaltrees)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript adaptive_cv_worker.R <n> <batch_start> <batch_size> <output_dir> <seed_offset>")
}

n <- as.integer(args[1])
batch_start <- as.integer(args[2])
batch_size <- as.integer(args[3])
output_dir <- args[4]
seed_offset <- as.integer(args[5])

# Load DGP
source('code/dgps.R')

# Configuration
dgp_name <- "complex"
true_att <- 0.15
batch_end <- batch_start + batch_size - 1

cat("\n")
cat("================================================================\n")
cat("Adaptive CV Validation - Worker\n")
cat("================================================================\n\n")
cat("Configuration:\n")
cat("  DGP:", dgp_name, "\n")
cat("  Sample size:", n, "\n")
cat("  True ATT:", true_att, "\n")
cat("  Batch:", batch_start, "-", batch_end, "(", batch_size, "reps )\n")
cat("  Seed offset:", seed_offset, "\n")
cat("  Output:", output_dir, "\n\n")

# Storage
results <- data.frame(
  rep = integer(),
  theta = numeric(),
  sigma = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  covers = logical(),
  bias = numeric(),
  time_sec = numeric()
)

# Run batch
cat("Running replications...\n")

for (rep in batch_start:batch_end) {
  if ((rep - batch_start + 1) %% 10 == 0) {
    cat("  Rep", rep, "/", batch_end, "\n")
  }

  set.seed(seed_offset + rep)

  # Generate data
  data <- generate_dgp_complex(n = n)

  # Time the estimation
  start_time <- Sys.time()

  result <- tryCatch({
    estimate_att(
      X = data$X,
      A = data$A,
      Y = data$Y,
      K = 5,
      outcome_type = "binary",
      cv_regularization = TRUE,  # Uses adaptive CV
      cv_K = 5,
      verbose = FALSE,
      seed = seed_offset + rep
    )
  }, error = function(e) {
    message("Rep ", rep, " failed: ", e$message)
    return(NULL)
  })

  end_time <- Sys.time()
  time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (!is.null(result)) {
    covers <- result$ci_95[1] <= true_att && true_att <= result$ci_95[2]
    bias <- result$theta - true_att

    results <- rbind(results, data.frame(
      rep = rep,
      theta = result$theta,
      sigma = result$sigma,
      ci_lower = result$ci_95[1],
      ci_upper = result$ci_95[2],
      covers = covers,
      bias = bias,
      time_sec = time_sec
    ))
  }
}

cat("\n")
cat("Batch complete: ", nrow(results), "/", batch_size, "successful\n\n")

# Save results
output_file <- file.path(output_dir, sprintf("adaptive_cv_n%d_batch%d-%d.rds", n, batch_start, batch_end))
saveRDS(results, output_file)
cat("Results saved to:", output_file, "\n\n")

# Return success if at least 90% completed
success_rate <- nrow(results) / batch_size
if (success_rate < 0.9) {
  stop("Batch failed: only ", round(100 * success_rate, 1), "% successful (< 90%)")
}

cat("================================================================\n\n")
