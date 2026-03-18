# Quick Timing Test: 1 rep per config to get accurate time estimates
# Goal: Determine realistic runtime for full simulations

library(optimaltrees)
library(dplyr)

cat("\n=== QUICK TIMING TEST ===\n\n")

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source DGPs
source("dgps/dgps_smooth.R")

# Source baseline methods
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

K_FOLDS <- 5
SEED <- 99999

# Test each method on n=400 with both DGPs
configs <- expand.grid(
  dgp = c("binary", "continuous"),
  n = c(400),
  method = c("tree", "rashomon", "forest", "linear"),
  stringsAsFactors = FALSE
)

cat("Testing:", nrow(configs), "configurations (1 rep each)\n")
cat("Sample size: n=400 only\n\n")

timing_results <- data.frame()

for (i in 1:nrow(configs)) {
  config <- configs[i,]

  cat(sprintf("[%d/%d] %s, n=%d, %s... ",
              i, nrow(configs), config$dgp, config$n, config$method))

  # Generate data
  if (config$dgp == "binary") {
    d <- generate_dgp_binary_att(n = config$n, tau = 0.10, seed = SEED)
  } else {
    d <- generate_dgp_continuous_att(n = config$n, tau = 0.10, seed = SEED)
  }

  # Time the fit
  start_time <- Sys.time()

  result <- tryCatch({
    if (config$method == "tree") {
      fit <- estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        K = K_FOLDS,
        outcome_type = "binary",
        regularization = log(config$n) / config$n,
        cv_regularization = FALSE,
        use_rashomon = FALSE,
        worker_limit = 4,
        verbose = FALSE
      )
      list(success = TRUE, theta = fit$theta)

    } else if (config$method == "rashomon") {
      fit <- estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        K = K_FOLDS,
        outcome_type = "binary",
        regularization = log(config$n) / config$n,
        cv_regularization = FALSE,
        use_rashomon = TRUE,
        worker_limit = 4,
        verbose = FALSE
      )
      list(success = TRUE, theta = fit$theta)

    } else if (config$method == "forest") {
      fit <- att_forest(
        X = d$X, A = d$A, Y = d$Y,
        K = K_FOLDS,
        num.trees = 500,
        seed = SEED
      )
      list(success = TRUE, theta = fit$theta)

    } else {  # linear
      fit <- att_linear(
        X = d$X, A = d$A, Y = d$Y,
        K = K_FOLDS,
        interactions = FALSE,
        seed = SEED
      )
      list(success = TRUE, theta = fit$theta)
    }
  }, error = function(e) {
    list(success = FALSE, theta = NA)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  cat(sprintf("%.2fs ", elapsed))
  if (result$success) {
    cat("✓\n")
  } else {
    cat("✗ FAILED\n")
  }

  timing_results <- rbind(timing_results, data.frame(
    dgp = config$dgp,
    n = config$n,
    method = config$method,
    time_sec = elapsed,
    success = result$success,
    stringsAsFactors = FALSE
  ))
}

cat("\n=== TIMING RESULTS ===\n\n")
print(timing_results, row.names = FALSE)

# Calculate extrapolations
cat("\n=== EXTRAPOLATIONS ===\n\n")

# By method
method_times <- timing_results %>%
  group_by(method) %>%
  summarize(avg_time = mean(time_sec), .groups = "drop") %>%
  arrange(desc(avg_time))

cat("Average time per rep by method (n=400):\n")
for (i in 1:nrow(method_times)) {
  cat(sprintf("  %s: %.2f sec/rep\n", method_times$method[i], method_times$avg_time[i]))
}

# Full simulation estimates (binary + continuous only)
# 24 configs (2 DGPs × 3 n × 4 methods) × 500 reps = 12,000 runs
overall_avg <- mean(timing_results$time_sec)
cat(sprintf("\nOverall average: %.2f sec/rep\n", overall_avg))
cat("\nFull simulation estimates (binary + continuous):\n")
cat(sprintf("  24 configs × 500 reps = 12,000 runs\n"))
cat(sprintf("  Estimated time: %.1f hours (%.1f days)\n",
            (12000 * overall_avg) / 3600,
            (12000 * overall_avg) / 3600 / 24))

# Breakdown by n scaling (rough estimate)
cat("\nRough estimates by sample size:\n")
cat("  n=400: 1.0x (baseline)\n")
cat("  n=800: ~1.5-2x slower\n")
cat("  n=1600: ~2-3x slower\n")

cat("\n=== DONE ===\n")
