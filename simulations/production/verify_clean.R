# Clean Verification: Small test with real-time progress
# Tests representative sample: 2 DGPs × 1 n × 4 methods × 3 reps = 24 runs

library(optimaltrees)
library(dplyr)

cat("\n=== CLEAN VERIFICATION RUN ===\n\n")
flush.console()

# Source files
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

N_REPS <- 3
K_FOLDS <- 5
N_VAL <- 400  # Single sample size for quick test

# Configs: 2 DGPs × 4 methods = 8 configs
configs <- expand.grid(
  dgp = c("binary", "continuous"),
  method = c("tree", "rashomon", "forest", "linear"),
  stringsAsFactors = FALSE
)

cat("Configuration:\n")
cat("  DGPs: binary, continuous\n")
cat("  Sample size: n=400\n")
cat("  Methods: tree, rashomon, forest, linear\n")
cat("  Reps:", N_REPS, "\n")
cat("  Total:", nrow(configs) * N_REPS, "runs\n\n")
flush.console()

start_total <- Sys.time()
all_results <- list()

for (i in 1:nrow(configs)) {
  config <- configs[i,]

  cat(sprintf("\n[Config %d/%d] %s + %s\n",
              i, nrow(configs), config$dgp, config$method))
  flush.console()

  config_results <- list()
  times <- numeric(N_REPS)

  for (rep in 1:N_REPS) {
    cat(sprintf("  Rep %d/%d... ", rep, N_REPS))
    flush.console()

    seed <- 90000 + i*1000 + rep

    # Generate data
    if (config$dgp == "binary") {
      d <- generate_dgp_binary_att(n = N_VAL, tau = 0.10, seed = seed)
    } else {
      d <- generate_dgp_continuous_att(n = N_VAL, tau = 0.10, seed = seed)
    }

    # Time the fit
    start_rep <- Sys.time()

    result <- tryCatch({
      if (config$method == "tree") {
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                      outcome_type = "binary", regularization = log(N_VAL)/N_VAL,
                      cv_regularization = FALSE, use_rashomon = FALSE,
                      worker_limit = 4, verbose = FALSE)
        list(theta = fit$theta, sigma = fit$sigma, success = TRUE)
      } else if (config$method == "rashomon") {
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                      outcome_type = "binary", regularization = log(N_VAL)/N_VAL,
                      cv_regularization = FALSE, use_rashomon = TRUE,
                      worker_limit = 4, verbose = FALSE)
        list(theta = fit$theta, sigma = fit$sigma, success = TRUE)
      } else if (config$method == "forest") {
        fit <- att_forest(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                            num.trees = 500, seed = seed)
        list(theta = fit$theta, sigma = fit$sigma, success = TRUE)
      } else {  # linear
        fit <- att_linear(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                            interactions = FALSE, seed = seed)
        list(theta = fit$theta, sigma = fit$sigma, success = TRUE)
      }
    }, error = function(e) {
      list(theta = NA, sigma = NA, success = FALSE)
    })

    elapsed <- as.numeric(difftime(Sys.time(), start_rep, units = "secs"))
    times[rep] <- elapsed

    if (result$success) {
      cat(sprintf("%.2fs ✓\n", elapsed))
    } else {
      cat(sprintf("%.2fs ✗ FAILED\n", elapsed))
    }
    flush.console()

    config_results[[rep]] <- data.frame(
      dgp = config$dgp,
      method = config$method,
      rep = rep,
      true_att = d$true_att,
      theta = result$theta,
      sigma = result$sigma,
      time_sec = elapsed,
      success = result$success,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("  Config average: %.2f sec/rep\n", mean(times)))
  flush.console()

  all_results[[i]] <- do.call(rbind, config_results)
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
results_df <- do.call(rbind, all_results)

cat(sprintf("\n=== COMPLETED in %.1f sec (%.1f min) ===\n\n", total_time, total_time/60))

# Summary
n_total <- nrow(results_df)
n_success <- sum(results_df$success)

cat(sprintf("Success rate: %d/%d (%.1f%%)\n\n", n_success, n_total, 100*n_success/n_total))

# Timing by method
timing_summary <- results_df %>%
  filter(success) %>%
  group_by(method) %>%
  summarize(
    n = n(),
    mean_time = mean(time_sec),
    sd_time = sd(time_sec),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_time))

cat("Timing by method:\n")
print(as.data.frame(timing_summary), row.names = FALSE)

# Extrapolation
overall_avg <- mean(results_df$time_sec[results_df$success])
cat(sprintf("\nOverall average: %.2f sec/rep\n", overall_avg))

# Estimate full run with scaling
cat("\n=== FULL SIMULATION ESTIMATES ===\n\n")
cat("Assuming scaling factors:\n")
cat("  n=400: 1.0x (baseline)\n")
cat("  n=800: 1.5x\n")
cat("  n=1600: 2.5x\n")

# 24 configs (2 DGPs × 3 n × 4 methods) × 500 reps
n400_configs <- 8  # 2 DGPs × 4 methods
n800_configs <- 8
n1600_configs <- 8

time_n400 <- n400_configs * 500 * overall_avg
time_n800 <- n800_configs * 500 * overall_avg * 1.5
time_n1600 <- n1600_configs * 500 * overall_avg * 2.5

total_estimated <- (time_n400 + time_n800 + time_n1600) / 3600

cat(sprintf("\nEstimated time: %.1f hours\n", total_estimated))

if (total_estimated < 2) {
  cat("✓ FEASIBLE: < 2 hours for overnight run\n")
} else if (total_estimated < 5) {
  cat("✓ REASONABLE: 2-5 hours\n")
} else {
  cat("⚠️  SLOW: > 5 hours\n")
}

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_file <- sprintf("verify_clean_%s.rds", timestamp)
saveRDS(results_df, save_file)
cat(sprintf("\nResults saved: %s\n", save_file))

cat("\n")
