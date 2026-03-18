# Production Run: Primary Simulations
# Binary + Continuous DGPs × 3 sample sizes × 4 methods × 500 reps = 12,000 runs
# Estimated time: 1-2 hours

library(optimaltrees)
library(dplyr)

cat("\n=== PRIMARY SIMULATIONS: PRODUCTION RUN ===\n\n")
flush.console()

# Source dmltree
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

# Configuration
N_REPS <- 500
K_FOLDS <- 5
SEED_OFFSET <- 100000  # Different from verification

# Build grid: 2 DGPs × 3 n × 4 methods = 24 configs
configs <- expand.grid(
  dgp = c("binary", "continuous"),
  n = c(400, 800, 1600),
  method = c("tree", "rashomon", "forest", "linear"),
  stringsAsFactors = FALSE
)

cat("Configuration:\n")
cat("  DGPs: binary, continuous\n")
cat("  Sample sizes: 400, 800, 1600\n")
cat("  Methods: tree, rashomon, forest, linear\n")
cat("  Replications per config:", N_REPS, "\n")
cat("  Total simulations:", nrow(configs) * N_REPS, "\n")
cat("  K-fold CV:", K_FOLDS, "\n")
cat("  Worker limit: 4 (parallel)\n\n")
flush.console()

# Create results directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- sprintf("results/primary_%s", timestamp)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Results directory: %s/\n\n", results_dir))
flush.console()

start_total <- Sys.time()
all_results <- list()
config_times <- numeric(nrow(configs))

for (i in 1:nrow(configs)) {
  config <- configs[i,]

  cat(sprintf("\n=== Config %d/%d: %s, n=%d, %s ===\n",
              i, nrow(configs), config$dgp, config$n, config$method))
  flush.console()

  start_config <- Sys.time()
  config_results <- list()
  n_converged <- 0
  n_failed <- 0

  # Progress reporting every 50 reps
  for (rep in 1:N_REPS) {
    if (rep %% 50 == 0 || rep == 1) {
      cat(sprintf("  Rep %d/%d... ", rep, N_REPS))
      flush.console()
    }

    seed <- SEED_OFFSET + i*10000 + rep

    # Generate data
    if (config$dgp == "binary") {
      d <- generate_dgp_binary_att(n = config$n, tau = 0.10, seed = seed)
    } else {
      d <- generate_dgp_continuous_att(n = config$n, tau = 0.10, seed = seed)
    }

    # Fit model
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
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)

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
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)

      } else if (config$method == "forest") {
        fit <- att_forest(
          X = d$X, A = d$A, Y = d$Y,
          K = K_FOLDS,
          num.trees = 500,
          seed = seed
        )
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)

      } else {  # linear
        fit <- att_linear(
          X = d$X, A = d$A, Y = d$Y,
          K = K_FOLDS,
          interactions = FALSE,
          seed = seed
        )
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)
      }
    }, error = function(e) {
      list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
           converged = FALSE, error = conditionMessage(e))
    })

    if (result$converged) {
      n_converged <- n_converged + 1
    } else {
      n_failed <- n_failed + 1
    }

    config_results[[rep]] <- data.frame(
      dgp = config$dgp,
      n = config$n,
      method = config$method,
      rep = rep,
      true_att = d$true_att,
      theta = result$theta,
      sigma = result$sigma,
      ci_lower = result$ci_lower,
      ci_upper = result$ci_upper,
      converged = result$converged,
      error = as.character(result$error),
      stringsAsFactors = FALSE
    )

    if (rep %% 50 == 0) {
      cat(sprintf("%.0f%% converged\n", 100*n_converged/rep))
      flush.console()
    }
  }

  time_config <- as.numeric(difftime(Sys.time(), start_config, units = "secs"))
  config_times[i] <- time_config

  cat(sprintf("  Completed: %d/%d converged (%.1f%%), %.1f sec (%.2f sec/rep)\n",
              n_converged, N_REPS, 100*n_converged/N_REPS,
              time_config, time_config/N_REPS))

  if (n_failed > 0) {
    cat(sprintf("  WARNING: %d failures\n", n_failed))
  }
  flush.console()

  # Save config results incrementally
  config_df <- do.call(rbind, config_results)
  all_results[[i]] <- config_df

  config_file <- sprintf("%s/config_%02d_%s_n%d_%s.rds",
                        results_dir, i, config$dgp, config$n, config$method)
  saveRDS(config_df, config_file)

  # Estimate remaining time
  if (i > 1) {
    avg_time_per_config <- mean(config_times[1:i])
    remaining_configs <- nrow(configs) - i
    est_remaining <- remaining_configs * avg_time_per_config
    cat(sprintf("  Estimated remaining: %.1f min (%.1f hours)\n",
                est_remaining/60, est_remaining/3600))
    flush.console()
  }

  # Garbage collection every 5 configs
  if (i %% 5 == 0) {
    gc(verbose = FALSE)
  }
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
results_df <- do.call(rbind, all_results)

cat(sprintf("\n=== COMPLETED in %.1f minutes (%.2f hours) ===\n\n",
            total_time/60, total_time/3600))
flush.console()

# Overall summary
converged <- results_df[results_df$converged, ]
n_total <- nrow(results_df)
n_conv <- nrow(converged)

cat(sprintf("Overall Results:\n"))
cat(sprintf("  Total simulations: %d\n", n_total))
cat(sprintf("  Converged: %d (%.1f%%)\n", n_conv, 100*n_conv/n_total))
cat(sprintf("  Failed: %d (%.1f%%)\n", n_total-n_conv, 100*(n_total-n_conv)/n_total))
cat(sprintf("  Average time: %.2f sec/rep\n\n", total_time/n_total))
flush.console()

# Convergence by method
conv_by_method <- results_df %>%
  group_by(method) %>%
  summarize(
    total = n(),
    converged = sum(converged),
    pct = 100 * mean(converged),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))

cat("Convergence by method:\n")
print(as.data.frame(conv_by_method), row.names = FALSE)
cat("\n")
flush.console()

# Performance by DGP and method
if (n_conv > 0) {
  stats <- converged %>%
    group_by(dgp, method) %>%
    summarize(
      n = n(),
      bias = mean(theta - true_att),
      rmse = sqrt(mean((theta - true_att)^2)),
      coverage = 100 * mean(ci_lower <= true_att & ci_upper >= true_att),
      .groups = "drop"
    )

  cat("Performance summary:\n")
  print(as.data.frame(stats), row.names = FALSE)
  cat("\n")
  flush.console()
}

# Save complete results
saveRDS(results_df, file.path(results_dir, "primary_simulations_complete.rds"))
write.csv(results_df, file.path(results_dir, "primary_simulations_complete.csv"), row.names = FALSE)

# Save metadata
metadata <- list(
  timestamp = timestamp,
  n_configs = nrow(configs),
  n_reps = N_REPS,
  k_folds = K_FOLDS,
  total_runs = n_total,
  converged = n_conv,
  total_time_sec = total_time,
  avg_time_per_rep = total_time/n_total,
  configs = configs
)
saveRDS(metadata, file.path(results_dir, "metadata.rds"))

cat(sprintf("Results saved to: %s/\n", results_dir))
cat(sprintf("  - primary_simulations_complete.rds (main results)\n"))
cat(sprintf("  - primary_simulations_complete.csv (for inspection)\n"))
cat(sprintf("  - config_*.rds (24 individual config files)\n"))
cat(sprintf("  - metadata.rds (run information)\n\n"))
flush.console()

cat("=== PRODUCTION RUN COMPLETE ===\n\n")
