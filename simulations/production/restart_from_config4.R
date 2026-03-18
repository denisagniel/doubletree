# Restart Production Run from Config 4
# Configs 1-3 already complete, continue from Config 4/24

library(optimaltrees)
library(dplyr)

cat("\n=== RESTARTING PRODUCTION RUN FROM CONFIG 4 ===\n\n")
flush.console()

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

N_REPS <- 500
K_FOLDS <- 5
SEED_OFFSET <- 100000

# Build full grid
configs <- expand.grid(
  dgp = c("binary", "continuous"),
  n = c(400, 800, 1600),
  method = c("tree", "rashomon", "forest", "linear"),
  stringsAsFactors = FALSE
)

# Use same results directory
results_dir <- "results/primary_20260313_133335"

cat("Resuming from Config 4/24\n")
cat("Configs 1-3 already completed\n")
cat("Remaining:", nrow(configs) - 3, "configs\n\n")
flush.console()

# Load existing results
all_results <- list()
for (i in 1:3) {
  file <- sprintf("%s/config_%02d_%s_n%d_%s.rds",
                  results_dir, i, configs$dgp[i], configs$n[i], configs$method[i])
  all_results[[i]] <- readRDS(file)
  cat(sprintf("Loaded: Config %d ✓\n", i))
}
flush.console()

start_total <- Sys.time()
config_times <- numeric(nrow(configs))

# Resume from Config 4
for (i in 4:nrow(configs)) {
  config <- configs[i,]

  cat(sprintf("\n=== Config %d/%d: %s, n=%d, %s ===\n",
              i, nrow(configs), config$dgp, config$n, config$method))
  flush.console()

  start_config <- Sys.time()
  config_results <- list()
  n_converged <- 0
  n_failed <- 0

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
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                      outcome_type = "binary", regularization = log(config$n)/config$n,
                      cv_regularization = FALSE, use_rashomon = FALSE,
                      worker_limit = 4, verbose = FALSE)
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)
      } else if (config$method == "rashomon") {
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                      outcome_type = "binary", regularization = log(config$n)/config$n,
                      cv_regularization = FALSE, use_rashomon = TRUE,
                      worker_limit = 4, verbose = FALSE)
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)
      } else if (config$method == "forest") {
        fit <- att_forest(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                            num.trees = 500, seed = seed)
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)
      } else {  # linear
        fit <- att_linear(X = d$X, A = d$A, Y = d$Y, K = K_FOLDS,
                            interactions = FALSE, seed = seed)
        list(theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
             converged = TRUE, error = NA)
      }
    }, error = function(e) {
      cat(sprintf("\n  ERROR in rep %d: %s\n", rep, conditionMessage(e)))
      flush.console()
      list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
           converged = FALSE, error = conditionMessage(e))
    })

    if (result$converged) {
      n_converged <- n_converged + 1
    } else {
      n_failed <- n_failed + 1
    }

    config_results[[rep]] <- data.frame(
      dgp = config$dgp, n = config$n, method = config$method, rep = rep,
      true_att = d$true_att, theta = result$theta, sigma = result$sigma,
      ci_lower = result$ci_lower, ci_upper = result$ci_upper,
      converged = result$converged, error = as.character(result$error),
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
              n_converged, N_REPS, 100*n_converged/N_REPS, time_config, time_config/N_REPS))
  if (n_failed > 0) cat(sprintf("  WARNING: %d failures\n", n_failed))
  flush.console()

  # Save config results
  config_df <- do.call(rbind, config_results)
  all_results[[i]] <- config_df
  config_file <- sprintf("%s/config_%02d_%s_n%d_%s.rds",
                        results_dir, i, config$dgp, config$n, config$method)
  saveRDS(config_df, config_file)

  # Estimate remaining
  completed_configs <- sum(!is.na(config_times))
  if (completed_configs > 3) {
    avg_time <- mean(config_times[4:i], na.rm = TRUE)
    remaining <- (nrow(configs) - i) * avg_time
    cat(sprintf("  Estimated remaining: %.1f min (%.1f hours)\n",
                remaining/60, remaining/3600))
    flush.console()
  }

  if (i %% 5 == 0) gc(verbose = FALSE)
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
results_df <- do.call(rbind, all_results)

cat(sprintf("\n=== COMPLETED in %.1f minutes (%.2f hours) ===\n\n",
            total_time/60, total_time/3600))

# Summary
n_total <- nrow(results_df)
n_conv <- sum(results_df$converged)
cat(sprintf("Total: %d runs, %d converged (%.1f%%)\n\n", n_total, n_conv, 100*n_conv/n_total))

# Save complete results
saveRDS(results_df, file.path(results_dir, "primary_simulations_complete.rds"))
write.csv(results_df, file.path(results_dir, "primary_simulations_complete.csv"), row.names = FALSE)

cat(sprintf("Results saved to: %s/\n\n", results_dir))
cat("=== PRODUCTION RUN COMPLETE ===\n")
