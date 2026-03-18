# Verification Run: 10 reps per config to test all DGPs/methods
# Before committing to full 60-hour run

library(optimaltrees)
library(dplyr)

cat("\n=== VERIFICATION RUN: All Configs ===\n\n")

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

N_REPS <- 10  # Quick verification (instead of 500)
K_FOLDS <- 5
SEED_OFFSET <- 90000

# Build full grid
configs <- expand.grid(
  dgp = c("binary", "continuous", "moderate"),
  n = c(400, 800, 1600),
  method = c("tree", "rashomon", "forest", "linear"),
  stringsAsFactors = FALSE
)

cat("Configuration:\n")
cat("  DGPs:", paste(unique(configs$dgp), collapse=", "), "\n")
cat("  Sample sizes:", paste(unique(configs$n), collapse=", "), "\n")
cat("  Methods:", paste(unique(configs$method), collapse=", "), "\n")
cat("  Replications per config:", N_REPS, "\n")
cat("  Total simulations:", nrow(configs) * N_REPS, "\n\n")

start_total <- Sys.time()
all_results <- list()
failed_configs <- list()

for (i in 1:nrow(configs)) {
  config <- configs[i,]
  
  cat(sprintf("[%d/%d] DGP=%s, n=%d, method=%s...", 
              i, nrow(configs), config$dgp, config$n, config$method))
  
  start_config <- Sys.time()
  config_results <- list()
  n_converged <- 0
  n_failed <- 0
  
  for (rep in 1:N_REPS) {
    seed <- SEED_OFFSET + i*1000 + rep
    
    # Generate data
    if (config$dgp == "binary") {
      d <- generate_dgp_binary_att(n = config$n, tau = 0.10, seed = seed)
    } else if (config$dgp == "continuous") {
      d <- generate_dgp_continuous_att(n = config$n, tau = 0.10, seed = seed)
    } else {  # moderate
      d <- generate_dgp_moderate_att(n = config$n, tau = 0.10, seed = seed)
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
  }
  
  time_config <- as.numeric(difftime(Sys.time(), start_config, units = "secs"))
  
  cat(sprintf(" %.1fs, %d/%d converged", time_config, n_converged, N_REPS))
  
  if (n_failed > 0) {
    cat(sprintf(" [%d FAILED]", n_failed))
    failed_configs[[length(failed_configs) + 1]] <- list(
      config = config,
      n_failed = n_failed,
      errors = sapply(config_results, function(x) if(!x$converged) x$error else NA)
    )
  }
  cat("\n")
  
  all_results[[i]] <- do.call(rbind, config_results)
  
  # Periodic garbage collection
  if (i %% 5 == 0) gc(verbose = FALSE)
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
results_df <- do.call(rbind, all_results)

cat(sprintf("\n=== COMPLETED in %.1f minutes (%.1f seconds) ===\n\n", 
            total_time/60, total_time))

# Summary statistics
converged <- results_df[results_df$converged, ]
n_total <- nrow(results_df)
n_conv <- nrow(converged)

cat(sprintf("Overall Convergence: %d/%d (%.1f%%)\n\n", 
            n_conv, n_total, 100*n_conv/n_total))

# Convergence by method
conv_by_method <- results_df %>%
  group_by(method) %>%
  summarize(
    n = n(),
    converged = sum(converged),
    pct = 100 * mean(converged)
  )

cat("Convergence by Method:\n")
print(as.data.frame(conv_by_method))
cat("\n")

# Convergence by DGP
conv_by_dgp <- results_df %>%
  group_by(dgp) %>%
  summarize(
    n = n(),
    converged = sum(converged),
    pct = 100 * mean(converged)
  )

cat("Convergence by DGP:\n")
print(as.data.frame(conv_by_dgp))
cat("\n")

# Basic statistics (only converged runs)
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
  
  cat("Performance Summary (converged runs):\n")
  print(as.data.frame(stats), row.names = FALSE)
  cat("\n")
}

# Failed configurations
if (length(failed_configs) > 0) {
  cat("=== FAILED CONFIGURATIONS ===\n\n")
  for (i in 1:length(failed_configs)) {
    fc <- failed_configs[[i]]
    cat(sprintf("%d. DGP=%s, n=%d, method=%s: %d failures\n",
                i, fc$config$dgp, fc$config$n, fc$config$method, fc$n_failed))
    errors <- fc$errors[!is.na(fc$errors)]
    if (length(errors) > 0) {
      cat(sprintf("   Error: %s\n", substr(errors[1], 1, 100)))
    }
  }
  cat("\n")
}

# Extrapolate to full run
avg_time_per_rep <- total_time / n_total
cat("=== FULL SIMULATION ESTIMATES ===\n\n")
cat(sprintf("Average time per rep: %.3f seconds\n", avg_time_per_rep))
cat("Full primary simulations (500 reps per config):\n")
cat(sprintf("  - %d configs × 500 reps = %d runs\n", nrow(configs), nrow(configs) * 500))
cat(sprintf("  - Estimated time: %.1f hours\n", (nrow(configs) * 500 * avg_time_per_rep) / 3600))
cat("\n")

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_dir <- sprintf("results/verification_%s", timestamp)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(results_df, file.path(save_dir, "verification_results.rds"))
write.csv(results_df, file.path(save_dir, "verification_results.csv"), row.names = FALSE)

cat(sprintf("Results saved to: %s/\n\n", save_dir))

# Final verdict
cat("=== VERDICT ===\n\n")

if (n_conv == n_total) {
  cat("✅ PASS: All configurations completed successfully!\n")
  cat("   - 100% convergence across all DGPs and methods\n")
  cat("   - Ready for full production run\n")
} else if (n_conv >= 0.95 * n_total) {
  cat("✅ MOSTLY PASS: High convergence rate\n")
  cat(sprintf("   - %.1f%% convergence\n", 100*n_conv/n_total))
  cat("   - A few failures acceptable for production\n")
} else if (n_conv >= 0.8 * n_total) {
  cat("⚠️  CAUTION: Some failures detected\n")
  cat(sprintf("   - %.1f%% convergence\n", 100*n_conv/n_total))
  cat("   - Review failed configurations before full run\n")
} else {
  cat("❌ FAIL: Too many failures\n")
  cat(sprintf("   - Only %.1f%% convergence\n", 100*n_conv/n_total))
  cat("   - Do not proceed with full run\n")
  cat("   - Debug failed configurations first\n")
}

cat("\n")
