# Estimate full simulation runtime
library(optimaltrees)

source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

cat("\n=== SIMULATION SPEED ESTIMATION ===\n\n")

# Test one full configuration (like in production simulations)
N_TEST <- 30  # 30 reps instead of 500
configs <- expand.grid(
  dgp = c("binary", "continuous", "moderate"),
  n = c(400, 800, 1600),
  stringsAsFactors = FALSE
)

cat("Testing ", nrow(configs), " configurations with ", N_TEST, " reps each\n")
cat("Total: ", nrow(configs) * N_TEST, " simulations\n\n")

start_total <- Sys.time()
all_results <- list()

for (i in 1:nrow(configs)) {
  config <- configs[i,]
  cat(sprintf("Config %d/%d: DGP=%s, n=%d...", 
              i, nrow(configs), config$dgp, config$n))
  
  start_config <- Sys.time()
  
  # Run N_TEST replications
  for (rep in 1:N_TEST) {
    if (config$dgp == "binary") {
      d <- generate_dgp_binary_att(n = config$n, tau = 0.10, seed = 80000 + i*1000 + rep)
    } else if (config$dgp == "continuous") {
      d <- generate_dgp_continuous_att(n = config$n, tau = 0.10, seed = 80000 + i*1000 + rep)
    } else {
      d <- generate_dgp_moderate_att(n = config$n, tau = 0.10, seed = 80000 + i*1000 + rep)
    }
    
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = 5,
      outcome_type = "binary",
      regularization = log(config$n) / config$n,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      worker_limit = 4,  # Use parallelization
      verbose = FALSE
    )
  }
  
  time_config <- as.numeric(difftime(Sys.time(), start_config, units = "secs"))
  cat(sprintf(" %.2f sec (%.3f sec/rep)\n", time_config, time_config/N_TEST))
  
  all_results[[i]] <- data.frame(
    dgp = config$dgp,
    n = config$n,
    time_total = time_config,
    time_per_rep = time_config/N_TEST
  )
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
results_df <- do.call(rbind, all_results)

cat(sprintf("\nCompleted in %.1f seconds (%.2f min)\n\n", total_time, total_time/60))

cat("=== TIMING BREAKDOWN ===\n\n")
print(results_df)

# Extrapolate to full simulation
cat("\n=== FULL SIMULATION ESTIMATES ===\n\n")
cat("Primary simulations (Table 1):\n")
cat("  - 3 DGPs × 3 sample sizes × 4 methods × 500 reps = 18,000 runs\n")

avg_time_per_rep <- mean(results_df$time_per_rep)
cat(sprintf("  - Average time per rep: %.3f seconds\n", avg_time_per_rep))

# 4 methods means 4x the runs
estimated_primary <- avg_time_per_rep * 18000
cat(sprintf("  - Estimated time (worker_limit=4): %.1f hours\n", estimated_primary/3600))

cat("\nStress simulations (Table 2):\n")
cat("  - 3 DGPs × 2 sample sizes × 2 methods × 200 reps = 2,400 runs\n")
estimated_stress <- avg_time_per_rep * 2400
cat(sprintf("  - Estimated time (worker_limit=4): %.1f hours\n", estimated_stress/3600))

cat(sprintf("\nTotal estimated time: %.1f hours\n", (estimated_primary + estimated_stress)/3600))

cat("\n=== VERDICT ===\n\n")
if (avg_time_per_rep < 0.2) {
  cat("EXCELLENT: Simulations are fast!\n")
  cat(sprintf("  - Average %.3f sec/rep\n", avg_time_per_rep))
  cat(sprintf("  - Full primary simulations: ~%.1f hours\n", estimated_primary/3600))
  cat(sprintf("  - Full stress tests: ~%.1f hours\n", estimated_stress/3600))
  cat("\n  Simulations are ready for production runs.\n")
} else if (avg_time_per_rep < 0.5) {
  cat("GOOD: Reasonable simulation speed\n")
  cat(sprintf("  - Full simulations will take ~%.1f hours\n", (estimated_primary + estimated_stress)/3600))
} else {
  cat("SLOW: May want to optimize before full runs\n")
}

cat("\n")
