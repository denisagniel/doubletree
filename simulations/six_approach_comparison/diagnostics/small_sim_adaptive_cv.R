#!/usr/bin/env Rscript
# Small simulation to validate adaptive CV impact
# Focus: n=2000, complex DGP (where coverage was worst at 88%)
# Goal: Measure if adaptive CV improves bias and coverage

library(doubletree)
library(optimaltrees)
library(parallel)
source('code/dgps.R')

cat("\n")
cat("================================================================\n")
cat("Adaptive CV Validation Simulation\n")
cat("================================================================\n\n")

# Simulation parameters
n_reps <- 50
n <- 2000
dgp_name <- "complex"
true_att <- 0.15
seed_start <- 1000

cat("Design:\n")
cat("  Replications:", n_reps, "\n")
cat("  Sample size:", n, "\n")
cat("  DGP:", dgp_name, "\n")
cat("  True ATT:", true_att, "\n")
cat("  Estimator: estimate_att (full-sample crossfit)\n")
cat("  CV: Adaptive (max_iterations=10)\n\n")

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

cat("Running simulation...\n\n")

for (rep in 1:n_reps) {
  if (rep %% 10 == 0) cat("  Rep", rep, "/", n_reps, "\n")

  set.seed(seed_start + rep)

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
      cv_regularization = TRUE,  # Uses adaptive CV now
      cv_K = 5,
      verbose = FALSE,
      seed = seed_start + rep
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
cat("================================================================\n")
cat("RESULTS\n")
cat("================================================================\n\n")

# Compute summary statistics
n_success <- nrow(results)
mean_theta <- mean(results$theta)
se_theta <- sd(results$theta)
mean_bias <- mean(results$bias)
rmse <- sqrt(mean(results$bias^2))
coverage <- mean(results$covers) * 100
mean_se <- mean(results$sigma)
mean_time <- mean(results$time_sec)

cat("Success rate:", n_success, "/", n_reps,
    "(", round(100 * n_success / n_reps, 1), "%)\n\n")

cat("Point estimation:\n")
cat("  Mean estimate:", round(mean_theta, 4), "\n")
cat("  True ATT:", true_att, "\n")
cat("  Bias:", round(mean_bias, 4), "\n")
cat("  RMSE:", round(rmse, 4), "\n")
cat("  Std dev:", round(se_theta, 4), "\n\n")

cat("Inference:\n")
cat("  Coverage:", round(coverage, 1), "%\n")
cat("  Mean SE:", round(mean_se, 4), "\n")
cat("  SE/SD ratio:", round(mean_se / se_theta, 3), "\n\n")

cat("Computation:\n")
cat("  Mean time:", round(mean_time, 1), "sec\n")
cat("  Total time:", round(sum(results$time_sec) / 60, 1), "min\n\n")

cat("----------------------------------------------------------------\n")
cat("Comparison to Baseline (from original results)\n")
cat("----------------------------------------------------------------\n\n")

cat("Original results (n=2000, complex DGP, standard CV):\n")
cat("  Bias: -0.020 (underestimation)\n")
cat("  Coverage: 88% (severe undercoverage)\n\n")

cat("Current results (n=2000, complex DGP, adaptive CV):\n")
cat("  Bias:", round(mean_bias, 4), "\n")
cat("  Coverage:", round(coverage, 1), "%\n\n")

# Improvement assessment
bias_improved <- abs(mean_bias) < 0.020
coverage_improved <- coverage > 88

if (bias_improved && coverage_improved) {
  cat("✓ SUCCESS: Adaptive CV improved both bias and coverage!\n")
} else if (bias_improved) {
  cat("⚠ PARTIAL: Bias improved but coverage still low\n")
} else if (coverage_improved) {
  cat("⚠ PARTIAL: Coverage improved but bias still present\n")
} else {
  cat("✗ NO IMPROVEMENT: Problem persists\n")
}

cat("\n")
cat("================================================================\n\n")

# Save results
saveRDS(results, "diagnostics/adaptive_cv_validation_results.rds")
cat("Results saved to: diagnostics/adaptive_cv_validation_results.rds\n\n")
