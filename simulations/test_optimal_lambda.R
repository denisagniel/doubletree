# test_optimal_lambda.R
# Quick test: Does optimal lambda (0.01) improve coverage vs default (0.1)?

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}
devtools::load_all()

library(cli)

# Load DGP3
source("simulations/run_simulations_extended.R", local = TRUE)

# Configuration ----------------------------------------------------------------

n_val <- 3200  # Largest n where we saw worst coverage (75-82%)
n_reps <- 50   # Quick test
lambdas <- c(0.01, 0.1)  # Optimal vs default

cli_h1("Optimal Lambda Test")
cli_text("Goal: Test if λ=0.01 improves coverage vs λ=0.1")
cli_text("Sample size: n = {n_val}")
cli_text("Replications: {n_reps}")
cli_text("")

# Run --------------------------------------------------------------------------

results <- expand.grid(
  lambda = lambdas,
  rep = 1:n_reps,
  stringsAsFactors = FALSE
)
results$theta <- NA_real_
results$sigma <- NA_real_
results$true_tau <- NA_real_
results$ci_lower <- NA_real_
results$ci_upper <- NA_real_
results$coverage <- NA_integer_

cli_progress_bar("Test runs", total = nrow(results))

for (i in 1:nrow(results)) {
  lambda <- results$lambda[i]
  rep <- results$rep[i]
  seed <- 5000 + i  # Different seeds from other tests

  # Generate data
  data <- generate_data_dgp3(n = n_val, seed = seed)

  # Fit with specified lambda
  tryCatch({
    fit <- dml_att(
      X = data$X,
      A = data$A,
      Y = data$Y,
      use_rashomon = FALSE,
      K = 5,
      regularization = lambda,
      seed = seed
    )

    results$theta[i] <- fit$theta
    results$sigma[i] <- fit$sigma
    results$true_tau[i] <- data$tau
    results$ci_lower[i] <- fit$ci_95[1]
    results$ci_upper[i] <- fit$ci_95[2]
    results$coverage[i] <- as.integer(fit$ci_95[1] <= data$tau && data$tau <= fit$ci_95[2])
  }, error = function(e) {
    cli_alert_warning("Failed (rep {rep}, λ={lambda}): {e$message}")
  })

  cli_progress_update()
}

cli_progress_done()

# Analyze ----------------------------------------------------------------------

cli_h1("Results")

# Aggregate by lambda
agg <- aggregate(cbind(theta, sigma, coverage) ~ lambda, results, mean, na.rm = TRUE)

cli_h2("Coverage by Lambda (n={n_val})")
for (i in 1:nrow(agg)) {
  lambda_val <- agg$lambda[i]
  cov <- agg$coverage[i]
  status <- if (cov >= 0.93 && cov <= 0.97) "✓" else "⚠"
  cli_text("{status} λ={lambda_val}: coverage = {round(cov, 3)} ({round(cov * 100, 1)}%)")
}

# Statistical test
lambda_01 <- results$coverage[results$lambda == 0.01]
lambda_10 <- results$coverage[results$lambda == 0.1]
diff_mean <- mean(lambda_01, na.rm = TRUE) - mean(lambda_10, na.rm = TRUE)
n_01 <- sum(!is.na(lambda_01))
n_10 <- sum(!is.na(lambda_10))

cli_h2("Comparison")
cli_text("Coverage improvement (λ=0.01 vs λ=0.1): {round(diff_mean, 3)} ({round(diff_mean * 100, 1)} percentage points)")

if (diff_mean > 0.05) {
  cli_alert_success("✓ Optimal lambda IMPROVES coverage substantially")
  cli_text("  → Tree complexity WAS the main problem")
  cli_text("  → Recommendation: Re-run focused test with λ=0.01")
} else if (diff_mean > 0.02) {
  cli_alert_info("Optimal lambda helps modestly ({round(diff_mean * 100, 1)}% improvement)")
  cli_text("  → Tree complexity is part of the problem")
  cli_text("  → May still need to investigate variance estimation")
} else {
  cli_alert_warning("⚠ Optimal lambda does NOT fix coverage problem")
  cli_text("  → Coverage improvement: only {round(diff_mean * 100, 1)}%")
  cli_text("  → Problem is primarily about variance estimation, not tree complexity")
  cli_text("  → Need to investigate:")
  cli_text("    - Score function centering")
  cli_text("    - Cross-fitting with K=5 and binary features")
  cli_text("    - Variance formula validity for 4-leaf trees")
}

# Save
saveRDS(results, "simulations/optimal_lambda_test.rds")
write.csv(results, "simulations/optimal_lambda_test.csv", row.names = FALSE)

cli_h2("Saved")
cli_text("- simulations/optimal_lambda_test.rds")
cli_text("- simulations/optimal_lambda_test.csv")
