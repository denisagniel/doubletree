# diagnose_variance_estimation.R
# Check if variance estimation is correct by examining:
# 1. Is mean(score) ≈ 0? (should be by definition of theta_hat)
# 2. Does empirical SE match estimated SE?
# 3. Does coverage improve with lower lambda?

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}
devtools::load_all()

library(cli)

# Load DGP3
source("simulations/run_simulations_extended.R", local = TRUE)

# Configuration ----------------------------------------------------------------

sample_sizes <- c(800, 1600, 3200)
lambda_values <- c(0.01, 0.05, 0.1)  # Test optimal vs default
n_reps <- 100
test_lambdas <- TRUE

cli_h1("Variance Estimation Diagnostic")
cli_text("Goal: Check if variance formula is working correctly")
cli_text("Sample sizes: {paste(sample_sizes, collapse=', ')}")
cli_text("Lambda values: {paste(lambda_values, collapse=', ')}")
cli_text("Replications: {n_reps}")
cli_text("")

# Run diagnostic ---------------------------------------------------------------

results <- expand.grid(
  n = sample_sizes,
  lambda = lambda_values,
  rep = 1:n_reps,
  stringsAsFactors = FALSE
)
results$theta <- NA_real_
results$sigma <- NA_real_
results$true_tau <- NA_real_
results$bias <- NA_real_
results$ci_lower <- NA_real_
results$ci_upper <- NA_real_
results$ci_width <- NA_real_
results$coverage <- NA_integer_
results$mean_score <- NA_real_
results$sd_score <- NA_real_

cli_progress_bar("Diagnostic runs", total = nrow(results))

for (i in 1:nrow(results)) {
  n <- results$n[i]
  lambda <- results$lambda[i]
  rep <- results$rep[i]
  seed <- 1000 * i + rep

  # Generate data
  data <- generate_data_dgp3(n = n, seed = seed)

  # Fit with specified lambda
  tryCatch({
    fit <- estimate_att(
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
    results$bias[i] <- fit$theta - data$tau
    results$ci_lower[i] <- fit$ci_95[1]
    results$ci_upper[i] <- fit$ci_95[2]
    results$ci_width[i] <- fit$ci_95[2] - fit$ci_95[1]
    results$coverage[i] <- as.integer(fit$ci_95[1] <= data$tau && data$tau <= fit$ci_95[2])

    # Check score properties
    results$mean_score[i] <- mean(fit$score_values)
    results$sd_score[i] <- sd(fit$score_values)
  }, error = function(e) {
    # Leave as NA
  })

  cli_progress_update()
}

cli_progress_done()

# Save
saveRDS(results, "simulations/variance_diagnostic.rds")
write.csv(results, "simulations/variance_diagnostic.csv", row.names = FALSE)

# Analyze ----------------------------------------------------------------------

cli_h1("Diagnostic Results")

# Aggregate by n and lambda
agg <- aggregate(cbind(bias, sigma, coverage, ci_width, mean_score, sd_score) ~ n + lambda,
                 data = results, FUN = mean, na.rm = TRUE)

cli_h2("Coverage by Sample Size and Lambda")
for (n_val in unique(agg$n)) {
  cli_text("n = {n_val}")
  subset_agg <- agg[agg$n == n_val, ]
  for (j in 1:nrow(subset_agg)) {
    lambda_val <- subset_agg$lambda[j]
    cov <- subset_agg$coverage[j]
    status <- if (cov >= 0.93 && cov <= 0.97) "✓" else "⚠"
    cli_text("  {status} λ={lambda_val}: coverage={round(cov, 3)}, CI width={round(subset_agg$ci_width[j], 4)}")
  }
}

cli_h2("Score Function Properties")
cli_text("Checking if mean(score) ≈ 0 (should be by definition of theta_hat)")
for (n_val in unique(agg$n)) {
  subset_agg <- agg[agg$n == n_val, ]
  for (j in 1:nrow(subset_agg)) {
    lambda_val <- subset_agg$lambda[j]
    mean_score <- subset_agg$mean_score[j]
    sd_score <- subset_agg$sd_score[j]
    status <- if (abs(mean_score) < 0.001) "✓" else "⚠"
    cli_text("{status} n={n_val}, λ={lambda_val}: mean(score)={round(mean_score, 6)}, sd(score)={round(sd_score, 4)}")
  }
}

cli_h2("Empirical vs Estimated SE")
cli_text("Compare estimated σ to empirical SD of bias")
for (n_val in unique(agg$n)) {
  subset_res <- results[results$n == n_val, ]
  for (lambda_val in unique(agg$lambda)) {
    subset_lambda <- subset_res[subset_res$lambda == lambda_val, ]
    empirical_sd <- sd(subset_lambda$bias, na.rm = TRUE)
    estimated_se <- mean(subset_lambda$sigma, na.rm = TRUE) / sqrt(n_val)
    ratio <- empirical_sd / estimated_se
    status <- if (ratio >= 0.8 && ratio <= 1.2) "✓" else "⚠"
    cli_text("{status} n={n_val}, λ={lambda_val}: empirical SD={round(empirical_sd, 4)}, estimated SE={round(estimated_se, 4)}, ratio={round(ratio, 2)}")
  }
}

cli_h2("Coverage Improvement with Lower Lambda?")
cov_by_lambda <- aggregate(coverage ~ lambda, agg, mean)
cli_text("Average coverage by lambda (across all n):")
for (i in 1:nrow(cov_by_lambda)) {
  cli_text("  λ={cov_by_lambda$lambda[i]}: {round(cov_by_lambda$coverage[i], 3)}")
}

if (cov_by_lambda$coverage[1] > cov_by_lambda$coverage[nrow(cov_by_lambda)]) {
  cli_alert_success("✓ Lower lambda improves coverage")
  cli_text("  Recommendation: Use λ=0.01 or λ=0.05")
} else {
  cli_alert_warning("⚠ Lower lambda does NOT fix coverage problem")
  cli_text("  → Problem is not primarily about tree complexity")
  cli_text("  → Need to investigate variance estimation formula or score function")
}

cli_h2("Saved")
cli_text("- simulations/variance_diagnostic.rds")
cli_text("- simulations/variance_diagnostic.csv")
