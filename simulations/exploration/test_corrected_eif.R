# test_corrected_eif.R
# Clean test with corrected EIF to verify coverage

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}

# Reload package to get corrected EIF
devtools::load_all()

library(cli)

# Load DGP3
source("simulations/run_simulations_extended.R", local = TRUE)

# Configuration ----------------------------------------------------------------

n_val <- 3200  # Large n where theory should work
lambda <- 0.01  # Optimal lambda from diagnostic
n_reps <- 200  # Good balance of precision and speed

cli_h1("Test with Corrected EIF")
cli_text("Goal: Verify that corrected EIF fixes coverage")
cli_text("Sample size: n = {n_val}")
cli_text("Lambda: {lambda} (optimal)")
cli_text("Replications: {n_reps}")
cli_text("")

# Run --------------------------------------------------------------------------

results <- data.frame(
  rep = 1:n_reps,
  theta = NA_real_,
  sigma = NA_real_,
  true_tau = NA_real_,
  bias = NA_real_,
  ci_lower = NA_real_,
  ci_upper = NA_real_,
  coverage = NA_integer_,
  mean_score = NA_real_,
  sd_score = NA_real_
)

cli_progress_bar("Simulations", total = n_reps)

for (i in 1:n_reps) {
  seed <- 7000 + i  # Different seeds from previous tests

  # Generate data
  data <- generate_data_dgp3(n = n_val, seed = seed)

  # Fit with optimal lambda
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
    results$coverage[i] <- as.integer(fit$ci_95[1] <= data$tau && data$tau <= fit$ci_95[2])
    results$mean_score[i] <- mean(fit$score_values)
    results$sd_score[i] <- sd(fit$score_values)
  }, error = function(e) {
    cli_alert_warning("Failed (rep {i}): {e$message}")
  })

  cli_progress_update()
}

cli_progress_done()

# Save
saveRDS(results, "simulations/corrected_eif_test.rds")
write.csv(results, "simulations/corrected_eif_test.csv", row.names = FALSE)

# Analyze ----------------------------------------------------------------------

cli_h1("Results with Corrected EIF")

# Basic stats
cov <- mean(results$coverage, na.rm = TRUE)
cov_se <- sqrt(cov * (1 - cov) / sum(!is.na(results$coverage)))
cov_ci_lower <- cov - 1.96 * cov_se
cov_ci_upper <- cov + 1.96 * cov_se

cli_h2("Coverage")
status <- if (cov >= 0.93 && cov <= 0.97) "✓" else "⚠"
cli_text("{status} Coverage: {round(cov, 3)} ({round(cov * 100, 1)}%)")
cli_text("   95% CI: [{round(cov_ci_lower, 3)}, {round(cov_ci_upper, 3)}]")

if (cov >= 0.93 && cov <= 0.97) {
  cli_alert_success("✓ Corrected EIF achieves nominal coverage!")
} else if (cov >= 0.90) {
  cli_alert_info("Corrected EIF improves coverage (90-93%)")
  cli_text("  → Close to nominal, may need slight adjustment")
} else {
  cli_alert_warning("⚠ Coverage still below target")
  cli_text("  → May need further investigation")
}

# Variance check
empirical_sd <- sd(results$bias, na.rm = TRUE)
estimated_se <- mean(results$sigma, na.rm = TRUE) / sqrt(n_val)
ratio <- empirical_sd / estimated_se

cli_h2("Variance Estimation")
cli_text("Empirical SD of bias: {round(empirical_sd, 4)}")
cli_text("Estimated SE: {round(estimated_se, 4)}")
cli_text("Ratio (empirical / estimated): {round(ratio, 2)}")

if (ratio >= 0.95 && ratio <= 1.05) {
  cli_alert_success("✓ Variance estimation is accurate!")
} else if (ratio >= 0.90 && ratio <= 1.10) {
  cli_alert_info("Variance estimation is close (ratio={round(ratio, 2)})")
} else if (ratio > 1.10) {
  cli_alert_warning("⚠ Variance still underestimated (ratio={round(ratio, 2)})")
} else {
  cli_alert_warning("⚠ Variance overestimated (ratio={round(ratio, 2)})")
}

# Score centering
mean_score_overall <- mean(results$mean_score, na.rm = TRUE)
cli_h2("Score Function")
cli_text("Mean of score values: {format(mean_score_overall, scientific=TRUE)}")
status <- if (abs(mean_score_overall) < 1e-10) "✓" else "⚠"
cli_text("{status} Score is properly centered")

# Comparison to previous
cli_h2("Comparison to Previous Results")
cli_text("Previous optimal lambda test (λ=0.1, n=3200, 50 reps): 94%")
cli_text("Previous variance diagnostic (λ=0.1, n=3200, 100 reps): 82%")
cli_text("Previous focused test (λ=0.1, n=3200, 100 reps): 75-82%")
cli_text("Current test (λ=0.01, n=3200, {n_reps} reps): {round(cov * 100, 1)}%")
cli_text("")

diff_optimal <- cov - 0.94
diff_variance_diag <- cov - 0.82

if (abs(diff_optimal) < 0.05) {
  cli_text("→ Similar to optimal lambda test (may have been statistical noise)")
} else if (cov > 0.94) {
  cli_alert_success("→ IMPROVEMENT over previous tests!")
} else {
  cli_text("→ Different from previous tests (may reflect corrected EIF)")
}

cli_h2("Bias and RMSE")
bias_mean <- mean(results$bias, na.rm = TRUE)
rmse <- sqrt(mean(results$bias^2, na.rm = TRUE))
cli_text("Mean bias: {round(bias_mean, 5)}")
cli_text("RMSE: {round(rmse, 4)}")

cli_h2("Saved")
cli_text("- simulations/corrected_eif_test.rds")
cli_text("- simulations/corrected_eif_test.csv")
