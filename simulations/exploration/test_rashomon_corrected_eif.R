# test_rashomon_corrected_eif.R
# Test Rashomon-DML with corrected EIF

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}

# Reload package to get corrected EIF
devtools::load_all()

library(cli)

# Load DGP3 and oracle
source("simulations/run_simulations_extended.R", local = TRUE)

# Configuration ----------------------------------------------------------------

n_val <- 3200  # Large n where theory works
lambda <- 0.01  # Optimal lambda
n_reps <- 100  # Good balance for all three methods

# Theory-guided epsilon: ε_n = c * sqrt(log(n)/n)
epsilon_theory <- function(n, c = 2) {
  c * sqrt(log(n) / n)
}

eps <- epsilon_theory(n_val, c = 2)  # Moderate choice

cli_h1("Test Rashomon-DML with Corrected EIF")
cli_text("Goal: Verify Rashomon implementation works with corrected EIF")
cli_text("Sample size: n = {n_val}")
cli_text("Lambda: {lambda} (optimal)")
cli_text("Epsilon: {round(eps, 4)} (theory-guided, c=2)")
cli_text("Replications: {n_reps}")
cli_text("")

# Run --------------------------------------------------------------------------

results <- expand.grid(
  method = c("fold_specific", "rashomon", "oracle"),
  rep = 1:n_reps,
  stringsAsFactors = FALSE
)
results$theta <- NA_real_
results$sigma <- NA_real_
results$true_tau <- NA_real_
results$bias <- NA_real_
results$ci_lower <- NA_real_
results$ci_upper <- NA_real_
results$coverage <- NA_integer_
results$mean_score <- NA_real_
results$pct_nonempty_e <- NA_real_
results$pct_nonempty_m0 <- NA_real_
results$pct_nonempty_m1 <- NA_real_
results$n_intersecting_e <- NA_integer_

cli_progress_bar("Simulations", total = nrow(results))

for (i in 1:nrow(results)) {
  method <- results$method[i]
  rep <- results$rep[i]
  seed <- 8000 + i

  # Generate data
  data <- generate_data_dgp3(n = n_val, seed = seed)

  if (method == "oracle") {
    # Oracle with true nuisances
    tryCatch({
      fit <- dml_att_oracle(data, K = 5, seed = seed)

      results$theta[i] <- fit$theta
      results$sigma[i] <- fit$sigma
      results$true_tau[i] <- data$tau
      results$bias[i] <- fit$theta - data$tau
      results$ci_lower[i] <- fit$ci_95[1]
      results$ci_upper[i] <- fit$ci_95[2]
      results$coverage[i] <- as.integer(fit$ci_95[1] <= data$tau && data$tau <= fit$ci_95[2])
      results$mean_score[i] <- mean(fit$score_values)
    }, error = function(e) {
      cli_alert_warning("Oracle failed (rep {rep}): {e$message}")
    })

  } else if (method == "fold_specific") {
    # Fold-specific trees
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
    }, error = function(e) {
      cli_alert_warning("Fold-specific failed (rep {rep}): {e$message}")
    })

  } else if (method == "rashomon") {
    # Rashomon-DML
    tryCatch({
      fit <- estimate_att(
        X = data$X,
        A = data$A,
        Y = data$Y,
        use_rashomon = TRUE,
        rashomon_bound_multiplier = eps,
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

      # Rashomon-specific metadata
      if (!is.null(fit$pct_nonempty_e)) results$pct_nonempty_e[i] <- fit$pct_nonempty_e
      if (!is.null(fit$pct_nonempty_m0)) results$pct_nonempty_m0[i] <- fit$pct_nonempty_m0
      if (!is.null(fit$pct_nonempty_m1)) results$pct_nonempty_m1[i] <- fit$pct_nonempty_m1
      if (!is.null(fit$n_intersecting_e)) results$n_intersecting_e[i] <- fit$n_intersecting_e
    }, error = function(e) {
      cli_alert_warning("Rashomon failed (rep {rep}): {e$message}")
    })
  }

  cli_progress_update()
}

cli_progress_done()

# Save
saveRDS(results, "simulations/rashomon_corrected_eif_test.rds")
write.csv(results, "simulations/rashomon_corrected_eif_test.csv", row.names = FALSE)

# Analyze ----------------------------------------------------------------------

cli_h1("Results with Corrected EIF")

# Aggregate by method
agg <- aggregate(cbind(theta, sigma, bias, coverage, mean_score,
                       pct_nonempty_e, pct_nonempty_m0, pct_nonempty_m1) ~ method,
                 results, mean, na.rm = TRUE)

cli_h2("Coverage by Method")
for (i in 1:nrow(agg)) {
  method_name <- agg$method[i]
  cov <- agg$coverage[i]
  n_obs <- sum(results$method == method_name & !is.na(results$coverage))
  cov_se <- sqrt(cov * (1 - cov) / n_obs)
  cov_ci <- c(cov - 1.96 * cov_se, cov + 1.96 * cov_se)

  status <- if (cov >= 0.93 && cov <= 0.97) "✓" else "⚠"
  cli_text("{status} {method_name}: {round(cov, 3)} ({round(cov * 100, 1)}%) [{round(cov_ci[1], 3)}, {round(cov_ci[2], 3)}]")
}

# Compare methods
oracle_cov <- agg$coverage[agg$method == "oracle"]
fold_cov <- agg$coverage[agg$method == "fold_specific"]
rashomon_cov <- agg$coverage[agg$method == "rashomon"]

cli_h2("Method Comparison")
cli_text("Oracle: {round(oracle_cov * 100, 1)}%")
cli_text("Fold-specific: {round(fold_cov * 100, 1)}%")
cli_text("Rashomon: {round(rashomon_cov * 100, 1)}%")

if (oracle_cov >= 0.93 && oracle_cov <= 0.97) {
  cli_alert_success("✓ Oracle achieves nominal coverage (validates DML framework)")
} else {
  cli_alert_warning("⚠ Oracle coverage off target: {round(oracle_cov * 100, 1)}%")
}

if (fold_cov >= 0.93 && fold_cov <= 0.97) {
  cli_alert_success("✓ Fold-specific achieves nominal coverage")
} else if (fold_cov >= 0.90) {
  cli_alert_info("Fold-specific close to nominal: {round(fold_cov * 100, 1)}%")
} else {
  cli_alert_warning("⚠ Fold-specific coverage low: {round(fold_cov * 100, 1)}%")
}

if (rashomon_cov >= 0.93 && rashomon_cov <= 0.97) {
  cli_alert_success("✓ Rashomon achieves nominal coverage!")
  cli_text("  → Theory validated!")
} else if (rashomon_cov >= 0.90) {
  cli_alert_info("Rashomon close to nominal: {round(rashomon_cov * 100, 1)}%")
  cli_text("  → Good, but slightly conservative")
} else {
  cli_alert_warning("⚠ Rashomon coverage low: {round(rashomon_cov * 100, 1)}%")
  cli_text("  → May need further investigation")
}

# Rashomon intersection
cli_h2("Rashomon Intersection Success")
rashomon_data <- results[results$method == "rashomon", ]
if (any(!is.na(rashomon_data$pct_nonempty_e))) {
  pct_e <- mean(rashomon_data$pct_nonempty_e, na.rm = TRUE)
  pct_m0 <- mean(rashomon_data$pct_nonempty_m0, na.rm = TRUE)
  pct_m1 <- mean(rashomon_data$pct_nonempty_m1, na.rm = TRUE)
  pct_avg <- mean(c(pct_e, pct_m0, pct_m1), na.rm = TRUE)

  cli_text("Non-empty intersection rate:")
  cli_text("  e(X): {round(pct_e * 100, 1)}%")
  cli_text("  m0(X): {round(pct_m0 * 100, 1)}%")
  cli_text("  m1(X): {round(pct_m1 * 100, 1)}%")
  cli_text("  Average: {round(pct_avg * 100, 1)}%")

  if (pct_avg >= 0.85) {
    cli_alert_success("✓ High intersection success rate")
  } else if (pct_avg >= 0.70) {
    cli_alert_info("Moderate intersection success rate")
  } else {
    cli_alert_warning("⚠ Low intersection success rate")
  }
} else {
  cli_alert_warning("⚠ Intersection metadata still missing (all NA)")
  cli_text("  → Check Rashomon implementation")
}

# Variance estimation
cli_h2("Variance Estimation by Method")
for (method_name in c("fold_specific", "rashomon", "oracle")) {
  subset_data <- results[results$method == method_name, ]
  empirical_sd <- sd(subset_data$bias, na.rm = TRUE)
  estimated_se <- mean(subset_data$sigma, na.rm = TRUE) / sqrt(n_val)
  ratio <- empirical_sd / estimated_se

  status <- if (ratio >= 0.95 && ratio <= 1.05) "✓" else if (ratio >= 0.90 && ratio <= 1.10) "~" else "⚠"
  cli_text("{status} {method_name}: ratio={round(ratio, 2)} (empirical SD={round(empirical_sd, 4)}, estimated SE={round(estimated_se, 4)})")
}

# Bias comparison
cli_h2("Bias and RMSE by Method")
for (i in 1:nrow(agg)) {
  method_name <- agg$method[i]
  bias_mean <- agg$bias[i]
  subset_data <- results[results$method == method_name, ]
  rmse <- sqrt(mean(subset_data$bias^2, na.rm = TRUE))

  cli_text("{method_name}: bias={round(bias_mean, 5)}, RMSE={round(rmse, 4)}")
}

cli_h2("Previous Results Comparison")
cli_text("Previous focused test (λ=0.1, wrong EIF):")
cli_text("  Oracle: 94.3% ✓")
cli_text("  Fold-specific: 82.7%")
cli_text("  Rashomon: 82.0%")
cli_text("")
cli_text("Current test (λ=0.01, corrected EIF):")
cli_text("  Oracle: {round(oracle_cov * 100, 1)}%")
cli_text("  Fold-specific: {round(fold_cov * 100, 1)}%")
cli_text("  Rashomon: {round(rashomon_cov * 100, 1)}%")

if (rashomon_cov >= 0.90 && fold_cov >= 0.90) {
  cli_alert_success("✓ IMPROVEMENT: Both methods now achieve good coverage!")
  cli_text("  → Corrected EIF fixed the problem")
} else {
  cli_text("  → Further investigation needed")
}

cli_h2("Saved")
cli_text("- simulations/rashomon_corrected_eif_test.rds")
cli_text("- simulations/rashomon_corrected_eif_test.csv")
