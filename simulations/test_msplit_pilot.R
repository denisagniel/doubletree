#!/usr/bin/env Rscript
# Pilot M-Split Simulations
# Test M-split implementation with simple DGPs

library(optimaltrees)
library(doubletree)
library(cli)

cli_h1("M-Split Pilot Simulations")

# Simple DGP: Binary outcome with treatment effect (BINARY COVARIATES)
generate_simple_att_binary <- function(n, att = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Covariates
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = rbinom(n, 1, 0.5)
  )

  # Treatment (confounded by x1)
  ps <- plogis(-0.5 + 1.0 * X$x1 + 0.5 * X$x2)
  A <- rbinom(n, 1, ps)

  # Outcome (treatment effect + covariate effect)
  y0 <- plogis(-1.0 + 0.8 * X$x1 + 0.6 * X$x2)
  y1 <- plogis(-1.0 + 0.8 * X$x1 + 0.6 * X$x2 + att)
  Y <- ifelse(A == 1, rbinom(n, 1, y1), rbinom(n, 1, y0))

  # True ATT (average over treated)
  true_att <- mean(y1[A == 1] - y0[A == 1])

  list(X = X, A = A, Y = Y, true_att = true_att, n = n, n_treated = sum(A))
}

# DGP with continuous covariates
generate_simple_att <- function(n, att = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Continuous covariates (8 features)
  X <- data.frame(
    x1 = rnorm(n, 0, 1),
    x2 = rnorm(n, 0, 1),
    x3 = rnorm(n, 0, 1),
    x4 = rnorm(n, 0, 1),
    x5 = rnorm(n, 0, 1),
    x6 = runif(n, 0, 1),
    x7 = runif(n, 0, 1),
    x8 = runif(n, 0, 1)
  )

  # Treatment (confounded by x1, x2, x3)
  ps <- plogis(-0.5 + 0.6 * X$x1 + 0.4 * X$x2 - 0.3 * X$x3)
  A <- rbinom(n, 1, ps)

  # Outcome (nonlinear effects)
  # y0 depends on x1, x2, x4, x6 with some interactions
  mu0 <- plogis(-1.0 +
                0.5 * X$x1 +
                0.4 * X$x2 +
                0.3 * X$x4 +
                0.5 * X$x6 +
                0.2 * X$x1 * X$x2)  # interaction

  y0 <- mu0
  y1 <- plogis(qlogis(mu0) + att)  # additive on logit scale

  Y <- ifelse(A == 1, rbinom(n, 1, y1), rbinom(n, 1, y0))

  # True ATT (average over treated)
  true_att <- mean(y1[A == 1] - y0[A == 1])

  list(X = X, A = A, Y = Y, true_att = true_att, n = n, n_treated = sum(A))
}

# Test 1: Single dataset, varying M
cli_h2("Test 1: Effect of M on variance")
cli_text("Fixed dataset (n=300), M = 1, 5, 10, 20")
cli_text("Using BINARY covariates (3 features) for consistency check")

data <- generate_simple_att_binary(n = 300, att = 0.2, seed = 123)
cli_alert_info("True ATT: {round(data$true_att, 3)}")
cli_alert_info("n = {data$n}, n_treated = {data$n_treated} ({round(100*data$n_treated/data$n, 1)}%)")

M_values <- c(1, 5, 10, 20)
results_test1 <- vector("list", length(M_values))

for (i in seq_along(M_values)) {
  M <- M_values[i]
  cli_alert("Running M = {M}...")

  result <- estimate_att_msplit(
    data$X, data$A, data$Y,
    M = M,
    K = 3,
    seed_base = 100,
    verbose = FALSE
  )

  results_test1[[i]] <- list(
    M = M,
    theta = result$theta,
    sigma = result$sigma,
    ci = result$ci_95,
    coverage = (result$ci_95[1] <= data$true_att && data$true_att <= result$ci_95[2]),
    bias = result$theta - data$true_att,
    struct_freq_e = result$diagnostics$structure_frequency_e,
    struct_freq_m0 = result$diagnostics$structure_frequency_m0,
    pred_var_e = result$diagnostics$mean_prediction_variance_e,
    pred_var_m0 = result$diagnostics$mean_prediction_variance_m0
  )

  cli_alert_success("M={M}: θ̂={round(result$theta, 3)}, SE={round(result$sigma, 3)}, bias={round(result$theta - data$true_att, 3)}")
}

# Summary table
cli_h3("Summary Table")
cat("\n")
cat("M  | θ̂     | SE    | Bias  | Coverage | Struct Freq (e/m0) | Pred Var (e/m0)\n")
cat("---|-------|-------|-------|----------|--------------------|-----------------\n")
for (res in results_test1) {
  cat(sprintf(
    "%-2d | %5.3f | %.3f | %5.3f | %-8s | %.2f / %.2f       | %.4f / %.4f\n",
    res$M, res$theta, res$sigma, res$bias,
    ifelse(res$coverage, "Yes", "No"),
    res$struct_freq_e, res$struct_freq_m0,
    res$pred_var_e, res$pred_var_m0
  ))
}
cat("\n")

# Check variance reduction
var_e_m1 <- results_test1[[1]]$pred_var_e
if (!is.na(var_e_m1) && var_e_m1 > 1e-10) {
  cli_alert_info("Prediction variance reduction (e):")
  for (i in 2:length(results_test1)) {
    ratio <- results_test1[[i]]$pred_var_e / var_e_m1
    theoretical <- 1 / results_test1[[i]]$M
    cli_text("  M={results_test1[[i]]$M}: {round(ratio, 3)} (theoretical: {round(theoretical, 3)})")
  }
} else {
  if (is.na(var_e_m1)) {
    cli_alert_warning("Prediction variance is NA for M=1 (single split, no variance to compute)")
  } else {
    cli_alert_warning("Prediction variance is zero (perfect consistency across splits)")
  }
}

# Test 2: Replication (check coverage by sample size)
cli_h2("Test 2: Coverage by sample size (M=10 vs M=1)")
cli_text("n = 200, 500, 1000 | M=10 and M=1, K=3, 50 replications each")
cli_text("Using CONTINUOUS covariates (8 features) with interactions")

n_values <- c(200, 500, 1000)
n_reps <- 50
results_by_n_m10 <- vector("list", length(n_values))
results_by_n_m1 <- vector("list", length(n_values))

for (n_idx in seq_along(n_values)) {
  n <- n_values[n_idx]
  cli_alert("Testing n = {n}...")

  results_m10 <- vector("list", n_reps)
  results_m1 <- vector("list", n_reps)

  for (rep in 1:n_reps) {
    data <- generate_simple_att(n = n, att = 0.2, seed = 1000 + n_idx * 1000 + rep)

    # M=10
    result_m10 <- estimate_att_msplit(
      data$X, data$A, data$Y,
      M = 10,
      K = 3,
      seed_base = 100,
      verbose = FALSE
    )

    results_m10[[rep]] <- list(
      theta = result_m10$theta,
      sigma = result_m10$sigma,
      bias = result_m10$theta - data$true_att,
      coverage = (result_m10$ci_95[1] <= data$true_att && data$true_att <= result_m10$ci_95[2]),
      ci_width = result_m10$ci_95[2] - result_m10$ci_95[1]
    )

    # M=1 (single-split)
    result_m1 <- estimate_att_msplit(
      data$X, data$A, data$Y,
      M = 1,
      K = 3,
      seed_base = 100,
      verbose = FALSE
    )

    results_m1[[rep]] <- list(
      theta = result_m1$theta,
      sigma = result_m1$sigma,
      bias = result_m1$theta - data$true_att,
      coverage = (result_m1$ci_95[1] <= data$true_att && data$true_att <= result_m1$ci_95[2]),
      ci_width = result_m1$ci_95[2] - result_m1$ci_95[1]
    )

    if (rep %% 10 == 0) cli_alert("  Completed {rep}/{n_reps}")
  }

  results_by_n_m10[[n_idx]] <- results_m10
  results_by_n_m1[[n_idx]] <- results_m1

  # Summary for this n
  thetas_m10 <- sapply(results_m10, function(r) r$theta)
  biases_m10 <- sapply(results_m10, function(r) r$bias)
  coverages_m10 <- sapply(results_m10, function(r) r$coverage)
  ci_widths_m10 <- sapply(results_m10, function(r) r$ci_width)

  thetas_m1 <- sapply(results_m1, function(r) r$theta)
  biases_m1 <- sapply(results_m1, function(r) r$bias)
  coverages_m1 <- sapply(results_m1, function(r) r$coverage)
  ci_widths_m1 <- sapply(results_m1, function(r) r$ci_width)

  cli_h3("n = {n} Summary")
  cli_alert_info("M=10: Coverage {round(100 * mean(coverages_m10), 1)}%, Mean θ̂={round(mean(thetas_m10), 3)}, RMSE={round(sqrt(mean(biases_m10^2)), 3)}")
  cli_alert_info("M=1:  Coverage {round(100 * mean(coverages_m1), 1)}%, Mean θ̂={round(mean(thetas_m1), 3)}, RMSE={round(sqrt(mean(biases_m1^2)), 3)}")
  cat("\n")
}

# Overall summary table
cli_h3("Coverage Summary Table: M=10")
cat("\n")
cat("n    | Mean θ̂ | SD(θ̂) | Bias  | RMSE  | Coverage | Mean CI Width\n")
cat("-----|--------|--------|-------|-------|----------|---------------\n")
for (n_idx in seq_along(n_values)) {
  n <- n_values[n_idx]
  results <- results_by_n_m10[[n_idx]]
  thetas <- sapply(results, function(r) r$theta)
  biases <- sapply(results, function(r) r$bias)
  coverages <- sapply(results, function(r) r$coverage)
  ci_widths <- sapply(results, function(r) r$ci_width)

  cat(sprintf(
    "%-4d | %6.3f | %6.3f | %5.3f | %5.3f | %6.1f%%  | %6.3f\n",
    n, mean(thetas), sd(thetas), mean(biases), sqrt(mean(biases^2)),
    100 * mean(coverages), mean(ci_widths)
  ))
}
cat("\n")

cli_h3("Coverage Summary Table: M=1 (Single-Split)")
cat("\n")
cat("n    | Mean θ̂ | SD(θ̂) | Bias  | RMSE  | Coverage | Mean CI Width\n")
cat("-----|--------|--------|-------|-------|----------|---------------\n")
for (n_idx in seq_along(n_values)) {
  n <- n_values[n_idx]
  results <- results_by_n_m1[[n_idx]]
  thetas <- sapply(results, function(r) r$theta)
  biases <- sapply(results, function(r) r$bias)
  coverages <- sapply(results, function(r) r$coverage)
  ci_widths <- sapply(results, function(r) r$ci_width)

  cat(sprintf(
    "%-4d | %6.3f | %6.3f | %5.3f | %5.3f | %6.1f%%  | %6.3f\n",
    n, mean(thetas), sd(thetas), mean(biases), sqrt(mean(biases^2)),
    100 * mean(coverages), mean(ci_widths)
  ))
}
cat("\n")

# Store for later
results_test2 <- list(m10 = results_by_n_m10, m1 = results_by_n_m1)

# Test 3: Continuous outcome
cli_h2("Test 3: Continuous outcome")

generate_continuous_att <- function(n, att = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5)
  )

  ps <- plogis(-0.3 + 0.8 * X$x1)
  A <- rbinom(n, 1, ps)

  # Continuous outcome
  y0 <- 5 + 2 * X$x1 + 1 * X$x2 + rnorm(n, 0, 1)
  y1 <- y0 + att
  Y <- ifelse(A == 1, y1, y0)

  true_att <- att  # Constant treatment effect

  list(X = X, A = A, Y = Y, true_att = true_att, n = n)
}

data_cont <- generate_continuous_att(n = 200, att = 1.5, seed = 456)
cli_alert_info("True ATT (continuous): {data_cont$true_att}")

result_cont <- estimate_att_msplit(
  data_cont$X, data_cont$A, data_cont$Y,
  M = 10,
  K = 3,
  outcome_type = "continuous",
  seed_base = 200,
  verbose = FALSE
)

cli_alert_success("Continuous outcome: θ̂={round(result_cont$theta, 3)}, SE={round(result_cont$sigma, 3)}")
cli_alert_info("Bias: {round(result_cont$theta - data_cont$true_att, 3)}")
cli_alert_info("95% CI: [{round(result_cont$ci_95[1], 3)}, {round(result_cont$ci_95[2], 3)}]")
cli_alert_info("Coverage: {ifelse(result_cont$ci_95[1] <= data_cont$true_att && data_cont$true_att <= result_cont$ci_95[2], 'Yes', 'No')}")

cli_h1("Pilot Simulations Complete")
cli_alert_success("All tests completed successfully!")

# Save results
results_all <- list(
  test1_varying_M = results_test1,
  test2_coverage_by_n = list(
    n_values = n_values,
    n_reps = n_reps,
    results_m10 = results_by_n_m10,
    results_m1 = results_by_n_m1
  ),
  test3_continuous = result_cont
)

saveRDS(results_all, "simulations/results/msplit_pilot_results.rds")
cli_alert_info("Results saved to: simulations/results/msplit_pilot_results.rds")

# Quick comparison: does coverage improve with n?
coverages_by_n_m10 <- sapply(results_by_n_m10, function(res) {
  mean(sapply(res, function(r) r$coverage))
})

coverages_by_n_m1 <- sapply(results_by_n_m1, function(res) {
  mean(sapply(res, function(r) r$coverage))
})

cli_h3("Key Finding: Coverage vs. Sample Size (M=10 vs M=1)")
cat("\n")
cat("n    | M=10 Coverage | M=1 Coverage | Difference\n")
cat("-----|---------------|--------------|------------\n")
for (i in seq_along(n_values)) {
  diff <- coverages_by_n_m10[i] - coverages_by_n_m1[i]
  cat(sprintf("%-4d | %6.1f%%      | %5.1f%%      | %+.1f%%\n",
              n_values[i], 100*coverages_by_n_m10[i], 100*coverages_by_n_m1[i], 100*diff))
}
cat("\n")

# Diagnosis
m10_declining <- coverages_by_n_m10[length(coverages_by_n_m10)] < coverages_by_n_m10[1]
m1_declining <- coverages_by_n_m1[length(coverages_by_n_m1)] < coverages_by_n_m1[1]

if (m10_declining && !m1_declining) {
  cli_alert_warning("M=10 coverage declines with n, but M=1 does not → Issue specific to M-split")
} else if (m10_declining && m1_declining) {
  cli_alert_warning("Both M=10 and M=1 coverage decline with n → Issue in base estimator")
} else if (!m10_declining && !m1_declining) {
  cli_alert_success("Coverage stable or improving for both M=10 and M=1")
} else {
  cli_alert_info("M=1 declining but M=10 not → Unexpected pattern")
}
