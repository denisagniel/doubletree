# run_simulations.R
# Main script for running simulations for DML causal estimation with interpretable trees
#
# Uses dmltree::dml_att() with treefarmr for nuisance trees.
# Supports multiple DGPs, configurable replications, and manuscript-style metrics
# (bias, MSE, empirical SE, 95% CI coverage, mean CI width).

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. Install with: install.packages('devtools')")
}
devtools::load_all()

if (!requireNamespace("treefarmr", quietly = TRUE)) {
  stop("treefarmr is required. Install from source or set up .Rprofile to load from ../treefarmr")
}

results_dir <- "simulations/results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# DGPs: binary X, A, Y --------------------------------------------------------
# DGP 1: X ~ Bernoulli(0.5), P(A=1|X) = plogis(0.5*X1 - 0.2), Y(a) ~ Bernoulli(0.3 + 0.2*X1 + tau*a)
generate_data_dgp1 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.5 * X1 - 0.2)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X1))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X1 + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp1")
}

# DGP 2: stronger confounding (X1 and X2 both in propensity and outcome)
generate_data_dgp2 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.4 * X1 + 0.3 * X2 - 0.35)
  A <- as.integer(runif(n) < e)
  m0 <- 0.25 + 0.3 * X1 + 0.2 * X2
  m1 <- m0 + tau
  Y0 <- as.integer(runif(n) < m0)
  Y1 <- as.integer(runif(n) < pmin(1, m1))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp2")
}

# Backward compatibility
generate_data <- generate_data_dgp1

# Compute manuscript-style metrics from a list of dml_att results
compute_metrics <- function(results, true_att) {
  theta_hats <- vapply(results, function(r) r$theta, numeric(1))
  ci_low <- vapply(results, function(r) r$ci_95[1L], numeric(1))
  ci_high <- vapply(results, function(r) r$ci_95[2L], numeric(1))
  n_sim <- length(results)
  bias <- mean(theta_hats) - true_att
  mse <- mean((theta_hats - true_att)^2)
  empirical_se <- sd(theta_hats)
  coverage_95 <- mean(ci_low <= true_att & true_att <= ci_high)
  mean_ci_width <- mean(ci_high - ci_low)
  data.frame(
    n_sim = n_sim,
    true_att = true_att,
    mean_theta_hat = mean(theta_hats),
    bias = bias,
    mse = mse,
    empirical_se = empirical_se,
    coverage_95 = coverage_95,
    mean_ci_width = mean_ci_width
  )
}

# Run one replication and multiple replications ---------------------------------
n_obs <- 400
K <- 5
n_sim <- 50L   # replications per DGP for manuscript-style summary (set to 2–5 for quick runs)
true_tau <- 0.15

set.seed(12345)
data <- generate_data_dgp1(n_obs, tau = true_tau, seed = 12345)
result <- dmltree::dml_att(data$X, data$A, data$Y, K = K)

message("Single run (DGP 1):")
message("  Point estimate (theta): ", round(result$theta, 4))
message("  95% CI: ", paste(round(result$ci_95, 4), collapse = ", "))
message("  True ATT: ", data$tau)

saveRDS(result, file.path(results_dir, "dml_att_result.rds"))

# Multiple replications, DGP 1
message("\nRunning ", n_sim, " replications (DGP 1)...")
results_dgp1 <- list()
for (i in seq_len(n_sim)) {
  d <- generate_data_dgp1(n_obs, tau = true_tau, seed = 10000L + i)
  results_dgp1[[i]] <- dmltree::dml_att(d$X, d$A, d$Y, K = K)
}
metrics_dgp1 <- compute_metrics(results_dgp1, true_tau)
message("DGP 1 metrics: bias = ", round(metrics_dgp1$bias, 4), ", MSE = ", round(metrics_dgp1$mse, 5), ", coverage_95 = ", round(metrics_dgp1$coverage_95, 3))

# Optional: DGP 2 replications
message("\nRunning ", n_sim, " replications (DGP 2)...")
results_dgp2 <- list()
for (i in seq_len(n_sim)) {
  d <- generate_data_dgp2(n_obs, tau = true_tau, seed = 20000L + i)
  results_dgp2[[i]] <- dmltree::dml_att(d$X, d$A, d$Y, K = K)
}
metrics_dgp2 <- compute_metrics(results_dgp2, true_tau)
message("DGP 2 metrics: bias = ", round(metrics_dgp2$bias, 4), ", MSE = ", round(metrics_dgp2$mse, 5), ", coverage_95 = ", round(metrics_dgp2$coverage_95, 3))

# Save replication results and summary table for manuscript
saveRDS(list(theta_hats = vapply(results_dgp1, function(r) r$theta, numeric(1)), results = results_dgp1),
        file.path(results_dir, "simulation_results_dgp1.rds"))
saveRDS(list(theta_hats = vapply(results_dgp2, function(r) r$theta, numeric(1)), results = results_dgp2),
        file.path(results_dir, "simulation_results_dgp2.rds"))

summary_table <- rbind(
  cbind(dgp = "dgp1", metrics_dgp1),
  cbind(dgp = "dgp2", metrics_dgp2)
)
saveRDS(summary_table, file.path(results_dir, "simulation_summary.rds"))
utils::write.csv(summary_table, file.path(results_dir, "simulation_summary.csv"), row.names = FALSE)

message("\nResults and summary table saved to ", results_dir)
