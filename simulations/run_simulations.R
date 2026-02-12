# run_simulations.R
# Main script for running simulations for DML causal estimation with interpretable trees
#
# Uses dmltree::dml_att() with treefarmr for nuisance trees.

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

# Simple DGP: binary X, A, Y ----------------------------------------------------
# X ~ Bernoulli(0.5), P(A=1|X) = plogis(0.5*X - 0.2), Y(a) ~ Bernoulli(0.3 + 0.2*X + tau*a), tau = ATT
generate_data <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.5 * X1 - 0.2)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X1))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X1 + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau)
}

# Run one replication and optionally multiple ----------------------------------
n_obs <- 400
K <- 5

set.seed(12345)
data <- generate_data(n_obs, tau = 0.15, seed = 12345)
result <- dmltree::dml_att(data$X, data$A, data$Y, K = K)

message("Point estimate (theta): ", round(result$theta, 4))
message("95% CI: ", paste(round(result$ci_95, 4), collapse = ", "))
message("True ATT: ", data$tau)

# Save single run
saveRDS(result, file.path(results_dir, "dml_att_result.rds"))

# Optional: multiple replications
n_sim <- 10
results <- list()
for (i in seq_len(n_sim)) {
  d <- generate_data(n_obs, tau = 0.15, seed = 10000 + i)
  results[[i]] <- dmltree::dml_att(d$X, d$A, d$Y, K = K)
}
theta_hats <- vapply(results, function(r) r$theta, numeric(1))
message("\nOver ", n_sim, " replications: mean(theta_hat) = ", round(mean(theta_hats), 4), ", true ATT = 0.15")
saveRDS(list(theta_hats = theta_hats, results = results), file.path(results_dir, "simulation_results.rds"))

message("Results saved to ", results_dir)
