# test_repeated_debug.R
# Quick debug test to see variance components

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/run_simulations_extended.R", local = TRUE)

# Single run with verbose output
n <- 400
tau <- 0.15
K <- 5
seed <- 123

dgp_fn <- generate_data_dgp1
d <- dgp_fn(n, tau, seed = seed)

message("=== Single Cross-Fit ===")
result_single <- estimate_att(
  d$X, d$A, d$Y, K = K,
  use_rashomon = FALSE,
  regularization = log(n)/n,
  verbose = TRUE,
  seed = seed
)

message("\nSingle result:")
message(sprintf("  theta: %.4f", result_single$theta))
message(sprintf("  sigma: %.4f (sqrt(n) scale)", result_single$sigma))
message(sprintf("  SE: %.4f (theta scale)", result_single$sigma / sqrt(n)))
message(sprintf("  CI: [%.4f, %.4f]", result_single$ci_95[1], result_single$ci_95[2]))

message("\n=== Repeated Cross-Fit (M=10) ===")
result_repeated <- dml_att_repeated(
  d$X, d$A, d$Y, K = K,
  use_rashomon = FALSE,
  regularization = log(n)/n,
  verbose = TRUE,
  seed = seed,
  n_splits = 10,
  aggregation = "median"
)

message("\nRepeated result:")
message(sprintf("  theta: %.4f", result_repeated$theta))
message(sprintf("  sigma: %.4f (sqrt(n) scale)", result_repeated$sigma))
message(sprintf("  SE: %.4f (theta scale)", result_repeated$sigma / sqrt(n)))
message(sprintf("  CI: [%.4f, %.4f]", result_repeated$ci_95[1], result_repeated$ci_95[2]))
message(sprintf("  theta_splits range: [%.4f, %.4f]",
                min(result_repeated$theta_splits), max(result_repeated$theta_splits)))
message(sprintf("  sigma_splits range: [%.4f, %.4f]",
                min(result_repeated$sigma_splits), max(result_repeated$sigma_splits)))
