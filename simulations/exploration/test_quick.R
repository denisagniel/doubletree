# test_quick.R
# Quick functionality test (no full simulation run)

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

message("=== Quick Functionality Test ===\n")

# Source functions without running full grid
source("simulations/run_simulations_extended.R", local = TRUE)

# Test 1: DGP generation
message("Test 1: All DGPs generate valid data...")
for (dgp_name in c("dgp1", "dgp2", "dgp3", "dgp4")) {
  dgp_fn <- get(paste0("generate_data_", dgp_name))
  d <- dgp_fn(n = 100, tau = 0.15, seed = 123)
  stopifnot(
    nrow(d$X) == 100,
    all(d$A %in% c(0, 1)),
    all(d$Y %in% c(0, 1)),
    "true_e" %in% names(d)
  )
  message("  ", dgp_name, ": OK")
}

# Test 2: Oracle DML
message("\nTest 2: Oracle DML...")
d <- generate_data_dgp1(n = 100, tau = 0.15, seed = 456)
oracle_result <- dml_att_oracle(d, K = 5, seed = 456)
stopifnot("theta" %in% names(oracle_result), !is.na(oracle_result$theta))
message("  OK (theta = ", round(oracle_result$theta, 4), ")")

# Test 3: Fold-specific DML
message("\nTest 3: Fold-specific DML...")
fold_result <- estimate_att(d$X, d$A, d$Y, K = 5, use_rashomon = FALSE,
                      regularization = 0.1, verbose = FALSE, seed = 789)
stopifnot("theta" %in% names(fold_result), !is.na(fold_result$theta))
message("  OK (theta = ", round(fold_result$theta, 4), ")")

# Test 4: Rashomon DML
message("\nTest 4: Rashomon DML...")
rashomon_result <- estimate_att(d$X, d$A, d$Y, K = 5, use_rashomon = TRUE,
                          rashomon_bound_multiplier = 0.05,
                          regularization = 0.1, verbose = FALSE, seed = 789)
stopifnot("theta" %in% names(rashomon_result), !is.na(rashomon_result$theta))
message("  OK (theta = ", round(rashomon_result$theta, 4), ")")

# Test 5: Metrics
message("\nTest 5: Metrics computation...")
small_results <- list(fold_result, rashomon_result)
metrics_standard <- compute_metrics(small_results, true_att = 0.15)
stopifnot("bias" %in% names(metrics_standard))
message("  Standard metrics: OK")

small_rashomon <- list(rashomon_result, rashomon_result)
metrics_rashomon <- compute_metrics_rashomon(small_rashomon, true_att = 0.15)
message("  Rashomon metrics: ", if ("pct_nonempty_e" %in% names(metrics_rashomon)) "OK" else "WARNING")

# Test 6: Mini three-way comparison (2 reps)
message("\nTest 6: Three-way comparison (2 reps)...")
comparison_result <- run_comparison(
  dgp_fn = generate_data_dgp1,
  n = 100,
  K = 5,
  tau = 0.15,
  n_reps = 2,
  epsilon = 0.05,
  seed_start = 999
)
stopifnot(
  "fold_specific" %in% names(comparison_result),
  "rashomon" %in% names(comparison_result),
  "oracle" %in% names(comparison_result)
)
message("  OK")
message("    Biases: fold=", round(comparison_result$fold_specific$bias, 4),
        ", rash=", round(comparison_result$rashomon$bias, 4),
        ", oracle=", round(comparison_result$oracle$bias, 4))

message("\n=== All Tests Passed ===")
message("Infrastructure is ready for full simulation run.")
message("\nTo run full simulations:")
message("  Rscript simulations/run_simulations_extended.R")
message("\nExpected time: ~4-6 hours (64 configurations × 100 reps)")
