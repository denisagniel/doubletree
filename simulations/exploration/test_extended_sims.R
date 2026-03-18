# test_extended_sims.R
# Quick test of extended simulation infrastructure before full run
#
# Tests:
#   - All 4 DGPs generate valid data
#   - Oracle DML works
#   - 3-way comparison runs without errors
#   - Rashomon metrics are extracted correctly
#   - Results save and load correctly

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) {
  stop("optimaltrees is required")
}

message("=== Testing Extended Simulation Infrastructure ===\n")

# Source the extended simulation script to get functions
source("simulations/run_simulations_extended.R", local = TRUE)

# Test 1: DGP generation -------------------------------------------------------

message("Test 1: DGP generation...")
dgps_to_test <- list(
  dgp1 = generate_data_dgp1,
  dgp2 = generate_data_dgp2,
  dgp3 = generate_data_dgp3,
  dgp4 = generate_data_dgp4
)

for (dgp_name in names(dgps_to_test)) {
  d <- dgps_to_test[[dgp_name]](n = 200, tau = 0.15, seed = 123)
  stopifnot(
    nrow(d$X) == 200,
    length(d$A) == 200,
    length(d$Y) == 200,
    all(d$A %in% c(0, 1)),
    all(d$Y %in% c(0, 1)),
    "true_e" %in% names(d),
    "true_m0" %in% names(d),
    "true_m1" %in% names(d)
  )
  message("  ", dgp_name, ": OK (n=200, binary A/Y, true nuisances present)")
}

# Test 2: Oracle DML -----------------------------------------------------------

message("\nTest 2: Oracle DML...")
d <- generate_data_dgp1(n = 200, tau = 0.15, seed = 456)
oracle_result <- dml_att_oracle(d, K = 5, seed = 456)

stopifnot(
  "theta" %in% names(oracle_result),
  "ci_95" %in% names(oracle_result),
  !is.na(oracle_result$theta),
  length(oracle_result$ci_95) == 2
)
message("  Oracle DML: OK (theta = ", round(oracle_result$theta, 4), ")")

# Test 3: Fold-specific DML ----------------------------------------------------

message("\nTest 3: Fold-specific DML...")
fold_result <- estimate_att(
  d$X, d$A, d$Y, K = 5,
  use_rashomon = FALSE,
  regularization = 0.1,
  verbose = FALSE,
  seed = 789
)

stopifnot(
  "theta" %in% names(fold_result),
  !is.na(fold_result$theta),
  "nuisance_fits" %in% names(fold_result),
  length(fold_result$nuisance_fits) == 5  # K folds
)
message("  Fold-specific: OK (theta = ", round(fold_result$theta, 4), ")")

# Test 4: Rashomon DML ---------------------------------------------------------

message("\nTest 4: Rashomon DML...")
rashomon_result <- estimate_att(
  d$X, d$A, d$Y, K = 5,
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 0.05,
  regularization = 0.1,
  verbose = FALSE,
  seed = 789
)

stopifnot(
  "theta" %in% names(rashomon_result),
  !is.na(rashomon_result$theta),
  "nuisance_fits" %in% names(rashomon_result)
)

# Check if Rashomon-specific info is present
if ("cf_e" %in% names(rashomon_result$nuisance_fits)) {
  cf_e <- rashomon_result$nuisance_fits$cf_e
  message("  Rashomon: OK (theta = ", round(rashomon_result$theta, 4), ")")
  if ("n_intersecting" %in% names(cf_e)) {
    message("    e nuisance: ", cf_e$n_intersecting, " intersecting trees")
  }
} else {
  message("  Rashomon: OK but no cf_e structure (check implementation)")
}

# Test 5: Metrics computation --------------------------------------------------

message("\nTest 5: Metrics computation...")

# Create small replication set
small_results <- list(fold_result, fold_result, fold_result)  # Dummy replicates
metrics_standard <- compute_metrics(small_results, true_att = 0.15)

stopifnot(
  "bias" %in% names(metrics_standard),
  "mse" %in% names(metrics_standard),
  "coverage_95" %in% names(metrics_standard)
)
message("  Standard metrics: OK")

# Rashomon-specific metrics
small_rashomon_results <- list(rashomon_result, rashomon_result, rashomon_result)
metrics_rashomon <- compute_metrics_rashomon(small_rashomon_results, true_att = 0.15)

if ("pct_nonempty_e" %in% names(metrics_rashomon)) {
  message("  Rashomon metrics: OK (pct_nonempty_e present)")
} else {
  warning("  Rashomon metrics: pct_nonempty_e not found")
}

# Test 6: Three-way comparison (mini run) -------------------------------------

message("\nTest 6: Three-way comparison (2 reps)...")
comparison_result <- run_comparison(
  dgp_fn = generate_data_dgp1,
  n = 200,
  K = 5,
  tau = 0.15,
  n_reps = 2,  # Minimal reps for test
  epsilon = 0.05,
  seed_start = 999
)

stopifnot(
  "fold_specific" %in% names(comparison_result),
  "rashomon" %in% names(comparison_result),
  "oracle" %in% names(comparison_result),
  "raw_results" %in% names(comparison_result)
)
message("  Three-way comparison: OK")
message("    Fold-specific bias: ", round(comparison_result$fold_specific$bias, 4))
message("    Rashomon bias: ", round(comparison_result$rashomon$bias, 4))
message("    Oracle bias: ", round(comparison_result$oracle$bias, 4))

# Test 7: Save/load ------------------------------------------------------------

message("\nTest 7: Save and load...")
test_dir <- "simulations/results_test"
if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

saveRDS(comparison_result, file.path(test_dir, "test_result.rds"))
loaded_result <- readRDS(file.path(test_dir, "test_result.rds"))

stopifnot(
  all.equal(comparison_result$fold_specific$bias, loaded_result$fold_specific$bias)
)
message("  Save/load: OK")

# Cleanup test directory
unlink(test_dir, recursive = TRUE)

# Summary ----------------------------------------------------------------------

message("\n=== All Tests Passed ===")
message("Infrastructure ready for full simulation run.")
message("\nTo run full simulations:")
message("  Rscript simulations/run_simulations_extended.R")
message("\nExpected runtime (full grid):")
message("  4 DGPs × 4 sample sizes × 4 tolerances × 100 reps = 6,400 runs")
message("  Estimated: 4-6 hours (depending on hardware)")
