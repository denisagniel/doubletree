# Test parallel speed improvements from thread-safety fixes
library(optimaltrees)
library(dplyr)

cat("\n=== PARALLEL SPEED TEST ===\n\n")

# Source dmltree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

N <- 800
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 10
SEED_OFFSET <- 70000

cat("Configuration:\n")
cat("  n =", N, "\n")
cat("  Replications:", N_REPS, "\n")
cat("  Testing: worker_limit 1 vs 4\n\n")

# Test 1: Single-threaded (worker_limit=1)
cat("Test 1: worker_limit = 1\n")
start_time1 <- Sys.time()

results1 <- lapply(1:N_REPS, function(rep) {
  d <- generate_dgp_binary_att(n = N, tau = TAU, seed = SEED_OFFSET + rep)
  
  fit <- estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    K = K_FOLDS,
    outcome_type = "binary",
    regularization = log(N) / N,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    worker_limit = 1,  # Single-threaded
    verbose = FALSE
  )
  
  fit$theta
})

time1 <- as.numeric(difftime(Sys.time(), start_time1, units = "secs"))
cat(sprintf("  Completed in %.2f seconds (%.2f sec/rep)\n\n", time1, time1/N_REPS))

# Test 2: Multi-threaded (worker_limit=4)
cat("Test 2: worker_limit = 4\n")
start_time2 <- Sys.time()

results2 <- lapply(1:N_REPS, function(rep) {
  d <- generate_dgp_binary_att(n = N, tau = TAU, seed = SEED_OFFSET + rep)
  
  fit <- estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    K = K_FOLDS,
    outcome_type = "binary",
    regularization = log(N) / N,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    worker_limit = 4,  # Multi-threaded
    verbose = FALSE
  )
  
  fit$theta
})

time2 <- as.numeric(difftime(Sys.time(), start_time2, units = "secs"))
cat(sprintf("  Completed in %.2f seconds (%.2f sec/rep)\n\n", time2, time2/N_REPS))

# Results comparison
cat("=== RESULTS ===\n\n")
cat(sprintf("Single-threaded: %.2f seconds\n", time1))
cat(sprintf("Multi-threaded:  %.2f seconds\n", time2))
cat(sprintf("Speedup:         %.2fx\n\n", time1/time2))

# Verify consistency
results1_vec <- unlist(results1)
results2_vec <- unlist(results2)
max_diff <- max(abs(results1_vec - results2_vec))

cat("Consistency check:\n")
cat(sprintf("  Max difference: %.6f\n", max_diff))
if (max_diff < 1e-6) {
  cat("  Status: PASS - Results identical\n")
} else {
  cat("  Status: WARNING - Results differ slightly\n")
}

cat("\n=== VERDICT ===\n\n")
if (max_diff < 1e-6 && time2 < time1) {
  cat(sprintf("PASS: Thread-safety working correctly!\n"))
  cat(sprintf("  - Results identical between worker_limit=1 and worker_limit=4\n"))
  cat(sprintf("  - Speedup: %.2fx with 4 workers\n", time1/time2))
  cat(sprintf("  - Simulations can now use parallelization\n"))
} else if (max_diff < 1e-6) {
  cat("PASS: Results consistent but no speedup\n")
  cat("  (Dataset may be too small to benefit from parallelization)\n")
} else {
  cat("FAIL: Results differ between single and multi-threaded\n")
  cat("  Thread-safety issue may remain\n")
}

cat("\n")
