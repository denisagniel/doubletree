# Intensive DML integration stress test
# Verifies dmltree works correctly with updated optimaltrees thread-safety fixes

library(dmltree)

cat("═══════════════════════════════════════════════════════\n")
cat("DML INTEGRATION STRESS TEST\n")
cat("Testing dmltree with optimaltrees thread-safety fixes\n")
cat("═══════════════════════════════════════════════════════\n\n")

# Setup DGP
set.seed(2026)
n <- 200
X <- data.frame(
  X1 = rbinom(n, 1, 0.5),
  X2 = rbinom(n, 1, 0.5),
  X3 = rbinom(n, 1, 0.5)
)

# True treatment effect
tau_true <- 0.4

# Propensity
e <- plogis(-0.5 + 0.8 * X$X1 - 0.3 * X$X2 + 0.2 * X$X3)
A <- rbinom(n, 1, e)

# Outcome with treatment effect
Y <- rbinom(n, 1, plogis(-0.3 + 0.5 * X$X1 + 0.4 * X$X2 - 0.2 * X$X3 + tau_true * A))

cat("DGP Setup:\n")
cat("  n:", n, "\n")
cat("  Features: 3 binary\n")
cat("  True ATT:", tau_true, "\n")
cat("  Pr(A=1):", mean(A), "\n")
cat("  Pr(Y=1):", mean(Y), "\n\n")

# ─────────────────────────────────────────
# Test 1: Consistency across worker_limit
# ─────────────────────────────────────────
cat("Test 1: Consistency across worker_limit settings\n")
cat("─────────────────────────────────────────────────\n")

results_w1 <- list()
results_w4 <- list()
n_runs <- 5

for (i in 1:n_runs) {
  cat(sprintf("  Run %d/%d...", i, n_runs))

  # worker_limit = 1
  fit_w1 <- dml_att(
    X, A, Y,
    K = 3,
    regularization = 0.1,
    seed = 1000 + i,
    verbose = FALSE,
    worker_limit = 1
  )
  results_w1[[i]] <- fit_w1$theta

  # worker_limit = 4
  fit_w4 <- dml_att(
    X, A, Y,
    K = 3,
    regularization = 0.1,
    seed = 1000 + i,
    verbose = FALSE,
    worker_limit = 4
  )
  results_w4[[i]] <- fit_w4$theta

  cat(sprintf(" θ(w=1)=%.4f, θ(w=4)=%.4f, diff=%.2e\n",
              fit_w1$theta, fit_w4$theta, abs(fit_w1$theta - fit_w4$theta)))
}

theta_w1 <- unlist(results_w1)
theta_w4 <- unlist(results_w4)
diffs <- abs(theta_w1 - theta_w4)

cat("\nResults:\n")
cat("  worker_limit=1: mean=", mean(theta_w1), " sd=", sd(theta_w1), "\n")
cat("  worker_limit=4: mean=", mean(theta_w4), " sd=", sd(theta_w4), "\n")
cat("  Max difference:", max(diffs), "\n")
cat("  All differences < 1e-10:", all(diffs < 1e-10), "\n")

if (all(diffs < 1e-10)) {
  cat("  ✓ PASS: Results identical across worker_limit\n\n")
} else {
  cat("  ✗ FAIL: Results differ across worker_limit\n\n")
  stop("Integration test failed: worker_limit consistency")
}

# ─────────────────────────────────────────
# Test 2: Continuous outcomes
# ─────────────────────────────────────────
cat("Test 2: Continuous outcomes\n")
cat("────────────────────────────\n")

Y_cont <- rnorm(n, mean = 1 + 0.6 * X$X1 + 0.4 * X$X2 + tau_true * A, sd = 0.8)

fit_cont_w1 <- dml_att(
  X, A, Y_cont,
  K = 3,
  outcome_type = "continuous",
  regularization = 0.1,
  seed = 2000,
  verbose = FALSE,
  worker_limit = 1
)

fit_cont_w4 <- dml_att(
  X, A, Y_cont,
  K = 3,
  outcome_type = "continuous",
  regularization = 0.1,
  seed = 2000,
  verbose = FALSE,
  worker_limit = 4
)

cat(sprintf("  worker_limit=1: θ=%.4f, σ=%.4f\n", fit_cont_w1$theta, fit_cont_w1$sigma))
cat(sprintf("  worker_limit=4: θ=%.4f, σ=%.4f\n", fit_cont_w4$theta, fit_cont_w4$sigma))
cat(sprintf("  Difference: %.2e\n", abs(fit_cont_w1$theta - fit_cont_w4$theta)))

if (abs(fit_cont_w1$theta - fit_cont_w4$theta) < 1e-10) {
  cat("  ✓ PASS: Continuous outcomes consistent\n\n")
} else {
  cat("  ✗ FAIL: Continuous outcomes differ\n\n")
  stop("Integration test failed: continuous outcomes")
}

# ─────────────────────────────────────────
# Test 3: Rashomon-DML
# ─────────────────────────────────────────
cat("Test 3: Rashomon-DML integration\n")
cat("─────────────────────────────────\n")

fit_rash_w1 <- dml_att(
  X, A, Y,
  K = 3,
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 0.1,
  regularization = 0.1,
  seed = 3000,
  verbose = FALSE,
  worker_limit = 1
)

fit_rash_w4 <- dml_att(
  X, A, Y,
  K = 3,
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 0.1,
  regularization = 0.1,
  seed = 3000,
  verbose = FALSE,
  worker_limit = 4
)

cat(sprintf("  worker_limit=1: θ=%.4f, σ=%.4f\n", fit_rash_w1$theta, fit_rash_w1$sigma))
cat(sprintf("  worker_limit=4: θ=%.4f, σ=%.4f\n", fit_rash_w4$theta, fit_rash_w4$sigma))
cat(sprintf("  Difference: %.2e\n", abs(fit_rash_w1$theta - fit_rash_w4$theta)))

if (abs(fit_rash_w1$theta - fit_rash_w4$theta) < 1e-10) {
  cat("  ✓ PASS: Rashomon-DML consistent\n\n")
} else {
  cat("  ✗ FAIL: Rashomon-DML differs\n\n")
  stop("Integration test failed: Rashomon-DML")
}

# ─────────────────────────────────────────
# Test 4: Stability with repeated runs
# ─────────────────────────────────────────
cat("Test 4: Stability over 20 runs (worker_limit=4)\n")
cat("────────────────────────────────────────────────\n")

results_stability <- replicate(20, {
  fit <- dml_att(
    X, A, Y,
    K = 3,
    regularization = 0.1,
    seed = 4000,
    verbose = FALSE,
    worker_limit = 4
  )
  fit$theta
})

cat(sprintf("  Mean: %.4f\n", mean(results_stability)))
cat(sprintf("  SD: %.2e\n", sd(results_stability)))
cat(sprintf("  Range: %.4f - %.4f\n", min(results_stability), max(results_stability)))
cat(sprintf("  Unique values: %d\n", length(unique(results_stability))))

if (sd(results_stability) < 1e-10) {
  cat("  ✓ PASS: Perfect stability across 20 runs\n\n")
} else {
  cat("  ✗ FAIL: Estimates vary across runs\n\n")
  stop("Integration test failed: stability")
}

# ─────────────────────────────────────────
# Summary
# ─────────────────────────────────────────
cat("═══════════════════════════════════════════════════════\n")
cat("ALL TESTS PASSED ✓\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("Summary:\n")
cat("  ✓ worker_limit consistency (5 runs each)\n")
cat("  ✓ Continuous outcomes\n")
cat("  ✓ Rashomon-DML integration\n")
cat("  ✓ Stability over 20 runs\n\n")

cat("Conclusion:\n")
cat("  dmltree works correctly with optimaltrees v0.4.0\n")
cat("  Thread-safety fixes do not break DML workflows\n")
cat("  Multi-threading produces identical results to single-threading\n\n")

cat("Test completed at:", format(Sys.time()), "\n")
