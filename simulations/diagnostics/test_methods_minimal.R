# Minimal test - just test if baseline method functions work
# Assumes packages are already compiled/loaded

cat("Sourcing DGPs and methods...\n")
source("../production/dgps/dgps_smooth.R")
source("../production/methods/method_forest_dml.R")
source("../production/methods/method_linear_dml.R")

cat("✓ Functions loaded\n\n")

# Generate test data
set.seed(123)
n <- 100  # Small for speed
tau <- 0.10
d <- generate_dgp_binary_att(n, tau = tau, seed = 123)

cat("Test Data Generated:\n")
cat(sprintf("  n = %d, n_treated = %d\n", n, sum(d$A)))
cat(sprintf("  True ATT = %.4f\n\n", d$true_att))

# Test Forest-DML
cat("Testing Forest-DML...\n")
start_time <- Sys.time()
result_forest <- tryCatch({
  att_forest(
    X = d$X, A = d$A, Y = d$Y,
    K = 3,  # Fewer folds for speed
    seed = 123,
    num.trees = 100,  # Fewer trees for speed
    verbose = FALSE
  )
}, error = function(e) {
  list(error = e$message)
})

if (is.null(result_forest$error)) {
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("  ✓ SUCCESS (%.1f sec)\n", elapsed))
  cat(sprintf("    θ̂ = %.4f, SE = %.4f\n", result_forest$theta, result_forest$sigma))
  cat(sprintf("    95%% CI: [%.4f, %.4f]\n", result_forest$ci[1], result_forest$ci[2]))
  cat(sprintf("    Covers truth: %s\n\n",
              if(result_forest$ci[1] <= d$true_att && result_forest$ci[2] >= d$true_att) "YES" else "NO"))
} else {
  cat(sprintf("  ✗ FAILED: %s\n\n", result_forest$error))
}

# Test Linear-DML
cat("Testing Linear-DML...\n")
start_time <- Sys.time()
result_linear <- tryCatch({
  att_linear(
    X = d$X, A = d$A, Y = d$Y,
    K = 3,
    seed = 123,
    verbose = FALSE
  )
}, error = function(e) {
  list(error = e$message)
})

if (is.null(result_linear$error)) {
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("  ✓ SUCCESS (%.1f sec)\n", elapsed))
  cat(sprintf("    θ̂ = %.4f, SE = %.4f\n", result_linear$theta, result_linear$sigma))
  cat(sprintf("    95%% CI: [%.4f, %.4f]\n", result_linear$ci[1], result_linear$ci[2]))
  cat(sprintf("    Covers truth: %s\n\n",
              if(result_linear$ci[1] <= d$true_att && result_linear$ci[2] >= d$true_att) "YES" else "NO"))
} else {
  cat(sprintf("  ✗ FAILED: %s\n\n", result_linear$error))
}

# Summary
cat(strrep("=", 60), "\n")
if (is.null(result_forest$error) && is.null(result_linear$error)) {
  cat("✓ Both baseline methods work correctly\n")
  cat(sprintf("  Forest estimate: %.4f\n", result_forest$theta))
  cat(sprintf("  Linear estimate: %.4f\n", result_linear$theta))
  cat(sprintf("  Difference:      %.4f\n", abs(result_forest$theta - result_linear$theta)))
} else {
  cat("✗ One or more methods failed\n")
}
cat(strrep("=", 60), "\n")
