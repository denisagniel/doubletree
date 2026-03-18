# Quick verification that reinstalled package works
# Test just the core optimaltrees functionality without verbose output

library(optimaltrees)

cat("=== Quick Verification Test ===\n\n")

# Test 1: Simple fit_tree with model_limit safeguard
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
y <- rbinom(n, 1, 0.5)

cat("Test 1: fit_tree with n=100, p=3...\n")
m1 <- fit_tree(
  X, y,
  loss_function = "log_loss",
  regularization = 0.1,
  discretize_bins = "adaptive",
  verbose = FALSE
)

cat("  n_trees:", m1$n_trees, "\n")
cat("  tree_json exists:", !is.null(m1$model$tree_json), "\n")
cat("  result_data exists:", !is.null(m1$model$result_data), "\n")

if (!is.null(m1$model$tree_json) || !is.null(m1$model$result_data)) {
  cat("  ✓ Model extraction successful\n")
} else {
  cat("  ✗ Model extraction failed\n")
  quit(status = 1)
}

# Test 2: Check discretization metadata exists
cat("  Discretization metadata exists:", !is.null(m1$discretization_metadata), "\n")

# Test 3: Larger sample (n=800) like in actual simulations
cat("\nTest 2: fit_tree with n=800, p=4 (simulation scale)...\n")
set.seed(123)
n <- 800
X <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
y <- rbinom(n, 1, 0.5)

m2 <- fit_tree(
  X, y,
  loss_function = "log_loss",
  regularization = log(n) / n,
  discretize_bins = "adaptive",
  verbose = FALSE
)

cat("  n_trees:", m2$n_trees, "\n")
cat("  tree_json exists:", !is.null(m2$model$tree_json), "\n")
cat("  result_data exists:", !is.null(m2$model$result_data), "\n")

if (!is.null(m2$model$tree_json) || !is.null(m2$model$result_data)) {
  cat("  ✓ Model extraction successful\n")
} else {
  cat("  ✗ Model extraction failed\n")
  quit(status = 1)
}

cat("  Discretization metadata exists:", !is.null(m2$discretization_metadata), "\n")

cat("\n=== All Verification Tests Passed ===\n")
cat("Ready for full simulations!\n")
