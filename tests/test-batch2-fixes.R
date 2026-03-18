# Quick integration test for Batch 2 fixes (dmltree)
# Run with: Rscript tests/test-batch2-fixes.R

# Load development version
devtools::load_all()

cat("=== Testing Batch 2 Fixes (dmltree) ===\n\n")

# Set up basic test data
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)

# Test Issue #16: logical parameter validation
cat("Test 1: logical parameter validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, verbose = "yes", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, stratified = 1, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, cv_regularization = "TRUE", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, use_rashomon = 1, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, auto_tune_intersecting = "false", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: logical parameter validation works\n")

# Test Issue #17: seed validation
cat("Test 2: seed validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, seed = "abc", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, seed = TRUE, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: seed validation works\n")

# Test Issue #18: max_leaves validation
cat("Test 3: max_leaves validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, max_leaves = 0, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, max_leaves = 2.5, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, max_leaves = -1, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: max_leaves validation works\n")

# Test Issue #19: K must be integer
cat("Test 4: K must be integer\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3.5, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: K integer validation works\n")

# Test Issue #20: discretize_method validation
cat("Test 5: discretize_method validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, discretize_method = "invalid", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, discretize_method = "mean", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: discretize_method validation works\n")

# Test Issue #21: rashomon_bound_adder validation
cat("Test 6: rashomon_bound_adder validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, rashomon_bound_adder = -0.1, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, rashomon_bound_adder = "0.1", regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: rashomon_bound_adder validation works\n")

# Test valid inputs still work
cat("Test 7: Valid inputs still work\n")
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)
fit <- estimate_att(X, A, Y, K = 3, regularization = 0.1, verbose = FALSE,
                    seed = 42, max_leaves = NULL, discretize_method = "quantiles")
stopifnot(!is.null(fit))
stopifnot(!is.null(fit$theta))
stopifnot(is.numeric(fit$theta))
cat("  PASS: Valid inputs work correctly\n")

cat("\n=== All Batch 2 (dmltree) tests passed! ===\n")
