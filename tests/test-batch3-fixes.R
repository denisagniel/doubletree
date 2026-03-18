# Quick integration test for Batch 3 fixes (dmltree)
# Run with: Rscript tests/test-batch3-fixes.R

# Load development versions
devtools::load_all("../optimaltrees")
devtools::load_all()

cat("=== Testing Batch 3 Fixes (dmltree) ===\n\n")

# Set up basic test data
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)

# Test Issue #27: cv_K must be integer
cat("Test 1: cv_K must be integer\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, cv_regularization = TRUE, cv_K = 3.5, regularization = 0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: Non-integer cv_K caught\n")

# Test Issue #28: regularization validation (already existed, just verify it works)
cat("Test 2: regularization validation\n")
result <- tryCatch({
  estimate_att(X, A, Y, K = 3, regularization = -0.1)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")

result <- tryCatch({
  estimate_att(X, A, Y, K = 3, regularization = 0)
  "NO_ERROR"
}, error = function(e) "ERROR")
stopifnot(result == "ERROR")
cat("  PASS: Invalid regularization caught\n")

# Test Issue #29: outcome_type validation
cat("Test 3: outcome_type validation (internal)\n")
# This is validated inside predict_nuisances_fold, which is internal
# We test it indirectly by ensuring valid inputs work
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y_binary <- rbinom(n, 1, 0.5)
Y_continuous <- rnorm(n)

fit_binary <- estimate_att(X, A, Y_binary, K = 3, outcome_type = "binary",
                           regularization = 0.1, verbose = FALSE)
stopifnot(!is.null(fit_binary$theta))

fit_continuous <- estimate_att(X, A, Y_continuous, K = 3, outcome_type = "continuous",
                               regularization = 0.1, verbose = FALSE)
stopifnot(!is.null(fit_continuous$theta))
cat("  PASS: Both outcome types work correctly\n")

# Test Issue #30: Propensity clamping parameters exist (infrastructure already in place)
cat("Test 4: Propensity clamping infrastructure\n")
# Constants are defined in score_att.R and used in get_fold_specific_eta
# Verify they exist
propensity_lower <- get(".PROPENSITY_LOWER_BOUND", envir = asNamespace("doubletree"))
propensity_upper <- get(".PROPENSITY_UPPER_BOUND", envir = asNamespace("doubletree"))
stopifnot(propensity_lower == 1e-6)
stopifnot(propensity_upper == (1 - 1e-6))
cat("  PASS: Propensity bounds defined correctly\n")

# Test valid inputs still work
cat("Test 5: Valid inputs still work\n")
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)
fit <- estimate_att(X, A, Y, K = 3, regularization = 0.1, verbose = FALSE,
                    cv_regularization = FALSE, discretize_method = "quantiles")
stopifnot(!is.null(fit))
stopifnot(!is.null(fit$theta))
stopifnot(is.numeric(fit$theta))
cat("  PASS: Valid inputs work correctly\n")

cat("\n=== All Batch 3 (dmltree) tests passed! ===\n")
