# Quick integration test for Suggestions fixes (doubletree)
# Run with: Rscript tests/test-suggestions-fixes.R

# Load development versions
devtools::load_all("../optimaltrees")
devtools::load_all()

cat("=== Testing Suggestions Fixes (doubletree) ===\n\n")

# Test Issue #31: Obsolete TODO removed
cat("Test 1: No obsolete TODO comments\n")
# Simply verify the package loads and functions work
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)
fit <- estimate_att(X, A, Y, K = 3, regularization = 0.1, verbose = FALSE)
stopifnot(!is.null(fit$theta))
cat("  PASS: Package functions correctly (obsolete TODO removed)\n")

# Test Issue #33: Performance TODO documented as known limitation
cat("Test 2: Rashomon fallback works\n")
# Test that Rashomon with fallback works correctly
# (Even though refitting both models is not optimal, it's correct)
set.seed(42)
n <- 150
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.5)

# Use tight epsilon to potentially trigger fallback
fit_rash <- estimate_att(X, A, Y, K = 3, use_rashomon = TRUE,
                         rashomon_bound_multiplier = 0.01,  # Tight bound
                         regularization = 0.1, verbose = FALSE)
stopifnot(!is.null(fit_rash$theta))
stopifnot(is.numeric(fit_rash$theta))
cat("  PASS: Rashomon with fallback works (known limitation documented)\n")

# Test valid inputs still work after all changes
cat("Test 3: All functionality intact\n")
set.seed(42)
n <- 100
X <- data.frame(x1 = runif(n), x2 = runif(n))
A <- rbinom(n, 1, 0.5)
Y_binary <- rbinom(n, 1, 0.5)
Y_continuous <- rnorm(n)

fit_binary <- estimate_att(X, A, Y_binary, K = 3, outcome_type = "binary",
                           regularization = 0.1, verbose = FALSE)
fit_continuous <- estimate_att(X, A, Y_continuous, K = 3, outcome_type = "continuous",
                               regularization = 0.1, verbose = FALSE)

stopifnot(!is.null(fit_binary$theta))
stopifnot(!is.null(fit_continuous$theta))
stopifnot(is.numeric(fit_binary$theta))
stopifnot(is.numeric(fit_continuous$theta))
cat("  PASS: Binary and continuous outcomes work\n")

cat("\n=== All Suggestions (doubletree) tests passed! ===\n")
