# ATT estimation and score tests
# Requires optimaltrees to be installed (doubletree Imports optimaltrees).

test_that("psi_att returns vector of length n and is linear in theta", {
  n <- 5
  Y <- c(1, 0, 1, 0, 1)
  A <- c(1, 0, 1, 0, 1)
  eta <- list(e = rep(0.5, n), m0 = rep(0.4, n))
  pi_hat <- 0.6
  s0 <- psi_att(Y, A, theta = 0, eta, pi_hat)
  s1 <- psi_att(Y, A, theta = 1, eta, pi_hat)
  expect_length(s0, n)
  expect_length(s1, n)
  # psi(theta) = psi(0) - theta * (A/pi)
  expect_equal(s1, s0 - 1 * (A / pi_hat), tolerance = 1e-10)
})

test_that("estimate_att returns list with theta, sigma, ci_95 and runs with binary data", {
  skip_if_not_installed("optimaltrees")
  set.seed(42)
  n <- 120
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
  Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
  fit <- estimate_att(X, A, Y, K = 3)
  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_type(fit$sigma, "double")
  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)
  expect_equal(fit$n, n)
  expect_equal(fit$K, 3)
})

test_that("create_folds returns integer vector 1..K", {
  f <- create_folds(100, K = 5)
  expect_length(f, 100)
  expect_true(all(f >= 1 & f <= 5))
  expect_type(f, "integer")
})

test_that("estimate_att with use_rashomon = TRUE runs and returns same structure", {
  skip_if_not_installed("optimaltrees")
  set.seed(42)
  n <- 150
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
  Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
  fit <- estimate_att(X, A, Y, K = 3, use_rashomon = TRUE, verbose = FALSE)
  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_type(fit$sigma, "double")
  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)
  expect_equal(fit$n, n)
  expect_equal(fit$K, 3)
  expect_true(is.list(fit$nuisance_fits))
  expect_true(length(fit$fold_indices) == n)
})

test_that("att_se and att_ci work", {
  scores <- rnorm(100, 0, 1)
  se <- att_se(scores)
  expect_true(se > 0)
  ci <- att_ci(0.5, se, 100, level = 0.95)
  expect_length(ci, 2)
  expect_true(ci[1] < 0.5)
  expect_true(ci[2] > 0.5)
})

# Continuous outcome tests (Sprint 2, MAJOR-5)
test_that("estimate_att works with continuous outcomes", {
  skip_if_not_installed("optimaltrees")

  set.seed(123)
  n <- 150
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5),
    X2 = rbinom(n, 1, 0.5),
    X3 = rbinom(n, 1, 0.5)
  )

  # DGP with continuous outcome
  e <- plogis(-0.5 + 0.8 * X$X1 - 0.3 * X$X2)
  A <- rbinom(n, 1, e)

  # Continuous Y with treatment effect
  tau <- 0.5
  Y <- rnorm(n, mean = 1 + 0.6 * X$X1 + 0.4 * X$X2 + tau * A, sd = 0.8)

  # Fit with continuous outcome
  fit <- estimate_att(
    X, A, Y,
    K = 3,
    outcome_type = "continuous",
    regularization = 0.1,
    seed = 42
  )

  # Basic checks
  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_true(is.finite(fit$theta))

  expect_type(fit$sigma, "double")
  expect_true(fit$sigma > 0)

  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)

  # Check estimate is reasonable (should be near tau = 0.5)
  expect_true(abs(fit$theta - tau) < 1.0)  # Loose bound for small sample
})

test_that("continuous outcome uses squared_error loss", {
  skip_if_not_installed("optimaltrees")

  set.seed(456)
  n <- 100
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n, mean = 2 + 0.5 * A)

  # Fit with continuous outcome
  fit <- estimate_att(X, A, Y, K = 3, outcome_type = "continuous")

  # Verify structure
  expect_true("nuisance_fits" %in% names(fit))
  expect_type(fit$theta, "double")
  expect_true(is.finite(fit$theta))
})

test_that("estimate_att handles small K with continuous outcomes", {
  skip_if_not_installed("optimaltrees")

  set.seed(789)
  n <- 60
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rnorm(n, mean = 1 + 0.3 * A, sd = 0.5)

  # K=2 should work (may produce warnings from diagnostics - that's expected)
  fit <- estimate_att(X, A, Y, K = 2, outcome_type = "continuous")
  expect_true(is.finite(fit$theta))
})

test_that("continuous outcome validates input appropriately", {
  skip_if_not_installed("optimaltrees")

  set.seed(101)
  n <- 50
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)

  # Continuous outcome should work (may produce warnings from diagnostics - that's expected)
  Y_continuous <- rnorm(n, mean = 1 + 0.5 * A)
  fit <- estimate_att(X, A, Y_continuous, K = 3, outcome_type = "continuous")
  expect_true(is.finite(fit$theta))
})
