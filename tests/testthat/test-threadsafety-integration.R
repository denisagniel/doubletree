# Integration tests for thread-safety with optimaltrees worker_limit > 1
# Verifies that the thread-safety fixes in optimaltrees (2026-03-13) work
# correctly in DML workflows.

test_that("dml_att produces identical results with worker_limit=1 vs worker_limit=4", {
  skip_if_not_installed("optimaltrees")

  set.seed(42)
  n <- 150
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5),
    X2 = rbinom(n, 1, 0.5),
    X3 = rbinom(n, 1, 0.5)
  )

  # DGP with treatment effect
  e <- plogis(-0.5 + 0.8 * X$X1 - 0.3 * X$X2)
  A <- rbinom(n, 1, e)
  tau <- 0.3
  Y <- rbinom(n, 1, plogis(-0.2 + 0.5 * X$X1 + 0.4 * X$X2 + tau * A))

  # Fit with worker_limit=1 (single-threaded)
  fit1 <- dml_att(
    X, A, Y,
    K = 3,
    regularization = 0.1,
    seed = 123,
    verbose = FALSE,
    worker_limit = 1
  )

  # Fit with worker_limit=4 (multi-threaded)
  fit4 <- dml_att(
    X, A, Y,
    K = 3,
    regularization = 0.1,
    seed = 123,
    verbose = FALSE,
    worker_limit = 4
  )

  # Results should be identical (within numerical precision)
  expect_equal(fit1$theta, fit4$theta, tolerance = 1e-10)
  expect_equal(fit1$sigma, fit4$sigma, tolerance = 1e-10)
  expect_equal(fit1$ci_95, fit4$ci_95, tolerance = 1e-10)
  expect_equal(fit1$score_values, fit4$score_values, tolerance = 1e-10)

  # Nuisance predictions should also match
  expect_equal(fit1$nuisance_fits$propensity, fit4$nuisance_fits$propensity, tolerance = 1e-10)
  expect_equal(fit1$nuisance_fits$outcome_control, fit4$nuisance_fits$outcome_control, tolerance = 1e-10)
})

test_that("dml_att with continuous outcome and worker_limit=4 is stable", {
  skip_if_not_installed("optimaltrees")

  set.seed(456)
  n <- 120
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5),
    X2 = rbinom(n, 1, 0.5)
  )

  # Continuous outcome DGP
  e <- plogis(0.5 * X$X1 - 0.2 * X$X2)
  A <- rbinom(n, 1, e)
  tau <- 0.5
  Y <- rnorm(n, mean = 1 + 0.6 * X$X1 + 0.4 * X$X2 + tau * A, sd = 0.8)

  # Fit with worker_limit=1
  fit1 <- dml_att(
    X, A, Y,
    K = 3,
    outcome_type = "continuous",
    regularization = 0.1,
    seed = 789,
    verbose = FALSE,
    worker_limit = 1
  )

  # Fit with worker_limit=4
  fit4 <- dml_att(
    X, A, Y,
    K = 3,
    outcome_type = "continuous",
    regularization = 0.1,
    seed = 789,
    verbose = FALSE,
    worker_limit = 4
  )

  # Results should be identical
  expect_equal(fit1$theta, fit4$theta, tolerance = 1e-10)
  expect_equal(fit1$sigma, fit4$sigma, tolerance = 1e-10)
  expect_equal(fit1$ci_95, fit4$ci_95, tolerance = 1e-10)
})

test_that("dml_att with use_rashomon=TRUE and worker_limit=4 works correctly", {
  skip_if_not_installed("optimaltrees")

  set.seed(101)
  n <- 150
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5),
    X2 = rbinom(n, 1, 0.5)
  )

  A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
  Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)

  # Fit Rashomon-DML with worker_limit=1
  fit1 <- dml_att(
    X, A, Y,
    K = 3,
    use_rashomon = TRUE,
    rashomon_bound_multiplier = 0.1,  # Loose bound to ensure intersection
    regularization = 0.1,
    seed = 202,
    verbose = FALSE,
    worker_limit = 1
  )

  # Fit Rashomon-DML with worker_limit=4
  fit4 <- dml_att(
    X, A, Y,
    K = 3,
    use_rashomon = TRUE,
    rashomon_bound_multiplier = 0.1,
    regularization = 0.1,
    seed = 202,
    verbose = FALSE,
    worker_limit = 4
  )

  # Results should be identical
  expect_equal(fit1$theta, fit4$theta, tolerance = 1e-10)
  expect_equal(fit1$sigma, fit4$sigma, tolerance = 1e-10)
  expect_equal(fit1$ci_95, fit4$ci_95, tolerance = 1e-10)

  # Both should have valid estimates
  expect_true(is.finite(fit1$theta))
  expect_true(is.finite(fit4$theta))
})

test_that("stress test: multiple dml_att runs with worker_limit=4 are stable", {
  skip_if_not_installed("optimaltrees")

  set.seed(303)
  n <- 100
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5),
    X2 = rbinom(n, 1, 0.5)
  )

  A <- rbinom(n, 1, plogis(0.3 * X$X1))
  Y <- rbinom(n, 1, 0.4 + 0.2 * X$X1 + 0.15 * A)

  # Run 10 times with same seed - should get identical results
  results <- replicate(10, {
    fit <- dml_att(
      X, A, Y,
      K = 3,
      regularization = 0.1,
      seed = 404,
      verbose = FALSE,
      worker_limit = 4
    )
    fit$theta
  })

  # All estimates should be identical (no race conditions)
  expect_equal(sd(results), 0, tolerance = 1e-10)
  expect_true(length(unique(results)) == 1)
})
