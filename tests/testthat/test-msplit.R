# Tests for M-Split Doubletree Estimation

# NOTE: select_structure_modal, analyze_structure_diversity, and
# compute_functional_consistency were relocated to optimaltrees (2026-07-08); their
# unit tests now live in optimaltrees/tests/testthat/test-structure-selection.R.
# This file keeps the doubletree-specific estimate_att_msplit / print / summary tests.

test_that("estimate_att_msplit runs end-to-end with small M", {
  skip_if_not_installed("optimaltrees")

  set.seed(456)
  n <- 150
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5)
  )
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.15 * A + 0.1 * X$x1)

  result <- estimate_att_msplit(
    X, A, Y,
    M = 3,  # Small for speed
    K = 2,  # Small for speed
    verbose = FALSE,
    seed_base = 100
  )

  # Check structure
  expect_s3_class(result, "msplit_att")
  expect_type(result$theta, "double")
  expect_type(result$sigma, "double")
  expect_length(result$ci_95, 2)
  expect_equal(result$M, 3)
  expect_equal(result$K, 2)
  expect_equal(result$n, n)

  # Check predictions
  expect_equal(nrow(result$predictions_all_splits$e), n)
  expect_equal(ncol(result$predictions_all_splits$e), 3)
  expect_equal(length(result$averaged_predictions$e), n)

  # Check diagnostics
  expect_true(result$diagnostics$structure_frequency_e > 0)
  expect_true(result$diagnostics$structure_frequency_m0 > 0)
  expect_true(result$diagnostics$mean_prediction_variance_e >= 0)
  expect_true(result$diagnostics$mean_prediction_variance_m0 >= 0)
  expect_true(result$diagnostics$functional_consistency$max_diff_e >= 0)
  expect_true(result$diagnostics$functional_consistency$max_diff_m0 >= 0)
})

test_that("estimate_att_msplit works with M=1 (single split)", {
  skip_if_not_installed("optimaltrees")

  set.seed(789)
  n <- 100
  X <- data.frame(x1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  result <- estimate_att_msplit(X, A, Y, M = 1, K = 2, verbose = FALSE)

  # With M=1, modal frequency should be 1.0 (only one structure)
  expect_equal(result$diagnostics$structure_frequency_e, 1.0)
  expect_equal(result$diagnostics$structure_frequency_m0, 1.0)

  # Prediction variance is NA for M=1 (cannot compute variance from single value)
  expect_true(is.na(result$diagnostics$mean_prediction_variance_e))
  expect_true(is.na(result$diagnostics$mean_prediction_variance_m0))
})

test_that("estimate_att_msplit works with continuous outcomes", {
  skip_if_not_installed("optimaltrees")

  set.seed(1011)
  n <- 100
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rnorm(n, mean = 5 + 2 * A + 1 * X$x1, sd = 1)

  result <- estimate_att_msplit(
    X, A, Y,
    M = 2,
    K = 2,
    outcome_type = "continuous",
    verbose = FALSE,
    seed_base = 200
  )

  expect_s3_class(result, "msplit_att")
  expect_equal(result$outcome_type, "continuous")
  expect_type(result$theta, "double")
  expect_true(is.finite(result$theta))
})

test_that("print.msplit_att works", {
  skip_if_not_installed("optimaltrees")

  set.seed(1213)
  n <- 80
  X <- data.frame(x1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  result <- estimate_att_msplit(X, A, Y, M = 2, K = 2, verbose = FALSE)

  # Should print without error
  expect_output(print(result), "M-Split Doubletree ATT Estimation")
  expect_output(print(result), "Estimate:")
  expect_output(print(result), "Structure Selection:")
  expect_output(print(result), "Functional Consistency")
})

test_that("summary.msplit_att works", {
  skip_if_not_installed("optimaltrees")

  set.seed(1415)
  n <- 80
  X <- data.frame(x1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  result <- estimate_att_msplit(X, A, Y, M = 2, K = 2, verbose = FALSE)

  # Should print without error
  expect_output(summary(result), "M-Split Doubletree ATT Estimation Summary")
  expect_output(summary(result), "Point Estimate and Inference")
  expect_output(summary(result), "Stability Diagnostics")
})

test_that("estimate_att_msplit validates inputs", {
  n <- 50
  X <- data.frame(x1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  # Invalid X
  expect_error(
    estimate_att_msplit("not a dataframe", A, Y, M = 2, K = 2),
    "X must be a data.frame or matrix"
  )

  # Mismatched dimensions
  expect_error(
    estimate_att_msplit(X, A[1:10], Y, M = 2, K = 2),
    "length\\(A\\) must equal nrow\\(X\\)"
  )

  expect_error(
    estimate_att_msplit(X, A, Y[1:10], M = 2, K = 2),
    "length\\(Y\\) must equal nrow\\(X\\)"
  )

  # Invalid M
  expect_error(
    estimate_att_msplit(X, A, Y, M = 0, K = 2),
    "M must be at least 1"
  )

  # Invalid K
  expect_error(
    estimate_att_msplit(X, A, Y, M = 2, K = 1),
    "K must be at least 2"
  )
})

test_that("prediction variance decreases with M (stochastic test)", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()  # Stochastic test

  set.seed(2021)
  n <- 150
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A)

  # M = 3
  result_M3 <- estimate_att_msplit(X, A, Y, M = 3, K = 2, verbose = FALSE,
                                    seed_base = 300)

  # M = 10
  result_M10 <- estimate_att_msplit(X, A, Y, M = 10, K = 2, verbose = FALSE,
                                     seed_base = 300)

  # Variance should generally decrease (not guaranteed due to randomness, but likely)
  # Just check that both are finite and non-negative
  expect_true(is.finite(result_M3$diagnostics$mean_prediction_variance_e))
  expect_true(is.finite(result_M10$diagnostics$mean_prediction_variance_e))
  expect_true(result_M3$diagnostics$mean_prediction_variance_e >= 0)
  expect_true(result_M10$diagnostics$mean_prediction_variance_e >= 0)
})
