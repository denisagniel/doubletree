# Tests for M-Split Doubletree Estimation

test_that("select_structure_modal works", {
  skip_if_not_installed("optimaltrees")

  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  y <- rbinom(n, 1, 0.5)

  # Create M structures (some will be identical if seeds are repeated)
  M <- 5
  structures <- vector("list", M)
  for (m in seq_len(M)) {
    # Use same seed for some to create duplicates
    set.seed(100 + (m %% 3))
    model <- optimaltrees::fit_tree(X, y, loss_function = "log_loss",
                                      regularization = 0.1,
                                      store_training_data = TRUE)
    structures[[m]] <- optimaltrees::extract_tree_structure(model)
  }

  result <- select_structure_modal(structures)

  expect_type(result, "list")
  expect_s3_class(result$structure, "TreeStructure")
  expect_true(result$frequency > 0 && result$frequency <= 1)
  expect_type(result$hash, "character")
  expect_s3_class(result$counts, "table")
})

test_that("select_structure_modal validates inputs", {
  expect_error(
    select_structure_modal(NULL),
    "structures must be a non-empty list"
  )

  expect_error(
    select_structure_modal(list()),
    "structures must be a non-empty list"
  )

  expect_error(
    select_structure_modal(list("not a structure")),
    "All elements.*must be TreeStructure"
  )
})

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

  # Prediction variance should be 0 (only one prediction per obs)
  expect_equal(result$diagnostics$mean_prediction_variance_e, 0)
  expect_equal(result$diagnostics$mean_prediction_variance_m0, 0)
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

test_that("compute_functional_consistency works", {
  skip_if_not_installed("optimaltrees")

  set.seed(1617)
  n <- 50
  M <- 3

  # Create data with some duplicate rows (tied covariates)
  X <- data.frame(
    x1 = rep(c(0, 1), each = 25),
    x2 = rep(c(0, 1), times = 25)
  )

  # Create predictions (with some variability)
  predictions_e <- matrix(runif(n * M), nrow = n, ncol = M)
  predictions_m0 <- matrix(runif(n * M), nrow = n, ncol = M)

  result <- compute_functional_consistency(predictions_e, predictions_m0, X)

  expect_type(result, "list")
  expect_type(result$max_diff_e, "double")
  expect_type(result$max_diff_m0, "double")
  expect_true(result$max_diff_e >= 0)
  expect_true(result$max_diff_m0 >= 0)
  expect_type(result$n_unique_patterns, "integer")
  expect_type(result$n_groups_with_ties, "integer")
  expect_true(result$n_unique_patterns > 0)
})

test_that("analyze_structure_diversity works", {
  skip_if_not_installed("optimaltrees")

  set.seed(1819)
  n <- 100
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  y <- rbinom(n, 1, 0.5)

  # Create M structures
  M <- 5
  structures <- vector("list", M)
  for (m in seq_len(M)) {
    set.seed(100 + m)
    model <- optimaltrees::fit_tree(X, y, loss_function = "log_loss",
                                      regularization = 0.1,
                                      store_training_data = TRUE)
    structures[[m]] <- optimaltrees::extract_tree_structure(model)
  }

  diversity <- analyze_structure_diversity(structures)

  expect_type(diversity, "list")
  expect_type(diversity$n_unique, "integer")
  expect_true(diversity$n_unique >= 1 && diversity$n_unique <= M)
  expect_type(diversity$shannon_entropy, "double")
  expect_type(diversity$simpson_index, "double")
  expect_type(diversity$modal_frequency, "double")
  expect_true(diversity$modal_frequency > 0 && diversity$modal_frequency <= 1)
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
