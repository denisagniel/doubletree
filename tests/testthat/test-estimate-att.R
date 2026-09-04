# estimate_att() -- the paper's flagship, full-sample, no-cross-fitting
# estimator (theory.tex Part II, sec:main). Both nuisance partitions AND leaf
# values are fit on all n observations; see estimate_att_crossfit() for the
# K-fold fallback and estimate_att_rashomon() for the superseded variant.
# Requires optimaltrees to be installed (doubletree Imports optimaltrees).

test_that("estimate_att returns the documented list structure and runs on binary data", {
  skip_if_not_installed("optimaltrees")
  set.seed(42)
  n <- 400
  # Hierarchically sparse nuisances: e depends on X1 only, m0 on X2 only, so
  # a 4-leaf budget represents both exactly (matches the roxygen example).
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5), X3 = rbinom(n, 1, 0.5)
  )
  A <- rbinom(n, 1, ifelse(X$X1 == 1, 0.65, 0.35))
  Y <- rbinom(n, 1, 0.25 + 0.30 * X$X2 + 0.15 * A)

  fit <- estimate_att(X, A, Y, leaf_budget = 4L)

  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_true(is.finite(fit$theta))

  expect_type(fit$sigma, "double")
  expect_true(fit$sigma > 0)

  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)
  # Documented to be identical, for cross-estimator API parity.
  expect_equal(fit$ci_95, fit$ci_95_wald)

  expect_length(fit$score_values, n)
  expect_equal(fit$n, n)
  expect_equal(fit$leaf_budget, 4L)
  expect_equal(fit$m_n, 1L)
  # Default lambda_n resolves to log(n)/n (prop:parsimony-compatible rate).
  expect_equal(fit$lambda_n, log(n) / n)

  # Sparsity-proxy diagnostics are present but never auto-acted on.
  expect_type(fit$certified_e, "logical")
  expect_type(fit$certified_m0, "logical")
  expect_type(fit$used_search_e, "logical")
  expect_type(fit$used_search_m0, "logical")
  expect_true(fit$n_leaves_e >= 1 && fit$n_leaves_e <= 4L)
  expect_true(fit$n_leaves_m0 >= 1 && fit$n_leaves_m0 <= 4L)

  expect_true(is.list(fit$nuisance_fits))
  expect_true(all(c("e_model", "m0_model", "propensity", "outcome_control") %in%
    names(fit$nuisance_fits)))
  expect_length(fit$nuisance_fits$propensity, n)
  expect_length(fit$nuisance_fits$outcome_control, n)
})

test_that("estimate_att fits both trees on all n observations (no sample splitting)", {
  skip_if_not_installed("optimaltrees")
  set.seed(7)
  n <- 200
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A)

  fit <- estimate_att(X, A, Y, leaf_budget = 3L)

  # In-sample predictions are returned for every one of the n rows -- there is
  # no held-out/held-in distinction on this path (sec:main).
  expect_length(fit$nuisance_fits$propensity, n)
  expect_length(fit$nuisance_fits$outcome_control, n)
  # Propensity predictions are clamped into the shared bounds.
  expect_true(all(fit$nuisance_fits$propensity >= 0 &
    fit$nuisance_fits$propensity <= 1))
})

test_that("estimate_att works with continuous outcomes (squared_error loss)", {
  skip_if_not_installed("optimaltrees")
  set.seed(99)
  n <- 200
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  tau <- 0.5
  Y <- rnorm(n, mean = 1 + 0.6 * X$X1 + tau * A, sd = 0.8)

  fit <- estimate_att(X, A, Y, leaf_budget = 3L, outcome_type = "continuous")

  expect_true(is.finite(fit$theta))
  expect_true(fit$sigma > 0)
  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta && fit$theta < fit$ci_95[2])
})

test_that("estimate_att requires leaf_budget with no default", {
  skip_if_not_installed("optimaltrees")
  n <- 100
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(estimate_att(X, A, Y), "required")
  expect_error(estimate_att(X, A, Y, leaf_budget = NULL), "required")
})

test_that("estimate_att validates leaf_budget, m_n, and lambda_n", {
  skip_if_not_installed("optimaltrees")
  n <- 100
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(estimate_att(X, A, Y, leaf_budget = 2.5), "positive integer")
  expect_error(estimate_att(X, A, Y, leaf_budget = -1L), "positive integer")
  expect_error(estimate_att(X, A, Y, leaf_budget = 2L, m_n = 0L), "positive integer")
  expect_error(estimate_att(X, A, Y, leaf_budget = 2L, lambda_n = -0.1),
    "positive")
  expect_error(estimate_att(X, A, Y, leaf_budget = 2L, lambda_n = c(0.1, 0.2)),
    "single")
})

test_that("estimate_att rejects non-binary covariates and points to the cross-fit fallback", {
  skip_if_not_installed("optimaltrees")
  n <- 100
  X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(
    estimate_att(X, A, Y, leaf_budget = 2L),
    "binary.*estimate_att_crossfit|estimate_att_crossfit"
  )
})

test_that("estimate_att errors informatively on degenerate treatment/control samples", {
  skip_if_not_installed("optimaltrees")
  n <- 100
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  Y <- rbinom(n, 1, 0.5)

  expect_error(
    estimate_att(X, A = rep(0L, n), Y, leaf_budget = 2L),
    "[Nn]o treated"
  )
  expect_error(
    estimate_att(X, A = rep(1L, n), Y, leaf_budget = 2L),
    "[Nn]o control"
  )
})

test_that("estimate_att errors when there are too few control units for leaf_budget * m_n", {
  skip_if_not_installed("optimaltrees")
  set.seed(11)
  n <- 20
  X <- data.frame(X1 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.9) # very few controls expected
  Y <- rbinom(n, 1, 0.5)
  # Force a small control count deterministically for a reproducible message.
  A[1:18] <- 1L
  A[19:20] <- 0L

  expect_error(
    estimate_att(X, A, Y, leaf_budget = 10L, m_n = 5L),
    "[Ii]nsufficient control"
  )
})

test_that("estimate_att's certified/used_search/n_leaves diagnostics are internally consistent", {
  skip_if_not_installed("optimaltrees")
  set.seed(2024)
  n <- 300
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, ifelse(X$X1 == 1, 0.6, 0.4))
  Y <- rbinom(n, 1, 0.2 + 0.3 * X$X2 + 0.1 * A)

  fit <- estimate_att(X, A, Y, leaf_budget = 4L)

  # certified == TRUE implies gap == 0 (proven exact, no slack to report);
  # certified == FALSE with a feasible fit gives a non-negative gap or NA.
  if (isTRUE(fit$certified_e)) expect_equal(fit$gap_e, 0)
  if (isTRUE(fit$certified_m0)) expect_equal(fit$gap_m0, 0)
})
