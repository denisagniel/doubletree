# estimate_att_twostage() -- Instantiation 2 (continuum sparsity, two-stage
# off-grid threshold refinement). Requires optimaltrees >= the 2026-09-04
# fit_twostage()/refine_tree() revision (commit 1cecfb4).
# Requires optimaltrees to be installed (doubletree Imports optimaltrees).

test_that("estimate_att_twostage returns the documented list structure and runs", {
  skip_if_not_installed("optimaltrees")
  set.seed(42)
  n <- 600
  # Continuum sparsity: e and m0 each depend on ONE continuous covariate
  # through a single off-grid threshold, matching the roxygen example.
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  A <- rbinom(n, 1, ifelse(X$X1 > 0.37, 0.65, 0.35))
  Y <- rbinom(n, 1, 0.25 + 0.30 * (X$X2 > 0.62) + 0.15 * A)

  fit <- estimate_att_twostage(X, A, Y, leaf_budget = 2L, verbose = FALSE)

  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_true(is.finite(fit$theta))

  expect_type(fit$sigma, "double")
  expect_true(fit$sigma > 0)

  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)
  expect_equal(fit$ci_95, fit$ci_95_wald)

  expect_length(fit$score_values, n)
  expect_equal(fit$n, n)
  expect_equal(fit$leaf_budget, 2L)

  # Per-nuisance diagnostics, not shared, unlike leaf_budget/m_ladder/M_n.
  expect_type(fit$certified_e, "logical")
  expect_type(fit$certified_m0, "logical")
  expect_type(fit$certified_full_class_e, "logical")
  expect_type(fit$certified_full_class_m0, "logical")
  expect_type(fit$stop_reason_e, "character")
  expect_type(fit$stop_reason_m0, "character")
  expect_true(fit$n_leaves_e >= 1 && fit$n_leaves_e <= 2L)
  expect_true(fit$n_leaves_m0 >= 1 && fit$n_leaves_m0 <= 2L)

  # lambda_n/m_n differ by nuisance (variance/sample-size-scaled defaults);
  # leaf_budget/depth_budget/m_ladder/M_n are shared.
  expect_true(is.numeric(fit$lambda_n_e) && fit$lambda_n_e > 0)
  expect_true(is.numeric(fit$lambda_n_m0) && fit$lambda_n_m0 > 0)
  expect_true(is.numeric(fit$m_n_e) && fit$m_n_e >= 1)
  expect_true(is.numeric(fit$m_n_m0) && fit$m_n_m0 >= 1)

  expect_true(is.list(fit$nuisance_fits))
  expect_true(all(c("e_fit", "m0_fit", "propensity", "outcome_control") %in%
    names(fit$nuisance_fits)))
  expect_s3_class(fit$nuisance_fits$e_fit, "optimaltrees_twostage_fit")
  expect_s3_class(fit$nuisance_fits$m0_fit, "optimaltrees_twostage_fit")
  expect_length(fit$nuisance_fits$propensity, n)
  expect_length(fit$nuisance_fits$outcome_control, n)
  expect_true(all(fit$nuisance_fits$propensity >= 0 &
    fit$nuisance_fits$propensity <= 1))
})

test_that("estimate_att_twostage's m_n/lambda_n defaults satisfy the documented rate conditions", {
  skip_if_not_installed("optimaltrees")
  set.seed(11)
  n <- 800
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.4 + 0.1 * A)

  fit <- estimate_att_twostage(X, A, Y, leaf_budget = 2L, verbose = FALSE)

  # (R5): m_n / sqrt(n) should exceed 1 well before n is very large -- a loose
  # but real check that the default is not the trivial estimate_att() m_n = 1.
  expect_true(fit$m_n_e > 1L)
  expect_true(fit$m_n_m0 > 1L)
  expect_true(fit$m_n_e / sqrt(n) > 1)

  # lambda_n is scaled to each response's own variance, so a fit on A (a
  # 0/1 variable, Var <= 0.25) and a fit on continuous Y0 should generally
  # differ, not share estimate_att()'s single log(n)/n value.
  expect_false(isTRUE(all.equal(fit$lambda_n_e, fit$lambda_n_m0)))
})

test_that("estimate_att_twostage works with continuous outcomes", {
  skip_if_not_installed("optimaltrees")
  set.seed(99)
  n <- 500
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  A <- rbinom(n, 1, 0.5)
  tau <- 0.4
  Y <- rnorm(n, mean = 1 + 0.5 * (X$X1 > 0.4) + tau * A, sd = 0.8)

  fit <- estimate_att_twostage(X, A, Y, leaf_budget = 2L, verbose = FALSE)

  expect_true(is.finite(fit$theta))
  expect_true(fit$sigma > 0)
  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta && fit$theta < fit$ci_95[2])
})

test_that("estimate_att_twostage requires leaf_budget with no default", {
  skip_if_not_installed("optimaltrees")
  n <- 200
  X <- data.frame(X1 = runif(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(estimate_att_twostage(X, A, Y), "required")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = NULL), "required")
})

test_that("estimate_att_twostage rejects purely binary covariates, pointing to estimate_att", {
  skip_if_not_installed("optimaltrees")
  n <- 200
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(
    estimate_att_twostage(X, A, Y, leaf_budget = 2L),
    "continuous.*estimate_att\\(\\)|estimate_att\\(\\)"
  )
})

test_that("estimate_att_twostage validates leaf_budget, depth_budget, m_ladder, M_n, m_n, lambda_n", {
  skip_if_not_installed("optimaltrees")
  n <- 200
  X <- data.frame(X1 = runif(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2.5), "positive integer")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = -1L), "positive integer")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, depth_budget = "partial"),
    "depth_budget")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, m_ladder = c(16)),
    "at least two")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, m_ladder = c(-1, 16)),
    "at least two")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, M_n = -1), "M_n")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, m_n = 0L), "positive integer")
  expect_error(estimate_att_twostage(X, A, Y, leaf_budget = 2L, lambda_n = -0.1), "positive")
})

test_that("estimate_att_twostage rejects reserved ... arguments owned by named parameters", {
  skip_if_not_installed("optimaltrees")
  n <- 200
  X <- data.frame(X1 = runif(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  for (nm in c("discretize_bins", "loss_function", "max_depth", "min_leaf_n",
               "group", "group_value", "m_n_group")) {
    args <- list(X = X, A = A, Y = Y, leaf_budget = 2L)
    args[[nm]] <- 1L
    expect_error(do.call(estimate_att_twostage, args), "reserves",
      info = paste("reserved name:", nm))
  }
})

test_that("estimate_att_twostage errors informatively on degenerate treatment samples", {
  skip_if_not_installed("optimaltrees")
  n <- 200
  X <- data.frame(X1 = runif(n))
  Y <- rbinom(n, 1, 0.5)

  expect_error(
    estimate_att_twostage(X, A = rep(0L, n), Y, leaf_budget = 2L),
    "[Nn]o treated"
  )
  expect_error(
    estimate_att_twostage(X, A = rep(1L, n), Y, leaf_budget = 2L),
    "[Nn]o control"
  )
})

test_that("estimate_att_twostage raises a classed, catchable condition when a fit produces no model", {
  skip_if_not_installed("optimaltrees")
  set.seed(7)
  n <- 300
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n, mean = 1 + 0.3 * (X$X1 > 0.4) + 0.2 * A)

  caught <- FALSE
  withCallingHandlers(
    tryCatch(
      estimate_att_twostage(X, A, Y, leaf_budget = 2L,
        time_budget = 1e-4, fit_time_limit = 1e-4, verbose = FALSE),
      doubletree_twostage_no_model = function(e) {
        caught <<- TRUE
        expect_true(!is.null(e$stop_reason))
        expect_true("ladder_topologies" %in% names(e))
      }
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_true(caught)
})

test_that("estimate_att_twostage warns (not errors) on off-unit-scale covariates", {
  skip_if_not_installed("optimaltrees")
  set.seed(5)
  n <- 300
  X <- data.frame(X1 = runif(n, 0, 1000), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  expect_warning(
    estimate_att_twostage(X, A, Y, leaf_budget = 2L, verbose = FALSE),
    "unit scale"
  )
})
