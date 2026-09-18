## ============================================================================
## application/tests/test-estimate-att-smoke.R
##
## Tier-0 smoke test: estimate_att() and estimate_att_crossfit() called at the
## REAL shape -- 15 binary covariates in the confirmed order,
## outcome_type = "continuous", leaf_budget = config_leaf_budget -- on a
## hand-built known-truth DGP.
##
## WHAT THIS PROVES: the application's design-matrix contract and the package's
## estimator contract fit together, and grid-exact sparsity is SATISFIABLE at
## this covariate shape in principle. In helper-toy.R's toy_att_data(), the
## propensity depends on one binary covariate and the control outcome on another,
## so a 2-leaf tree represents each exactly and a 4-leaf budget is more than
## enough. 15 binary columns at leaf_budget = 4 is therefore not a
## self-contradictory specification.
##
## WHAT THIS DOES NOT PROVE, and cannot: that sparsity holds on the REAL joint
## distribution of the 15 confirmed covariates. That is not checkable from data
## at all (theory.tex ass:sparsity is an identifying assumption), which is exactly
## why 07_estimate_att.R always reports the cross-fit companion and why
## OPEN_DECISIONS$diagnostics is open. A passing test here says the pipe fits
## together; it says nothing about the applied answer.
##
## Not a re-test of the estimators themselves -- tests/testthat/ (the package
## suite) covers those.
## ============================================================================

test_that("estimate_att() runs at the real 15-column shape and returns the documented structure", {
  skip_if_not_installed("doubletree")
  skip_if_not_installed("optimaltrees")

  d <- toy_att_data(n = 400L, tau = 250, sigma = 25, seed = 11L)

  ## Shape assertions first: if these fail, the test below is testing the wrong
  ## thing and its result is uninformative.
  expect_equal(ncol(d$X), 15L)
  expect_equal(names(d$X), design_matrix_columns())
  expect_true(isTRUE(assert_binary_design_matrix(d$X)))
  expect_true(all(d$A %in% c(0L, 1L)))
  expect_type(d$Y, "double")

  fit <- doubletree::estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    leaf_budget = config_leaf_budget,
    outcome_type = "continuous"
  )

  ## Return value is a named LIST (R/estimate_att.R's final expression), not a
  ## tibble. Named elements, verified against that function's own @return block.
  expect_type(fit, "list")
  expect_true(all(c(
    "theta", "sigma", "ci_95", "ci_95_wald", "score_values", "nuisance_fits",
    "n", "leaf_budget", "lambda_n", "m_n",
    "certified_e", "certified_m0", "used_search_e", "used_search_m0",
    "n_leaves_e", "n_leaves_m0", "gap_e", "gap_m0"
  ) %in% names(fit)))

  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$sigma) && fit$sigma > 0)
  expect_length(fit$ci_95, 2L)
  expect_true(fit$ci_95[1] < fit$ci_95[2])
  expect_length(fit$score_values, 400L)

  expect_equal(fit$n, 400L)
  expect_equal(fit$leaf_budget, config_leaf_budget)
  ## lambda_n defaults to log(n)/n -- hand-computed, not read off the fit.
  expect_equal(fit$lambda_n, log(400) / 400)

  ## The realised leaf counts must respect the budget. This is the budget
  ## mechanism working, not a statement about sparsity.
  expect_lte(fit$n_leaves_e, config_leaf_budget)
  expect_lte(fit$n_leaves_m0, config_leaf_budget)
})

test_that("estimate_att() recovers the known constant effect on the toy DGP", {
  skip_if_not_installed("doubletree")
  skip_if_not_installed("optimaltrees")

  ## tau = 250 with residual SD 25 at n = 400: a wide tolerance is deliberate.
  ## This is a mechanics check -- "the estimate lands in the right neighbourhood
  ## rather than being an order of magnitude off or the wrong sign" -- not a
  ## coverage or consistency study, which belongs in simulations/.
  d <- toy_att_data(n = 400L, tau = 250, sigma = 25, seed = 11L)
  fit <- doubletree::estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    leaf_budget = config_leaf_budget, outcome_type = "continuous"
  )

  expect_gt(fit$theta, 150)
  expect_lt(fit$theta, 350)
})

test_that("estimate_att_crossfit() runs at the same shape (07's companion path)", {
  skip_if_not_installed("doubletree")
  skip_if_not_installed("optimaltrees")

  d <- toy_att_data(n = 400L, tau = 250, sigma = 25, seed = 11L)

  ## K = 2 rather than config_crossfit_k (5): this is a shape check, and 2 folds
  ## fit 4 trees instead of 10.
  fit <- doubletree::estimate_att_crossfit(
    X = d$X, A = d$A, Y = d$Y,
    K = 2L, outcome_type = "continuous", seed = config_seed
  )

  expect_type(fit, "list")
  expect_true(all(c("theta", "sigma", "ci_95", "score_values") %in% names(fit)))
  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$sigma) && fit$sigma > 0)
})

test_that("estimate_att() rejects a continuous prior_cost -- the reason it is pre-discretized", {
  skip_if_not_installed("doubletree")

  d <- toy_att_data(n = 200L, seed = 3L)

  ## Replace the three quartile dummies with the raw dollar amount, i.e. exactly
  ## what happens if discretize_prior_cost() is skipped. estimate_att() rejects
  ## it; estimate_att_crossfit() would not (see test-design-matrix.R), which is
  ## why the discretization is a pipeline step rather than a caller's option.
  X_cont <- d$X
  X_cont$prior_cost_q2 <- NULL
  X_cont$prior_cost_q3 <- NULL
  X_cont$prior_cost_q4 <- NULL
  X_cont$prior_cost <- as.numeric(seq_len(nrow(X_cont))) * 100

  expect_error(
    doubletree::estimate_att(
      X = X_cont, A = d$A, Y = d$Y,
      leaf_budget = config_leaf_budget, outcome_type = "continuous"
    ),
    "binary"
  )
})

test_that("leaf_budget is required, with no default (ass:budget)", {
  skip_if_not_installed("doubletree")

  d <- toy_att_data(n = 120L, seed = 5L)
  expect_error(
    doubletree::estimate_att(X = d$X, A = d$A, Y = d$Y, outcome_type = "continuous"),
    "leaf_budget"
  )
})
