## ============================================================================
## application/tests/test-estimate-att-smoke.R
##
## Tier-0 smoke test: estimate_att() and estimate_att_crossfit() called at the
## REAL shape -- 15 binary covariates in the confirmed order,
## outcome_type = "continuous" -- on a hand-built known-truth DGP.
##
## LEAF_BUDGET IS DELIBERATELY DECOUPLED FROM config_leaf_budget, AS OF
## 2026-09-28. Benchmarked directly (not guessed) on this same 400-row toy
## fixture before this test suite ever ran at the raised value: estimate_att()
## took 4.8s at leaf_budget=4, 94.5s at 6, and was still running past 5+
## minutes at 15 (config_leaf_budget's confirmed value, raised from 4, PI
## 2026-09-28) when killed -- optimaltrees' GOSDT-style search cost grows
## steeply with leaf_budget, essentially independent of n. Using
## config_leaf_budget here directly would turn this file's job -- proving the
## design-matrix contract and the package's estimator contract fit together,
## fast, with no data and no server (see application/README.md's own
## "165 assertions... runs today" claim) -- into a multi-minute-or-worse local
## test run. smoke_leaf_budget below is a FIXED, small value chosen only to
## be feasible on 400 toy rows; it proves nothing about whether
## config_leaf_budget itself is feasible on the real cohort -- that is
## 91_leaf_budget_feasibility.R's job, on real data, on the server.
##
## WHAT THIS PROVES: the application's design-matrix contract and the package's
## estimator contract fit together, and grid-exact sparsity is SATISFIABLE at
## this covariate shape in principle. In helper-toy.R's toy_att_data(), the
## propensity depends on one binary covariate and the control outcome on another,
## so a 2-leaf tree represents each exactly and any smoke_leaf_budget >= 2 is
## more than enough. 15 binary columns at a small leaf_budget is therefore not
## a self-contradictory specification -- a claim that does not depend on which
## exact small value is used.
##
## WHAT THIS DOES NOT PROVE, and cannot: that sparsity holds on the REAL joint
## distribution of the 15 confirmed covariates, OR that config_leaf_budget is
## computationally feasible at real scale. Neither is checkable here -- the
## first is theory.tex's ass:sparsity, an identifying assumption never
## checkable from data; the second is exactly why
## 91_leaf_budget_feasibility.R exists as a separate, real-data-only script.
## This file's passing tests say the pipe fits together; they say nothing
## about the applied answer OR about feasibility at config_leaf_budget.
##
## Not a re-test of the estimators themselves -- tests/testthat/ (the package
## suite) covers those.
## ============================================================================

## Fixed, small, and independent of config_leaf_budget -- see header. Not
## read from _config.R on purpose: a future config_leaf_budget change must
## not silently change this file's own runtime again.
smoke_leaf_budget <- 4L

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
    leaf_budget = smoke_leaf_budget,
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
  expect_equal(fit$leaf_budget, smoke_leaf_budget)
  ## lambda_n defaults to log(n)/n -- hand-computed, not read off the fit.
  expect_equal(fit$lambda_n, log(400) / 400)

  ## The realised leaf counts must respect the budget. This is the budget
  ## mechanism working, not a statement about sparsity.
  expect_lte(fit$n_leaves_e, smoke_leaf_budget)
  expect_lte(fit$n_leaves_m0, smoke_leaf_budget)
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
    leaf_budget = smoke_leaf_budget, outcome_type = "continuous"
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
      leaf_budget = smoke_leaf_budget, outcome_type = "continuous"
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
