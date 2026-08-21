# Tests for the `propensity_loss` argument on estimate_att() (the full-sample
# Part II estimator) and its threading into predict_nuisances_fold().
# See quality_reports/specs/2026-08-20_propensity-loss-choice.md (decision:
# keep "log_loss" default; "squared_error" is a supported sensitivity option)
# and the Oracle review that closed this file's design (2026-08-21 session).

make_binary_dgp <- function(n, seed = 20260821) {
  set.seed(seed)
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                   X3 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(-0.5 + 1.2 * X$X1 + 0.8 * X$X2))
  Y <- rbinom(n, 1, 0.2 + 0.15 * X$X1 + 0.15 * X$X3 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

# -- A. Argument validation (no fitting; fast) --------------------------------

test_that("propensity_loss defaults to log_loss and is a no-op vs. the explicit default", {
  d <- make_binary_dgp(150)
  fit_default <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L, verbose = FALSE)
  fit_explicit <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                               propensity_loss = "log_loss", verbose = FALSE)

  expect_identical(eval(formals(estimate_att)$propensity_loss),
                    c("log_loss", "squared_error"))
  expect_equal(fit_default$theta, fit_explicit$theta)
  expect_equal(fit_default$sigma, fit_explicit$sigma)
  expect_identical(fit_default$n_leaves_e, fit_explicit$n_leaves_e)
})

test_that("propensity_loss rejects unsupported values", {
  d <- make_binary_dgp(150)
  expect_error(
    estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                 propensity_loss = "misclassification", verbose = FALSE),
    "should be one of"
  )
  expect_error(
    estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                 propensity_loss = NA_character_, verbose = FALSE)
  )
})

# -- B. Threading and non-interference ----------------------------------------

test_that("propensity_loss reaches the fitted e_model; the m0/mu tree is unaffected", {
  d <- make_binary_dgp(200)

  for (outcome_type in c("binary", "continuous")) {
    Y <- if (outcome_type == "continuous") d$Y + rnorm(length(d$Y), sd = 0.1) else d$Y
    expected_m0_loss <- if (outcome_type == "continuous") "squared_error" else "log_loss"

    for (ploss in c("log_loss", "squared_error")) {
      fit <- estimate_att(d$X, d$A, Y, leaf_budget = 4L,
                          outcome_type = outcome_type,
                          propensity_loss = ploss, verbose = FALSE)
      expect_identical(fit$nuisance_fits$e_model@loss_function, ploss)
      expect_identical(fit$nuisance_fits$m0_model@loss_function, expected_m0_loss)
    }
  }
})

test_that("estimate_att() with squared_error propensity returns the full documented structure", {
  d <- make_binary_dgp(200)
  fit <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                      propensity_loss = "squared_error", verbose = FALSE)

  expect_true(is.numeric(fit$theta) && length(fit$theta) == 1)
  expect_true(is.numeric(fit$sigma) && fit$sigma > 0)
  expect_length(fit$ci_95, 2)
  expect_identical(fit$ci_95, fit$ci_95_wald)
  expect_length(fit$score_values, 200)
  expect_true(is.numeric(fit$n_leaves_e) && fit$n_leaves_e >= 1)
  expect_true(is.logical(fit$certified_e))
})

test_that("squared_error propensities are non-degenerate (canary for a defeatable predict() guard)", {
  # Regression test locking in the Oracle-reviewed fix to predict_nuisances_fold():
  # before adding type = "response", a misrouted model could silently pass this
  # branch's format check while returning hard {0, 1} class predictions instead
  # of leaf-mean propensities. This test would have failed on that defect.
  d <- make_binary_dgp(200)
  fit <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                      propensity_loss = "squared_error", verbose = FALSE)
  e_hat <- fit$nuisance_fits$propensity

  expect_true(all(is.finite(e_hat)))
  expect_true(all(e_hat >= 0.01 & e_hat <= 0.99))
  expect_false(all(e_hat %in% c(0.01, 0.99)))
  expect_lte(length(unique(e_hat)), fit$n_leaves_e)
})

test_that("log_loss and squared_error agree within MC noise on a well-separated DGP", {
  d <- make_binary_dgp(400)
  fit_ll <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                         propensity_loss = "log_loss", verbose = FALSE)
  fit_se <- estimate_att(d$X, d$A, d$Y, leaf_budget = 4L,
                         propensity_loss = "squared_error", verbose = FALSE)
  expect_equal(fit_se$theta, fit_ll$theta, tolerance = 0.05)
})

# -- C. predict_nuisances_fold() in isolation ---------------------------------

test_that("predict_nuisances_fold() with a squared_error e_model returns a valid propensity vector", {
  d <- make_binary_dgp(150)
  e_fit <- optimaltrees::bisect_lambda_to_budget(d$X, d$A, leaf_budget = 4L,
                                                 loss_function = "squared_error")
  # outcome_type = "continuous" so m0's predict() path (plain numeric, no
  # type =) is compatible with reusing the same squared_error model for
  # m0_model -- this test is only exercising the e-branch.
  models <- list(e_model = e_fit$fit, m0_model = e_fit$fit, outcome_type = "continuous")
  pred <- predict_nuisances_fold(models, d$X, fold_rows = seq_len(150))

  expect_true(is.numeric(pred$e))
  expect_length(pred$e, 150)
  expect_true(all(pred$e >= 0 & pred$e <= 1))
})

test_that("predict_nuisances_fold() with a log_loss e_model matches predict(type = 'prob')[, 2]", {
  d <- make_binary_dgp(150)
  e_fit <- optimaltrees::bisect_lambda_to_budget(d$X, d$A, leaf_budget = 4L,
                                                 loss_function = "log_loss")
  models <- list(e_model = e_fit$fit, m0_model = e_fit$fit, outcome_type = "binary")
  pred <- predict_nuisances_fold(models, d$X, fold_rows = seq_len(150))

  expect_equal(pred$e, predict(e_fit$fit, d$X, type = "prob")[, 2])
})

test_that("predict_nuisances_fold() errors rather than silently misreading a misrouted model", {
  # TDD anchor for the Oracle-flagged Required #1 fix: without type = "response"
  # in the squared_error branch, a log_loss model whose @loss_function slot is
  # (incorrectly) "squared_error" would silently return its {0, 1} CLASS
  # predictions as propensities, passing the format check. With the fix, the
  # log_loss model's predict(type = "response") still returns its 2-column
  # matrix, so the length check fires -- loud, not silent.
  d <- make_binary_dgp(150)
  e_fit_ll <- optimaltrees::bisect_lambda_to_budget(d$X, d$A, leaf_budget = 4L,
                                                    loss_function = "log_loss")
  misrouted <- e_fit_ll$fit
  misrouted@loss_function <- "squared_error"
  models <- list(e_model = misrouted, m0_model = misrouted, outcome_type = "binary")

  expect_error(
    predict_nuisances_fold(models, d$X, fold_rows = seq_len(150)),
    "unexpected format"
  )
})

test_that("KNOWN GAP: the log_loss branch's format guard does not catch a squared_error model mis-flagged as log_loss", {
  # Discovered while testing the Required #1 fix, NOT introduced by it, and
  # NOT fixed here (out of scope for the propensity_loss diff; this is a
  # property of optimaltrees::predict.optimaltrees_model()'s dispatch, which
  # trusts the model's OWN @loss_function slot with no cross-check against
  # the tree's actual fitted content). Symmetric to Required #1's scenario,
  # but the log_loss branch's `!is.matrix(pe) || ncol(pe) != 2` guard is a
  # SHAPE check only: optimaltrees::get_probabilities_from_tree() still
  # returns a valid-shaped 2-column matrix for a squared_error-fitted tree
  # (verified empirically), it just silently thresholds the regression leaf
  # means to hard {0, 1} "probabilities" -- exactly the degenerate-propensity
  # failure mode the leaf-size-subgroup-enforcement feature exists to prevent,
  # reappearing through a different door. Like Required #1's scenario, this
  # requires manually desyncing @loss_function from the model's actual fit;
  # it is not reachable through any shipped call path today. Flagged as a
  # follow-up for optimaltrees, not fixed in this diff.
  d <- make_binary_dgp(150)
  e_fit_se <- optimaltrees::bisect_lambda_to_budget(d$X, d$A, leaf_budget = 4L,
                                                    loss_function = "squared_error")
  misrouted <- e_fit_se$fit
  misrouted@loss_function <- "log_loss"
  models <- list(e_model = misrouted, m0_model = misrouted, outcome_type = "continuous")

  pred <- predict_nuisances_fold(models, d$X, fold_rows = seq_len(150))
  expect_true(all(pred$e %in% c(0, 1)))
})

test_that("predict_nuisances_fold() with empty fold_rows returns empty vectors regardless of e_loss", {
  d <- make_binary_dgp(150)
  e_fit <- optimaltrees::bisect_lambda_to_budget(d$X, d$A, leaf_budget = 4L,
                                                 loss_function = "squared_error")
  models <- list(e_model = e_fit$fit, m0_model = e_fit$fit, outcome_type = "binary")
  pred <- predict_nuisances_fold(models, d$X, fold_rows = integer(0))

  expect_identical(pred, list(e = numeric(0), m0 = numeric(0)))
})

# -- D. Scope boundary: propensity_loss is full-sample estimate_att() ONLY ----

test_that("propensity_loss does not leak into estimate_att_crossfit()", {
  d <- make_binary_dgp(150)
  expect_false("propensity_loss" %in% names(formals(estimate_att_crossfit)))

  fit <- estimate_att_crossfit(d$X, d$A, d$Y, K = 3, verbose = FALSE)
  # nuisance_fits is a length-K list of per-fold fits with $propensity/
  # $outcome_control tacked on afterward (numeric vectors, not fold fits) --
  # index the K fold elements explicitly rather than iterating the whole list.
  for (k in seq_len(3)) {
    expect_identical(fit$nuisance_fits[[k]]$e_model@loss_function, "log_loss")
  }
})
