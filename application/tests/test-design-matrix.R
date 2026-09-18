## ============================================================================
## application/tests/test-design-matrix.R
##
## assemble_design_matrix() and assert_binary_design_matrix().
##
## The NEGATIVE test in this file is the point of the file. See its own comment
## and helpers/design_matrix.R's header: `check_att_data()` (R/utils.R) does not
## examine X's coding at all, and `estimate_att_crossfit()` adds no further
## check, so a 1/2-coded column reaches optimaltrees and returns a number.
## `estimate_att()` does reject it -- but only at estimation time, after the
## analytic data has been assembled. This guard fires at assembly time.
## ============================================================================

test_that("assemble_design_matrix() produces ID + 15 columns in declared order", {
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)
  out <- assemble_design_matrix(cov, dummies)

  expect_equal(names(out), c("ID", design_matrix_columns()))
  expect_equal(nrow(out), 64L)
  expect_equal(out$ID, cov$ID)
  ## The realized cutpoints travel with the matrix they define.
  expect_equal(attr(out, "cutpoints"), attr(dummies, "cutpoints"))
})

test_that("a valid 15-column binary matrix passes the guard", {
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)
  X <- assemble_design_matrix(cov, dummies)

  ## With and without the ID column: ID is a key, not a covariate, and is ignored.
  expect_true(isTRUE(assert_binary_design_matrix(X)))
  expect_true(isTRUE(assert_binary_design_matrix(
    X[, design_matrix_columns(), drop = FALSE]
  )))
})

test_that("assert_binary_design_matrix() catches a 1/2-coded column that the package would not", {
  ## THE CASE THIS GUARD EXISTS FOR. A SAS `_YN` flag arriving as 1/2 rather than
  ## 0/1 has two distinct values, so a naive "is it binary" check passes it. It
  ## then:
  ##   * passes check_att_data(), which validates NAs, A's coding, and Y against
  ##     outcome_type -- and never looks at X's coding;
  ##   * is ACCEPTED by estimate_att_crossfit(), whose only X validation is
  ##     check_att_data(), so it runs to completion and returns a number computed
  ##     on a grid the analyst never declared;
  ##   * is rejected by estimate_att(), but only after 03's multi-hour cost
  ##     aggregation and 06's write have already happened.
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)
  X <- assemble_design_matrix(cov, dummies)[, design_matrix_columns(), drop = FALSE]

  X_bad <- X
  X_bad$HOUSING_YN <- X_bad$HOUSING_YN + 1L          # now 1/2, not 0/1

  expect_error(assert_binary_design_matrix(X_bad), "not coded")
  expect_error(assert_binary_design_matrix(X_bad), "HOUSING_YN")

  ## And confirm the premise rather than asserting it: check_att_data() is not
  ## exported, so reach it through the estimator that performs no further check.
  ## The call is expected to SUCCEED on the mis-coded matrix -- that success is
  ## the failure mode. Kept small (K = 2) so it is cheap.
  skip_if_not_installed("doubletree")
  d <- toy_att_data(n = 120L, seed = 7L)
  d$X$HOUSING_YN <- d$X$HOUSING_YN + 1L
  fit <- doubletree::estimate_att_crossfit(
    X = d$X, A = d$A, Y = d$Y, K = 2L, outcome_type = "continuous", seed = 1L
  )
  expect_true(is.finite(fit$theta))   # a number, on an undeclared grid
})

test_that("assert_binary_design_matrix() catches a Y/N character column", {
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)
  X <- assemble_design_matrix(cov, dummies)[, design_matrix_columns(), drop = FALSE]

  X_bad <- X
  X_bad$FOOD_YN <- ifelse(X_bad$FOOD_YN == 1L, "Y", "N")

  expect_error(assert_binary_design_matrix(X_bad), "not coded")
})

test_that("assert_binary_design_matrix() catches NA, constant, and wrong-shape input", {
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)
  X <- assemble_design_matrix(cov, dummies)[, design_matrix_columns(), drop = FALSE]

  X_na <- X
  X_na$MH_IP_YN[1L] <- NA_integer_
  expect_error(assert_binary_design_matrix(X_na), "NA")

  X_const <- X
  X_const$POVERTY_YN <- 0L
  expect_error(assert_binary_design_matrix(X_const), "constant")

  ## Missing a column entirely.
  expect_error(
    assert_binary_design_matrix(X[, -1L, drop = FALSE]),
    "Missing"
  )
  ## An extra column.
  X_extra <- X
  X_extra$AGE_BINARY <- 1L
  expect_error(assert_binary_design_matrix(X_extra), "Unexpected")

  ## Right columns, WRONG ORDER -- the failure a column-set check would miss.
  X_reordered <- X[, rev(design_matrix_columns()), drop = FALSE]
  expect_error(assert_binary_design_matrix(X_reordered), "WRONG ORDER")
})

test_that("assemble_design_matrix() refuses misaligned or incomplete inputs", {
  cov <- toy_covariates(n = 64L)
  dummies <- discretize_prior_cost(cov$prior_cost)

  expect_error(
    assemble_design_matrix(cov[, -1L, drop = FALSE], dummies),
    "ID"
  )
  expect_error(
    assemble_design_matrix(cov[, setdiff(names(cov), "FOOD_YN"), drop = FALSE], dummies),
    "missing confirmed covariate"
  )
  expect_error(
    assemble_design_matrix(cov, dummies[seq_len(10L), , drop = FALSE]),
    "Row counts differ"
  )
  expect_error(
    assemble_design_matrix(cov, dummies[, c("prior_cost_q2"), drop = FALSE]),
    "prior_cost_q3"
  )
  ## Duplicated patients: the design matrix must be one row per patient.
  cov_dup <- dplyr::bind_rows(cov, cov[1L, ])
  expect_error(
    assemble_design_matrix(cov_dup, discretize_prior_cost(cov_dup$prior_cost)),
    "one row per patient"
  )
})
