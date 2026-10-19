## ============================================================================
## application/tests/test-winsorize-outcome.R
##
## winsorize_upper(): hand-computed cap, counts, and the no-op / invalid cases.
## Type-7 quantile of 1:100 at p = 0.99 is 1 + 0.99 * 99 = 99.01.
## ============================================================================

test_that("winsorize_upper caps at the hand-computed type-7 quantile", {
  y <- as.numeric(1:100)
  out <- winsorize_upper(y, 0.99)
  expect_equal(attr(out, "cap"), 99.01)
  ## Only y = 100 exceeds 99.01.
  expect_equal(attr(out, "n_capped"), 1L)
  expect_equal(as.numeric(out)[100], 99.01)
  expect_equal(as.numeric(out)[1:99], y[1:99])
  expect_equal(attr(out, "q"), 0.99)
})

test_that("q = 1 is the identity with an infinite cap", {
  y <- c(0, 5, 1e6)
  out <- winsorize_upper(y, 1)
  expect_equal(as.numeric(out), y)
  expect_equal(attr(out, "cap"), Inf)
  expect_equal(attr(out, "n_capped"), 0L)
})

test_that("winsorize_upper never lowers a value below the cap or changes zeros", {
  y <- c(0, 0, 1, 2, 3, 1000)
  out <- as.numeric(winsorize_upper(y, 0.8))
  expect_true(all(out <= y))
  expect_equal(out[y <= attr(winsorize_upper(y, 0.8), "cap")],
               y[y <= attr(winsorize_upper(y, 0.8), "cap")])
})

test_that("winsorize_upper rejects NA and bad quantiles", {
  expect_error(winsorize_upper(c(1, NA), 0.9), "no NA")
  expect_error(winsorize_upper(1:3, 0), "in \\(0, 1\\]")
  expect_error(winsorize_upper(1:3, 1.5), "in \\(0, 1\\]")
  expect_error(winsorize_upper("a", 0.9), "numeric")
})

test_that("the winsorization level is an open decision that moves the fingerprint", {
  expect_equal(OPEN_DECISIONS$outcome_winsor_quantile$status, "open_decision")
  expect_false(OPEN_DECISIONS$outcome_winsor_quantile$blocking_final)
  f1 <- config_fingerprint()
  old <- OPEN_DECISIONS
  withr::defer(assign("OPEN_DECISIONS", old, envir = globalenv()))
  od <- OPEN_DECISIONS
  od$outcome_winsor_quantile$value <- 0.995
  assign("OPEN_DECISIONS", od, envir = globalenv())
  expect_false(identical(config_fingerprint(), f1))
})
