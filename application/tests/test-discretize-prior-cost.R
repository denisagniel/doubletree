## ============================================================================
## application/tests/test-discretize-prior-cost.R
##
## discretize_prior_cost(): quartile assignment, the zero-inflation case that
## still works, and the collapsed-boundary abort.
##
## Every expected value is HAND-CALCULATED and the calculation is stated. R's
## type-7 quantile on 1:100 at p = 0.25 is 1 + 0.25 * 99 = 25.75, and the bins
## are (-Inf, c1], (c1, c2], ... because findInterval(left.open = TRUE) puts a
## value exactly AT a cutpoint in the lower bin.
## ============================================================================

test_that("quartiles on 1:100 split 25/25/25/25 at the hand-computed cutpoints", {
  out <- discretize_prior_cost(as.numeric(1:100))

  ## Type-7 quantiles of 1:100: p25 = 1 + 0.25*99 = 25.75, p50 = 50.5,
  ## p75 = 1 + 0.75*99 = 75.25.
  expect_equal(unname(attr(out, "cutpoints")), c(25.75, 50.5, 75.25))
  expect_equal(attr(out, "probs"), c(0.25, 0.5, 0.75))

  expect_named(out, c("prior_cost_q2", "prior_cost_q3", "prior_cost_q4"))
  expect_equal(nrow(out), 100L)

  ## Q1 = (-Inf, 25.75] = 1..25 -> 25 patients, so all three dummies are 0 there.
  ## Q2 = (25.75, 50.5]  = 26..50 -> 25
  ## Q3 = (50.5, 75.25]  = 51..75 -> 25
  ## Q4 = (75.25, Inf)   = 76..100 -> 25
  expect_equal(sum(out$prior_cost_q2), 25L)
  expect_equal(sum(out$prior_cost_q3), 25L)
  expect_equal(sum(out$prior_cost_q4), 25L)

  ## Q1's count is implied rather than stored: all three dummies 0.
  expect_equal(sum(rowSums(out) == 0L), 25L)

  ## Reference coding: exactly one dummy is 1 outside Q1, never two.
  expect_true(all(rowSums(out) %in% c(0L, 1L)))
})

test_that("the boundary value itself lands in the LOWER bin", {
  ## 8 sorted values 0,10,15,20,30,40,50,60. Type-7 quantiles:
  ##   p25: h = 1 + 0.25*7 = 2.75 -> 10 + 0.75*(15 - 10) = 13.75
  ##   p50: h = 1 + 0.50*7 = 4.5  -> 20 + 0.50*(30 - 20) = 25
  ##   p75: h = 1 + 0.75*7 = 6.25 -> 40 + 0.25*(50 - 40) = 42.5
  x <- c(0, 10, 15, 20, 30, 40, 50, 60)
  out <- discretize_prior_cost(x)
  cuts <- unname(attr(out, "cutpoints"))
  expect_equal(cuts, c(13.75, 25, 42.5))

  ## Hand-assigned against bins (-Inf, 13.75], (13.75, 25], (25, 42.5],
  ## (42.5, Inf): 0->Q1, 10->Q1, 15->Q2, 20->Q2, 30->Q3, 40->Q3, 50->Q4, 60->Q4.
  expect_equal(out$prior_cost_q2, c(0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L))
  expect_equal(out$prior_cost_q3, c(0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L))
  expect_equal(out$prior_cost_q4, c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L))

  ## The convention itself: a value exactly AT a cutpoint goes to the LOWER bin.
  ## 25 is the p50 cutpoint above, and it is not in x, so state it directly.
  expect_equal(1L + findInterval(25, cuts, left.open = TRUE), 2L)   # Q2, not Q3
  expect_equal(1L + findInterval(13.75, cuts, left.open = TRUE), 1L) # Q1, not Q2
})

test_that("zero-inflation BELOW 25% still discretizes, with zeros in Q1", {
  ## 20 zeros out of 100 = 20% at zero, below the 25th percentile, so p25 > 0 and
  ## no boundary collapses. Values 1..80 fill the rest.
  x <- c(rep(0, 20), as.numeric(1:80))
  out <- discretize_prior_cost(x)
  cuts <- unname(attr(out, "cutpoints"))

  expect_true(all(diff(cuts) > 0))
  expect_true(cuts[1] > 0)
  ## Every zero is in Q1, i.e. all dummies 0.
  expect_true(all(rowSums(out[seq_len(20), ]) == 0L))
  ## And no dummy is all-zero -- that is the property the abort protects.
  expect_true(all(vapply(out, function(col) sum(col) > 0L, logical(1))))
})

test_that("zero-inflation ABOVE 50% collapses the boundary and aborts", {
  ## 60 zeros out of 100 = 60% at zero, so p25 = p50 = 0 and the Q1/Q2 boundary
  ## collapses. "Q2" would be empty and prior_cost_q2 all-zero.
  x <- c(rep(0, 60), as.numeric(1:40))

  expect_error(discretize_prior_cost(x), "boundaries collapsed")
  expect_error(discretize_prior_cost(x), "60% of the sample")
})

test_that("an all-zero prior_cost aborts rather than returning three dead dummies", {
  expect_error(discretize_prior_cost(rep(0, 50)), "boundaries collapsed")
})

test_that("ties at a cutpoint that empty a bin abort too", {
  ## Strictly increasing cutpoints are NOT sufficient. c(1, 1, 2, 3): type-7
  ## p25 = 1 + 0.25*3 = 1, p50 = 1.5, p75 = 2.25 -- strictly increasing -- but
  ## both 1s land in (-Inf, 1] = Q1, so Q2 = (1, 1.5] is EMPTY.
  expect_error(discretize_prior_cost(c(1, 1, 2, 3)), "Empty bin")
})

test_that("NA is refused rather than binned", {
  expect_error(
    discretize_prior_cost(c(1, 2, NA, 4, 5, 6, 7, 8)),
    "NA"
  )
})

test_that("negative amounts are refused (claim reversals are a decision)", {
  expect_error(discretize_prior_cost(c(-5, 1, 2, 3, 4, 5)), "negative")
})

test_that("input validation covers type, emptiness, and probs", {
  expect_error(discretize_prior_cost("100"), "must be numeric")
  expect_error(discretize_prior_cost(numeric(0)), "empty")
  expect_error(discretize_prior_cost(as.numeric(1:100), probs = c(0.5, 0.25)),
               "strictly increasing")
  expect_error(discretize_prior_cost(as.numeric(1:100), probs = c(0, 0.5)),
               "strictly inside")
  expect_error(discretize_prior_cost(as.numeric(1:100), probs = numeric(0)),
               "strictly inside")
})

test_that("a non-quartile probs vector gives (k - 1) dummies", {
  ## Tertiles: 2 interior cutpoints -> 2 dummies, named q2 and q3.
  out <- discretize_prior_cost(as.numeric(1:99), probs = c(1 / 3, 2 / 3))

  expect_named(out, c("prior_cost_q2", "prior_cost_q3"))
  expect_length(attr(out, "cutpoints"), 2L)
  ## Type-7 on 1:99: p33.3 = 1 + (1/3)*98 = 33.667, p66.7 = 1 + (2/3)*98 = 66.333
  expect_equal(unname(attr(out, "cutpoints")), c(33 + 2 / 3, 66 + 1 / 3))
})
