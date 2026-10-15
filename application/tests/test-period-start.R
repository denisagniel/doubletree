## ============================================================================
## application/tests/test-period-start.R
##
## msr_yr_period_start(): T0 = first day of the MSR_YR month. Hand-calculated.
## ============================================================================

test_that("T0 is the first day of the MSR_YR month", {
  skip_if_not_installed("smidata")
  skip_if_not("smi_parse_msr_yr" %in% getNamespaceExports("smidata"))
  expect_equal(
    msr_yr_period_start(c("RY2019-06", "RY2020-12", "RY2017-01")),
    as.Date(c("2019-06-01", "2020-12-01", "2017-01-01"))
  )
})

test_that("an unparseable MSR_YR aborts rather than becoming NA", {
  skip_if_not_installed("smidata")
  skip_if_not("smi_parse_msr_yr" %in% getNamespaceExports("smidata"))
  expect_error(msr_yr_period_start("not a period"), "not match")
  expect_error(msr_yr_period_start(NA_character_))
})

test_that("T0 sits at signed gap -1, 0, +1 from an INDEX_DT in mid-month", {
  skip_if_not_installed("smidata")
  skip_if_not("smi_parse_msr_yr" %in% getNamespaceExports("smidata"))
  ## INDEX_DT = 2019-06-15 (month 6). Hand-calculated: the period starting in
  ## May / June / July 2019 gives T0 = 2019-05-01 / 06-01 / 07-01.
  expect_equal(
    msr_yr_period_start(c("RY2019-05", "RY2019-06", "RY2019-07")),
    as.Date(c("2019-05-01", "2019-06-01", "2019-07-01"))
  )
})

test_that("windows built on T0 abut exactly: prior | exposure | outcome", {
  skip_if_not_installed("smidata")
  ## T0 = 2019-06-01. Half-open windows from DOUBLETREE_WINDOWS, hand-calculated:
  ##   prior_cost [2018-06-01, 2019-06-01)   outcome [2020-06-01, 2021-06-01)
  ## and the exposure [T0, T0 + 12mo) = [2019-06-01, 2020-06-01) sits between them.
  w <- smidata::smi_patient_windows(
    tibble::tibble(id = "p1", index_dt = as.Date("2019-06-01")),
    DOUBLETREE_WINDOWS
  )
  prior <- w[w$window == "prior_cost", ]
  outcome <- w[w$window == "outcome", ]
  t0 <- as.Date("2019-06-01")
  exposure_end <- seq(t0, by = "12 months", length.out = 2L)[[2]]
  expect_equal(prior$window_start, as.Date("2018-06-01"))
  expect_equal(prior$window_end, t0)
  expect_equal(exposure_end, outcome$window_start)
  expect_equal(outcome$window_end, as.Date("2021-06-01"))
})
