## ============================================================================
## application/tests/test-windows-two-anchors.R
##
## The riskiest boundary of the 2026-10-15 rules: outcome is anchored on
## s (= T0), prior_cost on INDEX_DT, and the two differ whenever INDEX_DT is not
## the 1st. Every date below is worked out by hand.
## ============================================================================

mk <- function(index_dt) {
  idx <- as.Date(index_dt)
  doubletree_patient_windows(tibble::tibble(
    ID = "a", INDEX_DT = idx, T0 = exposure_start(idx)
  ))
}

test_that("INDEX_DT on the 1st: s = INDEX_DT and both anchors coincide", {
  w <- mk("2019-06-01")
  ## s = 2019-06-01. outcome [s+12mo, s+24mo) = [2020-06-01, 2021-06-01);
  ## prior_cost [INDEX_DT-12mo, INDEX_DT) = [2018-06-01, 2019-06-01).
  expect_equal(w$window_start, as.Date(c("2020-06-01", "2018-06-01")))
  expect_equal(w$window_end, as.Date(c("2021-06-01", "2019-06-01")))
})

test_that("mid-month INDEX_DT: outcome on s, prior_cost on INDEX_DT itself", {
  w <- mk("2019-06-15")
  ## s = first of the month on or after 2019-06-15 = 2019-07-01.
  ## outcome [2020-07-01, 2021-07-01); prior_cost [2018-06-15, 2019-06-15).
  ## The days 2019-06-15 .. 2019-06-30 (between INDEX_DT and s) are in neither.
  expect_equal(w$window_start, as.Date(c("2020-07-01", "2018-06-15")))
  expect_equal(w$window_end, as.Date(c("2021-07-01", "2019-06-15")))
  expect_true(all(w$window_end[w$window == "prior_cost"] <= as.Date("2019-06-15")))
  expect_true(w$window_start[w$window == "outcome"] > as.Date("2019-06-30"))
})

test_that("INDEX_DT on the last day of a month still moves s to the next month", {
  w <- mk("2019-12-31")
  ## s = 2020-01-01. outcome [2021-01-01, 2022-01-01); prior_cost [2018-12-31, 2019-12-31).
  expect_equal(w$window_start, as.Date(c("2021-01-01", "2018-12-31")))
  expect_equal(w$window_end, as.Date(c("2022-01-01", "2019-12-31")))
})

test_that("claims at the window edges are in or out as the half-open rule says", {
  skip_if_not_installed("smidata")
  cohort <- tibble::tibble(ID = "a", INDEX_DT = as.Date("2019-06-15"),
                           T0 = as.Date("2019-07-01"))
  claims <- tibble::tibble(
    ID = "a", AMOUNT_PAID = c(1, 10, 100, 1000, 10000, 100000),
    SS_DESC = "Inpatient",
    SRV_DT = as.Date(c(
      "2019-06-14",   # day before INDEX_DT: prior_cost (1)
      "2019-06-15",   # INDEX_DT itself: NOT prior_cost (end exclusive); before s: NOT outcome (10)
      "2020-06-30",   # day before s+12mo: NOT outcome (100)
      "2020-07-01",   # s+12mo: outcome (1000)
      "2021-06-30",   # last outcome day (10000)
      "2021-07-01"    # s+24mo: NOT outcome (100000)
    ))
  )
  src <- tibble::tibble(dataset_key = "toy", family = "f", suffix_year = 20L,
                        amount_col = "AMOUNT_PAID", date_col = "SRV_DT",
                        setting_col = "SS_DESC")
  out <- suppressMessages(suppressWarnings(smidata::smi_sum_cost_in_windows(
    doubletree_patient_windows(cohort), src,
    negative_amount_policy = config_negative_amount_policy,
    exclude_setting_regex = config_exclude_setting_regex,
    read_fn = function(dataset_key) claims, progress = FALSE
  )))
  by_window <- stats::setNames(out$cost, out$window)
  expect_equal(by_window[["prior_cost"]], 1)
  expect_equal(by_window[["outcome"]], 1000 + 10000)
})
