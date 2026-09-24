## ============================================================================
## application/tests/test-complete-case.R
##
## parse_year_month(), ym_month_index(), and compute_complete_case(): the real
## implementation behind OPEN_DECISIONS$enrollment_source's PI resolution
## (application/_config.R, confirmed 2026-09-24).
##
## Every expected value is HAND-CALCULATED and the calculation is stated, per
## this directory's convention (see test-discretize-prior-cost.R's header).
## window_months is kept small (3) throughout so every required-month set is
## checkable by inspection; config_complete_case_months (24) is exercised
## separately in the "uses the real config value" test at the bottom.
## ============================================================================

## ---- parse_year_month() -----------------------------------------------------

test_that("YYYY-MM format parses to the first of the month", {
  out <- parse_year_month(c("2019-06", "2019-07", "2020-01"))
  expect_equal(out, as.Date(c("2019-06-01", "2019-07-01", "2020-01-01")))
})

test_that("YYYYMM (flat) format parses to the first of the month", {
  out <- parse_year_month(c("201906", "201907", "202001"))
  expect_equal(out, as.Date(c("2019-06-01", "2019-07-01", "2020-01-01")))
})

test_that("Date/POSIXct input is floored to month-start, not re-parsed", {
  expect_equal(parse_year_month(as.Date("2020-03-17")), as.Date("2020-03-01"))
  expect_equal(
    parse_year_month(as.POSIXct("2020-03-17 08:00:00", tz = "UTC")),
    as.Date("2020-03-01")
  )
})

test_that("an empty input returns an empty Date vector, not an error", {
  expect_equal(parse_year_month(character()), as.Date(character()))
})

test_that("mixing YYYY-MM and YYYYMM in the same vector aborts", {
  ## Neither format matches ALL values, so this is the "matches neither"
  ## abort path, not a silent per-row dispatch.
  expect_error(
    parse_year_month(c("2019-06", "201907")),
    "does not match either anticipated format"
  )
})

test_that("an out-of-range month (13) aborts rather than wrapping", {
  expect_error(parse_year_month("2019-13"), "month component out of range")
})

test_that("an implausible year aborts even though the regex matches", {
  ## "0099-06" matches ^\\d{4}-\\d{2}$ (4 digits) but as.integer("0099") = 99,
  ## outside the 1900-2100 band.
  expect_error(parse_year_month("0099-06"), "year component outside")
})

test_that("NA in the input aborts rather than being coerced", {
  expect_error(parse_year_month(c("2019-06", NA)), "NA")
})

test_that("non-character, non-Date input aborts", {
  expect_error(parse_year_month(201906L), "must be character")
  expect_error(parse_year_month(list("2019-06")), "must be character")
})

## ---- ym_month_index() --------------------------------------------------------

test_that("adjacent calendar months differ by exactly 1", {
  expect_equal(
    ym_month_index(as.Date("2020-01-01")) - ym_month_index(as.Date("2019-12-01")),
    1L
  )
})

test_that("day-of-month does not affect the index -- only year and month", {
  expect_equal(
    ym_month_index(as.Date("2019-06-01")),
    ym_month_index(as.Date("2019-06-28"))
  )
})

test_that("the index is monotone across a full year", {
  months <- as.Date(sprintf("2019-%02d-01", 1:12))
  idx <- ym_month_index(months)
  expect_equal(diff(idx), rep(1L, 11L))
})

## ---- compute_complete_case(): shared fixture ---------------------------------
##
## window_months = 3, anchor 2019-06-01 (or 2019-06-15 for the mid-month case),
## so the required month set is always {2019-06, 2019-07, 2019-08} regardless
## of which day of June INDEX_DT falls on.

cc_cohort <- function() {
  tibble::tibble(
    ID = c("P1", "P2", "P3"),
    INDEX_DT = as.Date(c("2019-06-01", "2019-06-01", "2019-06-01"))
  )
}

## P1: fully covered all 3 required months, single TYPE.
## P2: covered Jun/Jul, Aug entirely ABSENT (no row at all).
## P3: absent from monthly_flag entirely.
cc_monthly_flag <- function() {
  tibble::tibble(
    ID = c(rep("P1", 3L), rep("P2", 2L)),
    INDEX_DT = as.Date("2019-06-01"),
    TYPE = "IP",
    YEAR_MONTH = c("2019-06", "2019-07", "2019-08", "2019-06", "2019-07"),
    MEDICAID_FLAG = c(1L, 1L, 1L, 1L, 1L)
  )
}

test_that("a fully-covered patient is complete_case = TRUE, reason = covered", {
  out <- compute_complete_case(cc_monthly_flag(), cc_cohort(), window_months = 3L)
  p1 <- out[out$ID == "P1", ]
  expect_true(p1$complete_case)
  expect_equal(p1$reason, "covered")
})

test_that("a patient missing exactly one required month is insufficient_coverage", {
  out <- compute_complete_case(cc_monthly_flag(), cc_cohort(), window_months = 3L)
  p2 <- out[out$ID == "P2", ]
  expect_false(p2$complete_case)
  expect_equal(p2$reason, "insufficient_coverage")
})

test_that("a patient with zero rows in monthly_flag is absent_from_flag_table", {
  out <- compute_complete_case(cc_monthly_flag(), cc_cohort(), window_months = 3L)
  p3 <- out[out$ID == "P3", ]
  expect_false(p3$complete_case)
  expect_equal(p3$reason, "absent_from_flag_table")
})

test_that("output is exactly one row per cohort patient, no more, no fewer", {
  out <- compute_complete_case(cc_monthly_flag(), cc_cohort(), window_months = 3L)
  expect_equal(nrow(out), 3L)
  expect_setequal(out$ID, c("P1", "P2", "P3"))
})

test_that("the tally attribute matches the hand-counted fixture exactly", {
  out <- compute_complete_case(cc_monthly_flag(), cc_cohort(), window_months = 3L)
  tally <- attr(out, "tally")
  expect_equal(tally$n_cohort, 3L)
  expect_equal(tally$n_complete, 1L)                 # P1
  expect_equal(tally$n_insufficient_coverage, 1L)    # P2
  expect_equal(tally$n_absent, 1L)                   # P3
  expect_equal(tally$n_type_disagreement_pairs, 0L)
})

## ---- mid-month INDEX_DT: the bug this function's Details section names ------

test_that("a mid-month INDEX_DT still requires INDEX_DT's OWN calendar month", {
  ## INDEX_DT = 2019-06-15 (day 15, not day 1). If the window were built from
  ## the RAW (un-floored) INDEX_DT via Date arithmetic, the half-open interval
  ## [2019-06-15, 2019-09-15) would exclude June entirely and include a sliver
  ## of September instead -- exactly the shift this function's Details
  ## section describes. Required months must still be {2019-06, 07, 08}.
  cohort_mid <- tibble::tibble(ID = "PM", INDEX_DT = as.Date("2019-06-15"))
  mf_covered <- tibble::tibble(
    ID = "PM", INDEX_DT = as.Date("2019-06-15"), TYPE = "IP",
    YEAR_MONTH = c("2019-06", "2019-07", "2019-08"), MEDICAID_FLAG = c(1L, 1L, 1L)
  )
  out_covered <- compute_complete_case(mf_covered, cohort_mid, window_months = 3L)
  expect_true(out_covered$complete_case)
  expect_equal(out_covered$reason, "covered")

  ## Same patient, June's row REMOVED (only Jul/Aug present). If June were
  ## wrongly excluded from the required set, this would ALSO show "covered"
  ## (2 rows would look sufficient against a shifted 2-month requirement) --
  ## it must instead be insufficient, because June is required and absent.
  mf_missing_june <- mf_covered[mf_covered$YEAR_MONTH != "2019-06", ]
  out_missing <- compute_complete_case(mf_missing_june, cohort_mid, window_months = 3L)
  expect_false(out_missing$complete_case)
  expect_equal(out_missing$reason, "insufficient_coverage")
})

## ---- out-of-window rows must not "leak in" to satisfy coverage --------------

test_that("flag=1 rows OUTSIDE the required window do not count toward it", {
  ## 4 rows with MEDICAID_FLAG == 1 total (May, Jun, Jul, Sep), but the
  ## REQUIRED window for anchor 2019-06 at window_months=3 is {Jun,Jul,Aug}.
  ## Aug is missing, so despite 4 "covered" rows total, this must be
  ## insufficient -- a naive "count of flag==1 rows >= 3" check would wrongly
  ## pass.
  cohort_ow <- tibble::tibble(ID = "PO", INDEX_DT = as.Date("2019-06-01"))
  mf_ow <- tibble::tibble(
    ID = "PO", INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
    YEAR_MONTH = c("2019-05", "2019-06", "2019-07", "2019-09"),
    MEDICAID_FLAG = c(1L, 1L, 1L, 1L)
  )
  out <- compute_complete_case(mf_ow, cohort_ow, window_months = 3L)
  expect_false(out$complete_case)
  expect_equal(out$reason, "insufficient_coverage")
})

## ---- TYPE is ignored: agreement is silently fine, disagreement aborts -------

test_that("two TYPE rows for the same month that AGREE do not abort and count once", {
  cohort_ty <- tibble::tibble(ID = "PT", INDEX_DT = as.Date("2019-06-01"))
  mf_agree <- tibble::tibble(
    ID = "PT", INDEX_DT = as.Date("2019-06-01"),
    TYPE = c("IP", "OP", "IP", "IP"),
    YEAR_MONTH = c("2019-06", "2019-06", "2019-07", "2019-08"),
    MEDICAID_FLAG = c(1L, 1L, 1L, 1L)
  )
  out <- compute_complete_case(mf_agree, cohort_ty, window_months = 3L)
  expect_true(out$complete_case)
  tally <- attr(out, "tally")
  expect_equal(tally$n_type_disagreement_pairs, 0L)
})

test_that("TYPE rows that DISAGREE on MEDICAID_FLAG abort, naming CHECK 4", {
  cohort_ty <- tibble::tibble(ID = "PT", INDEX_DT = as.Date("2019-06-01"))
  mf_disagree <- tibble::tibble(
    ID = "PT", INDEX_DT = as.Date("2019-06-01"),
    TYPE = c("IP", "OP", "IP", "IP"),
    YEAR_MONTH = c("2019-06", "2019-06", "2019-07", "2019-08"),
    MEDICAID_FLAG = c(1L, 0L, 1L, 1L)   # IP says 1, OP says 0 for 2019-06
  )
  expect_error(
    compute_complete_case(mf_disagree, cohort_ty, window_months = 3L),
    "DISAGREE"
  )
  expect_error(
    compute_complete_case(mf_disagree, cohort_ty, window_months = 3L),
    "CHECK 4"
  )
})

test_that("a disagreement OUTSIDE the required window does not abort", {
  ## Same disagreement pattern as above, but in 2019-04 -- before the
  ## anchor's required window ({2019-06,07,08}). Disagreement is checked at
  ## cohort/window-RESTRICTED scope on purpose (this function's Details),
  ## so this must not fire the abort.
  cohort_ty <- tibble::tibble(ID = "PT", INDEX_DT = as.Date("2019-06-01"))
  mf_out_of_window <- tibble::tibble(
    ID = "PT", INDEX_DT = as.Date("2019-06-01"),
    TYPE = c("IP", "OP", "IP", "IP", "IP"),
    YEAR_MONTH = c("2019-04", "2019-04", "2019-06", "2019-07", "2019-08"),
    MEDICAID_FLAG = c(1L, 0L, 1L, 1L, 1L)  # disagreement in 2019-04 only
  )
  out <- compute_complete_case(mf_out_of_window, cohort_ty, window_months = 3L)
  expect_true(out$complete_case)
})

## ---- NA MEDICAID_FLAG does not count as covered ------------------------------

test_that("an NA MEDICAID_FLAG warns and does not count its month as covered", {
  cohort_na <- tibble::tibble(ID = "PN", INDEX_DT = as.Date("2019-06-01"))
  mf_na <- tibble::tibble(
    ID = "PN", INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
    YEAR_MONTH = c("2019-06", "2019-07", "2019-08"),
    MEDICAID_FLAG = c(NA_integer_, 1L, 1L)
  )
  expect_warning(
    out <- compute_complete_case(mf_na, cohort_na, window_months = 3L),
    "NA"
  )
  expect_false(out$complete_case)
  expect_equal(out$reason, "insufficient_coverage")
})

## ---- input validation ---------------------------------------------------------

test_that("missing required columns abort, naming which table and column", {
  expect_error(
    compute_complete_case(
      dplyr::select(cc_monthly_flag(), -"MEDICAID_FLAG"), cc_cohort(), 3L
    ),
    "monthly_flag.*missing column"
  )
  expect_error(
    compute_complete_case(
      cc_monthly_flag(), dplyr::select(cc_cohort(), -"INDEX_DT"), 3L
    ),
    "cohort.*missing column"
  )
})

test_that("window_months validation rejects zero, negative, non-integer, and vectors", {
  mf <- cc_monthly_flag(); co <- cc_cohort()
  expect_error(compute_complete_case(mf, co, 0L), "positive integer")
  expect_error(compute_complete_case(mf, co, -1L), "positive integer")
  expect_error(compute_complete_case(mf, co, 1.5), "positive integer")
  expect_error(compute_complete_case(mf, co, c(3L, 4L)), "positive integer")
  expect_error(compute_complete_case(mf, co, NA_integer_), "positive integer")
})

test_that("an empty cohort aborts rather than returning zero rows", {
  expect_error(
    compute_complete_case(cc_monthly_flag(), cc_cohort()[0, ], 3L),
    "empty"
  )
})

test_that("a duplicate ID in cohort aborts", {
  dup_cohort <- dplyr::bind_rows(cc_cohort(), cc_cohort()[1, ])
  expect_error(compute_complete_case(cc_monthly_flag(), dup_cohort, 3L), "duplicate")
})

test_that("an NA INDEX_DT in cohort aborts", {
  na_cohort <- cc_cohort()
  na_cohort$INDEX_DT[1] <- NA
  expect_error(compute_complete_case(cc_monthly_flag(), na_cohort, 3L), "NA")
})

test_that("an ID class mismatch between the two tables aborts", {
  mf_int_id <- cc_monthly_flag()
  mf_int_id$ID <- match(mf_int_id$ID, c("P1", "P2", "P3"))  # character -> integer
  expect_error(
    compute_complete_case(mf_int_id, cc_cohort(), 3L),
    "class"
  )
})

test_that("more than one distinct INDEX_DT for one ID within monthly_flag aborts", {
  mf_bad <- cc_monthly_flag()
  mf_bad$INDEX_DT[mf_bad$ID == "P1"][1] <- as.Date("2019-07-01")  # was all 06-01
  expect_error(
    compute_complete_case(mf_bad, cc_cohort(), 3L),
    "more than one distinct"
  )
})

test_that("monthly_flag's INDEX_DT disagreeing with cohort's aborts", {
  mismatched_cohort <- cc_cohort()
  mismatched_cohort$INDEX_DT[mismatched_cohort$ID == "P1"] <- as.Date("2019-07-01")
  expect_error(
    compute_complete_case(cc_monthly_flag(), mismatched_cohort, 3L),
    "different"
  )
})

test_that("a MEDICAID_FLAG value outside {0, 1} (and not NA) aborts", {
  mf_bad_flag <- cc_monthly_flag()
  mf_bad_flag$MEDICAID_FLAG[1] <- 2L
  expect_error(compute_complete_case(mf_bad_flag, cc_cohort(), 3L), "outside")
})

test_that("every cohort patient coming back absent aborts as a likely join-key bug", {
  ## Same class, but genuinely disjoint ID values -- the class check passes,
  ## but nobody in cohort matches anybody in monthly_flag.
  disjoint_mf <- cc_monthly_flag()
  disjoint_mf$ID <- paste0("OTHER_", disjoint_mf$ID)
  expect_error(
    compute_complete_case(disjoint_mf, cc_cohort(), 3L),
    "absent_from_flag_table"
  )
})

## ---- uses the real config value, end to end ----------------------------------

test_that("a 24-month window (config_complete_case_months) works end to end", {
  ## config_complete_case_months is sourced from _config.R by helper-toy.R's
  ## source() chain (via 05_complete_case.R's own convention). Confirm it is
  ## 24 (the PI-confirmed value) and that a patient covered for exactly those
  ## 24 months is complete.
  expect_equal(config_complete_case_months, 24L)

  months_24 <- format(
    seq(as.Date("2019-06-01"), by = "month", length.out = config_complete_case_months),
    "%Y-%m"
  )
  cohort_24 <- tibble::tibble(ID = "P24", INDEX_DT = as.Date("2019-06-01"))
  mf_24 <- tibble::tibble(
    ID = "P24", INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
    YEAR_MONTH = months_24, MEDICAID_FLAG = 1L
  )
  out <- compute_complete_case(mf_24, cohort_24, config_complete_case_months)
  expect_true(out$complete_case)

  ## Drop the LAST required month (the one most likely to be silently
  ## skipped by an off-by-one): must flip to insufficient.
  mf_23 <- mf_24[-length(months_24), ]
  out_23 <- compute_complete_case(mf_23, cohort_24, config_complete_case_months)
  expect_false(out_23$complete_case)
  expect_equal(out_23$reason, "insufficient_coverage")
})
