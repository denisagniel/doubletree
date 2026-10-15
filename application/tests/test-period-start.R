## ============================================================================
## application/tests/test-period-start.R
##
## The 2026-10-15 project-wide rules as this pipeline applies them. Every
## expected value is HAND-CALCULATED in a comment.
##   * s = first day of the month ON OR AFTER INDEX_DT (exposure_start())
##   * the msr_saa period must START EXACTLY at s (smidata ceiling_month anchor);
##     no match, or MSR_DEN < 1 -> ineligible (dropped)
##   * achievement is MSR_NUM / MSR_DEN >= 0.5
##   * outcome [s+12mo, s+24mo) anchors on s; prior_cost [INDEX_DT-12mo, INDEX_DT)
##     anchors on INDEX_DT itself
##   * costs: service files only, negative lines netted
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

test_that("exposure_start() is the first of the month on or after INDEX_DT", {
  ## Hand-calculated:
  ##   2019-06-15 -> next 1st is 2019-07-01
  ##   2019-06-01 -> already a 1st, stays 2019-06-01
  ##   2019-12-31 -> 2020-01-01 (year rolls)
  ##   2019-01-31 -> 2019-02-01
  ##   2020-02-29 (leap day) -> 2020-03-01
  ##   2019-06-02 -> 2019-07-01 (the day after a 1st is NOT a 1st)
  expect_equal(
    exposure_start(as.Date(c("2019-06-15", "2019-06-01", "2019-12-31",
                             "2019-01-31", "2020-02-29", "2019-06-02"))),
    as.Date(c("2019-07-01", "2019-06-01", "2020-01-01",
              "2019-02-01", "2020-03-01", "2019-07-01"))
  )
})

test_that("the ceiling_month rule: exact-start only, dropped otherwise, >= 0.5", {
  skip_if_not_installed("smidata")
  skip_if_not("anchor" %in% names(formals(smidata::smi_select_measure_row)))
  ## Hand-calculated (s = first of the month on or after INDEX_DT):
  ##  P1 INDEX 2019-06-15 -> s 2019-07-01. RY2019-05, -06, -07, -09 exist; only
  ##     RY2019-07 starts at s. NUM 1 / DEN 1 = 1 >= 0.5 -> achieved TRUE.
  ##  P2 INDEX 2019-06-01 (a 1st) -> s 2019-06-01. RY2019-06 starts at s
  ##     (NUM 0 / DEN 1 = 0 -> achieved FALSE). Under the superseded rule this
  ##     period "started on INDEX_DT, not after" and was dropped; now it is the
  ##     exposure period.
  ##  P3 INDEX 2019-06-15 -> s 2019-07-01 but rows are only RY2019-06 and
  ##     RY2019-08, one month either side: a nearest-period rule would match; the
  ##     exact-start rule does NOT -> unmatched, ineligible.
  ##  P4 in the cohort, no msr_saa rows: unmatched.
  ##  P5 INDEX 2019-12-31 -> s 2020-01-01; RY2020-01 has DEN 0: matched but
  ##     elig FALSE (MSR_DEN < 1) -> dropped; achieved NA.
  ##  P6 INDEX 2019-06-20 -> s 2019-07-01; NUM 1 / DEN 2 = 0.5 EXACTLY. With >=
  ##     achieved TRUE (a strict > would give FALSE).
  ##  P7 INDEX 2019-06-20 -> s 2019-07-01; NUM 1 / DEN 3 = 0.333 -> FALSE.
  ##  X  an ID outside the cohort: ignored.
  saa <- tibble::tibble(
    ID = c("P1", "P1", "P1", "P1", "P2", "P2", "P3", "P3", "P5", "P6", "P7", "X"),
    MSR = "SAA",
    MSR_YR = c("RY2019-05", "RY2019-06", "RY2019-07", "RY2019-09",
               "RY2019-06", "RY2019-08", "RY2019-06", "RY2019-08",
               "RY2020-01", "RY2019-07", "RY2019-07", "RY2019-07"),
    MSR_DEN = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3, 1),
    MSR_NUM = c(1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1)
  )
  idx <- tibble::tibble(
    ID = c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
    INDEX_DT = as.Date(c("2019-06-15", "2019-06-01", "2019-06-15", "2019-06-15",
                         "2019-12-31", "2019-06-20", "2019-06-20"))
  )
  out <- smidata::smi_select_measure_row(
    saa, idx, achievement_threshold = config_saa_achievement_threshold,
    anchor = "ceiling_month"
  )
  by_id <- function(id) out[out$ID == id, ]
  expect_equal(nrow(out), 7L)                       # one row per cohort patient
  expect_equal(by_id("P1")$MSR_YR_selected, "RY2019-07")
  expect_true(by_id("P1")$achieved)
  expect_equal(by_id("P2")$MSR_YR_selected, "RY2019-06")
  expect_false(by_id("P2")$achieved)
  expect_false(by_id("P3")$matched)
  expect_false(by_id("P3")$elig)
  expect_false(by_id("P4")$matched)
  expect_true(by_id("P5")$matched)
  expect_false(by_id("P5")$elig)
  expect_true(is.na(by_id("P5")$achieved))
  expect_true(by_id("P6")$achieved)                 # 0.5 >= 0.5
  expect_false(by_id("P7")$achieved)
  expect_true(all(out$month_gap[out$matched] == 0L))
  ## Eligible patients are exactly P1, P2, P6, P7: everyone else is dropped.
  expect_equal(sort(out$ID[out$elig %in% TRUE]), c("P1", "P2", "P6", "P7"))
  ## The selected period's start equals exposure_start(INDEX_DT) for every match.
  m <- out[out$matched, ]
  expect_equal(
    msr_yr_period_start(m$MSR_YR_selected),
    exposure_start(idx$INDEX_DT[match(m$ID, idx$ID)])
  )
  ## Loosening the gap with this anchor is refused by smidata, not tolerated.
  expect_error(
    smidata::smi_select_measure_row(saa, idx, anchor = "ceiling_month", max_month_gap = 1L)
  )
})

test_that("no pharmacy file is in any cost source and there are 18 service files", {
  ## 2 service families (aim3_svc, new_svc) x 9 years (16:24) = 18 keys.
  expect_equal(nrow(COST_SOURCES), 18L)
  expect_false(any(grepl("pharm", COST_SOURCES$dataset_key)))
  expect_setequal(unique(COST_SOURCES$family), c("aim3_svc_cost", "new_svc_cost"))
  expect_equal(length(cost_dataset_keys()), 18L)
})

test_that("windows: outcome anchors on s, prior_cost on INDEX_DT itself", {
  skip_if_not_installed("smidata")
  ## Hand-calculated half-open windows.
  ##  P1 INDEX_DT 2019-06-15, s = 2019-07-01:
  ##    prior_cost [2019-06-15 - 12mo, 2019-06-15)       = [2018-06-15, 2019-06-15)
  ##    outcome    [2019-07-01 + 12mo, 2019-07-01 + 24mo) = [2020-07-01, 2021-07-01)
  ##    exposure   [s, s+12mo) = [2019-07-01, 2020-07-01) abuts outcome exactly
  ##    prior_cost ends 16 days (2019-06-15 .. 2019-07-01) BEFORE the exposure starts
  ##  P2 INDEX_DT 2019-06-01, s = 2019-06-01 (the two anchors coincide):
  ##    prior_cost [2018-06-01, 2019-06-01); outcome [2020-06-01, 2021-06-01)
  cohort <- tibble::tibble(
    ID = c("P1", "P2"),
    INDEX_DT = as.Date(c("2019-06-15", "2019-06-01"))
  ) |>
    dplyr::mutate(T0 = exposure_start(.data$INDEX_DT))
  w <- doubletree_patient_windows(cohort)
  expect_equal(nrow(w), 4L)                         # 2 patients x 2 windows
  get_w <- function(id, win) w[w$id == id & w$window == win, ]
  expect_equal(get_w("P1", "prior_cost")$window_start, as.Date("2018-06-15"))
  expect_equal(get_w("P1", "prior_cost")$window_end, as.Date("2019-06-15"))
  expect_equal(get_w("P1", "outcome")$window_start, as.Date("2020-07-01"))
  expect_equal(get_w("P1", "outcome")$window_end, as.Date("2021-07-01"))
  expect_equal(get_w("P2", "prior_cost")$window_start, as.Date("2018-06-01"))
  expect_equal(get_w("P2", "prior_cost")$window_end, as.Date("2019-06-01"))
  expect_equal(get_w("P2", "outcome")$window_start, as.Date("2020-06-01"))
  expect_equal(get_w("P2", "outcome")$window_end, as.Date("2021-06-01"))
  ## Exposure end == outcome start, built independently from s.
  exposure_end <- seq(as.Date("2019-07-01"), by = "12 months", length.out = 2L)[[2]]
  expect_equal(exposure_end, get_w("P1", "outcome")$window_start)
  ## prior_cost ends 16 days before s for P1 (June 15 -> July 1).
  expect_equal(as.integer(cohort$T0[1] - get_w("P1", "prior_cost")$window_end), 16L)
})

test_that("service costs are netted over the INDEX_DT-anchored prior_cost window", {
  skip_if_not_installed("smidata")
  skip_if_not("negative_amount_policy" %in% names(formals(smidata::smi_sum_cost_in_windows)))
  ## One service file (suffix 19 = calendar 2019). P1: INDEX_DT 2019-06-15.
  ## prior_cost window [2018-06-15, 2019-06-15). Claims in the 2019 file:
  ##   2019-01-10   100            in window, counted
  ##   2019-02-10   -30            in window, NEGATIVE, kept (net)
  ##   2019-03-01    50 + "Managed Care Invoice" -> excluded
  ##   2019-06-15   500            ON INDEX_DT: half-open end, NOT in the window
  ##                               (the index-event cost stays out of prior_cost)
  ##   2019-06-20   999            after INDEX_DT, not in the window
  ## net prior_cost = 100 + (-30) = 70; gross_positive = 100; negative_amount = -30.
  ## P2: 100 + (-120) = -20 (net negative): gross 100, negative -120, n_negative 1.
  src <- COST_SOURCES[COST_SOURCES$dataset_key == "aim3_svc_cost_19", ]
  fake_read <- function(key) {
    tibble::tibble(
      ID = c("P1", "P1", "P1", "P1", "P1", "P2", "P2"),
      SRV_DT = as.Date(c("2019-01-10", "2019-02-10", "2019-03-01", "2019-06-15",
                         "2019-06-20", "2019-01-05", "2019-01-06")),
      AMOUNT_PAID = c(100, -30, 50, 500, 999, 100, -120),
      SEQ_ID = 1L,
      SS_DESC = c("Inpatient", "Inpatient", "Managed Care Invoice", "Inpatient",
                  "Inpatient", "Inpatient", "Inpatient")
    )
  }
  cohort <- tibble::tibble(
    ID = c("P1", "P2"),
    INDEX_DT = as.Date(c("2019-06-15", "2019-06-15"))
  ) |>
    dplyr::mutate(T0 = exposure_start(.data$INDEX_DT))
  w <- doubletree_patient_windows(cohort)
  w <- w[w$window == "prior_cost", ]
  res <- smidata::smi_sum_cost_in_windows(
    w, src, negative_amount_policy = "net", read_fn = fake_read, progress = FALSE
  )
  p1 <- res[res$id == "P1", ]
  p2 <- res[res$id == "P2", ]
  expect_equal(p1$cost, 70)
  expect_equal(p1$gross_positive, 100)
  expect_equal(p1$negative_amount, -30)
  expect_equal(p1$n_negative, 1L)
  expect_equal(p1$n_excluded_mci, 1L)
  expect_equal(p2$cost, -20)
  ## The floor is applied AFTER netting, to the window total: -20 -> 0, while a
  ## positive net (70) is untouched. (TODO: open decision, see _config.R.)
  floored <- floor_negative_prior_cost(c(p1$cost, p2$cost), floor_value = 0)
  expect_equal(as.numeric(floored), c(70, 0))
  expect_equal(attr(floored, "n_floored"), 1L)
})

test_that("enrollment decision records gap 0, gap_mode, require_agreement, 24 months", {
  d <- confirmed_value("enrollment_source")
  expect_equal(d$gap_months, 0L)
  expect_equal(d$gap_mode, "total")
  expect_equal(d$type_policy, "require_agreement")
  expect_equal(config_complete_case_months, 24L)
})
