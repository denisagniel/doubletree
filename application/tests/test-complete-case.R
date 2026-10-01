## ============================================================================
## application/tests/test-complete-case.R
##
## THINNED 2026-10-01 by the smidata migration. The month-by-month logic this
## file used to test -- parse_year_month(), ym_month_index(), and
## compute_complete_case() -- now lives in smidata as
## smi_parse_year_month() / smi_date_month_index() / smi_enrollment_panel() /
## smi_enrollment_coverage(), and is tested there in
## tests/testthat/test-enrollment.R (104 assertions).
##
## Those tests are NOT duplicated here. Same split the 0.2.0 msr_aap promotion
## used (dual-bounds' 92 local assertions became 14 local + 78 package): the
## package suite owns the dataset mechanics, this file owns the part that is
## THIS PAPER'S -- the translation from OPEN_DECISIONS$enrollment_source to
## smidata's explicit arguments, and the shape 06_assemble_analytic_data.R
## depends on.
##
## Specifically moved OUT of this file, because they test smidata now:
##   * YEAR_MONTH format parsing (both formats, mixed-format abort,
##     out-of-range month, implausible year, NA, wrong input type)
##   * month-index adjacency / day-of-month invariance / monotonicity
##   * membership-not-a-count, including the mid-month INDEX_DT window shift
##   * out-of-window rows not leaking in to satisfy coverage
##   * the absent-patient and all-absent cases, and the tally internals
##   * the TYPE-disagreement abort
##
## What stays: everything whose expected value depends on a doubletree
## DECISION rather than on the dataset.
## ============================================================================

skip_if_not_installed("smidata")

## The migration raised this paper's smidata floor: these functions do not
## exist before 0.5.0. Checked here rather than letting a call fail with
## "object not found", which reads as a typo rather than a stale install.
test_that("the installed smidata carries the promoted enrollment functions", {
  needed <- c("smi_enrollment_panel", "smi_enrollment_coverage",
              "smi_parse_year_month")
  expect_true(
    all(needed %in% getNamespaceExports("smidata")),
    label = paste0(
      "smidata ", as.character(utils::packageVersion("smidata")),
      " exports the promoted enrollment functions"
    )
  )
})

## ---- the translation --------------------------------------------------------
##
## 05_complete_case.R turns two recorded decision values into three smidata
## arguments. Each assertion below is about that mapping, not about what
## smidata then does with it.

test_that("enrollment_source is still confirmed, with the two values the translation reads", {
  ## If this regresses to open_decision, 05 aborts via confirmed_value() --
  ## but it would abort at run time, on the server. Catch it here.
  d <- OPEN_DECISIONS$enrollment_source
  expect_identical(d$status, "confirmed")
  expect_identical(d$value$gap_months, 0L)
  ## The TYPE value is prose, not a token, so match on its meaning rather
  ## than its exact wording.
  expect_match(d$value$type_filter, "none", ignore.case = TRUE)
})

test_that("gap_months = 0 makes gap_mode unambiguous, which is why 05 may hardcode it", {
  ## smidata separates a gap COUNT ("total") from a gap RUN ("consecutive").
  ## At 0 they coincide -- zero uncovered months total is the same condition
  ## as a longest run of zero -- and that equivalence is the ONLY reason 05
  ## passing a literal "total" is defensible. Asserted, because if it stopped
  ## being true the literal would silently become a choice nobody made.
  ##
  ## Fixture: patient 1 covers all 3 required months; patient 2 is missing
  ## 2019-07 (one uncovered month, longest run 1).
  cohort <- tibble::tibble(ID = c(1, 2), INDEX_DT = as.Date("2019-06-01"))
  flag <- tibble::tibble(
    ID = c(1, 1, 1, 2, 2),
    INDEX_DT = as.Date("2019-06-01"),
    TYPE = "IP",
    YEAR_MONTH = c("2019-06", "2019-07", "2019-08", "2019-06", "2019-08"),
    MEDICAID_FLAG = 1L
  )
  panel <- smidata::smi_enrollment_panel(flag, cohort, "require_agreement")
  tot <- smidata::smi_enrollment_coverage(panel, cohort, 3L, 0L, "total")
  con <- smidata::smi_enrollment_coverage(panel, cohort, 3L, 0L, "consecutive")
  expect_identical(tot$covered, con$covered)
  ## Hand-counted: only patient 1 is covered at zero tolerance.
  expect_identical(tot$covered, c(TRUE, FALSE))
})

test_that("require_agreement reproduces the deleted helper's TYPE behaviour", {
  ## THE translation's one substantive choice. The deleted
  ## compute_complete_case() implemented the PI's "ignore TYPE" AND aborted on
  ## any (ID, month) TYPE disagreement, because "ignore" is only a no-op where
  ## no disagreement exists. smidata's "any" WARNS instead; only
  ## "require_agreement" aborts.
  cohort <- tibble::tibble(ID = 1, INDEX_DT = as.Date("2019-06-01"))
  agreeing <- tibble::tibble(
    ID = 1, INDEX_DT = as.Date("2019-06-01"), TYPE = c("IP", "OP"),
    YEAR_MONTH = "2019-06", MEDICAID_FLAG = c(1L, 1L)
  )
  disagreeing <- agreeing
  disagreeing$MEDICAID_FLAG <- c(1L, 0L)

  ## Where TYPE agrees, the two policies are identical -- so the migration
  ## changes no result on data without disagreement.
  p_req <- smidata::smi_enrollment_panel(agreeing, cohort, "require_agreement")
  p_any <- smidata::smi_enrollment_panel(agreeing, cohort, "any")
  expect_identical(p_req$enrolled, p_any$enrolled)

  ## Where it disagrees, only "require_agreement" refuses -- which is the
  ## deleted helper's behaviour, and the reason 05 does not pass "any".
  expect_error(
    smidata::smi_enrollment_panel(disagreeing, cohort, "require_agreement"),
    "DISAGREE"
  )
  expect_warning(
    smidata::smi_enrollment_panel(disagreeing, cohort, "any"),
    "disagreeing"
  )
})

## ---- the shape 06 depends on ------------------------------------------------

test_that("the coverage result converts to 06's expected complete_case contract", {
  ## 06_assemble_analytic_data.R does
  ##   left_join(complete_case, by = "ID") |> filter(complete_case %in% TRUE)
  ## so it needs an ID column and a LOGICAL complete_case that is never NA.
  ## 05 renames smidata's `covered`; this asserts the renamed shape, because a
  ## silent NA there would drop every row via %in% TRUE and read as
  ## "complete_case not yet resolved" rather than as a bug.
  cohort <- tibble::tibble(ID = c(1, 2, 3), INDEX_DT = as.Date("2019-06-01"))
  flag <- tibble::tibble(
    ID = c(1, 1, 1, 2, 2),
    INDEX_DT = as.Date("2019-06-01"),
    TYPE = "IP",
    YEAR_MONTH = c("2019-06", "2019-07", "2019-08", "2019-06", "2019-07"),
    MEDICAID_FLAG = 1L
  )
  panel <- smidata::smi_enrollment_panel(flag, cohort, "require_agreement")
  cov <- smidata::smi_enrollment_coverage(panel, cohort, 3L, 0L, "total")

  complete_case <- cov |>
    dplyr::transmute(ID = .data$ID, complete_case = .data$covered,
                     reason = .data$reason)

  expect_named(complete_case, c("ID", "complete_case", "reason"))
  expect_type(complete_case$complete_case, "logical")
  expect_false(anyNA(complete_case$complete_case))
  expect_equal(nrow(complete_case), nrow(cohort))
  ## Hand-counted from the fixture: patient 1 covers all three required
  ## months; patient 2 is missing 2019-08; patient 3 has no rows at all.
  expect_identical(complete_case$complete_case, c(TRUE, FALSE, FALSE))
  expect_identical(
    complete_case$reason,
    c("covered", "insufficient_coverage", "absent_from_flag_table")
  )
})

## ---- this paper's own 24-month window --------------------------------------

test_that("the real config window separates a complete patient from a one-month-short one", {
  ## config_complete_case_months (24) is THIS PAPER'S decision, so it is
  ## tested here and not in smidata, whose own suite uses small windows that
  ## are checkable by inspection. These are the same two known-truth patients
  ## 05's demonstration block prints.
  expect_equal(config_complete_case_months, 24L)

  cohort <- tibble::tibble(ID = c(1, 2), INDEX_DT = as.Date("2019-06-01"))
  months <- format(
    seq(as.Date("2019-06-01"), by = "month",
        length.out = config_complete_case_months),
    "%Y-%m"
  )
  flag <- dplyr::bind_rows(
    tibble::tibble(ID = 1, INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
                   YEAR_MONTH = months, MEDICAID_FLAG = 1L),
    ## Drop the LAST required month -- the one an off-by-one is most likely to
    ## skip silently.
    tibble::tibble(ID = 2, INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
                   YEAR_MONTH = months[-length(months)], MEDICAID_FLAG = 1L)
  )
  panel <- smidata::smi_enrollment_panel(flag, cohort, "require_agreement")
  cov <- smidata::smi_enrollment_coverage(
    panel, cohort, config_complete_case_months, 0L, "total"
  )
  expect_identical(cov$covered, c(TRUE, FALSE))
  ## Hand-counted: 24 required months, 23 covered for patient 2.
  expect_identical(cov$n_covered, c(24L, 23L))
  expect_identical(cov$n_uncovered, c(0L, 1L))
  expect_identical(cov$reason[[2]], "insufficient_coverage")
})

test_that("05_complete_case.R no longer sources a local complete-case helper", {
  ## The helper is DELETED, not shimmed. A re-introduced local copy is the
  ## exact failure the promotion was meant to end -- dual-bounds' and
  ## doubletree's copies had already drifted 82 lines apart.
  ##
  ## Uses app_dir_for_tests (helper-toy.R's own convention) rather than
  ## probing the working directory: testthat runs a test file from the test
  ## file's own directory, so a cwd-relative path works under run_tests.R and
  ## breaks under test_file(). That difference is how this was caught.
  expect_false(
    file.exists(file.path(app_dir_for_tests, "helpers", "complete_case.R"))
  )
  txt <- paste(
    readLines(file.path(app_dir_for_tests, "05_complete_case.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl('source\\(.*"complete_case\\.R"', txt))
  expect_match(txt, "smi_enrollment_coverage")
})
