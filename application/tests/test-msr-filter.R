## ============================================================================
## application/tests/test-msr-filter.R
##
## filter_to_msr_code(): the real logic behind 02_population_and_eligibility.R's
## MSR-contamination guard (see that file's header, and OPEN_DECISIONS$msr_code's
## note in application/_config.R, for why this filter exists at all -- it closes
## a gap in smidata::smi_select_aap_row()'s abort-on-tie protection that a
## non-tie-forming non-AAP row would otherwise pass through silently).
##
## Every expected value is HAND-CALCULATED, per this directory's convention.
## ============================================================================

test_that("an all-AAP input passes through with zero exclusions", {
  df <- tibble::tibble(ID = c("P1", "P2", "P3"), MSR = "AAP", MSR_DEN = 10, MSR_NUM = 5)
  out <- filter_to_msr_code(df, "AAP")

  expect_equal(nrow(out), 3L)
  expect_equal(out$ID, df$ID)
  expect_equal(attr(out, "n_excluded"), 0L)
  expect_equal(attr(out, "excluded_msr_values"), character(0))
})

test_that("mixed MSR values are filtered, and the exclusion is reported correctly", {
  df <- tibble::tibble(
    ID = c("P1", "P2", "P3", "P4"),
    MSR = c("AAP", "OTHER", "AAP", "THIRD"),
    MSR_DEN = 10, MSR_NUM = 5
  )
  ## Suppress the warning here to inspect the return value; the warning
  ## itself is checked separately below.
  out <- suppressWarnings(filter_to_msr_code(df, "AAP"))

  expect_equal(nrow(out), 2L)
  expect_setequal(out$ID, c("P1", "P3"))
  expect_equal(attr(out, "n_excluded"), 2L)
  expect_setequal(attr(out, "excluded_msr_values"), c("OTHER", "THIRD"))
})

test_that("blank-padded MSR values are matched after trimws()", {
  df <- tibble::tibble(ID = c("P1", "P2"), MSR = c("AAP  ", "  AAP"), MSR_DEN = 10, MSR_NUM = 5)
  out <- filter_to_msr_code(df, "AAP")

  expect_equal(nrow(out), 2L)
  expect_equal(attr(out, "n_excluded"), 0L)
})

test_that("an all-non-AAP input excludes every row, not a subset", {
  df <- tibble::tibble(ID = c("P1", "P2"), MSR = c("OTHER", "OTHER"), MSR_DEN = 10, MSR_NUM = 5)
  out <- suppressWarnings(filter_to_msr_code(df, "AAP"))

  expect_equal(nrow(out), 0L)
  expect_equal(attr(out, "n_excluded"), 2L)
  expect_equal(attr(out, "excluded_msr_values"), "OTHER")
})

test_that("a warning fires on exclusion, naming the excluded codes", {
  df <- tibble::tibble(ID = c("P1", "P2"), MSR = c("AAP", "BAD_CODE"), MSR_DEN = 10, MSR_NUM = 5)
  expect_warning(filter_to_msr_code(df, "AAP"), "EXCLUDED")
  expect_warning(filter_to_msr_code(df, "AAP"), "BAD_CODE")
})

test_that("no warning fires when there is nothing to exclude", {
  df <- tibble::tibble(ID = "P1", MSR = "AAP", MSR_DEN = 10, MSR_NUM = 5)
  expect_no_warning(filter_to_msr_code(df, "AAP"))
})

test_that("quiet = TRUE suppresses the warning but not the attributes", {
  df <- tibble::tibble(ID = c("P1", "P2"), MSR = c("AAP", "BAD_CODE"), MSR_DEN = 10, MSR_NUM = 5)
  expect_no_warning(out <- filter_to_msr_code(df, "AAP", quiet = TRUE))
  expect_equal(attr(out, "n_excluded"), 1L)
  expect_equal(attr(out, "excluded_msr_values"), "BAD_CODE")
})

test_that("a missing MSR column aborts", {
  df <- tibble::tibble(ID = "P1", MSR_DEN = 10, MSR_NUM = 5)
  expect_error(filter_to_msr_code(df, "AAP"), "MSR")
})

test_that("an invalid msr_code argument aborts: length, type, NA", {
  df <- tibble::tibble(ID = "P1", MSR = "AAP", MSR_DEN = 10, MSR_NUM = 5)
  expect_error(filter_to_msr_code(df, c("AAP", "OTHER")), "single")
  expect_error(filter_to_msr_code(df, 1L), "character")
  expect_error(filter_to_msr_code(df, NA_character_), "single")
})

test_that("other columns and row content are preserved for kept rows", {
  df <- tibble::tibble(
    ID = c("P1", "P2"), MSR = c("AAP", "OTHER"),
    MSR_YR = c("RY2020-01", "RY2020-02"), MSR_DEN = c(10, 20), MSR_NUM = c(5, 15)
  )
  out <- suppressWarnings(filter_to_msr_code(df, "AAP"))
  expect_equal(out$MSR_YR, "RY2020-01")
  expect_equal(out$MSR_DEN, 10)
  expect_equal(out$MSR_NUM, 5)
})

## ---- integration: uses application/_config.R's real AAP_MSR_CODE ------------

test_that("AAP_MSR_CODE (application/_config.R) is the literal 'AAP'", {
  expect_equal(AAP_MSR_CODE, "AAP")
  ## And it is NOT OPEN_DECISIONS$msr_code$value, which is an annotation, not
  ## a filter value -- filtering on that whole string would match zero rows.
  expect_false(identical(AAP_MSR_CODE, OPEN_DECISIONS$msr_code$value))
})

test_that("filtering to OPEN_DECISIONS$msr_code$value directly would be a real bug", {
  ## Demonstrates why AAP_MSR_CODE exists as a separate literal: using the
  ## decision's own annotated value string as a filter matches nothing.
  df <- tibble::tibble(ID = c("P1", "P2"), MSR = c("AAP", "AAP"), MSR_DEN = 10, MSR_NUM = 5)
  wrong <- suppressWarnings(filter_to_msr_code(df, OPEN_DECISIONS$msr_code$value))
  expect_equal(nrow(wrong), 0L)

  right <- filter_to_msr_code(df, AAP_MSR_CODE)
  expect_equal(nrow(right), 2L)
})
