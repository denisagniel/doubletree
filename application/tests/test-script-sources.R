## ============================================================================
## application/tests/test-script-sources.R
##
## Each numbered script sources the helpers it calls. Moving a function between
## helper files breaks a script that is not exercised by the test suite (they
## need the real data), so check the source() lines themselves.
## ============================================================================

scripts <- list.files(app_dir_for_tests, pattern = "^[0-9]{2}[a-z_].*[.]R$", full.names = TRUE)
sourced <- function(f) {
  src <- readLines(f, warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  m <- regmatches(src, regexpr('"helpers", "[A-Za-z_]+[.]R"', src))
  sub('.*"([A-Za-z_]+[.]R)"$', "\\1", m)
}

test_that("a script that sources decision_report.R also sources disclosure.R first", {
  for (f in scripts) {
    s <- sourced(f)
    if ("decision_report.R" %in% s) {
      expect_true("disclosure.R" %in% s, info = basename(f))
      expect_lt(match("disclosure.R", s), match("decision_report.R", s), label = basename(f))
    }
  }
})

test_that("a script that sources design_matrix.R also sources covariate_blocks.R", {
  ## design_matrix_columns() lives in covariate_blocks.R.
  for (f in scripts) {
    s <- sourced(f)
    if ("design_matrix.R" %in% s) expect_true("covariate_blocks.R" %in% s, info = basename(f))
  }
})

test_that("a script that calls design_matrix_columns() sources covariate_blocks.R", {
  for (f in scripts) {
    body <- readLines(f, warn = FALSE)
    body <- body[!grepl("^\\s*#", body)]
    if (any(grepl("design_matrix_columns\\(", body))) {
      expect_true("covariate_blocks.R" %in% sourced(f), info = basename(f))
    }
  }
})
