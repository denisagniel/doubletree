## ============================================================================
## application/tests/test-output-layout.R
## output_file(): export/ and private/ are separate, created on demand, marked.
## ============================================================================

test_that("output_file separates export/ from private/ and creates folders", {
  root <- withr::local_tempdir()
  e <- output_file("export", "a.csv", app_dir = root)
  p <- output_file("private", "b.rds", app_dir = root)
  expect_equal(normalizePath(dirname(e)), normalizePath(file.path(root, "output", "export")))
  expect_equal(normalizePath(dirname(p)), normalizePath(file.path(root, "output", "private")))
  expect_true(dir.exists(dirname(e)) && dir.exists(dirname(p)))
  ## Files are not created, only their folders.
  expect_false(file.exists(e))
})

test_that("output_file handles subfolders and writes a marker in each area", {
  root <- withr::local_tempdir()
  r <- output_file("export", "decision_report", "decision_report.md", app_dir = root)
  expect_true(dir.exists(file.path(root, "output", "export", "decision_report")))
  output_file("private", "x.rds", app_dir = root)
  expect_match(readLines(file.path(root, "output", "export", "_READ_ME.txt"))[1], "^EXPORT")
  expect_match(readLines(file.path(root, "output", "private", "_READ_ME.txt"))[1], "^PRIVATE")
  expect_match(paste(readLines(file.path(root, "output", "private", "_READ_ME.txt")), collapse = " "),
               "NEVER copy")
})

test_that("output_file rejects an unknown area", {
  expect_error(output_file("elsewhere", "a.csv", app_dir = withr::local_tempdir()), "should be one of")
})

test_that("no script writes patient-level files into export/", {
  scripts <- list.files(app_dir_for_tests, pattern = "^[0-9]{2}[a-z_].*[.]R$", full.names = TRUE)
  patient_level <- "[.]rds|analytic_cohort(_PROVISIONAL)?[.]parquet|_manifest[.]csv"
  for (f in scripts) {
    src <- readLines(f, warn = FALSE)
    src <- src[!grepl("^\\s*#", src)]
    hits <- src[grepl('output_file\\("export"', src) & grepl(patient_level, src)]
    expect_length(hits, 0L)
  }
})
