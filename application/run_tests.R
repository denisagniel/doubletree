## ============================================================================
## application/run_tests.R
##
## Run the application pipeline's own test suite.
##
##   Rscript application/run_tests.R
##
## SEPARATE FROM tests/testthat/ ON PURPOSE. That directory is the PACKAGE's
## suite and runs under R CMD check, so everything it needs must appear in
## DESCRIPTION. This directory needs dplyr, tibble, cli and (for
## 01_declare_requirements.R) smidata, none of which the package itself depends
## on. Keeping the two apart is what lets the application pipeline use them
## without adding a single line to the package's Imports -- see
## application/README.md.
##
## application/ is .Rbuildignore'd, so nothing here is inside R CMD check's
## scope.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."

for (pkg in c("testthat", "tibble", "dplyr", "cli")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("application/tests needs package '", pkg, "', which is not installed.")
  }
}

## doubletree itself must be loadable: test-estimate-att-smoke.R and one negative
## test in test-design-matrix.R call the estimators. They skip rather than fail if
## it is absent, so report the state explicitly rather than letting a silent skip
## look like a pass.
if (!requireNamespace("doubletree", quietly = TRUE)) {
  cli::cli_alert_warning(
    "{.pkg doubletree} is not installed; the estimator tests will SKIP.
     Run {.code devtools::install()} (or {.code devtools::load_all()}) first for a
     meaningful run."
  )
}

results <- testthat::test_dir(
  file.path(app_dir, "tests"),
  reporter = "summary",
  stop_on_failure = FALSE
)

df <- as.data.frame(results)
n_fail <- sum(df$failed) + sum(df$error)
n_skip <- sum(df$skipped)
n_pass <- sum(df$passed)

cli::cli_h1("application/tests summary")
cli::cli_inform("Passed: {n_pass}   Failed/errored: {n_fail}   Skipped: {n_skip}")

## Non-zero exit only under Rscript. The literal guard, written inline: sourcing
## this file interactively must not kill the session (r-interactive-entry F8).
if (n_fail > 0L && !interactive() && sys.nframe() == 0L) {
  quit(status = 1L)
}
if (n_fail > 0L) {
  cli::cli_abort("{n_fail} application test{?s} failed or errored.")
}
