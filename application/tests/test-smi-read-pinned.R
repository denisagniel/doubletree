## ============================================================================
## application/tests/test-smi-read-pinned.R
##
## Behavioral tests for smi_read_pinned() (application/_config.R): the wrapper
## every pipeline read goes through instead of calling smidata::smi_read()
## directly (see test-no-unpinned-reads.R for the static enforcement side of
## that contract).
##
## THE BUG THIS HELPER EXISTS TO PREVENT. smidata::smi_read() ABORTS the
## moment it is given a non-NULL `ingest_id` while smi_env() == "local" (the
## ingest cache is server-only). Before smi_read_pinned() existed, setting
## config_ingest_id to a real value ahead of a server run would have broken
## every one of this suite's local, fixture-backed tests the next time
## anyone ran them on a laptop. smi_read_pinned() drops the pin on the local
## branch specifically so that setting config_ingest_id can never do that.
##
## Both tests below require smidata itself (skip_if_not_installed) and
## target the LOCAL branch specifically (skip_if_not on smi_env()) -- per
## AGENTS.md, SMI_ALLOW_PROVISIONAL_READS must never be set in checked-in
## code, so nothing here exercises the server branch at all.
## ============================================================================

test_that("smidata::smi_read() aborts on a pinned local read -- the bug smi_read_pinned() prevents", {
  skip_if_not_installed("smidata")
  skip_if_not(identical(smidata::smi_env(), "local"), "targets the local branch only")

  expect_error(
    smidata::smi_read(SMI_KEYS$covariates, ingest_id = "2026-09-18T142233Z_0001"),
    "session is local"
  )
})

test_that("smi_read_pinned() drops ingest_id on the local branch, so a set config_ingest_id does not abort", {
  skip_if_not_installed("smidata")
  skip_if_not(identical(smidata::smi_env(), "local"), "targets the local branch only")

  ## config_ingest_id lives in globalenv() (sourced from _config.R by
  ## helper-toy.R, `source()`'s default local = FALSE) -- assign()/on.exit()
  ## here, matching _config.R's own echo_open_decisions() idiom for mutating
  ## a global, rather than relying on `<<-`'s implicit environment search.
  old_ingest_id <- config_ingest_id
  on.exit(assign("config_ingest_id", old_ingest_id, envir = globalenv()), add = TRUE)
  assign("config_ingest_id", "2026-09-18T142233Z_0001", envir = globalenv())  # non-NULL placeholder

  result <- NULL
  expect_no_error(
    result <- smi_read_pinned(SMI_KEYS$covariates, columns = SMI_COLS$covariates)
  )

  ## "returns a fixture": a real data frame with the requested columns, not
  ## just "didn't error".
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)
  expect_named(result, SMI_COLS$covariates, ignore.order = TRUE)
})
