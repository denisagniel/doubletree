## ============================================================================
## application/tests/test-config-fingerprint.R
##
## config_fingerprint() is what lets 07 refuse an analytic parquet built under
## different rules, and effective_ingest_id() is the pin that is actually applied.
## ============================================================================

test_that("config_fingerprint() is stable and sensitive to a rule", {
  f1 <- config_fingerprint()
  expect_type(f1, "character")
  expect_length(f1, 1L)
  expect_identical(config_fingerprint(), f1)

  ## Changing a rule that decides which numbers 02-06 produce changes the hash.
  old <- config_negative_amount_policy
  withr::defer(assign("config_negative_amount_policy", old, envir = globalenv()))
  assign("config_negative_amount_policy", "ignore", envir = globalenv())
  expect_false(identical(config_fingerprint(), f1))
})

test_that("config_fingerprint() notices an open decision changing status", {
  f1 <- config_fingerprint()
  old <- OPEN_DECISIONS
  withr::defer(assign("OPEN_DECISIONS", old, envir = globalenv()))
  od <- OPEN_DECISIONS
  od$prior_cost_negative_floor$status <- "confirmed"
  assign("OPEN_DECISIONS", od, envir = globalenv())
  expect_false(identical(config_fingerprint(), f1))
})

test_that("the prior_cost floor is an open decision read with a warning", {
  expect_equal(OPEN_DECISIONS$prior_cost_negative_floor$status, "open_decision")
  expect_false(OPEN_DECISIONS$prior_cost_negative_floor$blocking_final)
  expect_warning(v <- open_value("prior_cost_negative_floor"), "Unconfirmed")
  expect_equal(v, 0)
  expect_error(confirmed_value("prior_cost_negative_floor"), "not confirmed")
})

test_that("effective_ingest_id() is NA locally, whatever is configured", {
  skip_if_not_installed("smidata")
  skip_if_not(identical(smidata::smi_env(), "local"))
  expect_true(is.na(effective_ingest_id()))
})
