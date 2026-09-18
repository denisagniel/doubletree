## ============================================================================
## application/tests/test-covariate-blocks.R
##
## The PI-confirmed 15-column selection: presence, count, order, and binary
## coding in the toy frame.
##
## Every expected value is HAND-COUNTED from
## smidata::inst/analyses/doubletree__application.yml's
## covariates.grid_exact_sparsity field, not from running the code.
## ============================================================================

test_that("the selection is 4 blocks of 3, totalling 12 table-resident columns", {
  blocks <- covariate_blocks()

  ## Hand-counted from the registry: demographics 3, prior_utilization 3,
  ## ap_psych_history 3, comorbidity 3.
  expect_length(blocks, 4L)
  expect_named(
    blocks,
    c("demographics", "prior_utilization", "ap_psych_history", "comorbidity")
  )
  expect_equal(vapply(blocks, length, integer(1)), c(
    demographics = 3L, prior_utilization = 3L,
    ap_psych_history = 3L, comorbidity = 3L
  ))
  expect_length(unlist(blocks, use.names = FALSE), 12L)
})

test_that("the 12 columns are exactly the registry-confirmed names", {
  ## Transcribed from doubletree__application.yml, covariates.grid_exact_sparsity,
  ## status: confirmed, PI 2026-09-17. If this test fails, either the registry
  ## changed or helpers/covariate_blocks.R drifted from it -- check which before
  ## editing either.
  expect_equal(covariate_blocks()$demographics,
               c("HOUSING_YN", "FOOD_YN", "POVERTY_YN"))
  expect_equal(covariate_blocks()$prior_utilization,
               c("MH_IP_YN", "MH_ER_YN", "MED_IP_YN"))
  expect_equal(covariate_blocks()$ap_psych_history,
               c("VAR5D_SZ_YN", "VAR9E_ANTIPSY_YN", "VAR9G_ANTIDEP_YN"))
  expect_equal(covariate_blocks()$comorbidity,
               c("VAR5G_DM2_YN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN"))
})

test_that("all 12 selected columns exist in the confirmed 72-column source list", {
  all_cols <- larger_smi_covariates_columns()

  expect_length(all_cols, 72L)                   # smidata snapshot 2026-09-15_2824080e
  expect_false(anyDuplicated(all_cols) > 0L)
  expect_true(all(unlist(covariate_blocks(), use.names = FALSE) %in% all_cols))
  expect_true(isTRUE(assert_blocks_confirmed()))
})

test_that("the design matrix is 15 columns: 12 _YN then 3 prior_cost dummies", {
  cols <- design_matrix_columns()

  ## Hand-counted: 12 + 3 = 15. Q1 is the implicit reference, so quartiles give
  ## 3 dummies, not 4.
  expect_length(cols, 15L)
  expect_equal(cols[1:12], unlist(covariate_blocks(), use.names = FALSE))
  expect_equal(cols[13:15],
               c("prior_cost_q2", "prior_cost_q3", "prior_cost_q4"))
  expect_equal(prior_cost_dummy_names(), cols[13:15])
})

test_that("assert_blocks_confirmed() catches a typo'd column name", {
  expect_error(
    assert_blocks_confirmed(blocks = list(
      demographics = c("HOUSING_YN", "FOOD_YN", "POVERTY_YN"),
      prior_utilization = c("MH_IP_YN", "MH_ER_YN", "MED_IP_YN"),
      ap_psych_history = c("VAR5D_SZ_YN", "VAR9E_ANTIPSY_YN", "VAR9G_ANTIDEP_YN"),
      ## VAR5G_DM2_YNN -- one character off, the realistic transcription error.
      comorbidity = c("VAR5G_DM2_YNN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN")
    )),
    "absent from the source table"
  )
})

test_that("assert_blocks_confirmed() catches a column selected into two blocks", {
  expect_error(
    assert_blocks_confirmed(blocks = list(
      demographics = c("HOUSING_YN", "FOOD_YN", "POVERTY_YN"),
      prior_utilization = c("MH_IP_YN", "MH_ER_YN", "MED_IP_YN"),
      ap_psych_history = c("VAR5D_SZ_YN", "VAR9E_ANTIPSY_YN", "VAR9G_ANTIDEP_YN"),
      ## HOUSING_YN repeated from demographics.
      comorbidity = c("VAR5G_DM2_YN", "VAR5G_HTN_YN", "HOUSING_YN")
    )),
    "selected into >1 block"
  )
})

test_that("assert_blocks_confirmed() catches a block of the wrong size", {
  expect_error(
    assert_blocks_confirmed(blocks = list(
      demographics = c("HOUSING_YN", "FOOD_YN"),        # 2, not 3
      prior_utilization = c("MH_IP_YN", "MH_ER_YN", "MED_IP_YN"),
      ap_psych_history = c("VAR5D_SZ_YN", "VAR9E_ANTIPSY_YN", "VAR9G_ANTIDEP_YN"),
      comorbidity = c("VAR5G_DM2_YN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN")
    )),
    "Expected 12 table-resident covariates"
  )
})

test_that("all 12 selected columns are binary and varying in the toy frame", {
  cov <- toy_covariates(n = 64L)
  yn <- cov[, unlist(covariate_blocks(), use.names = FALSE), drop = FALSE]

  expect_true(all(vapply(yn, function(col) all(col %in% c(0L, 1L)), logical(1))))
  ## Every column must have BOTH levels: the longest period is 64, so n = 64 is
  ## the minimum at which this holds. A constant column would pass a naive
  ## "is it binary" check and still break the grid.
  expect_true(all(vapply(yn, function(col) length(unique(col)) == 2L, logical(1))))
})
