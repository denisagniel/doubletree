## ============================================================================
## application/helpers/covariate_blocks.R
##
## AUTHORITATIVE, not a candidate. PI-CONFIRMED 2026-09-17
## (smidata inst/analyses/doubletree__application.yml,
##  covariates.grid_exact_sparsity, status: confirmed).
##
## Contrast dual-bounds' analysis/real_data/helpers/covariate_blocks.R, which is
## a NON-authoritative candidate partition of all 72 columns, registered as an
## open decision. This file is the opposite kind of object: a SELECTION of 15
## covariates, deliberately small, confirmed, and not open. The two projects want
## different things from the same table:
##
##   dual-bounds  needs every column assigned to some block, because its
##                leave-block-out calibration benchmarks omit whole blocks in
##                turn -- so a dropped column is a silently missing confounder.
##   doubletree   needs FEW covariates, because estimate_att()'s identifying
##                assumption is grid-exact sparsity: a tree of at most
##                `leaf_budget` leaves must represent both nuisances EXACTLY on
##                the analyst's pre-specified grid. Every extra binary column
##                doubles the grid. A 72-column binary grid has 2^72 atoms and
##                no 4-leaf tree represents anything on it exactly.
##
## So this is not dual-bounds' list with columns dropped by oversight. Selecting
## 3 covariates per clinical block IS the design decision that makes the
## flagship estimator's assumption plausible at all, and it is why the two
## papers are a methodological contrast rather than a duplicate.
##
## Column names are VERBATIM from the real fingerprint (smidata
## inst/contracts/2026-09-15_2824080e/manifest.json, larger_smi_covariates:
## 381,018 rows x 72 columns). 55 of those 72 columns are already binary (_YN
## suffix), which is why 12 of the 15 need no construction at all.
##
## prior_cost's 3 dummies are NOT in this file's blocks: they do not exist in
## larger_smi_covariates (the table has zero dollar fields, checked against the
## real 72-column list 2026-09-17). prior_cost is constructed from the
## cost-claims family over [INDEX_DT - 12mo, INDEX_DT) in 03_cost_windows.R and
## then discretized in helpers/discretize_prior_cost.R.
## ============================================================================

#' The 12 confirmed table-resident covariates, by clinical block.
#'
#' Four blocks, 3 columns each, all drawn from `larger_smi_covariates`'s own
#' already-binary `_YN` columns. PI-confirmed 2026-09-17.
#'
#' The blocks are a documentation and review device -- they make the selection
#' auditable as "3 per clinical domain" rather than "12 columns someone liked".
#' Unlike dual-bounds', they are NOT consumed by the estimator: neither
#' `doubletree::estimate_att()` nor `doubletree::estimate_att_crossfit()` takes
#' a block/subset argument.
#'
#' @return Named list of four character vectors, 3 elements each.
covariate_blocks <- function() {
  list(
    demographics = c(
      ## Z-code social-determinant flags: housing instability, food insecurity,
      ## poverty. Chosen over AGE / RACE_ETHNICITY / SEX precisely because these
      ## three are already binary; the others would each need their own
      ## analyst-declared cutpoints or dummy expansion, spending grid atoms this
      ## estimator cannot afford.
      "HOUSING_YN", "FOOD_YN", "POVERTY_YN"
    ),
    prior_utilization = c(
      ## Any mental-health inpatient stay, any mental-health ED visit, any
      ## medical inpatient stay in the pre-index window. The paired _LOS / _CNT
      ## columns are counts, deliberately not selected.
      "MH_IP_YN", "MH_ER_YN", "MED_IP_YN"
    ),
    ap_psych_history = c(
      ## Schizophrenia diagnosis, any antipsychotic, any antidepressant.
      ##
      ## VAR9E_ANTIPSY_YN is PRE-index antipsychotic use and the exposure is
      ## POST-index AAP achievement, so this is a confounder, not the exposure.
      ## That distinction rests on the measurement windows being disjoint;
      ## 90_checks_tier1.R specifies the Tier-1 check that confirms it rather
      ## than assuming it.
      "VAR5D_SZ_YN", "VAR9E_ANTIPSY_YN", "VAR9G_ANTIDEP_YN"
    ),
    comorbidity = c(
      ## Type 2 diabetes, hypertension, and any glucose-or-lipid screening.
      ## GLU_OR_LIPID_YN is preferred over its GLUCOSE_YN / LIPID_YN /
      ## GLU_AND_LIPID_YN siblings: the four are mechanically related, so taking
      ## more than one spends grid atoms on near-duplicate information.
      "VAR5G_DM2_YN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN"
    )
  )
}

#' The 3 constructed `prior_cost` dummy names, in declared order.
#'
#' Q1 is the implicit reference category, so quartiles become 3 dummies, not 4.
#' `doubletree::estimate_att()` requires each covariate to be individually
#' binary, so a k-level discretization must be (k-1) indicators rather than one
#' multi-valued column.
#'
#' @return Character vector of length 3.
prior_cost_dummy_names <- function() {
  c("prior_cost_q2", "prior_cost_q3", "prior_cost_q4")
}

#' The full 15-column design-matrix specification, in fixed declared order.
#'
#' The ORDER is part of the specification, not an implementation detail: the
#' fitted trees' split indices, and every table reporting them, are read against
#' this order. Reordering it silently reinterprets a saved fit.
#'
#' @return Character vector of length 15.
design_matrix_columns <- function() {
  c(unlist(covariate_blocks(), use.names = FALSE), prior_cost_dummy_names())
}

#' Verify the confirmed selection against the source table's column list.
#'
#' A membership check, not a data check: it runs on the column-name vector, so
#' it works with no data and no smidata present. It catches the two failure
#' modes a hand-transcribed selection is most likely to have -- a typo in a
#' column name, and the same column selected into two blocks.
#'
#' It deliberately does NOT require the selection to partition the source table:
#' 60 of the 72 columns are intentionally unselected (see this file's header).
#'
#' @param all_columns Character vector of the source table's columns. Default
#'   the confirmed 72.
#' @param blocks Named list of blocks. Default `covariate_blocks()`.
#' @return `invisible(TRUE)`, or aborts.
assert_blocks_confirmed <- function(all_columns = larger_smi_covariates_columns(),
                                    blocks = covariate_blocks()) {
  selected <- unlist(blocks, use.names = FALSE)

  dupes <- selected[duplicated(selected)]
  if (length(dupes) > 0L) {
    cli::cli_abort(
      "Column(s) selected into >1 block:
       {.field {paste(unique(dupes), collapse = ', ')}}."
    )
  }

  unknown <- setdiff(selected, all_columns)
  if (length(unknown) > 0L) {
    cli::cli_abort(c(
      "Selected column(s) absent from the source table:
       {.field {paste(unknown, collapse = ', ')}}.",
      "i" = "Column names come from smidata snapshot 2026-09-15_2824080e."
    ))
  }

  if (length(selected) != 12L) {
    cli::cli_abort(
      "Expected 12 table-resident covariates (4 blocks x 3); got {length(selected)}."
    )
  }

  n_per_block <- vapply(blocks, length, integer(1))
  if (any(n_per_block != 3L)) {
    cli::cli_abort(
      "Every block must have 3 columns; these do not:
       {.field {paste(names(blocks)[n_per_block != 3L], collapse = ', ')}}."
    )
  }

  invisible(TRUE)
}

#' The confirmed 72-column list of `larger_smi_covariates`.
#'
#' Verbatim from smidata snapshot 2026-09-15_2824080e. Held here so
#' `assert_blocks_confirmed()` can run with no smidata and no data present.
#'
#' @return Character vector of length 72.
larger_smi_covariates_columns <- function() {
  c(
    "ID", "COHORT", "INDEX_DT", "AID_CATEGORY_NEW", "RAND_NEW_RACE_ETHNICITY",
    "RECIPIENT_SEX_1210", "AGE", "MH_IP_YN", "MH_IP_LOS", "MH_ER_YN", "MH_ER_CNT",
    "SUICIDE_YN", "HOUSING_YN", "FOOD_YN", "POVERTY_YN", "OTHER_Z_YN", "OP_MH_CNT",
    "OP_MED_CNT", "MED_IP_YN", "MED_IP_LOS", "MED_ER_YN", "MED_ER_CNT",
    "GLUCOSE_YN", "LIPID_YN", "GLU_AND_LIPID_YN", "GLU_OR_LIPID_YN",
    "MBR_LVL_RES_POSTAL_CD", "MBR_LVL_RES_COUNTY_CD", "COUNTY_LABEL",
    "VAR5D_SZ_YN", "VAR5D_BD_YN", "VAR5D_OTH_PSY_YN", "VAR5D_NON_ACUTE_YN",
    "VAR5D_SUD_YN", "VAR5E_DEMENTIA_YN", "VAR5F_NA_ORG_YN", "VAR5G_DYSLIPID_YN",
    "VAR5G_HTN_YN", "VAR5G_DM2_YN", "VAR5G_CAD_YN", "VAR5G_HF_YN", "VAR5G_CEREB_YN",
    "VAR5G_OTH_VASC_YN", "VAR5G_REMOTE_CVD_YN", "VAR5G_META_SYN_YN", "VAR5G_OBES_YN",
    "VAR5G_PRE_DM_YN", "VAR5G_GES_DM_YN", "VAR5G_TOBACCO_YN", "VAR5G_PERSON_FAM_YN",
    "VAR5G_META_RISK_YN", "VAR5G_T1DM_YN", "VAR5G_OTH_DM_YN", "VAR5G_OTH_DM1_YN",
    "VAR5G_AP_YN", "VAR5G_CHRONIC_COM_YN", "VAR5G_NONDM_MORBID_YN",
    "VAR9A_ANTIDM_ALL_YN", "VAR9B_ANTILIPEMIC_YN", "VAR9C_CVD_HTN_YN",
    "VAR9D_WEIGHT_YN", "VAR9E_ANTIPSY_YN", "VAR9F_MOOD_YN", "VAR9G_ANTIDEP_YN",
    "VAR9H_OTH_PSY_YN", "VAR9I_OTHER_YN", "VAR10A_PSY_YN", "VAR10B_NON_PSY_YN",
    "PRE_START", "VAR5D_MDD_SEV_YN", "VAR5D_OTH_MDD_YN", "VAR5D_PSY_COM_YN"
  )
}
