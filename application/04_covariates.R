## ============================================================================
## 04_covariates.R
##
## Assemble X: the 15-column, all-binary design matrix.
##
##   12 columns  read straight off larger_smi_covariates, already binary (_YN),
##               no construction of any kind. PI-CONFIRMED 2026-09-17.
##    3 columns  prior_cost's quartile dummies, discretized here from the dollar
##               amounts 03_cost_windows.R computed.
##
## Design (smidata inst/analyses/doubletree__application.yml):
##   covariates.source              status: inherited (from dual-bounds)
##   covariates.grid_exact_sparsity status: CONFIRMED, PI 2026-09-17 -- this
##     paper's own decision, and the one place its covariate handling genuinely
##     departs from dual-bounds'. dual-bounds uses ~68 covariates across four
##     blocks because marbounds does not require binary X; this paper uses 15
##     because estimate_att()'s grid-exact sparsity assumption is stated over the
##     grid the covariates generate, and a 68-column binary grid has 2^68 atoms.
##
## THE SELECTION IS NOT AN OPEN DECISION HERE. Contrast dual-bounds'
## 04_covariates.R, which reads open_value("covariate_block_map") and falls back
## to a candidate. There is no fallback in this script because there is nothing
## to fall back from: the 15 columns are confirmed, so they are asserted.
##
## REAL AS OF 2026-09-24: covariates_raw/covariates/the preview design matrix
## are now live (cascade-guarded on 03, which cascades to 02), not commented
## pseudocode. Explicitly a PREVIEW, not the authoritative sample -- 06 owns
## the real, final-sample discretization; see the paragraph inline for why
## the two are not interchangeable.
##
## Run time: seconds locally against Tier-0 fixtures; on the server, seconds
## once 03's expensive cost-window pass has completed (this script joins
## already-built tables, spending no time of its own).
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "discretize_prior_cost.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("04 -- covariates")

## ---- the confirmed selection, verified against the source column list ------
##
## This runs for real: assert_blocks_confirmed() works on the column-name
## vector, so it needs no data and no smidata. It catches a typo in a
## transcribed column name and a column selected into two blocks.

assert_blocks_confirmed()
blocks <- covariate_blocks()
cli::cli_alert_success(
  "All {length(unlist(blocks, use.names = FALSE))} confirmed {.field _YN}
   covariates are present in the {length(larger_smi_covariates_columns())}-column
   source list."
)
for (b in names(blocks)) {
  cli::cli_alert_info("Block {.field {b}}: {.field {blocks[[b]]}}")
}
cli::cli_alert_info(
  "Plus {.field {prior_cost_dummy_names()}} (constructed; Q1 is the implicit
   reference)."
)

x_cols <- design_matrix_columns()
cli::cli_alert_info(
  "X has {length(x_cols)} binary column{?s}, so the declared grid has
   2^{length(x_cols)} atoms and {.code leaf_budget = {config_leaf_budget}} leaves
   must represent both nuisances exactly on it. That is the assumption, stated in
   the units it is actually about."
)

## ---- covariates + prior_cost: REAL as of 2026-09-24 -------------------------
##
## covariates_raw/covariates/prior_cost_dummies/X_with_id here are a PREVIEW,
## not the authoritative analytic sample -- see the paragraph below on why
## quartiles get recomputed in 06. Cascade-guarded on 03 (which itself
## cascades to 02), matching 02/03's `require_stage()` convention.

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; skipping the real covariate assembly
     entirely (same convention as 01/02/03)."
  )
} else {
  require_stage(c("analysis_cohort", "prior_cost"), "03_cost_windows.R")

  covariates_raw <- smi_read_pinned(
    SMI_KEYS$covariates,
    columns = c("ID", "INDEX_DT", unlist(blocks, use.names = FALSE))
  )
  ## smi_read_pinned() (i.e. smidata::smi_read()) already aborts on a missing
  ## requested column, so there is nothing to re-check here; a second probe
  ## would just be noise.

  ## The 12 confirmed _YN columns, ID-keyed, no INDEX_DT or prior_cost -- the
  ## explicit interface 06_assemble_analytic_data.R's own (pre-existing)
  ## comments already name and expect.
  covariates_yn <- dplyr::select(covariates_raw, "ID", dplyr::all_of(unlist(blocks, use.names = FALSE)))

  ## prior_cost from 03_cost_windows.R -- ALREADY one row per cohort patient,
  ## with the coverage-aware NA-vs-0 distinction resolved there (03's header
  ## explains why: absence of a claim under FULL file coverage is a real
  ## zero; absence/truncation under PARTIAL coverage is unknowable and left
  ## NA, not coalesced). No coalesce() here -- that would silently undo 03's
  ## fix and reintroduce the exact truncation bug this file's own comments
  ## originally warned about.
  prior_cost_x <- dplyr::transmute(prior_cost, ID = .data$ID, prior_cost = .data$cost)

  covariates <- covariates_raw |>
    dplyr::semi_join(analysis_cohort, by = "ID") |>   # eligible patients only
    dplyr::left_join(prior_cost_x, by = "ID")

  n_prior_cost_na <- sum(is.na(covariates$prior_cost))
  if (n_prior_cost_na > 0L) {
    cli::cli_alert_warning(
      "{n_prior_cost_na} of {nrow(covariates)} patients have {.val NA}
       prior_cost (file-coverage truncation, per 03_cost_windows.R) --
       discretize_prior_cost() below will ABORT on this, by design: an
       {.val NA} silently binned into a quartile is a fabricated covariate
       value. This is the decision 04's own header names -- do not coalesce
       past it here."
    )
  }

  ## WHICH SAMPLE THE QUARTILES ARE COMPUTED WITHIN IS A DECISION, MADE IN 06.
  ## The confirmed design says "computed within the analytic sample". The
  ## analytic sample is the AAP-eligible cohort AFTER the complete-case
  ## restriction of 05 -- so strictly, discretization belongs after that
  ## filter, not before it. Run here anyway, for a PREVIEW only (labelled as
  ## such below): discretizing on the pre-filter cohort and then subsetting
  ## gives different boundaries than discretizing the filtered sample, and
  ## the two are not interchangeable -- 06 MUST re-discretize on its own,
  ## final sample; this call does not gate or feed anything downstream.
  if (n_prior_cost_na == 0L) {
    prior_cost_dummies_preview <- discretize_prior_cost(covariates$prior_cost)
    cli::cli_alert_info(
      "PREVIEW quartile cutpoints, pre-filter cohort (NOT the reported
       cutpoints -- 06 recomputes these on the final, complete-case-filtered
       sample): {.val {round(attr(prior_cost_dummies_preview, 'cutpoints'), 2)}}"
    )

    x_with_id_preview <- assemble_design_matrix(covariates, prior_cost_dummies_preview)
    ## The guard the package does not give you. See helpers/design_matrix.R's
    ## header: check_att_data() does not examine X's coding at all, and
    ## estimate_att_crossfit() adds no further check -- so a "Y"/"N" or 1/2
    ## coded _YN column reaches optimaltrees and returns a number. Checked
    ## HERE, at preview time, not at estimation time hours later -- and again
    ## in 06 on the real sample, since a preview pass does not guarantee the
    ## final one.
    assert_binary_design_matrix(dplyr::select(x_with_id_preview, -"ID"))
    cli::cli_alert_success(
      "Preview design matrix ({nrow(x_with_id_preview)} pre-filter patients)
       passes assert_binary_design_matrix()."
    )
  } else {
    cli::cli_alert_info(
      "Skipping the preview discretize/assemble/assert steps -- see the
       {.val NA} prior_cost warning above."
    )
  }
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info(
    "04_covariates.R complete
     ({if (config_has_smidata) 'ran the real preview assembly' else 'structure only, no smidata'})."
  )
}
