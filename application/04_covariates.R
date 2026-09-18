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
## Run time: seconds. This script joins two already-built tables; the expensive
## work happened in 03.
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

## ---- INTENDED PIPELINE (needs real data) -----------------------------------
##
## covariates_raw <- smidata::smi_read(
##   SMI_KEYS$covariates,
##   columns = c("ID", "INDEX_DT", unlist(blocks, use.names = FALSE))
## )
## ## smidata::smi_read() already aborts on a missing requested column, so there
## ## is nothing to re-check here; a second probe would just be noise.
##
## ## prior_cost from 03_cost_windows.R's cost_by_window, "prior_cost" slice.
## ## A patient with NO pre-index claims has prior_cost 0, not NA -- absence of a
## ## claim is a real zero here, unlike absence of an enrollment record. Coding it
## ## as 0 is a substantive statement, made explicitly:
## prior_cost_x <- prior_cost |>
##   dplyr::transmute(ID = .data$id, prior_cost = .data$cost)
##
## covariates <- covariates_raw |>
##   dplyr::semi_join(analysis_cohort, by = "ID") |>   # eligible patients only
##   dplyr::left_join(prior_cost_x, by = "ID") |>
##   dplyr::mutate(prior_cost = dplyr::coalesce(.data$prior_cost, 0))
##
## ## WHICH SAMPLE THE QUARTILES ARE COMPUTED WITHIN IS A DECISION, AND THIS IS
## ## WHERE IT IS MADE. The confirmed design says "computed within the analytic
## ## sample". The analytic sample is the AAP-eligible cohort AFTER the
## ## complete-case restriction of 05 -- so strictly, discretization belongs after
## ## that filter, not before it. It is written here for readability and MUST be
## ## re-applied in 06 on the final sample; discretizing on the pre-filter cohort
## ## and then subsetting gives different boundaries than discretizing the filtered
## ## sample, and the two are not interchangeable.
## prior_cost_dummies <- discretize_prior_cost(covariates$prior_cost)
## cli::cli_alert_info(
##   "Realized quartile cutpoints (REPORT THESE in the methods section):
##    {.val {attr(prior_cost_dummies, 'cutpoints')}}"
## )
##
## X_with_id <- assemble_design_matrix(covariates, prior_cost_dummies)
##
## ## The guard the package does not give you. See helpers/design_matrix.R's
## ## header: check_att_data() does not examine X's coding at all, and
## ## estimate_att_crossfit() adds no further check -- so a "Y"/"N" or 1/2 coded
## ## _YN column reaches optimaltrees and returns a number. Checked HERE, at
## ## assembly time, not at estimation time hours later.
## assert_binary_design_matrix(dplyr::select(X_with_id, -"ID"))

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("04_covariates.R complete (structure only; no data read).")
}
