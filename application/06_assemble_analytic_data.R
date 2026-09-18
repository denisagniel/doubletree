## ============================================================================
## 06_assemble_analytic_data.R
##
## Join stages 02-05 into the one-row-per-patient analytic dataset, validate the
## design matrix, then gate the write on the open decisions.
##
## THIS SCRIPT IS EXPECTED TO ABORT TODAY. That is the design, not a bug:
## assert_no_blocking_open() refuses to write a final artifact while four
## blocking decisions are unresolved, and complete_case is NA for every patient
## (see 05_complete_case.R). An analytic dataset that wrote successfully right now
## would be one whose exposure is nondeterministic (msr_code), whose cost may be
## double-counted in BOTH Y and prior_cost's quartile boundaries
## (cost_family_scope), and whose inclusion criterion is fabricated
## (enrollment_source). The abort is the only correct outcome.
##
## THE RELEASE VALVE, if a number is wanted before the decisions land, is the
## _PROVISIONAL artifact named below -- NOT relaxing the gate. A provisional file
## is self-labelling and cannot be mistaken for a result; a relaxed gate is
## invisible in the output.
##
## Run time: seconds, once its inputs exist.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "discretize_prior_cost.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("06 -- assemble analytic data")

output_dir <- file.path(app_dir, "output")
output_stem <- "analytic_cohort"

cli::cli_alert_info(
  "This script is EXPECTED to stop at the gate below while blocking decisions
   remain."
)
cli::cli_alert_info("Blocking now: {.val {blocking_open_ids()}}")

## ---- INTENDED ASSEMBLY (needs 02-05's real outputs) ------------------------
##
## ## Step 1: the complete-case restriction, applied as a FILTER. This is the one
## ## structural difference from dual-bounds' 06, which keeps every patient and
## ## carries R as a column for its estimator to model. doubletree's estimators
## ## have no censoring argument, so the restriction happens here and the excluded
## ## patients are gone from the analytic data entirely -- which is exactly why the
## ## excluded COUNT must be reported.
## ##
## ## With complete_case all NA (today), this filter yields ZERO rows. That is the
## ## correct behaviour and not something to coalesce away: an unknown inclusion
## ## status is not an inclusion.
## analytic_pre <- analysis_cohort |>                         # 02
##   dplyr::select("ID", "INDEX_DT", "A", "elig_aap", "MSR_YR_selected", "month_gap") |>
##   dplyr::left_join(                                        # 03: Y
##     dplyr::transmute(outcome_cost, ID = .data$id, Y = .data$cost),
##     by = "ID"
##   ) |>
##   dplyr::left_join(                                        # 03: prior_cost ($)
##     dplyr::transmute(prior_cost, ID = .data$id, prior_cost = .data$cost),
##     by = "ID"
##   ) |>
##   dplyr::mutate(prior_cost = dplyr::coalesce(.data$prior_cost, 0)) |>
##   dplyr::left_join(covariates_yn, by = "ID") |>            # 04: the 12 _YN cols
##   dplyr::left_join(complete_case, by = "ID")               # 05: all NA today
##
## n_pre <- nrow(analytic_pre)
## analytic_cc <- dplyr::filter(analytic_pre, isTRUE(.data$complete_case))
## cli::cli_alert_info(
##   "Complete-case restriction: kept {nrow(analytic_cc)} of {n_pre}
##    ({round(100 * nrow(analytic_cc) / n_pre, 1)}%); excluded
##    {n_pre - nrow(analytic_cc)}."
## )
## if (nrow(analytic_cc) == 0L) {
##   cli::cli_abort(
##     "Complete-case restriction removed every patient. With
##      {.field complete_case} unresolved (05_complete_case.R) this is the expected
##      outcome, not a data problem."
##   )
## }
##
## ## Step 2: discretize prior_cost ON THE FINAL SAMPLE. Order matters and this is
## ## the reason the call lives here rather than only in 04: the confirmed design
## ## says quartiles are computed "within the analytic sample", and the analytic
## ## sample is the post-restriction one. Discretizing before the filter and then
## ## subsetting produces different boundaries -- boundaries partly determined by
## ## patients who are not in the analysis.
## prior_cost_dummies <- discretize_prior_cost(analytic_cc$prior_cost)
## cutpoints <- attr(prior_cost_dummies, "cutpoints")
## cli::cli_alert_info(
##   "Realized quartile cutpoints on the analytic sample (REPORT THESE):
##    {.val {cutpoints}}"
## )
##
## X_with_id <- assemble_design_matrix(analytic_cc, prior_cost_dummies)
## assert_binary_design_matrix(dplyr::select(X_with_id, -"ID"))
##
## analytic_data <- X_with_id |>
##   dplyr::left_join(
##     dplyr::select(analytic_cc, "ID", "A", "Y", "prior_cost", "MSR_YR_selected",
##                   "month_gap"),
##     by = "ID"
##   ) |>
##   dplyr::select("ID", "A", "Y", dplyr::all_of(design_matrix_columns()),
##                 "prior_cost", "MSR_YR_selected", "month_gap")
## ## prior_cost (dollars) is retained ALONGSIDE its dummies deliberately: it is
## ## not passed to the estimator, but the methods section needs its distribution
## ## and the cutpoints have to be reproducible from the saved data.
## attr(analytic_data, "prior_cost_cutpoints") <- cutpoints
##
## ## Assertions to run immediately after assembly:
## ##   * one row per patient
## ##   * A is 0/1 with no NA
## ##   * Y is numeric with no NA (guaranteed by the complete-case filter -- if it
## ##     is not, the filter and the outcome disagree and THAT is the finding)
## ##   * both arms non-empty, and n_control >= config_leaf_budget * m_n, since
## ##     estimate_att() fits its mu tree on controls alone
## ##   * nrow() matches the kept count reported above

## ---- the gate ---------------------------------------------------------------
##
## Reachable, and will not execute today. Written this way rather than as a
## comment so that resolving the decisions makes the write happen with no
## further edit.

if (has_blocking_open()) {
  cli::cli_alert_warning(
    "Any artifact written now would carry the {.field _PROVISIONAL} suffix:
     {.file {file.path(output_dir, paste0(output_stem, '_PROVISIONAL.parquet'))}}"
  )
  cli::cli_alert_warning(
    "That suffix -- not a loosened gate -- is the release valve if a number is
     needed before the decisions land."
  )
}

## Aborts while anything blocks. Everything after this line is unreachable
## today, by design.
assert_no_blocking_open()

fs::dir_create(output_dir)
out_path <- file.path(output_dir, paste0(output_stem, ".parquet"))
arrow::write_parquet(analytic_data, out_path)
cli::cli_alert_success("Wrote {.file {out_path}} ({nrow(analytic_data)} rows).")

## The cutpoints are an ATTRIBUTE, and parquet does not carry R attributes.
## Written as a sidecar so the grid the sparsity assumption is stated over is
## recoverable from disk rather than only from this session.
cutpoint_path <- file.path(output_dir, paste0(output_stem, "_prior_cost_cutpoints.csv"))
utils::write.csv(
  data.frame(
    prob = names(attr(analytic_data, "prior_cost_cutpoints")),
    cutpoint = as.numeric(attr(analytic_data, "prior_cost_cutpoints"))
  ),
  cutpoint_path, row.names = FALSE
)
cli::cli_alert_success("Wrote {.file {cutpoint_path}}.")

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("06_assemble_analytic_data.R complete.")
}
