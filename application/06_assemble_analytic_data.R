## ============================================================================
## 06_assemble_analytic_data.R
##
## Join stages 02-05 into the one-row-per-patient analytic dataset, validate
## the design matrix, then write it -- with the PROVISIONAL suffix while any
## decision blocks, or without it once nothing does.
##
## REAL AS OF 2026-09-24: the "intended assembly" is now live code, cascade-
## guarded on 02-05 via `require_stage()`. Locally, the cascade correctly
## stops at 05_complete_case.R's own Tier-0 limitation (larger_smi_medicaid_
## monthly_flag's YEAR_MONTH has no real-data census levels, so the fixture
## cannot exercise compute_complete_case() -- see that file's header): this
## is an honest Tier-0 boundary, not a defect introduced here. On the server,
## the full chain reaches this script for real.
##
## THE PROVISIONAL RELEASE VALVE IS NOW REAL, NOT JUST DOCUMENTED. Until
## today it was pure aspiration: the write step below was preceded by an
## unconditional `assert_no_blocking_open()`, which always aborted while
## anything blocked -- so a `_PROVISIONAL.parquet` could never actually be
## written, contradicting this file's own header ("the release valve... is
## the _PROVISIONAL artifact"). Design-reviewed by `oracle` (session
## ses_f2b12eccdffel35nxc8t7UppvR) before fixing: `has_blocking_open()`'s own
## roxygen already said it exists "to choose a _PROVISIONAL output suffix" --
## the helper was built for this and never used. `assert_no_blocking_open()`
## moved OUT of this script entirely (kept only in 07_estimate_att.R, where a
## REPORTABLE estimate is produced): `diagnostics`, the one decision still
## `blocking_final = TRUE` as of this writing, is specifically about which
## ESTIMATOR result to trust in 07, not about whether the assembled DATA is
## correct -- so it should not block an inspectable intermediate artifact
## here, only the reportable one there. Every STRUCTURAL check (binary
## coding, no NA, both arms populated, nonzero rows) still aborts
## unconditionally, regardless of decision status -- a _PROVISIONAL file
## must not be a way to ship a broken matrix.
##
## Y AND prior_cost ARE NOW COVERAGE-AWARE, NOT BLANKET-COALESCED. Reading
## 03's real streaming output against a Tier-0 fixture (oracle, same
## session) found a real gap in this file's own original comment: it claimed
## "Y is numeric with no NA (guaranteed by the complete-case filter)" -- that
## guarantee does not hold. A patient can be complete-case (enrolled every
## month) yet have an UNKNOWABLE Y because their outcome window extends past
## the cost files' 2016-2024 coverage -- a DIFFERENT reason than genuine
## zero-utilization (Y = 0, a real value 03 already resolves). Coalescing
## the unknowable case to 0 would silently conflate them. This script
## excludes patients with unknowable Y (or prior_cost) from the analytic
## sample, reporting that count SEPARATELY from the complete-case exclusion
## count -- three different attrition reasons, not one.
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
  "Blocking now: {.val {blocking_open_ids()}} -- any write below will carry
   the {.field _PROVISIONAL} suffix while this is non-empty (see this file's
   header for why {.val diagnostics} does not block the write outright)."
)

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; skipping the real assembly entirely
     (same convention as 01-05)."
  )
} else {
  ## Each require_stage() call cascades to the script that creates its
  ## objects, which itself cascades further upstream as needed -- see
  ## _config.R's require_stage() for why this is safe to call redundantly.
  require_stage(c("outcome_cost", "prior_cost"), "03_cost_windows.R")
  require_stage("covariates_yn", "04_covariates.R")
  require_stage("complete_case", "05_complete_case.R")

  ## ---- Step 1: join, THEN the complete-case restriction as a FILTER -------
  ##
  ## This is the one structural difference from dual-bounds' 06, which keeps
  ## every patient and carries R as a column for its estimator to model.
  ## doubletree's estimators have no censoring argument, so the restriction
  ## happens here and the excluded patients are gone from the analytic data
  ## entirely -- which is exactly why the excluded COUNT must be reported.
  analytic_pre <- analysis_cohort |>                          # 02
    dplyr::select("ID", "INDEX_DT", "A", "elig_aap", "MSR_YR_selected", "month_gap") |>
    dplyr::left_join(                                         # 03: Y
      dplyr::transmute(outcome_cost, ID = .data$ID, Y = .data$cost),
      by = "ID"
    ) |>
    dplyr::left_join(                                         # 03: prior_cost ($)
      dplyr::transmute(prior_cost, ID = .data$ID, prior_cost = .data$cost),
      by = "ID"
    ) |>
    dplyr::left_join(covariates_yn, by = "ID") |>             # 04: the 12 _YN cols
    dplyr::left_join(complete_case, by = "ID")                # 05

  n_pre <- nrow(analytic_pre)

  ## Three DISTINCT exclusion reasons, reported separately -- see this
  ## file's header on why conflating them would be wrong. Order: complete
  ## case first (the largest, most substantive restriction), then the two
  ## coverage-truncation exclusions on what remains.
  ##
  ## `.data$complete_case %in% TRUE`, NOT `isTRUE(.data$complete_case)` --
  ## found while verifying this script 2026-09-24 (hand-built inputs,
  ## bypassing the local cascade's own Tier-0 boundary): isTRUE() is NOT
  ## vectorized. It collapses the whole COLUMN to a single logical (almost
  ## always a single FALSE, since a vector is never `identical(x, TRUE)`),
  ## and dplyr::filter() then applies that ONE value to every row --
  ## silently dropping ALL rows regardless of complete_case's real values.
  ## This bug predates this pass (inherited from the original pseudocode's
  ## same line) and would have made analytic_cc EMPTY on every real run,
  ## for a reason unrelated to complete_case actually being unresolved.
  ## `%in%` is vectorized and treats NA as non-matching, i.e. excluded --
  ## the correct behaviour for an unresolved complete-case status.
  analytic_cc <- dplyr::filter(analytic_pre, .data$complete_case %in% TRUE)
  n_not_complete_case <- n_pre - nrow(analytic_cc)

  analytic_known_y <- dplyr::filter(analytic_cc, !is.na(.data$Y))
  n_unknown_y <- nrow(analytic_cc) - nrow(analytic_known_y)

  analytic_known_x <- dplyr::filter(analytic_known_y, !is.na(.data$prior_cost))
  n_unknown_prior_cost <- nrow(analytic_known_y) - nrow(analytic_known_x)

  cli::cli_alert_info(c(
    "Attrition, in order: {n_pre} eligible",
    " -> {nrow(analytic_cc)} complete-case ({n_not_complete_case} excluded)",
    " -> {nrow(analytic_known_y)} with known Y ({n_unknown_y} excluded,
       file-coverage truncation)",
    " -> {nrow(analytic_known_x)} with known prior_cost ({n_unknown_prior_cost}
       excluded, file-coverage truncation)."
  ))
  if (nrow(analytic_known_x) == 0L) {
    cli::cli_abort(c(
      "Every patient was excluded; nothing to assemble.",
      "i" = "{n_not_complete_case} were not complete-case; of the rest,
             {n_unknown_y} had unknowable Y and {n_unknown_prior_cost} had
             unknowable prior_cost. If {n_not_complete_case} == {n_pre}, this
             is the expected outcome of complete_case not yet being
             computable for real (see 05's own Tier-0 limitation) -- not a
             data problem to work around here."
    ))
  }

  ## ---- Step 2: discretize prior_cost ON THE FINAL SAMPLE ------------------
  ##
  ## Order matters and this is why the call lives here, not only in 04: the
  ## confirmed design says quartiles are computed "within the analytic
  ## sample", and the analytic sample is the post-restriction,
  ## known-Y-and-prior_cost one -- not 04's pre-filter preview. Discretizing
  ## before the filter and then subsetting produces different boundaries.
  prior_cost_dummies <- discretize_prior_cost(analytic_known_x$prior_cost)
  cutpoints <- attr(prior_cost_dummies, "cutpoints")
  cli::cli_alert_info(
    "Realized quartile cutpoints on the analytic sample (REPORT THESE):
     {.val {round(cutpoints, 2)}}"
  )

  x_with_id <- assemble_design_matrix(analytic_known_x, prior_cost_dummies)
  ## The guard the package does not give you -- unconditional, regardless of
  ## OPEN_DECISIONS status. A _PROVISIONAL file must not be a way to ship a
  ## broken matrix.
  assert_binary_design_matrix(dplyr::select(x_with_id, -"ID"))

  analytic_data <- x_with_id |>
    dplyr::left_join(
      dplyr::select(analytic_known_x, "ID", "A", "Y", "prior_cost",
                    "MSR_YR_selected", "month_gap"),
      by = "ID"
    ) |>
    dplyr::select("ID", "A", "Y", dplyr::all_of(design_matrix_columns()),
                  "prior_cost", "MSR_YR_selected", "month_gap")
  ## prior_cost (dollars) is retained ALONGSIDE its dummies deliberately: it
  ## is not passed to the estimator, but the methods section needs its
  ## distribution and the cutpoints have to be reproducible from the saved
  ## data.
  attr(analytic_data, "prior_cost_cutpoints") <- cutpoints

  ## ---- Structural checks: UNCONDITIONAL, regardless of decision status ----
  if (nrow(analytic_data) == 0L) {
    cli::cli_abort("analytic_data has zero rows after assembly.")
  }
  if (anyDuplicated(analytic_data$ID) > 0L) {
    cli::cli_abort("analytic_data is not one row per patient.")
  }
  if (anyNA(analytic_data$A) || !all(analytic_data$A %in% c(0L, 1L))) {
    cli::cli_abort("{.field A} must be 0/1 with no NA.")
  }
  if (anyNA(analytic_data$Y)) {
    ## Should be structurally impossible given the known-Y filter above --
    ## asserted anyway, matching this codebase's "assert defensively even
    ## when guaranteed" convention (see helpers/design_matrix.R's header).
    cli::cli_abort(
      "{.field Y} has {.val NA} despite the known-Y filter -- the filter and
       the outcome disagree, and THAT is the finding, not a bug to patch
       around here."
    )
  }
  n_treated <- sum(analytic_data$A == 1L)
  n_control <- sum(analytic_data$A == 0L)
  if (n_treated == 0L || n_control == 0L) {
    cli::cli_abort(
      "One arm is empty (A=1: {n_treated}, A=0: {n_control}) -- neither
       estimator can fit on a single arm."
    )
  }
  ## Informational only, not an abort: estimate_att()'s OWN internal check
  ## (R/estimate_att.R, n_control < leaf_budget * m_n) uses the REAL m_n
  ## passed at call time in 07, which this script does not know yet.
  ## Reported here so a downstream estimator failure is traceable to the
  ## cohort rather than to the estimator.
  cli::cli_alert_info(
    "Arms: A=1 (treated) {n_treated}; A=0 (control) {n_control}. At
     {.code m_n = 1L} (optimaltrees' own default), estimate_att() needs
     n_control >= leaf_budget * m_n = {config_leaf_budget}; got {n_control}."
  )

  ## ---- the write: PROVISIONAL while anything blocks, real once nothing does
  suffix <- if (has_blocking_open()) "_PROVISIONAL" else ""
  if (nzchar(suffix)) {
    cli::cli_alert_warning(
      "Writing {.field {suffix}}: {.val {blocking_open_ids()}} still block
       the REPORTABLE artifact (07_estimate_att.R still refuses to run) --
       this file is an inspectable intermediate, not a result."
    )
  } else {
    cli::cli_alert_success("Nothing blocks -- writing the REPORTABLE artifact.")
  }

  fs::dir_create(output_dir)
  out_path <- file.path(output_dir, paste0(output_stem, suffix, ".parquet"))
  arrow::write_parquet(analytic_data, out_path)
  cli::cli_alert_success("Wrote {.file {out_path}} ({nrow(analytic_data)} rows).")

  ## The cutpoints are an ATTRIBUTE, and parquet does not carry R attributes.
  ## Written as a sidecar so the grid the sparsity assumption is stated over
  ## is recoverable from disk rather than only from this session.
  cutpoint_path <- file.path(output_dir, paste0(output_stem, suffix, "_prior_cost_cutpoints.csv"))
  utils::write.csv(
    data.frame(
      prob = names(attr(analytic_data, "prior_cost_cutpoints")),
      cutpoint = as.numeric(attr(analytic_data, "prior_cost_cutpoints"))
    ),
    cutpoint_path, row.names = FALSE
  )
  cli::cli_alert_success("Wrote {.file {cutpoint_path}}.")

  ## A manifest, so a _PROVISIONAL file that outlives the session still says
  ## WHICH decision made it provisional -- a filename suffix alone does not.
  manifest_path <- file.path(output_dir, paste0(output_stem, suffix, "_manifest.csv"))
  utils::write.csv(
    data.frame(
      field = c("timestamp", "provisional", "blocking_open_ids",
                "n_eligible", "n_complete_case", "n_known_y", "n_final",
                "n_treated", "n_control"),
      value = c(
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        nzchar(suffix),
        paste(blocking_open_ids(), collapse = ";"),
        n_pre, nrow(analytic_cc), nrow(analytic_known_y), nrow(analytic_data),
        n_treated, n_control
      )
    ),
    manifest_path, row.names = FALSE
  )
  cli::cli_alert_success("Wrote {.file {manifest_path}}.")
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info(
    "06_assemble_analytic_data.R complete
     ({if (config_has_smidata) 'ran the real assembly (or stopped at a real Tier-0 boundary -- see above)' else 'skipped, no smidata'})."
  )
}
