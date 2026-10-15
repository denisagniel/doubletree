## ============================================================================
## 02_population_and_eligibility.R
##
## Base cohort -> SAA eligibility (elig_saa) -> exposure A (saa_achieved).
## EXPOSURE SWITCHED from AAP (msr_aap) to SAA, antipsychotic adherence
## (msr_saa), 2026-10-09. The registry text below that mentions AAP is the
## ORIGINAL dual-bounds-inherited design and is kept as provenance.
##
## EVERY DESIGN ELEMENT IN THIS SCRIPT IS INHERITED FROM dual-bounds, VERBATIM
## (smidata inst/analyses/doubletree__application.yml sections 2-4, all
## `status: inherited`), plus one criterion ported from dual-bounds (PI,
## 2026-09-30) that concerns the extract, not the patient:
##   population.index         = INDEX_DT (determines s; anchors prior_cost and the _YN covariates)
##   population.eligibility   = elig_saa == 1 iff an msr_saa period STARTS EXACTLY
##                              at s (first of the month on or after INDEX_DT)
##                              with MSR_DEN >= 1; all others are DROPPED
##                              (PI, 2026-10-15; no nearest period, no gap
##                              tolerance, no ties)
##   population.coverage      = elig_coverage: every config_coverage_required_windows
##                              window lies inside the cost files' date coverage
##   s (called T0 here)       = first day of the month on or after INDEX_DT
##                              = the selected period's start; outcome and
##                              complete case anchor on s, prior_cost on
##                              INDEX_DT itself; see _config.R
##   base_cohort.source_table = larger_smi_covariates (381,018 rows)
##   exposure.window          = [s, s + 12mo): the selected 12-month period
##   exposure.definition      = A = 1 iff MSR_NUM/MSR_DEN >= 0.5 (>=, PI 2026-10-15)
##
## Inheriting the design does NOT mean calling dual-bounds' code. The shared
## implementation is smidata::smi_select_measure_row(), which both projects call
## directly; nothing here reads anything from ~/RAND/rprojects/smi/dual-bounds.
##
## REAL AS OF 2026-09-24. OPEN_DECISIONS$msr_code was confirmed 2026-09-17 but
## this script's own gate was never updated to read that value -- found while
## implementing 05_complete_case.R the same day, fixed here. `msr_code <-
## open_value("msr_code")` always returned the confirmed decision's VALUE
## ("AAP (no filter needed)"), never NULL, so the old `if (is.null(msr_code))`
## abort had been silently dead code since 2026-09-17, not a live gate -- the
## script kept reporting "unknown MSR code" about a code that had, in fact,
## been known for a week.
##
## THE MSR CONTAMINATION GUARD BELOW IS NOT OPTIONAL POLISH. The confirmed
## resolution's own registry note originally claimed
## "smidata::smi_select_aap_row()'s abort-on-tie behavior remains the safety
## net regardless" if `msr_aap` ever contains a non-AAP row. Reading that
## function's actual source (~/RAND/tools/smidata/R/msr_aap.R) shows this is
## not true: its winner-selection groups and joins on `ID` only, and never
## reads `MSR` except inside the tie-diagnostic message -- the abort fires
## ONLY when multiple rows tie at the exact same winning month_gap for one
## patient. A non-AAP row that is merely CLOSER in time to the target month
## than any true AAP row, with no tie, would be silently selected as the
## winner, corrupting that patient's `elig_aap`/`aap_achieved` with no abort
## and no warning at all. `SAA_MSR_CODE` (application/_config.R) is filtered
## on explicitly here, before calling smi_select_measure_row(), specifically to
## close that gap -- see OPEN_DECISIONS$msr_code's note for the full
## correction (design-reviewed by `oracle`, 2026-09-24). NOT YET propagated
## to smidata's shared registry or to dual-bounds' identical script -- both
## still carry the overstated claim and the identical exposure.
##
## Run time when unblocked: msr_saa is ~3.5M rows (msr_aap was 19,012,819). Filtering it to one MSR
## code before the join is what keeps smidata::smi_select_measure_row() tractable;
## expect minutes, not seconds, on the server. Locally, smidata::smi_read()
## transparently returns a Tier-0 fixture instead (structurally valid,
## deliberately unrealistic values) -- see this file's Tier-0 section below.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "msr_filter.R"))
source(file.path(app_dir, "helpers", "period_start.R"))
source(file.path(app_dir, "helpers", "decision_report.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("02 -- population and eligibility")

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed on this machine; skipping population and
     eligibility assembly entirely (same convention as
     01_declare_requirements.R)."
  )
} else {

  ## `confirmed_value()`, not `open_value()`: this decision is no longer a
  ## placeholder. Aborts (rather than proceeding on a stale assumption) if
  ## OPEN_DECISIONS$msr_code ever regresses to open_decision.
  msr_code <- confirmed_value("msr_code")
  cli::cli_alert_success(
    "OPEN_DECISIONS$msr_code is confirmed: {.val {msr_code}}. Filtering
     {.val {SMI_KEYS$saa}} to the literal {.val {SAA_MSR_CODE}} below rather
     than trusting 'no filter needed' unverified -- see this file's header."
  )

  is_local <- identical(smidata::smi_env(), "local")
  if (is_local) {
    cli::cli_alert_warning(
      "Local environment: {.fn smidata::smi_read} returns a Tier-0 FIXTURE
       (structurally valid, deliberately unrealistic values -- see
       {.fn smi_fixture}'s own warning on read). Never report a number
       computed below as a real finding; this run proves the pipeline
       EXECUTES, nothing about the real data."
    )
  }

  ## ---- base cohort -----------------------------------------------------
  base_cohort <- smi_read_pinned(SMI_KEYS$covariates, columns = SMI_COLS$covariates)

  if (nrow(base_cohort) == 0L) {
    cli::cli_abort("{.val {SMI_KEYS$covariates}} returned zero rows.")
  }
  if (anyDuplicated(base_cohort$ID) > 0L) {
    cli::cli_abort("{.val {SMI_KEYS$covariates}} is not one row per patient.")
  }
  if (any(is.na(base_cohort$INDEX_DT))) {
    cli::cli_abort(
      "{.field INDEX_DT} determines s (and so T0); {.val NA} is not
       recoverable and must not be dropped silently."
    )
  }
  ## Confirmed row count, smidata snapshot 2026-09-15_2824080e -- real data
  ## only; a Tier-0 fixture's row count (default n=200) is expected to differ
  ## and is not itself a finding.
  if (!is_local && nrow(base_cohort) != 381018L) {
    cli::cli_warn(
      "Expected 381,018 rows; got {nrow(base_cohort)}. A mismatch means the
       source was replaced in place (it has no upstream archive) -- a
       finding, not a nuisance."
    )
  }

  ## ---- msr_saa: contamination guard, THEN smi_select_measure_row() ------
  saa_raw <- smi_read_pinned(SMI_KEYS$saa, columns = SMI_COLS$saa)

  ## filter_to_msr_code() (helpers/msr_filter.R) is pure and independently
  ## tested; `quiet` suppresses its own warning under a Tier-0 fixture, where
  ## an exclusion is expected and reported below in fixture-specific words,
  ## not a real finding.
  saa_raw <- filter_to_msr_code(saa_raw, SAA_MSR_CODE, quiet = is_local)
  n_excluded <- attr(saa_raw, "n_excluded")
  excluded_values <- attr(saa_raw, "excluded_msr_values")

  if (is_local && n_excluded > 0L) {
    cli::cli_alert_info(
      "{n_excluded} Tier-0 fixture row{?s} had {.field MSR} !=
       {.val {SAA_MSR_CODE}} (value{?s}: {.val {excluded_values}}) and
       {?was/were} excluded. Expected -- fixture MSR values are synthetic,
       not a real finding."
    )
  } else if (!is_local && n_excluded > 0L) {
    cli::cli_alert_warning(
      "Report this back to smidata's dual-bounds__application.yml registry
       (shared with dual-bounds) -- see OPEN_DECISIONS$msr_code's
       registry_ref. filter_to_msr_code()'s own warning above is the
       detail; this is the pointer to where it needs to land."
    )
  }
  if (nrow(saa_raw) == 0L) {
    if (is_local) {
      cli::cli_abort(
        "No {.val {SAA_MSR_CODE}} rows in the Tier-0 fixture for
         {.val {SMI_KEYS$saa}} this seed -- the fixture's synthetic MSR
         values did not happen to include {.val {SAA_MSR_CODE}}. Not a real
         finding; this pipeline cannot proceed further locally without one."
      )
    } else {
      cli::cli_abort(
        "No {.val {SMI_KEYS$saa}} rows with MSR == {.val {SAA_MSR_CODE}} --
         cannot proceed."
      )
    }
  }

  ## Grain change: one row per patient. PROJECT-WIDE RULE (PI, 2026-10-15):
  ## s = first day of the month on or after INDEX_DT, and the msr_saa period
  ## must START EXACTLY at s -- anchor = "ceiling_month": no nearest period, no
  ## gap tolerance (max_month_gap is left at Inf, which that anchor requires),
  ## no ties across periods. saa_raw was filtered to the single SAA code above,
  ## BEFORE this call. tie_break stays at smidata's "earlier" and is not an
  ## open decision here. MSR_YR is the START of the 12-month period (smidata
  ## 06dbbe2, 2026-10-09); that start convention is not yet confirmed against
  ## msr_saa itself.
  base_index <- dplyr::transmute(base_cohort, ID = .data$ID, INDEX_DT = as.Date(.data$INDEX_DT))
  saa_selected <- smidata::smi_select_measure_row(
    saa_raw,
    index_dt = base_index,
    achievement_threshold = config_saa_achievement_threshold,
    anchor = config_period_anchor
  )
  n_saa_rows_out_of_cohort <- sum(!saa_raw$ID %in% base_cohort$ID)
  cli::cli_alert_info(
    "{n_saa_rows_out_of_cohort} msr_saa row{?s} {?is/are} for IDs not in
     {.val {SMI_KEYS$covariates}} (ignored)."
  )

  ## ---- s (= T0): the exposure start -----------------------------------------
  ##
  ## For every matched patient the selected period's start IS s by construction
  ## of the anchor; recomputing s from INDEX_DT independently and comparing is
  ## the check that this script and smidata agree on the rule. EVERY window
  ## downstream is anchored on s or on INDEX_DT as DOUBLETREE_WINDOWS (_config.R)
  ## says: outcome [s+12mo, s+24mo) and complete case [s, s+24mo) on s; prior_cost
  ## [INDEX_DT-12mo, INDEX_DT) on INDEX_DT itself.
  matched_rows <- saa_selected[saa_selected$matched %in% TRUE, , drop = FALSE] |>
    dplyr::inner_join(base_index, by = "ID")
  t0_tbl <- tibble::tibble(
    ID = matched_rows$ID,
    T0 = msr_yr_period_start(matched_rows$MSR_YR_selected),
    INDEX_DT = matched_rows$INDEX_DT
  ) |>
    dplyr::mutate(
      t0_lag_months = smidata::smi_month_index(
        lubridate::year(.data$T0), lubridate::month(.data$T0)
      ) - smidata::smi_month_index(
        lubridate::year(.data$INDEX_DT), lubridate::month(.data$INDEX_DT)
      )
    )
  if (nrow(t0_tbl) == 0L) {
    cli::cli_abort(
      "Zero patients have an {.val {SAA_MSR_CODE}} period starting exactly at
       s (the first of the month on or after INDEX_DT): there is no T0 to
       anchor the windows on."
    )
  }
  if (any(t0_tbl$T0 != exposure_start(t0_tbl$INDEX_DT)) ||
      any(matched_rows$month_gap != 0L)) {
    cli::cli_abort(
      "A selected period does not start exactly at s = first of the month on
       or after INDEX_DT; smidata's ceiling_month anchor should have
       guaranteed it."
    )
  }
  ## s is INDEX_DT's own month (INDEX_DT on the 1st: lag 0) or the next (lag 1).
  if (!all(t0_tbl$t0_lag_months %in% c(0L, 1L))) {
    cli::cli_abort("t0_lag_months must be 0 (INDEX_DT on the 1st) or 1.")
  }
  lag_tab <- table(t0_tbl$t0_lag_months)
  cli::cli_alert_info(
    "s (T0) minus INDEX_DT's month, in months: {paste0(names(lag_tab), ': ',
     lag_tab, collapse = '; ')}."
  )

  ## ---- observability: file-date coverage (ported from dual-bounds) --------
  ##
  ## Exclude anyone whose study period extends beyond the observable cost
  ## files (PI, 2026-09-30, via dual-bounds). An eligibility criterion about
  ## the EXTRACT's vintage, not the patient. A pure function of s and INDEX_DT,
  ## computed here ahead of 03's expensive pass. Uses COST_SOURCES (_config.R),
  ## the same object 03 derives its own `coverage` from, so 03's
  ## contract-violation abort stays meaningful. Evaluated only for patients
  ## with an exposure period (a T0); the rest are already dropped as ineligible.
  coverage_wide <- smidata::smi_window_coverage(
    doubletree_patient_windows(t0_tbl),
    COST_SOURCES
  ) |>
    dplyr::select("id", "window", "coverage") |>
    tidyr::pivot_wider(names_from = "window", values_from = "coverage",
                       names_prefix = "coverage_")
  want_cov_cols <- paste0("coverage_", DOUBLETREE_WINDOWS$window)
  if (nrow(coverage_wide) != nrow(t0_tbl) || anyDuplicated(coverage_wide$id) > 0L ||
      !all(want_cov_cols %in% names(coverage_wide)) ||
      anyNA(coverage_wide[want_cov_cols])) {
    cli::cli_abort("smi_window_coverage() did not return one complete row per patient.")
  }
  t0_tbl <- dplyr::select(t0_tbl, -"INDEX_DT")

  ## Report EVERY window's rate, gating or not.
  for (w in DOUBLETREE_WINDOWS$window) {
    n_bad <- sum(coverage_wide[[paste0("coverage_", w)]] != "full")
    cli::cli_alert_info(
      "{.field {w}} window: {n_bad} of {nrow(coverage_wide)} SAA-matched
       patients not fully covered by the files{if (w %in%
       config_coverage_required_windows) ' -- EXCLUDED below' else ' -- NOT
       gating eligibility'}."
    )
  }

  coverage_wide$elig_coverage <- Reduce(`&`, lapply(
    config_coverage_required_windows,
    function(w) coverage_wide[[paste0("coverage_", w)]] == "full"
  ))

  population <- base_cohort |>
    dplyr::left_join(
      dplyr::rename(saa_selected, elig_saa = "elig", saa_achieved = "achieved"),
      by = "ID"
    ) |>
    dplyr::left_join(t0_tbl, by = "ID") |>
    dplyr::left_join(
      dplyr::select(coverage_wide, ID = "id",
                    dplyr::starts_with("coverage_"), "elig_coverage"),
      by = "ID"
    ) |>
    dplyr::mutate(A = dplyr::if_else(.data$saa_achieved, 1L, 0L))

  ## Every SAA-eligible patient has an exposure period, hence a T0 and a
  ## coverage verdict; an NA there means a join lost rows. Patients without an
  ## exposure period have no T0 and no verdict: they are ineligible via
  ## elig_saa already, so elig_coverage is set FALSE (not evaluated), not NA.
  if (anyNA(population$T0[population$elig_saa]) ||
      anyNA(population$elig_coverage[population$elig_saa])) {
    cli::cli_abort("An SAA-eligible patient has no T0 or coverage verdict: a join lost rows.")
  }
  population$elig_coverage <- population$elig_coverage %in% TRUE

  ## TWO independent criteria, reported marginally AND jointly: elig_saa is a
  ## property of the patient, elig_coverage a property of the extract. The
  ## analysis cohort satisfies BOTH; every other patient is DROPPED (PI,
  ## 2026-10-15: ineligible patients are dropped, never coded 0 or kept as NA).
  ## Report the attrition explicitly -- an eligibility filter that silently
  ## halves the cohort is the kind of thing that gets discovered in review.
  analysis_cohort <- dplyr::filter(population, .data$elig_saa, .data$elig_coverage)
  ## "Not matched" has two causes: no msr_saa row at all, or rows exist but none
  ## STARTS exactly at s. They must add up to the unmatched count.
  n_no_saa_rows <- sum(!base_cohort$ID %in% saa_raw$ID)
  n_no_row_at_s <- sum(!population$matched) - n_no_saa_rows
  if (n_no_row_at_s < 0L) {
    cli::cli_abort(
      "Unmatched-patient partition does not add up ({n_no_saa_rows} no rows,
       {n_no_row_at_s} none starting at s): the msr_saa and covariates tables
       disagree on patient IDs."
    )
  }
  n_with_rows <- sum(base_cohort$ID %in% saa_raw$ID)
  if (n_with_rows == 0L) {
    cli::cli_abort("No cohort ID has an msr_saa row: the two tables share no IDs.")
  }
  cli::cli_alert_info(
    "SAA-matched (a period starting exactly at s): {sum(population$matched)} of
     {n_with_rows} patients with msr_saa rows
     ({round(100 * sum(population$matched) / n_with_rows, 1)}%). A low rate
     means the table's period structure does not put a period start on the
     first of the month on or after INDEX_DT -- check 90 CHECK 6 before
     trusting the cohort."
  )
  cli::cli_alert_info(
    "Attrition: {nrow(population)} base cohort -> {sum(population$elig_saa)}
     elig_saa -> {sum(population$elig_saa & population$elig_coverage)} of those
     also elig_coverage (the analysis cohort: {nrow(analysis_cohort)}).
     DROPPED as not matched: {n_no_saa_rows} with no {SAA_MSR_CODE} rows,
     {n_no_row_at_s} with rows but none starting exactly at s. DROPPED as
     matched but MSR_DEN < 1: {sum(population$matched & !population$elig_saa)}.
     The per-window coverage counts above are over all SAA-matched patients,
     so they include those."
  )
  ## achieved is NA when MSR_NUM is NA although elig is TRUE. Stop here, not
  ## after 03's multi-hour pass.
  if (anyNA(analysis_cohort$A)) {
    cli::cli_abort(
      "{sum(is.na(analysis_cohort$A))} eligible patient{?s} ha{?s/ve} NA SAA
       achievement (NA MSR_NUM with MSR_DEN >= 1)."
    )
  }
  if (nrow(analysis_cohort) == 0L) {
    cli::cli_abort("Zero patients satisfy elig_saa AND elig_coverage.")
  }

  ## ---- record the attrition for 08_decision_report.R ----------------------
  ## The counts the report cites are RECORDED here, as the pipeline runs, not
  ## recomputed later. Raw (unsuppressed) and local: output/ is gitignored and
  ## the report suppresses everything it renders or exports.
  report_record("02_population", list(
    flow = data.frame(
      step = c("base cohort (larger_smi_covariates)",
               sprintf("has >= 1 msr_saa row with MSR == %s", SAA_MSR_CODE),
               "a period starts exactly at s (first of the month on or after INDEX_DT)",
               "MSR_DEN >= 1 (eligible for SAA)",
               "cost-file coverage of the outcome and prior_cost windows (analysis cohort)"),
      n = c(nrow(base_cohort), n_with_rows, sum(population$matched),
            sum(population$elig_saa), nrow(analysis_cohort)),
      reason = c("",
                 "no msr_saa rows for the SAA measure",
                 "rows exist but none starts exactly at s (no nearest period, no tolerance)",
                 "MSR_DEN < 1",
                 "a window in config_coverage_required_windows extends past the cost files"),
      stringsAsFactors = FALSE
    ),
    selection = saa_selected[c("ID", "matched", "elig", "achieved",
                               "MSR_YR_selected", "MSR_DEN")],
    index_dt = base_index,
    coverage_counts = tibble::tibble(
      window = DOUBLETREE_WINDOWS$window,
      gating = DOUBLETREE_WINDOWS$window %in% config_coverage_required_windows,
      n_matched = nrow(coverage_wide),
      n_not_full = vapply(DOUBLETREE_WINDOWS$window, function(w) {
        sum(coverage_wide[[paste0("coverage_", w)]] != "full")
      }, integer(1))
    ),
    n_msr_rows_excluded = n_excluded,
    n_saa_rows_out_of_cohort = n_saa_rows_out_of_cohort
  ), app_dir)

  ## Both ATT entry points need interior treatment proportion and enough mass
  ## in BOTH arms; estimate_att() additionally needs at least
  ## leaf_budget * m_n CONTROL units, because its mu tree is fit on controls
  ## alone. Report the two arm counts here so a downstream estimator failure
  ## is traceable to the cohort rather than to the estimator.
  cli::cli_alert_info(
    "A = 1 (SAA achieved): {sum(analysis_cohort$A == 1L)};
     A = 0: {sum(analysis_cohort$A == 0L)}."
  )
  if (is_local) {
    cli::cli_alert_warning(
      "The counts above are Tier-0 FIXTURE counts (n={nrow(base_cohort)}
       synthetic rows) -- proof this pipeline executes end to end, not a
       real eligibility or exposure distribution."
    )
  }
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info(
    "02_population_and_eligibility.R complete
     ({if (config_has_smidata) 'ran for real (fixture or server)' else 'skipped, no smidata'})."
  )
}
