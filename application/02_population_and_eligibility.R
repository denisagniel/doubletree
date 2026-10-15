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
## `status: inherited`). This paper adds no eligibility criterion of its own:
##   population.time_zero     = INDEX_DT
##   population.eligibility   = elig_saa == 1 iff MSR_DEN >= 1 in msr_saa
##   base_cohort.source_table = larger_smi_covariates (381,018 rows)
##   exposure.window          = [INDEX_DT, INDEX_DT + 12mo]
##   exposure.definition      = A = 1 iff MSR_NUM/MSR_DEN >= 0.5
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
      "{.field INDEX_DT} is the time zero for every window; {.val NA} is not
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

  ## Grain change: one row per patient. tie_break is left at smidata's own
  ## default, "earlier", and is NOT an open decision in this project (see
  ## OPEN_DECISIONS's header comment on what is deliberately absent):
  ## smi_select_measure_row() aborts on any other value, so no choice is being
  ## silently made here to flag.
  saa_selected <- smidata::smi_select_measure_row(
    saa_raw,
    index_dt = dplyr::select(base_cohort, "ID", "INDEX_DT"),
    ## MSR_YR is the START of the 12-month measurement period (smidata
    ## 06dbbe2, 2026-10-09), so the period starting closest to INDEX_DT
    ## (offset 0) IS the exposure window [INDEX_DT, INDEX_DT + 12mo]. The
    ## old 12L read MSR_YR as the period END. Passed explicitly rather than
    ## relying on smidata's default, which flipped 12L -> 0L in that commit.
    ## Start convention not yet confirmed against msr_saa itself.
    target_offset_months = 0L,
    achievement_threshold = config_saa_achievement_threshold,
    max_month_gap = Inf   # see 90_checks_tier1.R for the calibrating check
  )

  ## ---- observability: file-date coverage (ported from dual-bounds) --------
  ##
  ## Exclude anyone whose study period extends beyond the observable cost
  ## files (PI, 2026-09-30, via dual-bounds). An eligibility criterion about
  ## the EXTRACT's vintage, not the patient. Pure in INDEX_DT, so it sits here
  ## ahead of 03's expensive pass. Uses COST_SOURCES (_config.R), the same
  ## object 03 derives its own `coverage` from, so 03's contract-violation
  ## abort stays meaningful.
  coverage_wide <- smidata::smi_window_coverage(
    smidata::smi_patient_windows(
      dplyr::transmute(base_cohort, id = .data$ID,
                       index_dt = as.Date(.data$INDEX_DT)),
      windows = DOUBLETREE_WINDOWS
    ),
    COST_SOURCES
  ) |>
    dplyr::select("id", "window", "coverage") |>
    tidyr::pivot_wider(names_from = "window", values_from = "coverage",
                       names_prefix = "coverage_")

  ## Report EVERY window's rate, gating or not.
  for (w in DOUBLETREE_WINDOWS$window) {
    n_bad <- sum(coverage_wide[[paste0("coverage_", w)]] != "full")
    cli::cli_alert_info(
      "{.field {w}} window: {n_bad} of {nrow(coverage_wide)} patients not
       fully covered by the files{if (w %in% config_coverage_required_windows)
       ' -- EXCLUDED below' else ' -- NOT gating eligibility'}."
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
    dplyr::left_join(
      dplyr::select(coverage_wide, ID = "id",
                    dplyr::starts_with("coverage_"), "elig_coverage"),
      by = "ID"
    ) |>
    dplyr::mutate(A = dplyr::if_else(.data$saa_achieved, 1L, 0L))

  ## elig_coverage is a pure function of INDEX_DT; an NA here means the join
  ## lost rows, not that coverage is unknown.
  stopifnot(!anyNA(population$elig_coverage))

  ## TWO independent criteria, reported marginally AND jointly: elig_saa is a
  ## property of the patient, elig_coverage a property of the extract. The
  ## analysis cohort satisfies BOTH. Report the attrition explicitly -- an
  ## eligibility filter that silently halves the cohort is the kind of thing
  ## that gets discovered in review.
  analysis_cohort <- dplyr::filter(population, .data$elig_saa, .data$elig_coverage)
  cli::cli_alert_info(
    "Attrition: {nrow(population)} base cohort -> {sum(population$elig_saa)}
     elig_saa -> {sum(population$elig_coverage)} elig_coverage ->
     {nrow(analysis_cohort)} both (the analysis cohort).
     Unmatched to msr_saa: {sum(!population$matched)}."
  )
  if (nrow(analysis_cohort) == 0L) {
    cli::cli_abort("Zero patients satisfy elig_saa AND elig_coverage.")
  }

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
