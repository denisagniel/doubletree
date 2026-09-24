## ============================================================================
## 02_population_and_eligibility.R
##
## Base cohort -> AAP eligibility (elig_aap) -> exposure A (aap_achieved).
##
## EVERY DESIGN ELEMENT IN THIS SCRIPT IS INHERITED FROM dual-bounds, VERBATIM
## (smidata inst/analyses/doubletree__application.yml sections 2-4, all
## `status: inherited`). This paper adds no eligibility criterion of its own:
##   population.time_zero     = INDEX_DT
##   population.eligibility   = elig_aap == 1 iff MSR_DEN >= 1 in msr_aap
##   base_cohort.source_table = larger_smi_covariates (381,018 rows)
##   exposure.window          = [INDEX_DT, INDEX_DT + 12mo]
##   exposure.definition      = A = 1 iff MSR_NUM/MSR_DEN >= 0.5
##
## Inheriting the design does NOT mean calling dual-bounds' code. The shared
## implementation is smidata::smi_select_aap_row(), which both projects call
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
## and no warning at all. `AAP_MSR_CODE` (application/_config.R) is filtered
## on explicitly here, before calling smi_select_aap_row(), specifically to
## close that gap -- see OPEN_DECISIONS$msr_code's note for the full
## correction (design-reviewed by `oracle`, 2026-09-24). NOT YET propagated
## to smidata's shared registry or to dual-bounds' identical script -- both
## still carry the overstated claim and the identical exposure.
##
## Run time when unblocked: msr_aap is 19,012,819 rows. Filtering it to one MSR
## code before the join is what keeps smidata::smi_select_aap_row() tractable;
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
     {.val {SMI_KEYS$aap}} to the literal {.val {AAP_MSR_CODE}} below rather
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
  base_cohort <- smidata::smi_read(SMI_KEYS$covariates, columns = SMI_COLS$covariates)

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

  ## ---- msr_aap: contamination guard, THEN smi_select_aap_row() ----------
  aap_raw <- smidata::smi_read(SMI_KEYS$aap, columns = SMI_COLS$aap)

  ## filter_to_msr_code() (helpers/msr_filter.R) is pure and independently
  ## tested; `quiet` suppresses its own warning under a Tier-0 fixture, where
  ## an exclusion is expected and reported below in fixture-specific words,
  ## not a real finding.
  aap_raw <- filter_to_msr_code(aap_raw, AAP_MSR_CODE, quiet = is_local)
  n_excluded <- attr(aap_raw, "n_excluded")
  excluded_values <- attr(aap_raw, "excluded_msr_values")

  if (is_local && n_excluded > 0L) {
    cli::cli_alert_info(
      "{n_excluded} Tier-0 fixture row{?s} had {.field MSR} !=
       {.val {AAP_MSR_CODE}} (value{?s}: {.val {excluded_values}}) and
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
  if (nrow(aap_raw) == 0L) {
    if (is_local) {
      cli::cli_abort(
        "No {.val {AAP_MSR_CODE}} rows in the Tier-0 fixture for
         {.val {SMI_KEYS$aap}} this seed -- the fixture's synthetic MSR
         values did not happen to include {.val {AAP_MSR_CODE}}. Not a real
         finding; this pipeline cannot proceed further locally without one."
      )
    } else {
      cli::cli_abort(
        "No {.val {SMI_KEYS$aap}} rows with MSR == {.val {AAP_MSR_CODE}} --
         cannot proceed."
      )
    }
  }

  ## Grain change: one row per patient. tie_break is left at smidata's own
  ## default, "earlier", and is NOT an open decision in this project (see
  ## OPEN_DECISIONS's header comment on what is deliberately absent):
  ## smi_select_aap_row() aborts on any other value, so no choice is being
  ## silently made here to flag.
  aap_selected <- smidata::smi_select_aap_row(
    aap_raw,
    index_dt = dplyr::select(base_cohort, "ID", "INDEX_DT"),
    target_offset_months = 12L,
    achievement_threshold = config_aap_achievement_threshold,
    max_month_gap = Inf   # see 90_checks_tier1.R for the calibrating check
  )

  population <- base_cohort |>
    dplyr::left_join(aap_selected, by = "ID") |>
    dplyr::mutate(A = dplyr::if_else(.data$aap_achieved, 1L, 0L))

  ## The analysis cohort is the AAP-eligible subset. Report the attrition
  ## explicitly -- an eligibility filter that silently halves the cohort is
  ## the kind of thing that gets discovered in review.
  analysis_cohort <- dplyr::filter(population, .data$elig_aap)
  cli::cli_alert_info(
    "Eligible: {nrow(analysis_cohort)} of {nrow(population)};
     unmatched to msr_aap: {sum(!population$matched)}."
  )

  ## Both ATT entry points need interior treatment proportion and enough mass
  ## in BOTH arms; estimate_att() additionally needs at least
  ## leaf_budget * m_n CONTROL units, because its mu tree is fit on controls
  ## alone. Report the two arm counts here so a downstream estimator failure
  ## is traceable to the cohort rather than to the estimator.
  cli::cli_alert_info(
    "A = 1 (AAP achieved): {sum(analysis_cohort$A == 1L)};
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
