## ============================================================================
## 05_complete_case.R
##
## MIGRATED to smidata 2026-10-01. The month-by-month logic that used to live
## in helpers/complete_case.R::compute_complete_case() is now
## smidata::smi_enrollment_panel() + smi_enrollment_coverage(). That helper is
## DELETED, not kept as a shim: it had already been copied into dual-bounds
## (2026-09-28) and the two copies had drifted 82 lines apart, which is the
## promotion trigger smidata's own scope statement names -- "facts about the
## dataset promote; decisions about the analysis do not."
##
## THE TRANSLATION, stated because it is not a rename. The local function
## hardcoded two of OPEN_DECISIONS$enrollment_source's values; the promoted
## functions require them as explicit arguments with no defaults, so each is
## now passed from the decision rather than baked in:
##
##   type_filter = "none -- all TYPE values counted"  ->  type_policy
##   gap_months  = 0L                                 ->  gap_months + gap_mode
##
## type_policy = "require_agreement", NOT "any". This is the one place the
## translation needs care. The local function implemented "ignore TYPE" AND
## aborted if TYPE ever made MEDICAID_FLAG disagree for an (ID, month) pair --
## because the PI's "ignore" default is only a no-op if no such disagreement
## exists. smidata's "any" policy WARNS on disagreement instead of aborting;
## "require_agreement" aborts. The two policies are identical wherever there
## is no disagreement, so "require_agreement" reproduces the local function's
## behaviour EXACTLY, and "any" would silently weaken it. Keep it.
##
## gap_mode = "total" is passed explicitly and currently does not matter:
## smidata separates a gap COUNT ("total", matching the upstream CE rule's
## "consecutive or non-consecutive") from a gap RUN ("consecutive"), and at
## gap_months = 0 the two coincide -- zero uncovered months total is the same
## condition as a longest run of zero. IT WILL MATTER the moment the PI
## revises gap_months upward, which that decision's own note says is expected
## ("We can change this later"). At that point this paper must elect a mode;
## see NOTE-GAP-MODE below, which aborts rather than letting the literal
## stand.
##
## WHAT DID NOT CHANGE. Every correctness property the local function
## documented is in the promoted code and tested there: abort-on-ambiguous-
## YEAR_MONTH-format checking every value; the ID join-key guard; the
## all-absent tripwire; set membership against the exact required month set
## rather than a covered-month count (and the mid-month INDEX_DT shift that
## makes a count check wrong); the cohort's T0 (= s) as the window anchor (INDEX_DT only
## for the cross-check against the table's own copy). One guard was NARROWED on
## promotion, deliberately: the local ID check compared identical(class()) and
## so rejected integer-vs-double, a pairing dplyr joins correctly and the one
## that actually occurs here.
##
## WHAT STILL DOES NOT RUN THE SAME WAY ON A FIXTURE AS ON REAL DATA --
## unchanged by this migration, both limitations are in smidata's promoted
## code for the same reasons:
##
##   * The INDEX_DT cross-check fires SPURIOUSLY against Tier-0 fixtures:
##     smi_fixture() generates INDEX_DT independently per dataset with no
##     cross-table correlation for the same ID (verified 2026-09-24: 0 of 10
##     sampled IDs agreed). On real data that divergence is exactly the signal
##     the check exists to catch. The local patch below still applies, still
##     ONLY when is_local.
##   * YEAR_MONTH has no real-data levels in the census, so smi_fixture()
##     emits unparseable placeholder text and smi_enrollment_panel() correctly
##     aborts. Detected below so the script reports an expected Tier-0
##     limitation instead of crashing. The unblock is a real smi_census() run
##     on larger_smi_medicaid_monthly_flag -- smidata's own next-session
##     priority, not something to work around here.
##
## WHAT THIS SCRIPT MUST NOT DO, still: proxy enrollment from claims presence
## ("had a claim in month 23, therefore observed through month 24"). That
## proxy fails in two directional ways -- false negatives on zero-utilization
## complete patients, false positives on partially-enrolled patients with an
## early claim -- and MEDICAID_FLAG is not that proxy: it is an ENROLLMENT
## indicator recorded independent of claims activity, so neither failure mode
## applies. See smidata's inst/docs/dictionary/ENROLLMENT.md.
##
## SEPARATELY, AND STILL UNRESOLVED: death. There is no death source in the
## confirmed design, so a patient who dies at month 18 is indistinguishable
## from one who disenrolls at month 18, and both are excluded as incomplete.
## Excluding deaths is a substantive choice about the estimand -- the ATT is
## then defined among 24-month survivors -- and the manuscript must say so
## rather than treating it as a data-cleaning step.
##
## WINDOW (PI, 2026-10-15): this paper REQUIRES full enrollment through s + 24mo,
## i.e. every month of [s, s + 24mo) -- 24 months -- with type_policy
## "require_agreement", gap 0, gap_mode recorded. s (= T0 here) is the first day
## of the month on or after INDEX_DT. dual-bounds instead keeps everyone and
## models R.
##
## ACKNOWLEDGED LIMITATION. Complete-case status requires enrollment through
## s + 24mo, a window that CONTAINS the exposure period [s, s + 12mo), so it
## is determined partly AFTER treatment. The
## complete-case ATT equals the full-eligible-cohort ATT only if retention R is
## independent of (Y(0), Y(1)) given (A, X) AND the retained treated have the
## same X distribution as all treated (or effects do not vary with X); under
## the first condition alone, IPCW recovers the full-cohort ATT. Retention
## plausibly depends on adherence, i.e. on the treatment itself. No sensitivity
## analysis is run yet (IPCW or dual-bounds' modelled-R bounds are the
## candidates); the arm-specific complete-case rates are saved to
## output/attrition_by_arm.parquet (rate = n_complete / n, where n includes the
## n_unresolved patients) so the limitation can be reported with numbers.
##
## Run time: seconds against a Tier-0 fixture; hours against the real table
## (5.8GB / 24M rows -- see smidata's inst/server/07_enrollment_coverage.R
## header cost warning, which applies identically here).
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "disclosure.R"))
source(file.path(app_dir, "helpers", "decision_report.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("05 -- complete-case restriction")

## `confirmed_value()`, not `open_value()`: this decision is no longer a
## placeholder, and open_value()'s "Unconfirmed value" warning would be
## actively wrong here. Aborts (rather than proceeding on a stale assumption)
## if OPEN_DECISIONS$enrollment_source ever regresses to open_decision.
enrollment_decision <- confirmed_value("enrollment_source")

## The translation. gap_months is read FROM the decision rather than written
## as a literal, so a revision cannot leave this call stale.
cc_type_policy <- enrollment_decision$type_policy   # "require_agreement"; see header on why not "any"
stopifnot(identical(cc_type_policy, "require_agreement"))
cc_gap_months <- as.integer(enrollment_decision$gap_months)
## NOTE-GAP-MODE: "total" vs "consecutive" is a real election once
## cc_gap_months > 0. It is unambiguous only at 0. The mode is RECORDED in
## OPEN_DECISIONS$enrollment_source$value (PI project rule, 2026-10-15: gap 0,
## gap_mode recorded); the guard below aborts rather than letting it stand
## through an upward revision of the gap.
cc_gap_mode <- enrollment_decision$gap_mode   # recorded in OPEN_DECISIONS$enrollment_source (2026-10-15)
stopifnot(identical(cc_gap_mode, "total"))

cli::cli_alert_success(
  "OPEN_DECISIONS$enrollment_source is confirmed: TYPE filter =
   {.val {enrollment_decision$type_filter}}, gap tolerance =
   {.val {cc_gap_months}} month(s)."
)
cli::cli_alert_info(
  "Translated for smidata: type_policy = {.val {cc_type_policy}},
   gap_months = {cc_gap_months}, gap_mode = {.val {cc_gap_mode}}. See this
   file's header -- {.val require_agreement} is not a synonym for the PI's
   {.q ignore TYPE}, it is the faithful encoding of it."
)
if (cc_gap_months > 0L) {
  cli::cli_abort(c(
    "{.field gap_months} is {cc_gap_months}, so {.field gap_mode} is now a
     real election and the literal {.val {cc_gap_mode}} in this script is no
     longer defensible.",
    "i" = "{.val total} counts uncovered months anywhere in the window
           (matching the upstream CE rule); {.val consecutive} bounds the
           longest run. A patient with scattered single-month lapses passes
           one and fails the other.",
    "i" = "Record the mode in OPEN_DECISIONS$enrollment_source$value and read
           it here. See NOTE-GAP-MODE in this file."
  ))
}
cli::cli_alert_warning(
  "The TYPE default is provisional on application/90_checks_tier1.R CHECK 4
   (whole-table TYPE-disagreement check), NOT YET RUN as of this writing.
   {.fn smidata::smi_enrollment_panel} checks the identical question at
   cohort/window-restricted scope and ABORTS on disagreement under
   {.val require_agreement} -- but that is not a substitute for CHECK 4's
   whole-table run."
)

## ---- complete-case ----------------------------------------------------------

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; skipping the complete-case computation
     entirely (same convention as 01/02/03/04). NOTE: unlike before this
     file's 2026-10-01 migration there is no local fallback -- the logic now
     lives in smidata and this script cannot demonstrate it without the
     package. That is the intended consequence of promotion, not a gap."
  )
} else {
  require_stage("analysis_cohort", "02_population_and_eligibility.R")
  is_local <- identical(smidata::smi_env(), "local")

  smidata::smi_require(
    SMI_KEYS$medicaid_monthly_flag,
    columns = SMI_COLS$medicaid_monthly_flag
  )
  monthly_flag <- smi_read_pinned(
    SMI_KEYS$medicaid_monthly_flag,
    columns = SMI_COLS$medicaid_monthly_flag
  )

  ## Two cohort frames, deliberately. The panel is built on the true INDEX_DT
  ## so smi_enrollment_panel()'s cross-check against the table's own INDEX_DT
  ## keeps working; coverage is anchored on s (= T0, 02), the first day of the
  ## month on or after INDEX_DT, which is where the 24 required months start:
  ## window [s, s + 24mo). s is a 1st, so smidata's month floor leaves it alone.
  cohort_for_cc <- dplyr::select(analysis_cohort, "ID", "INDEX_DT")
  cohort_t0 <- dplyr::transmute(analysis_cohort, ID = .data$ID, INDEX_DT = .data$T0)

  if (is_local) {
    ## Fixture-only patch -- see this file's header. NEVER done against real
    ## data: there, INDEX_DT disagreement between the two tables is exactly
    ## the signal smi_enrollment_panel()'s cross-check exists to catch.
    cli::cli_alert_warning(
      "Local environment: overwriting the Tier-0 fixture's own
       {.field INDEX_DT} from {.field analysis_cohort} before building the
       panel -- smi_fixture() does not correlate {.field INDEX_DT} across
       tables for the same {.field ID}, so the (real, load-bearing)
       cross-check would abort on every patient. Tier-0-only."
    )
    monthly_flag <- monthly_flag |>
      dplyr::select(-"INDEX_DT") |>
      dplyr::left_join(cohort_for_cc, by = "ID")
  }

  ## smi_enrollment_panel() aborts on an unparseable YEAR_MONTH, which is the
  ## desired behaviour against real data with a genuinely wrong format and a
  ## guaranteed abort against the Tier-0 fixture's placeholder text. Detected
  ## here so the script reports a clear, expected limitation.
  ##
  ## Reads smidata's own "non_conforming" attribute (0.5.0+) rather than
  ## pattern-matching "FIXTURE_TEXT_". That pattern was duplicated here and in
  ## dual-bounds, which is why smi_fixture() now declares its own limitation;
  ## the placeholder's spelling is an implementation detail. The attribute
  ## empties itself once YEAR_MONTH gains real levels (census level-enumeration
  ## cap raised to 150, 2026-10-01 -- takes effect on the next capture), at
  ## which point this branch simply stops being taken, with no edit here.
  year_month_is_fixture_placeholder <-
    is_local && "YEAR_MONTH" %in% attr(monthly_flag, "non_conforming")

  if (year_month_is_fixture_placeholder) {
    cli::cli_alert_warning(
      "Local environment: {.field YEAR_MONTH} has no real-data levels in
       smidata's census, so the Tier-0 fixture emits unparseable placeholder
       text ({.val FIXTURE_TEXT_<n>}). {.fn smidata::smi_parse_year_month}
       correctly aborts on it, so the panel is NOT built here. This is a
       Tier-0 limitation, not a code defect -- the same computation is
       exercised for real, at the true {config_complete_case_months}-month
       window, on hand-built data in the demonstration below."
    )
  } else {
    enrollment_panel <- smidata::smi_enrollment_panel(
      monthly_flag = monthly_flag,
      cohort = cohort_for_cc,
      type_policy = cc_type_policy
    )
    ## A month absent from the table reads as NON-enrollment. The last required
    ## month is s's month + (window - 1); if the table ends before it, late-T0
    ## patients would be marked incomplete because of the table's range.
    obs_max <- attr(enrollment_panel, "observed_range")$max
    need_max <- max(smidata::smi_month_index(
      lubridate::year(cohort_t0$INDEX_DT), lubridate::month(cohort_t0$INDEX_DT)
    )) + config_complete_case_months - 1L
    if (need_max > obs_max) {
      cli::cli_abort(
        "The monthly-flag table ends {need_max - obs_max} month{?s} before the last
         required month; absent months would read as non-enrollment."
      )
    }
    panel_tally <- attr(enrollment_panel, "tally")

    coverage <- smidata::smi_enrollment_coverage(
      panel = enrollment_panel,
      cohort = cohort_t0,
      window_months = config_complete_case_months,
      gap_months = cc_gap_months,
      gap_mode = cc_gap_mode
    )
    cc_tally <- attr(coverage, "tally")

    ## 06_assemble_analytic_data.R joins on `ID` and filters
    ## `complete_case %in% TRUE`, so the column keeps its name here rather
    ## than renaming a downstream contract for a cosmetic reason.
    complete_case <- coverage |>
      dplyr::transmute(
        ID = .data$ID,
        complete_case = .data$covered,
        reason = .data$reason
      )

    cli::cli_alert_info(
      "Complete-case: {cc_tally$n_covered} of {cc_tally$n_cohort}
       ({round(100 * cc_tally$n_covered / cc_tally$n_cohort, 1)}%).
       Insufficient coverage: {cc_tally$n_insufficient}; absent from
       {.val {SMI_KEYS$medicaid_monthly_flag}} entirely: {cc_tally$n_absent}."
    )
    ## Two counts the local function did not report, both worth reading. A
    ## nonzero n_uncovered_at_index usually means a cohort/table mismatch
    ## rather than a real pattern, since INDEX_DT is by construction a date
    ## the patient had a qualifying claim.
    if (cc_tally$n_uncovered_at_index > 0L) {
      cli::cli_alert_warning(
        "{cc_tally$n_uncovered_at_index} patient{?s} {?is/are} uncovered in
         the s month (the first month of the exposure period, which is INDEX_DT's
         own month only when INDEX_DT is the 1st); compare it by lag in the
         table below."
      )
    }
    cli::cli_alert_info(
      "TYPE-disagreeing (ID, month) pairs within the cohort's windows:
       {panel_tally$n_type_disagreement_pairs} (zero by construction under
       {.val require_agreement} -- it aborts otherwise). Record this for
       CHECK 4."
    )
    if (is_local) {
      cli::cli_alert_warning(
        "The counts above are Tier-0 FIXTURE counts -- proof of execution,
         not a real complete-case rate."
      )
    }

    ## Differential attrition by arm is the substantive threat a 24-month
    ## continuous-enrollment filter poses (Oracle-flagged, 2026-09-24 design
    ## review): retention plausibly correlates with the exposure itself
    ## (medication adherence). Report it, don't just compute it.
    ## By lag (months from INDEX_DT's month to s: 0 or 1): the required enrollment
    ## months start at s, so a later s pushes the 24-month window later and
    ## nearer the end of the monthly-flag table. A whole stratum excluded would
    ## be a table-range artifact, not an enrollment pattern.
    cc_by_gap <- analysis_cohort |>
      dplyr::left_join(complete_case, by = "ID") |>
      dplyr::count(.data$t0_lag_months, .data$A, .data$reason)
    cli::cli_alert_info("Complete-case reason by T0 lag and arm:")
    cli::cli_verbatim(utils::capture.output(print(cc_by_gap, n = Inf)))
    fs::dir_create(file.path(app_dir, "output"))
    arrow::write_parquet(cc_by_gap, file.path(app_dir, "output", "complete_case_by_gap.parquet"))

    attrition_by_arm <- analysis_cohort |>
      dplyr::left_join(complete_case, by = "ID") |>
      dplyr::summarise(
        n = dplyr::n(),
        n_complete = sum(.data$complete_case %in% TRUE),
        n_unresolved = sum(is.na(.data$complete_case)),
        .by = "A"
      ) |>
      dplyr::mutate(
        rate = .data$n_complete / .data$n,
        env = smidata::smi_env(), ingest_id = effective_ingest_id(),
        gap_months = cc_gap_months, gap_mode = cc_gap_mode
      ) |>
      dplyr::arrange(.data$A)
    fs::dir_create(file.path(app_dir, "output"))
    arrow::write_parquet(
      attrition_by_arm, file.path(app_dir, "output", "attrition_by_arm.parquet")
    )
    cli::cli_alert_info(
      "Complete-case rate by arm (saved to output/attrition_by_arm.parquet):
       {paste0('A=', attrition_by_arm$A, ': ',
               round(100 * attrition_by_arm$rate, 1), '%', collapse = '; ')}."
    )

    ## Record for 08_decision_report.R: the per-patient result is reduced to
    ## n_covered and reason (no IDs) with smidata's tally attribute carried across.
    ## `reason` is needed: an absent patient has n_covered = window_months, so the
    ## months-covered table is split by reason (smi_report_enrollment()).
    coverage_for_report <- coverage[c("n_covered", "reason")]
    attr(coverage_for_report, "tally") <- cc_tally
    report_record("05_enrollment", list(
      coverage = coverage_for_report,
      window_months = config_complete_case_months,
      by_lag_arm = cc_by_gap,
      by_arm = attrition_by_arm
    ), app_dir)
  }
}

## ---- Tier-0 demonstration: the promoted functions, hand-built fixture ------
##
## Not a stand-in for real data (see application/README.md on why a
## hand-written fixture must never substitute for smi_fixture()'s
## census-backed one). Kept because these two patients exercise exact
## known-truth values -- fully covered vs. one required month short -- that a
## fixture's random MEDICAID_FLAG values do not guarantee to hit, and because
## it is the only proof of execution available while the YEAR_MONTH Tier-0
## limitation stands.

if (!config_has_smidata) {
  cli::cli_alert_info(
    "Skipping the demonstration too: it now calls smidata directly."
  )
} else {
  demo_cohort <- tibble::tibble(
    ID = c(1, 2),
    INDEX_DT = as.Date(c("2019-06-01", "2019-06-01"))
  )
  demo_months <- format(
    seq(as.Date("2019-06-01"), by = "month",
        length.out = config_complete_case_months),
    "%Y-%m"
  )
  demo_monthly_flag <- dplyr::bind_rows(
    tibble::tibble(
      ID = 1, INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
      YEAR_MONTH = demo_months, MEDICAID_FLAG = 1L
    ),
    ## demo 2: missing the LAST required month's row entirely, on purpose --
    ## demonstrates the computation distinguishes "complete" from "one
    ## required month short" rather than only exercising the covered path.
    tibble::tibble(
      ID = 2, INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
      YEAR_MONTH = demo_months[-length(demo_months)], MEDICAID_FLAG = 1L
    )
  )

  demo_cc <- smidata::smi_enrollment_coverage(
    panel = smidata::smi_enrollment_panel(
      demo_monthly_flag, demo_cohort, type_policy = cc_type_policy
    ),
    cohort = demo_cohort,
    window_months = config_complete_case_months,
    gap_months = cc_gap_months,
    gap_mode = cc_gap_mode
  )
  print(demo_cc)
  cli::cli_alert_info(
    "Demo (hand-built, {config_complete_case_months}-month window):
     ID 1 covered = {.val {demo_cc$covered[demo_cc$ID == 1]}} (fully covered);
     ID 2 covered = {.val {demo_cc$covered[demo_cc$ID == 2]}}
     (last required month absent)."
  )
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info(
    "05_complete_case.R complete
     ({if (!config_has_smidata) 'skipped, no smidata' else if (exists('year_month_is_fixture_placeholder') && isTRUE(year_month_is_fixture_placeholder)) 'Tier-0 YEAR_MONTH limitation, see above; demo ran instead' else 'ran the real complete-case computation'})."
  )
}
