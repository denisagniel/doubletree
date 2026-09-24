## ============================================================================
## 05_complete_case.R
##
## REAL IMPLEMENTATION as of 2026-09-24. Complete-case status IS now
## computable: OPEN_DECISIONS$enrollment_source resolved (PI, 2026-09-24,
## application/_config.R) -- MEDICAID_FLAG==1 in larger_smi_medicaid_
## monthly_flag for every one of the config_complete_case_months (24)
## consecutive months starting at INDEX_DT's own calendar month, no TYPE
## filtering, zero gap tolerance. The actual month-by-month logic lives in
## helpers/complete_case.R::compute_complete_case() -- a pure function, tested
## with hand-built fixtures in tests/test-complete-case.R, with no server and
## no smidata dependency (same split as helpers/discretize_prior_cost.R /
## 06_assemble_analytic_data.R, for the same reason: the WINDOWING/MEMBERSHIP
## logic must not be trusted to "it ran on the server and looked plausible").
##
## WHAT STILL DOES NOT RUN THE SAME WAY ON A FIXTURE AS ON REAL DATA. As of
## 2026-09-24, `02` produces a real `analysis_cohort` (fixture-backed
## locally; real on the server -- see that file's header), so this script's
## orchestration below IS live code now, cascade-guarded via
## `require_stage()`. But `compute_complete_case()`'s `INDEX_DT` cross-check
## (a real correctness feature -- see its own Details) fires SPURIOUSLY
## against Tier-0 fixtures: `smi_fixture()` generates `INDEX_DT`
## independently per dataset, with no cross-table correlation for the same
## `ID` (verified 2026-09-24: 0 of 10 sampled IDs agreed between the
## `larger_smi_covariates` and `larger_smi_medicaid_monthly_flag` fixtures).
## On real data this divergence is exactly the signal the check exists to
## catch (an ID collision or a stale copy); on a fixture it is a fixture
## artifact. The orchestration below overwrites the fixture's own
## `monthly_flag$INDEX_DT` from `analysis_cohort` ONLY when `is_local` --
## never on real data, where disagreement must still abort.
##
## WHAT THIS SCRIPT MUST NOT DO, still: proxy enrollment from claims presence
## ("had a claim in month 23, therefore observed through month 24"). That
## proxy fails in the two directional ways the previous version of this file
## documented at length -- false negatives on zero-utilization complete
## patients, false positives on partially-enrolled patients with an early
## claim -- and MEDICAID_FLAG (an ENROLLMENT indicator, recorded independent
## of claims activity) is not that proxy: see
## helpers/complete_case.R::compute_complete_case()'s Details for why neither
## failure mode applies to this input.
##
## SEPARATELY, AND STILL UNRESOLVED: death. There is no death source in the
## confirmed design, so a patient who dies at month 18 is indistinguishable
## from one who disenrolls at month 18, and both are excluded as incomplete.
## Excluding deaths is a substantive choice about the estimand -- the ATT is
## then defined among 24-month survivors -- and the manuscript must say so
## rather than treating it as a data-cleaning step.
##
## Run time: seconds against a Tier-0 fixture; hours against the real table
## (5.8GB / 24M rows -- see ~/RAND/tools/smidata/inst/server/
## 07_enrollment_coverage.R's header cost warning, which applies identically
## here).
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "complete_case.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("05 -- complete-case restriction")

## `confirmed_value()`, not `open_value()`: this decision is no longer a
## placeholder, and open_value()'s "Unconfirmed value" warning would be
## actively wrong here. Aborts (rather than proceeding on a stale assumption)
## if OPEN_DECISIONS$enrollment_source ever regresses to open_decision.
enrollment_decision <- confirmed_value("enrollment_source")
cli::cli_alert_success(
  "OPEN_DECISIONS$enrollment_source is confirmed: TYPE filter =
   {.val {enrollment_decision$type_filter}}, gap tolerance =
   {.val {enrollment_decision$gap_months}} month(s)."
)
cli::cli_alert_warning(
  "This default is provisional on application/90_checks_tier1.R CHECK 4
   (whole-table TYPE-disagreement check), NOT YET RUN as of this writing.
   compute_complete_case() self-checks the identical question at
   cohort/window-restricted scope and ABORTS on disagreement -- see its
   Details -- but that is not a substitute for CHECK 4's whole-table run."
)

## ---- complete-case: REAL as of 2026-09-24 -----------------------------------

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; skipping the real complete-case
     computation entirely (same convention as 01/02/03/04)."
  )
} else {
  require_stage("analysis_cohort", "02_population_and_eligibility.R")
  is_local <- identical(smidata::smi_env(), "local")

  smidata::smi_require(
    SMI_KEYS$medicaid_monthly_flag,
    columns = SMI_COLS$medicaid_monthly_flag
  )
  monthly_flag <- smidata::smi_read(
    SMI_KEYS$medicaid_monthly_flag,
    columns = SMI_COLS$medicaid_monthly_flag
  )

  cohort_for_cc <- dplyr::select(analysis_cohort, "ID", "INDEX_DT")

  if (is_local) {
    ## Fixture-only patch -- see this file's header. NEVER done against real
    ## data: there, INDEX_DT disagreement between the two tables is exactly
    ## the signal compute_complete_case()'s cross-check exists to catch.
    cli::cli_alert_warning(
      "Local environment: overwriting the Tier-0 fixture's own
       {.field INDEX_DT} from {.field analysis_cohort} before calling
       compute_complete_case() -- smi_fixture() does not correlate
       {.field INDEX_DT} across tables for the same {.field ID}, so the
       (real, load-bearing) cross-check would abort on every patient. This
       patch is Tier-0-only; see this file's header."
    )
    monthly_flag <- monthly_flag |>
      dplyr::select(-"INDEX_DT") |>
      dplyr::left_join(cohort_for_cc, by = "ID")
  }

  ## YEAR_MONTH is a text column with no real-data levels in the census, so
  ## smi_fixture() falls to its terminal placeholder branch
  ## ("FIXTURE_TEXT_<n>") -- same class of gap as 03's SS_DESC/MCI-exclusion
  ## limitation, found while testing this 2026-09-24. parse_year_month()
  ## correctly ABORTS on this (that is the desired behaviour against real
  ## data with a genuinely wrong format); detected here so the script
  ## reports a clear, expected Tier-0 limitation instead of an unexplained
  ## crash.
  year_month_is_fixture_placeholder <-
    is_local && all(grepl("^FIXTURE_TEXT_", monthly_flag$YEAR_MONTH))

  if (year_month_is_fixture_placeholder) {
    cli::cli_alert_warning(
      "Local environment: {.field YEAR_MONTH} is a text column with no
       real-data levels in smidata's census, so the Tier-0 fixture emits
       unparseable placeholder text ({.val FIXTURE_TEXT_<n>}) rather than a
       real {.val YYYY-MM}-shaped value. compute_complete_case() cannot be
       exercised against this fixture for that reason -- correctly aborts
       if asked to try (see its own error) -- so it is NOT called here.
       This is a Tier-0 limitation, not a code defect; the same computation
       is exercised for real, at the true {config_complete_case_months}-month
       window, on hand-built data in the demonstration below."
    )
  } else {
    complete_case <- compute_complete_case(
      monthly_flag = monthly_flag,
      cohort = cohort_for_cc,
      window_months = config_complete_case_months
    )

    tally <- attr(complete_case, "tally")
    cli::cli_alert_info(
      "Complete-case: {tally$n_complete} of {tally$n_cohort}
       ({round(100 * tally$n_complete / tally$n_cohort, 1)}%). Insufficient
       coverage: {tally$n_insufficient_coverage}; absent from
       {.val {SMI_KEYS$medicaid_monthly_flag}} entirely: {tally$n_absent}."
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
    attrition_by_arm <- analysis_cohort |>
      dplyr::left_join(complete_case, by = "ID") |>
      dplyr::summarise(
        n = dplyr::n(), n_complete = sum(.data$complete_case, na.rm = TRUE),
        .by = "A"
      )
    cli::cli_alert_info("Complete-case rate by arm (A): report this table verbatim.")
    print(attrition_by_arm)
  }
}

## ---- Tier-0 demonstration: the real function, on a tiny hand-built fixture -
##
## Not a stand-in for real data (see application/README.md on why a
## hand-written fixture must never substitute for smi_fixture()'s
## census-backed one). Redundant with the real pipeline above when smidata
## IS installed (both exercise the same function); kept unconditional
## because it's the ONLY proof available when smidata is absent entirely,
## and because these two hand-built patients exercise exact known-truth
## values (fully covered vs. one month short) that a fixture's random
## MEDICAID_FLAG values don't guarantee to hit.

demo_cohort <- tibble::tibble(
  ID = c("demo_1", "demo_2"),
  INDEX_DT = as.Date(c("2019-06-01", "2019-06-01"))
)
demo_months <- format(
  seq(as.Date("2019-06-01"), by = "month", length.out = config_complete_case_months),
  "%Y-%m"
)
demo_monthly_flag <- dplyr::bind_rows(
  tibble::tibble(
    ID = "demo_1", INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
    YEAR_MONTH = demo_months, MEDICAID_FLAG = 1L
  ),
  ## demo_2: missing the LAST required month's row entirely, on purpose --
  ## demonstrates the function distinguishes "complete" from "one required
  ## month short" rather than only exercising the all-covered path.
  tibble::tibble(
    ID = "demo_2", INDEX_DT = as.Date("2019-06-01"), TYPE = "IP",
    YEAR_MONTH = demo_months[-length(demo_months)], MEDICAID_FLAG = 1L
  )
)

demo_cc <- compute_complete_case(
  demo_monthly_flag, demo_cohort, config_complete_case_months
)
print(demo_cc)
cli::cli_alert_info(
  "Demo (hand-built, {config_complete_case_months}-month window):
   demo_1 complete_case = {.val {demo_cc$complete_case[demo_cc$ID == 'demo_1']}}
   (fully covered); demo_2 complete_case =
   {.val {demo_cc$complete_case[demo_cc$ID == 'demo_2']}}
   (last required month absent)."
)

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info(
    "05_complete_case.R complete
     ({if (!config_has_smidata) 'demo only, no smidata' else if (exists('year_month_is_fixture_placeholder') && isTRUE(year_month_is_fixture_placeholder)) 'Tier-0 YEAR_MONTH limitation, see above; demo ran instead' else 'ran the real complete-case computation'})."
  )
}
