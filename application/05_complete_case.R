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
## WHAT STILL DOES NOT RUN TODAY, AND WHY THAT IS A DIFFERENT GAP THAN BEFORE.
## The real orchestration below (Section "INTENDED PIPELINE") needs
## `analysis_cohort` (ID, INDEX_DT) from 02_population_and_eligibility.R --
## which, as of this writing, does not yet produce one for real: msr_code
## resolved (PI, 2026-09-17, OPEN_DECISIONS$msr_code) but 02's own gate check
## has not been updated to read that resolution, so 02 still aborts on its
## msr_code gate rather than reaching its (already-written) real assembly.
## THIS IS A SEPARATE, PRE-EXISTING GAP IN 02, DISCOVERED WHILE IMPLEMENTING
## THIS FILE, NOT SOMETHING THIS FILE CAUSES OR FIXES. Until it is closed,
## this script's own gate (below) reports its readiness honestly and stops
## before attempting a server read that has nothing real to read against.
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
## Run time: instant against a fixture; hours against the real table
## (5.8GB / 24M rows -- see ~/RAND/tools/smidata/inst/server/
## 07_enrollment_coverage.R's header cost warning, which applies identically
## here) once the analysis_cohort gap above is closed.
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

## ---- INTENDED PIPELINE (needs analysis_cohort from a real 02; not reachable
## today -- see this file's header for why that is a separate, pre-existing
## gap, not something to work around here) ------------------------------------
##
## smidata::smi_require(
##   SMI_KEYS$medicaid_monthly_flag,
##   columns = SMI_COLS$medicaid_monthly_flag
## )
## monthly_flag <- smidata::smi_read(SMI_KEYS$medicaid_monthly_flag,
##                                    columns = SMI_COLS$medicaid_monthly_flag)
##
## complete_case <- compute_complete_case(
##   monthly_flag = monthly_flag,
##   cohort = dplyr::select(analysis_cohort, "ID", "INDEX_DT"),   # from 02
##   window_months = config_complete_case_months
## )
##
## tally <- attr(complete_case, "tally")
## cli::cli_alert_info(
##   "Complete-case: {tally$n_complete} of {tally$n_cohort}
##    ({round(100 * tally$n_complete / tally$n_cohort, 1)}%). Insufficient
##    coverage: {tally$n_insufficient_coverage}; absent from
##    {.val {SMI_KEYS$medicaid_monthly_flag}} entirely: {tally$n_absent}."
## )
## ## Differential attrition by arm is the substantive threat a 24-month
## ## continuous-enrollment filter poses (Oracle-flagged, 2026-09-24 design
## ## review): retention plausibly correlates with the exposure itself
## ## (medication adherence). Report it, don't just compute it.
## attrition_by_arm <- analysis_cohort |>
##   dplyr::left_join(complete_case, by = "ID") |>
##   dplyr::summarise(
##     n = dplyr::n(), n_complete = sum(complete_case, na.rm = TRUE),
##     .by = "A"
##   )
## cli::cli_alert_info("Complete-case rate by arm (A): report this table verbatim.")
## print(attrition_by_arm)

## ---- Tier-0 demonstration: the real function, on a tiny hand-built fixture -
##
## Not a stand-in for real data (see application/README.md on why a
## hand-written fixture must never substitute for smi_fixture()'s
## census-backed one) -- this exists only to prove the function this script
## now calls is real, callable code, exercised at the real
## config_complete_case_months window length, the same way
## 04_covariates.R's checks run today without a server.

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
    "05_complete_case.R complete (real logic; server read still blocked --
     see this file's header)."
  )
}
