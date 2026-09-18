## ============================================================================
## 05_complete_case.R
##
## STUB, AND DELIBERATELY SO. COMPLETE-CASE STATUS IS NOT COMPUTABLE FROM THE
## CONFIRMED DESIGN.
##
## The design says (outcome.outcome_missingness, status: confirmed, PI
## 2026-09-16):
##   "Complete-case restriction: estimation sample restricted to patients with
##    Y observed through INDEX_DT+24mo."
##
## That is a well-formed research decision. It is also, right now, an
## unanswerable question, for the same reason dual-bounds cannot compute its
## response indicator R: THE ENROLLMENT SPANS ARE MISSING ENTIRELY. No table in
## the confirmed design carries enrollment spans reaching INDEX_DT+24mo:
##
##   larger_smi_covariates  a covariate SNAPSHOT -- 72 columns, zero span or
##                          enrollment-date fields.
##   msr_aap                MEASURE periods (RY<yyyy>-<mm>), not enrollment. A
##                          patient can be measure-eligible in a period without
##                          that telling you their enrollment span.
##   aim1_smi_charac        the sibling table, not elected by this paper, and in
##                          any case encodes only ~12mo post-index enrollment --
##                          not +24mo. (That filter is exactly what makes it a
##                          strict subset of larger_smi_covariates; see smidata
##                          D001, resolved 2026-09-17.)
##   the cost-claims family CLAIMS, not enrollment.
##
## NOTE WHAT IS AND IS NOT INHERITED HERE. dual-bounds MODELS missingness (it is
## the subject of its method); this paper EXCLUDES it. So the two papers differ
## in what they do with the answer -- but they need the identical missing input,
## which is why OPEN_DECISIONS$enrollment_source is worded for this paper's use
## and still describes the same absent table. Resolving it unblocks both.
##
## WHAT THIS SCRIPT MUST NOT DO: infer complete-case status from claims presence.
## "Had a claim in month 23, therefore observed through month 24" is a proxy, and
## it fails in two directions at once, both of which bias the ATT rather than
## merely adding noise:
##
##   (a) FALSE NEGATIVES -- it DROPS genuine zero-cost complete cases. A patient
##       enrolled the whole 24 months who simply used no services has no claim in
##       the window, so a claims-presence proxy calls them incomplete and removes
##       them. Y = 0 is a legitimate, informative outcome value; deleting the
##       low tail of a cost outcome shifts the estimated mean upward in whichever
##       arm has more of it. Since the exposure IS medication adherence, the arms
##       do not have equal shares of zero-utilization patients.
##
##   (b) FALSE POSITIVES -- it KEEPS partially-enrolled patients. A patient
##       enrolled for 6 of the 24 months, with claims in month 5, is called
##       complete, and their 6-month partial sum enters as if it were a 24-month
##       total. Their Y is not missing and not censored; it is a low number that
##       looks exactly like a genuinely low complete one. Nothing downstream can
##       distinguish them.
##
## (b) is the worse of the two: (a) at least removes rows visibly, while (b)
## silently substitutes a shorter accumulation window for the declared one. Both
## make the proxy a fabricated inclusion criterion, not a conservative
## approximation of the right one.
##
## Run time: instant. It computes nothing.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("05 -- complete-case restriction (STUB)")

cli::cli_alert_danger(
  "Complete-case status cannot be computed. This is a data-availability gap, not
   a threshold question."
)
cli::cli_alert_danger(OPEN_DECISIONS$enrollment_source$question)
cli::cli_alert_danger(OPEN_DECISIONS$enrollment_source$note)
cli::cli_alert_warning(
  "Required horizon: {.field INDEX_DT + {config_complete_case_months}mo}. No
   elected table reaches it."
)

#' Complete-case indicator for a cohort: `NA` for every patient, by design.
#'
#' Returns the correct SHAPE so `06_assemble_analytic_data.R` can be written and
#' read, while returning no VALUE, because no value is derivable. Every
#' downstream consumer therefore fails on `NA` rather than succeeding on a
#' fabrication -- and in this pipeline the failure is guaranteed, not merely
#' likely: `check_att_data()` rejects any `NA` in `X`, `A` or `Y`, and the
#' restriction is applied by FILTERING on this column, so an all-`NA` column
#' yields an empty analytic sample rather than a full one.
#'
#' @param ids Patient IDs, one element per patient.
#' @return Tibble with `ID` and `complete_case` (`NA` throughout).
complete_case_stub <- function(ids) {
  if (length(ids) == 0L) {
    cli::cli_abort("{.arg ids} is empty; nothing to build complete-case status for.")
  }
  cli::cli_alert_danger(
    "Returning complete_case = NA for all {length(ids)} patient{?s}. No estimate
     may be reported on this basis."
  )
  tibble::tibble(ID = ids, complete_case = NA)
}

## ---- WHAT THE REAL IMPLEMENTATION WOULD LOOK LIKE ---------------------------
##
## Written out so that resolving OPEN_DECISIONS$enrollment_source is a
## substitution rather than a design task. It needs a table this pipeline does
## not currently have: one row per (patient, enrollment span), with
## enroll_start / enroll_end covering at least
## [INDEX_DT, INDEX_DT + config_complete_case_months].
##
## complete_case <- enrollment_spans |>
##   dplyr::semi_join(analysis_cohort, by = "ID") |>
##   dplyr::inner_join(
##     dplyr::transmute(analysis_cohort,
##                      ID = .data$ID,
##                      cc_start = as.Date(.data$INDEX_DT),
##                      cc_end = lubridate::add_with_rollback(
##                        as.Date(.data$INDEX_DT),
##                        lubridate::period(months = config_complete_case_months)
##                      )),
##     by = "ID"
##   ) |>
##   ## Clip each span to the required interval, then check TOTAL uncovered days
##   ## rather than the largest single gap -- three 20-day gaps are not "briefly
##   ## disenrolled" even though no single one is long.
##   dplyr::mutate(
##     covered_start = pmax(.data$enroll_start, .data$cc_start),
##     covered_end   = pmin(.data$enroll_end, .data$cc_end)
##   ) |>
##   dplyr::filter(.data$covered_end > .data$covered_start) |>
##   dplyr::summarise(
##     covered_days = sum(as.integer(.data$covered_end - .data$covered_start)),
##     required_days = as.integer(dplyr::first(.data$cc_end - .data$cc_start)),
##     .by = "ID"
##   ) |>
##   dplyr::mutate(complete_case = .data$covered_days >= .data$required_days)
##
## NOTE: no allowed-gap tolerance appears above, and that absence is the point.
## dual-bounds registers admin_gap_days because its R is a MODELLED indicator
## whose definition can absorb a tolerance. Here the restriction is a hard
## inclusion criterion, so "complete" means complete; permitting a gap would make
## the retained patients' Y a sum over a window with holes in it, which is
## precisely failure mode (b) in this file's header, reintroduced by policy
## instead of by accident. If a tolerance is wanted it is a NEW PI decision with
## its own registry entry, not a default.
##
## SEPARATELY, AND ALSO UNRESOLVED: death. There is no death source in the
## confirmed design, so a patient who dies at month 18 is indistinguishable from
## one who disenrolls at month 18, and both are excluded as incomplete. Excluding
## deaths is a substantive choice about the estimand -- the ATT is then defined
## among 24-month survivors -- and the manuscript must say so rather than
## treating it as a data-cleaning step.

if (!interactive() && sys.nframe() == 0L) {
  demo_cc <- complete_case_stub(c("demo_1", "demo_2"))
  print(demo_cc)
  cli::cli_alert_info(
    "05_complete_case.R complete (stub; complete_case is NA by design)."
  )
}
