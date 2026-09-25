## ============================================================================
## application/_config.R
##
## Single source of gates, dataset keys, column contracts, analysis windows, and
## the OPEN_DECISIONS registry for doubletree's real-data application.
##
## Sourced by every numbered script (01-07, 90). Nothing here reads data.
##
## Configuration idiom (deliberate, matches smidata's inst/server/*.R
## convention as of 2026-09-17, and dual-bounds' analysis/real_data/_config.R):
## operational knobs are PLAIN TOP-LEVEL R VARIABLES, edited in place before a
## run -- NOT Sys.getenv(). Environment variables are reserved for smidata's own
## per-machine path overrides (SMI_DATA_DIR, SMI_CACHE_DIR, ANALYSIS_ENV, ...),
## which are a different kind of setting: they must resolve identically across
## five independent projects sharing one server. A per-run knob belongs to one
## invocation of one script and travels with the file, not the shell.
##
## RELATIONSHIP TO dual-bounds. This paper reuses dual-bounds' population,
## exposure, outcome and covariate SOURCE design verbatim (smidata
## inst/analyses/doubletree__application.yml, sections 2-6, mostly
## `status: inherited`). It differs in exactly two places: the ESTIMATOR
## (doubletree::estimate_att()/estimate_att_crossfit() instead of marbounds'
## partial-identification bounds), and OUTCOME MISSINGNESS (a complete-case
## restriction instead of a modelled response indicator R -- doubletree's
## estimator signature has no censoring argument for R to attach to).
##
## Inheriting a DESIGN is not the same as sharing CODE. Constants below that
## happen to equal dual-bounds' constants are still declared here, locally, and
## nothing in this directory reads dual-bounds' files. The shared IMPLEMENTATION
## lives in smidata >= 0.2.0 as exported smi_*() functions -- that is the one
## thing the two projects genuinely share, and it is shared by calling it, not
## by copying it.
## ============================================================================

## ---- 0 gates (edit these directly before a run) ----------------------------

## Emit the open-decision banner at the top of every numbered script.
config_echo_decisions <- TRUE

## Tiny smoke-check patient count for the pure-function calls that CAN run
## without real data (see 03_cost_windows.R's execution guard).
config_smoke_n_patients <- 2L

## Cost-claims file-suffix years present in the confirmed contract
## (aim3_svc_cost_16 .. _24 etc.). Verbatim from smidata snapshot
## 2026-09-15_2824080e: 4 families x 9 years = 36 dataset keys.
config_cost_years <- 16:24

## Which column of the cost-claims tables carries the service-setting
## description that the "Managed Care Invoice" exclusion matches on. Both
## SS_DESC and NEW_SERVICE_SETTINGS exist in all 36 files; SS_DESC is the
## declared default. Which one literally contains the string is a Tier-1 server
## check (90_checks_tier1.R), not something to paper over with column probing.
config_cost_setting_col <- "SS_DESC"

## Achievement threshold for AAP: MSR_NUM / MSR_DEN >= this counts as achieved.
## Confirmed PI 2026-09-17 (dual-bounds' exposure.definition, inherited here).
## Passed explicitly to smidata::smi_select_aap_row(), which exposes it because
## it means different things in different papers.
config_aap_achievement_threshold <- 0.5

## Literal msr_aap$MSR value for the AAP measure. DISTINCT from
## OPEN_DECISIONS$msr_code$value ("AAP (no filter needed)" -- a human-readable
## annotation of the resolution, not a filter value to plug into a query
## directly; filtering on that whole string would match zero real rows).
## 02_population_and_eligibility.R filters msr_aap to this literal before
## calling smi_select_aap_row() -- see that decision's note for why the
## registry's "abort-on-tie is the safety net" framing is not sufficient on
## its own, and why this filter (plus reporting any exclusion) is needed
## regardless of the "no filter needed" resolution's own intent.
AAP_MSR_CODE <- "AAP"

## Pinned ingest for a run against staged or published server data. NOT a
## per-run knob that gets cleared after use -- a pin persists even after the
## ingest it names PUBLISHES, because "newest" is a property of the frozen
## ingest cache at read time, not a property of this paper's design (see
## smidata::smi_read()'s own Details on pinned vs. unpinned reads).
##
## NA_character_, not NULL: smi_read_pinned() below tests it with is.na(), and
## smidata's own smi_run_start() normalizes a NULL ingest_id to NA_character_
## for the identical reason -- a length-1 character sentinel is type-stable
## (always character, never the length-0 NULL that would break a `||`
## boolean test below), and this file has no other convention for an unset
## SCALAR gate to follow instead (OPEN_DECISIONS's own `value = NULL`
## convention is for a different kind of registry -- a question with no
## defensible placeholder at all -- not an operational knob like this one).
config_ingest_id <- NA_character_   # <- edit before a run against staged/published data

## Every read in this pipeline goes through here -- never call
## smidata::smi_read() directly (enforced by
## application/tests/test-no-unpinned-reads.R). ingest_id is dropped whenever
## it would hit smi_read()'s own local-environment guard: that function
## ABORTS if given a non-NULL ingest_id while smi_env() == "local", because
## the ingest cache is server-only -- so an unset (NA) pin, or ANY pin at
## all while running locally, must never reach it. Tier-0 fixture runs must
## stay green with zero configuration.
smi_read_pinned <- function(dataset_key, columns = NULL, n = 200L) {
  smidata::smi_read(
    dataset_key, columns = columns, n = n,
    ingest_id = if (identical(smidata::smi_env(), "local") || is.na(config_ingest_id)) {
      NULL
    } else {
      config_ingest_id
    }
  )
}

## ---- estimator gates -------------------------------------------------------

## Which entry point 07_estimate_att.R treats as the FLAGSHIP result.
##   "att"      -> doubletree::estimate_att(), the paper's primary estimator.
##                 Valid only under grid-exact sparsity (ass:sparsity): a tree
##                 of at most leaf_budget leaves represents BOTH nuisances
##                 exactly on the analyst's pre-specified grid. Not checkable
##                 from data.
##   "crossfit" -> doubletree::estimate_att_crossfit(), the K-fold fallback,
##                 which does NOT require sparsity.
##
## Either way, 07 reports BOTH. See its own comments for why.
config_estimator <- "att"

## The leaf budget Lbar (estimate_att()'s `leaf_budget`, required, no default).
## A FIXED, analyst-chosen structural parameter, constant in n. 4L is the
## value doubletree's own README and test suite use as the worked example. It
## is an analysis decision, not a dataset fact, so it lives here.
config_leaf_budget <- 4L

## K for the cross-fitting companion estimate.
config_crossfit_k <- 5L

## Seed for the cross-fitting fold assignment. estimate_att() itself is
## deterministic given (X, A, Y); only the crossfit path randomizes.
config_seed <- 20260917L

## ---- smidata availability --------------------------------------------------

## This repo has no data and, on some machines, no smidata. Guard the version
## assertion so sourcing _config.R never hard-fails locally; the numbered
## scripts each report the absence explicitly rather than proceeding blind.
##
## The floor is 0.2.0, the release that PROMOTED the dataset-semantics helpers
## this pipeline needs: smi_parse_msr_yr(), smi_month_index(),
## smi_select_aap_row(), smi_classify_is_mci(), smi_patient_windows(),
## smi_window_coverage(), smi_stream_cost_files(), smi_sum_cost_in_windows(),
## and smi_read(). NOTHING in this directory reimplements or copies them; a
## local function with any of those names would be the drift the promotion
## exists to eliminate.
config_has_smidata <- requireNamespace("smidata", quietly = TRUE)
if (config_has_smidata) {
  stopifnot(utils::packageVersion("smidata") >= "0.2.0")
}

## ---- analysis windows (THIS PAPER's decision, not a dataset fact) ----------
##
## smidata::smi_patient_windows() deliberately has NO default `windows`
## argument: which months constitute the outcome ascertainment window and which
## the prior-covariate window is a per-paper research decision. This is that
## decision for this paper, and it is passed explicitly at every call site.
##
## Half-open: start INCLUSIVE, end EXCLUSIVE. prior_cost's end IS INDEX_DT, so
## the index day belongs to neither window and no claim is counted twice.
##
## Inherited from dual-bounds' confirmed design (outcome.ascertainment_window,
## covariates.measurement_window). The VALUES currently equal dual-bounds'
## DUAL_BOUNDS_WINDOWS. They are still declared here rather than read from that
## project: an inherited design decision that becomes a cross-project code
## dependency is a decision that changes silently in one repo and breaks another.
DOUBLETREE_WINDOWS <- tibble::tibble(
  window = c("outcome", "prior_cost"),
  start_months = c(12L, -12L),
  end_months = c(24L, 0L)
)

## The complete-case horizon (outcome.outcome_missingness, confirmed PI
## 2026-09-16): the estimation sample is restricted to patients with Y observed
## through INDEX_DT + this many months. It equals DOUBLETREE_WINDOWS's outcome
## end_months by construction -- observing Y means observing the whole outcome
## window -- and is named separately because 05_complete_case.R is about
## OBSERVABILITY, not about summing cost.
config_complete_case_months <- 24L

## ---- dataset keys and confirmed column contracts --------------------------

#' Dataset keys, as they appear in smidata's contract manifest.
#'
#' Verbatim from smidata snapshot 2026-09-15_2824080e. `cost_claims` holds
#' BASE keys; the real dataset keys append `_<yy>` for each of
#' `config_cost_years` (see `cost_dataset_keys()`).
#'
#' `medicaid_monthly_flag` was added 2026-09-24, when
#' `OPEN_DECISIONS$enrollment_source` resolved and 05_complete_case.R stopped
#' being a stub: it had never been declared here before, even though the
#' resolved decision's `note` already referenced it by name -- a real
#' registration gap, not an oversight in the decision itself.
#'
#' This paper requires NO datasets beyond the ones dual-bounds already
#' requires -- its only additional construction (prior_cost's quartile
#' discretization) is computed from data already being read, and
#' `medicaid_monthly_flag` is the SAME table dual-bounds' own (also still
#' pseudocode-only, as of 2026-09-24) response indicator `R` would read for
#' its identical, still-open registration gap -- see 05_complete_case.R's
#' header. This registers it here, not a claim that dual-bounds already has.
SMI_KEYS <- list(
  covariates            = "larger_smi_covariates",
  aap                   = "msr_aap",
  medicaid_monthly_flag = "larger_smi_medicaid_monthly_flag",
  cost_claims = list(
    aim3_svc   = "aim3_svc_cost",
    new_svc    = "new_svc_cost",
    aim3_pharm = "aim3_pharm_cost",
    new_pharm  = "new_pharm_cost"
  )
)

#' Columns this pipeline actually requires, by dataset.
#'
#' Confirmed against the real fingerprint, not inferred. `larger_smi_covariates`
#' has 72 columns and `msr_aap` has 5; these are the subsets the pipeline reads.
#' The 12 table-resident covariates are added on top -- see
#' `helpers/covariate_blocks.R`, which unlike dual-bounds' candidate map is
#' PI-CONFIRMED and therefore not an open decision.
#'
#' `medicaid_monthly_flag`'s five columns are the full confirmed schema (see
#' 90_checks_tier1.R CHECK 4 and smidata's inst/server/07_enrollment_coverage.R)
#' minus `COHORT`, which this pipeline's own use
#' (`helpers/complete_case.R::compute_complete_case()`) does not read.
SMI_COLS <- list(
  covariates            = c("ID", "INDEX_DT"),
  aap                   = c("ID", "MSR", "MSR_YR", "MSR_DEN", "MSR_NUM"),
  medicaid_monthly_flag = c("ID", "INDEX_DT", "TYPE", "YEAR_MONTH", "MEDICAID_FLAG"),
  cost_claims = c("ID", "SRV_DT", "AMOUNT_PAID", "SEQ_ID", config_cost_setting_col)
)

#' Expand the four cost-claims families across the confirmed suffix years.
#'
#' @param years Integer vector of two-digit file suffixes. Default
#'   `config_cost_years` (16:24), the confirmed span.
#' @return Character vector of 36 dataset keys, family-major.
cost_dataset_keys <- function(years = config_cost_years) {
  families <- unlist(SMI_KEYS$cost_claims, use.names = FALSE)
  as.vector(outer(families, years, function(f, y) paste0(f, "_", y)))
}

## ---- OPEN_DECISIONS registry ----------------------------------------------
##
## Every entry is a question that is NOT settled. `value` is either a
## PLACEHOLDER (loudly warned about on every access via open_value()) or NULL
## when no placeholder is even defensible. `blocking_final = TRUE` means no
## final artifact -- table, figure, or number in the manuscript -- may depend
## on it.
##
## SHARED ENTRIES. Two of these (cost_family_scope, msr_code) are the SAME
## question dual-bounds is asking, about the same files, because this paper
## inherits that part of the design verbatim. They carry `shared_with` and
## `registry_ref` fields rather than a re-derived question text: two projects
## independently wording the same open question is how the two answers end up
## differing. Resolving either resolves both, and the resolution belongs in
## smidata's registry, not in one of the two consuming repos.
##
## WHAT IS DELIBERATELY ABSENT, relative to dual-bounds' registry:
##   admin_gap_days  Superseded. There is no response indicator R here to apply
##                   an allowed-gap tolerance to; the complete-case restriction
##                   is a hard inclusion criterion, so the question collapses
##                   into enrollment_source itself.
##   msr_tie_break   Not applicable. This pipeline calls
##                   smidata::smi_select_aap_row() at its DEFAULT
##                   tie_break = "earlier" and passes no other value; smidata
##                   aborts on any other value anyway. There is no choice being
##                   made here to flag.
##   sl_lib          Not applicable. SuperLearner is marbounds' nuisance
##                   machinery. doubletree fits its nuisances with optimal trees
##                   via optimaltrees; there is no learner library to specify.
##
## LINE-LENGTH EXCEPTION (r-code-conventions §7). Several `question`/`note`
## strings below exceed 100 characters. They are single string LITERALS whose
## exact text is the registry's content: breaking one across source lines would
## embed a newline plus indentation into the recorded question.

OPEN_DECISIONS <- list(
  enrollment_source = list(
    question = "Table/columns providing continuous-enrollment spans through INDEX_DT+24mo, needed to determine COMPLETE-CASE STATUS (whether Y is observed for the whole outcome window). NO table in the confirmed design carries enrollment spans reaching that far: larger_smi_covariates is a covariate snapshot (72 columns, no span fields), msr_aap is measure periods (not enrollment), and the sibling aim1_smi_charac -- not elected by this paper -- encodes only ~12mo post-index enrollment.",
    status = "confirmed", blocking_final = FALSE,
    value = list(type_filter = "none -- all TYPE values counted", gap_months = 0L),
    shared_with = "dual-bounds",
    registry_ref = "smidata::inst/analyses/dual-bounds__application.yml -- IDENTICAL interim resolution, PI answered both papers in one sitting (2026-09-24)",
    note = "RESOLVED, INTERIM (PI, 2026-09-24, same finding as dual-bounds'
      identically-named entry). TYPE (IP/OP): PI explicitly deferred the
      IP-vs-OP choice -- 'note this as a question and just proceed ignoring
      this for now' -- so complete-case status is computed from
      larger_smi_medicaid_monthly_flag WITHOUT filtering on TYPE; every row
      counts as enrollment evidence regardless of IP/OP. PI separately
      flagged a prerequisite: verify no patient carries conflicting TYPE
      values before trusting that this is neutral -- see CHECK 4 in
      90_checks_tier1.R (server-only, not yet run; identical check to
      dual-bounds' CHECK 3). Gap-months: PI (2026-09-24) -- 'let's just not
      have a gap for now. We can change this later.' -- complete-case status
      now requires FULL continuity (0-month gap tolerance, MEDICAID_FLAG=1
      for every one of the ~24 months in [INDEX_DT, INDEX_DT+24mo]) as an
      explicit, PI-given, REVISABLE default, not an implementation guess.
      Both sub-decisions are provisional. 05_complete_case.R can now be
      implemented against these two values instead of remaining a stub."),
  cost_family_scope = list(
    question = "SHARED WITH dual-bounds -- see registry_ref. Whether all four cost-claims families are disjoint claim sources or aim3_*/new_* are overlapping extracts of the same claims (which would double-count cost if all four are summed).",
    status = "confirmed", blocking_final = FALSE, value = "all_four",
    shared_with = "dual-bounds",
    registry_ref = "smidata::inst/analyses/dual-bounds__application.yml#cost_family_scope",
    note = "RESOLVED (PI, 2026-09-17): 'They are disjoint subsets of patients - aim3* is the subset with 12-month enrollment after index and new* is without that enrollment guarantee.' A given patient's claims appear in exactly one family per claim type -- summing does NOT double-count, and does not shift prior_cost's quartile boundaries either, since no patient contributes to more than one family. IDENTICAL finding to dual-bounds' entry -- resolved once for both. Not yet independently verified against real claim-level data; the aim3_svc_cost_20-vs-new_svc_cost_20 comparison in dual-bounds' 90_checks_tier1.R remains worth running as corroboration."),
  msr_code = list(
    question = "SHARED WITH dual-bounds -- see registry_ref. Which msr_aap$MSR code value identifies the AAP measure specifically (msr_aap is not filtered to one measure, so 'the single nearest-month row' is ambiguous without it).",
    status = "confirmed", blocking_final = FALSE, value = "AAP (no filter needed)",
    shared_with = "dual-bounds",
    registry_ref = "smidata::inst/analyses/dual-bounds__application.yml#msr_code",
    note = "RESOLVED (PI, 2026-09-17): 'msr_code is always AAP in the msr_aap dataset' -- a measure-specific extract, not a combined all-measures file. Corroborated by contract evidence: msr_aap's MSR_NUM column carries the SAS label 'AAP_NUM'. IDENTICAL finding to dual-bounds' entry -- exposure A (inherited from dual-bounds) is now deterministic on this point. CORRECTION (2026-09-24, this repo, found while implementing 02_population_and_eligibility.R for real): the registry's claim that 'smidata::smi_select_aap_row()'s abort-on-tie behavior remains the safety net regardless' OVERSTATES that function's protection. Read smi_select_aap_row()'s actual source (~/RAND/tools/smidata/R/msr_aap.R): its winner-selection (candidates/winners) joins and groups on ID only and never references MSR at all except inside the tie-diagnostic message -- the abort fires ONLY when multiple rows tie at the exact same winning month_gap/row_month for one patient. A non-AAP row that is simply CLOSER in time to the target month than any true AAP row, with no tie, would be silently selected as the winner, corrupting elig_aap/aap_achieved for that patient with no abort and no warning. 02_population_and_eligibility.R now filters to AAP_MSR_CODE explicitly (below) and reports (does not silently trust) any exclusion, rather than relying on the abort-on-tie safety net this note previously (incorrectly) said was sufficient on its own. NOT YET propagated to smidata's actual dual-bounds__application.yml registry or to dual-bounds' own 02 -- both still carry the overstated claim and the identical silent-contamination exposure; see session_notes and the upstream item this should become."),
  diagnostics = list(
    question = "Which diagnostics gate a reportable ATT. estimate_att() returns sparsity-PROXY diagnostics (certified_e, certified_m0, used_search_e, used_search_m0, n_leaves_e, n_leaves_m0, gap_e, gap_m0) and never acts on them; which of them, at which thresholds, licenses reporting the flagship estimate rather than the cross-fit fallback is unspecified.",
    status = "open_decision", blocking_final = TRUE, value = NULL,
    note = "PARTIALLY RESOLVED (PI, 2026-09-24): the GENERAL diagnostics
      standard for this paper is now overlap/balance checks PLUS nuisance
      fit examined against OLS/GLM competitors (analysis_plan.diagnostics,
      smidata inst/analyses/doubletree__application.yml -- status: confirmed
      as of 2026-09-24). That answer does NOT resolve THIS entry's narrower
      question: which of estimate_att()'s own sparsity-proxy outputs
      (certified_e, gap_e, etc.), at which numeric thresholds, license
      reporting the flagship estimate rather than the cross-fit fallback.
      The underlying assumption -- grid-exact sparsity -- is not checkable
      from data, so no diagnostic can confirm it directly; the OLS/GLM
      nuisance-fit comparison may end up informing which proxy thresholds
      are trustworthy, but that is an inference to make once that
      comparison exists, not a PI statement resolving this entry now. Still
      open -- do not invent a threshold.")
)


#' Read an open decision's placeholder value, loudly.
#'
#' Every access warns. This is intentional: a placeholder that can be read
#' silently is a placeholder that ends up in a manuscript.
#'
#' @param id Name of an `OPEN_DECISIONS` entry.
#' @return The entry's `value` (possibly `NULL`).
open_value <- function(id) {
  d <- OPEN_DECISIONS[[id]]
  if (is.null(d)) cli::cli_abort("Unknown open decision {.val {id}}.")
  cli::cli_warn(c("!" = "Unconfirmed value for {.val {id}}: {.val {d$value}}", "i" = d$note))
  d$value
}

#' Read a RESOLVED (`status == "confirmed"`) open decision's value, quietly.
#'
#' [open_value()]'s "Unconfirmed value" warning is correct for a genuine
#' placeholder, but became misleading the moment `enrollment_source`,
#' `cost_family_scope`, and `msr_code` moved to `status = "confirmed"`
#' (2026-09-17/2026-09-24): calling [open_value()] on any of them would still
#' print "Unconfirmed" about a value the PI has actually given. This reads
#' the same `value` field without that warning, and aborts instead if the
#' entry is not actually confirmed -- so a status typo, or a decision that
#' later regresses to `open_decision`, cannot be read past silently either.
#'
#' @param id Name of an `OPEN_DECISIONS` entry.
#' @return The entry's `value`.
confirmed_value <- function(id) {
  d <- OPEN_DECISIONS[[id]]
  if (is.null(d)) cli::cli_abort("Unknown open decision {.val {id}}.")
  if (!identical(d$status, "confirmed")) {
    cli::cli_abort(c(
      "{.val {id}} is not confirmed (status: {.val {d$status}}).",
      "i" = "Use {.fn open_value} to read an unresolved placeholder -- it
             warns on every access, deliberately."
    ))
  }
  d$value
}

#' Print the full open-decision banner.
#'
#' Called at the top of every numbered script so no run is ever silent about
#' what is still unsettled. Prints ONCE per R session (a sentinel in
#' `globalenv()`), not once per call -- added 2026-09-24 alongside
#' `require_stage()`: a script cascading into several upstream scripts in one
#' session would otherwise print this banner once per cascade level, which is
#' noise, not new information, after the first time.
#'
#' @return `invisible(NULL)`, called for its side effect.
echo_open_decisions <- function() {
  if (isTRUE(get0(".doubletree_decisions_echoed", envir = globalenv(), ifnotfound = FALSE))) {
    return(invisible(NULL))
  }
  assign(".doubletree_decisions_echoed", TRUE, envir = globalenv())

  cli::cli_h2("Open decisions (unconfirmed -- see application/_config.R)")
  for (id in names(OPEN_DECISIONS)) {
    d <- OPEN_DECISIONS[[id]]
    ## cli has no cli_bullet(); cli_bullets() takes a named vector, "*" = plain
    ## bullet. Verified against cli 3.6.6.
    cli::cli_bullets(c("*" = paste0(
      id, " [", d$status,
      if (isTRUE(d$blocking_final)) ", BLOCKS FINAL ARTIFACT" else "",
      if (!is.null(d$shared_with)) paste0(", SHARED WITH ", d$shared_with) else "",
      "]: ", d$question
    )))
  }
  invisible(NULL)
}

#' Ensure upstream objects exist, sourcing the script that creates them if not.
#'
#' @description
#' Makes every numbered script genuinely runnable standalone
#' (`Rscript application/0N_....R`, no manual pre-sourcing) while ALSO
#' working efficiently when the whole pipeline is sourced in one session
#' (`run_pipeline.R`, or a human `source()`-ing 02, 03, 04, ... in order at
#' the console): if the named objects already exist in `globalenv()`, this
#' is a no-op; otherwise it sources `script` once and re-checks.
#'
#' @details
#' Re-sourcing `_config.R` and an upstream numbered script is safe here
#' because none of them read data unconditionally at the top level in a way
#' that depends on NOT having run before -- `_config.R` only defines
#' constants and functions, and each numbered script's real logic is
#' idempotent (re-running `02` twice produces the same `analysis_cohort`
#' from the same read, not a second copy appended to it).
#'
#' `source()`'s default `local = FALSE` evaluates in `globalenv()`, which is
#' where a script run via `Rscript` already sits at top level -- so an
#' object the cascade creates is visible to the caller with no extra
#' plumbing, the same way `_config.R`'s own constants already are.
#'
#' @param objects Character vector of object names required in `globalenv()`.
#' @param script Base filename (relative to `app_dir`, in scope where this is
#'   called) of the script that creates them.
#' @return `invisible(TRUE)`. Aborts, naming exactly which objects are STILL
#'   missing, if sourcing `script` did not produce them -- e.g. because
#'   `config_has_smidata` is `FALSE` and the creating script skipped its own
#'   real-data section entirely (see `02`'s own such guard).
require_stage <- function(objects, script) {
  missing <- objects[!vapply(
    objects, exists, logical(1), envir = globalenv(), inherits = FALSE
  )]
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }

  cli::cli_alert_info(
    "{.field {missing}} not yet in scope; sourcing {.file {script}} first."
  )
  ## app_dir is defined at the top of every numbered script, before that
  ## script's own source(_config.R) call -- it is in scope here because
  ## require_stage() is called from within that same script, after _config.R
  ## has been sourced into it.
  source(file.path(app_dir, script))

  still_missing <- objects[!vapply(
    objects, exists, logical(1), envir = globalenv(), inherits = FALSE
  )]
  if (length(still_missing) > 0L) {
    cli::cli_abort(c(
      "{.file {script}} did not produce {.field {still_missing}}.",
      "i" = "Two known causes: {.pkg smidata} is not installed (the script
             skips its real-data section entirely, via its own
             {.code config_has_smidata} guard); or a Tier-0 fixture-only
             limitation stopped it before producing this object (e.g.
             05_complete_case.R's YEAR_MONTH gap, 2026-09-24). Either way,
             {.file {script}}'s own console output above names which one --
             this abort does not re-derive it."
    ))
  }
  invisible(TRUE)
}

#' Refuse to proceed while any blocking open decision is unresolved.
#'
#' Called only where a REPORTABLE final result would otherwise be produced:
#' `07_estimate_att.R`, before any estimator call. NOT called from
#' `06_assemble_analytic_data.R` (moved out 2026-09-24, design-reviewed by
#' `oracle`): `06` assembles the analytic DATASET, and `diagnostics` -- the
#' one decision still `blocking_final = TRUE` as of this writing -- is
#' specifically about which ESTIMATOR result to trust in `07`, not about
#' whether the assembled data is correct. `06` instead branches on
#' [has_blocking_open()] to choose a `_PROVISIONAL` filename suffix, so an
#' inspectable intermediate artifact is reachable while the REPORTABLE
#' (non-suffixed) one still is not.
#'
#' @return `invisible(TRUE)` if nothing blocks; otherwise aborts.
assert_no_blocking_open <- function() {
  blocking <- vapply(OPEN_DECISIONS, function(d) isTRUE(d$blocking_final), logical(1))
  if (any(blocking)) {
    cli::cli_abort(c(
      "Cannot proceed: {sum(blocking)} blocking open decision(s) unresolved.",
      stats::setNames(
        paste0("{.val ", names(OPEN_DECISIONS)[blocking], "}: ",
               vapply(OPEN_DECISIONS[blocking], function(d) d$question, character(1))),
        rep("x", sum(blocking))
      )
    ))
  }
  invisible(TRUE)
}

#' TRUE iff at least one blocking open decision remains.
#'
#' Lets a script branch on the gate without triggering the abort (e.g. to
#' choose a `_PROVISIONAL` output suffix).
#'
#' @return Logical scalar.
has_blocking_open <- function() {
  any(vapply(OPEN_DECISIONS, function(d) isTRUE(d$blocking_final), logical(1)))
}

#' Names of the blocking open decisions.
#'
#' @return Character vector.
blocking_open_ids <- function() {
  names(OPEN_DECISIONS)[
    vapply(OPEN_DECISIONS, function(d) isTRUE(d$blocking_final), logical(1))
  ]
}
