## ============================================================================
## 07_estimate_att.R
##
## Tree-based ATT estimation via doubletree::estimate_att() (flagship, no sample
## splitting) and doubletree::estimate_att_crossfit() (K-fold fallback).
##
## Signatures below are REAL, taken from R/estimate_att.R and
## R/estimate_att_crossfit.R in this package. `leaf_budget` is required and has
## no default; `outcome_type = "continuous"` because Y is dollars.
##
## THE GATE IS NOW OPEN. As of 2026-09-30 all five OPEN_DECISIONS entries are
## `status = "confirmed"` with no `blocking_final`, so assert_no_blocking_open()
## below PASSES rather than aborting -- kept as a live regression guard, not as
## a known stop. Last to resolve was `diagnostics`: certified_e AND
## certified_m0 both TRUE means this proxy supports grid-exact sparsity for
## the flagship estimate (PI, 2026-09-30). This is a REPORTED FLAG, NOT AN
## ESTIMATOR SWITCH -- there is no fallback to estimate_att_crossfit(), silent
## or explicit; the flagship stays the flagship (per config_estimator) and
## crossfit stays its companion regardless of whether the licence holds. This
## script prints whether it holds and still never switches estimator on it.
##
## PREREQUISITE, newly binding: this reads the UN-SUFFIXED
## output/analytic_cohort.parquet, which 06 writes only once nothing blocks.
## Any analytic_cohort_PROVISIONAL.parquet on disk predates the resolution --
## re-run 06 first.
##
## WHY BOTH ESTIMATORS ARE ALWAYS REPORTED, regardless of config_estimator.
## estimate_att()'s validity rests on grid-exact sparsity (theory.tex
## ass:sparsity): a tree of at most `leaf_budget` leaves represents BOTH
## nuisances exactly on the 2^15-atom grid the 15 binary covariates generate.
## That assumption is NOT checkable from data -- estimate_att() reports proxies
## (certified_*, n_leaves_*, gap_*) and, by explicit design, never switches
## estimator on the basis of them. So if sparsity fails, estimate_att() does not
## fail; it returns a confidently-stated wrong number, because the assumption
## that licensed the plain Wald interval was false.
##
## estimate_att_crossfit() does not need sparsity. Running both and reporting the
## DISCREPANCY converts a silent bias into a visible disagreement. It is not a
## test of sparsity -- the two estimators differ for other reasons too (in-sample
## versus cross-fit nuisances, different regularization machinery) -- but a large
## gap is informative and a hidden gap is not. This is also why 07 does not
## auto-select: doing so would make estimator choice data-dependent
## post-selection inference the paper does not analyse.
##
## Run time when unblocked: both paths fit optimal trees via optimaltrees (GOSDT).
## config_leaf_budget/config_max_depth = 16L/4L (PI+Oracle, 2026-09-30,
## superseding 2026-09-28's 15L/default-14) guarantees the flagship's
## bisection needs exactly ONE fit per nuisance (16 = 2^4, so the budget
## cannot be overshot) -- see _config.R for the full rationale and the
## real-data-hang trace that motivated it. The crossfit path fits 2K trees
## rather than 2, and its cv_regularization sweep multiplies that again;
## it is NOT covered by this fix (a separate, confirmed root-loss-scale bug
## in optimaltrees::cv_regularization_adaptive() remains open there,
## crossfit explicitly de-prioritized 2026-09-30) -- budget its own timing
## as unknown, not a scaled-up version of anything already measured.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("07 -- tree-based ATT estimation")

## ---- the gate, before anything else -----------------------------------------
##
## Deliberately ahead of library(doubletree) and ahead of reading anything: the
## point is to refuse, not to get as far as possible before refusing.
assert_no_blocking_open()

## Fail before the (long) fits, not after them.
if (!(is.character(config_estimator) && length(config_estimator) == 1L &&
      config_estimator %in% c("att", "crossfit"))) {
  cli::cli_abort(
    "{.code config_estimator} must be {.val att} or {.val crossfit}, not
     {.val {config_estimator}}."
  )
}

source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))

library(doubletree)

set.seed(config_seed)

analytic_path <- file.path(app_dir, "output", "analytic_cohort.parquet")
if (!file.exists(analytic_path)) {
  cli::cli_abort(c(
    "No REPORTABLE analytic data at {.file {analytic_path}}.",
    "i" = "06_assemble_analytic_data.R writes this un-suffixed name only once
           nothing blocks. Any {.file analytic_cohort_PROVISIONAL.parquet}
           beside it predates {.val diagnostics}' 2026-09-30 resolution --
           re-run {.file 06_assemble_analytic_data.R}. Do not rename the
           provisional file."
  ))
}
analytic_data <- arrow::read_parquet(analytic_path)

x_cols <- design_matrix_columns()
X <- as.data.frame(analytic_data[, x_cols, drop = FALSE])
A <- as.integer(analytic_data$A)
Y <- as.numeric(analytic_data$Y)

## Re-assert at the estimator boundary. Not redundant with 06: the parquet
## round-trip is where a 0/1 integer column can come back as something else, and
## column ORDER is not guaranteed by a file format. See helpers/design_matrix.R
## for what the package itself does and does not check.
assert_binary_design_matrix(X)

cli::cli_alert_info(
  "n = {nrow(X)}; A = 1: {sum(A == 1L)}; A = 0: {sum(A == 0L)};
   Y in [{round(min(Y))}, {round(max(Y))}]"
)

## ---- flagship: estimate_att() (grid-exact sparsity required) ----------------

fit_att <- doubletree::estimate_att(
  X = X, A = A, Y = Y,
  leaf_budget = config_leaf_budget,
  outcome_type = "continuous",
  max_depth = config_max_depth,
  depth_restricted = TRUE,   # required whenever max_depth < leaf_budget - 1
  model_limit = 0L           # undoes fit_tree()'s worst-case-continuous
                              # feature-count heuristic (real X is already
                              # binary; see application/README.md's `07` row)
)

## Checkpoint BEFORE the crossfit companion: that path's timing is unknown and
## an open optimaltrees bug affects it (header). If it errors or is killed, the
## flagship fit, which may itself have taken hours, is not lost. Failures stay
## loud (no tryCatch); this only preserves work. Server-only, patient-level.
fs::dir_create(file.path(app_dir, "output"))
readr::write_rds(fit_att, file.path(app_dir, "output", "att_fit_flagship.rds"))

cli::cli_h2("estimate_att() -- flagship, requires grid-exact sparsity")
cli::cli_alert_info("theta = {signif(fit_att$theta, 4)}  (SE {signif(fit_att$sigma, 4)})")
cli::cli_alert_info("95% CI: [{signif(fit_att$ci_95[1], 4)}, {signif(fit_att$ci_95[2], 4)}]")

## Sparsity PROXIES, reported and never acted on automatically. Which of these
## licenses reporting this estimate is OPEN_DECISIONS$diagnostics, RESOLVED
## 2026-09-30 and evaluated explicitly below. Note what they can and cannot
## say: certified_* means the returned fit provably solves eq:select, i.e. it
## is the right tree for the budget; it does NOT mean the budget is large
## enough for the truth.
cli::cli_alert_info(
  "certified_e = {fit_att$certified_e}, certified_m0 = {fit_att$certified_m0};
   leaves: e = {fit_att$n_leaves_e}, m0 = {fit_att$n_leaves_m0} (budget
   {fit_att$leaf_budget})"
)
cli::cli_alert_info(
  "used_search_e = {fit_att$used_search_e}, used_search_m0 =
   {fit_att$used_search_m0}; gaps: {signif(fit_att$gap_e, 3)},
   {signif(fit_att$gap_m0, 3)}"
)
cli::cli_alert_info(
  "any_truncated_e = {fit_att$any_truncated_e}, any_truncated_m0 =
   {fit_att$any_truncated_m0}; clip_rate_e = {signif(fit_att$clip_rate_e, 3)}"
)
if (isTRUE(fit_att$any_truncated_e) || isTRUE(fit_att$any_truncated_m0)) {
  cli::cli_alert_warning(
    "A truncated fit's incumbent is not proven optimal -- this breaks the
     premise certified_* relies on for that tree even if certified_* itself
     still reads TRUE. Treat this as its own failure mode, distinct from
     certified_* == FALSE, not something the licence below already covers."
  )
}
## At config_leaf_budget >= 2^config_max_depth (16L/4L here, PI+Oracle,
## 2026-09-30) the budget cannot be overshot, so the FIRST fit is always the one
## returned and bisection never runs. certified_* then reduces to whether that
## first fit satisfies the m_n leaf-mass floor, within the depth-limited class:
## it still depends on the data, but it says NOTHING about whether the budget
## is large enough for the true nuisances. (If the first fit violates the
## floor it is still used, per optimaltrees' bisect_lambda_to_budget(); then
## certified_* is FALSE.) See _config.R's config_leaf_budget comment for the
## 15L/max_depth=4L alternative that keeps the budget check non-trivial, and
## the probe that should decide between them before this is a reportable run.
## Also worth stating in the write-up rather than discovering in review: this
## path's `sigma` is the in-sample EIF plug-in and is downward biased by a
## relative factor of order (2 * leaf_budget + 1) / n, so the nominal 95%
## interval can undercover at small n. estimate_att()'s own documentation says
## so; it is not implemented away because the fix would change the shared
## att_se() used by all three entry points.

## The CONFIRMED diagnostics rule (OPEN_DECISIONS$diagnostics, PI 2026-09-30),
## reported and still never acted on. Reporting the licence is not switching on
## it: the flagship below stays config_estimator's, because auto-switching is
## the post-selection inference this file's header rules out.
flagship_licensed <- isTRUE(fit_att$certified_e) && isTRUE(fit_att$certified_m0)
budget_check_vacuous <- config_leaf_budget >= 2L^config_max_depth
cli::cli_alert_info(
  "Flagship licensed under the confirmed diagnostics rule
   ({.val {confirmed_value('diagnostics')}}): {.val {flagship_licensed}}"
)
if (flagship_licensed) {
  cli::cli_alert_info(
    "Licence holds -- NECESSARY, not sufficient: certified_* says the search
     provably solved eq:select at leaf_budget = {fit_att$leaf_budget}, not that
     the budget is large enough for the true nuisances. The 2026-09-24 general
     standard (overlap/balance, nuisance fit versus OLS/GLM) still applies."
  )
  if (budget_check_vacuous) {
    cli::cli_alert_danger(
      "leaf_budget {config_leaf_budget} >= 2^max_depth
       ({2L^config_max_depth}): the budget cannot be overshot, so certified_*
       reduces to the first fit meeting the m_n leaf-mass floor within the
       depth-{config_max_depth} class. It carries no information about budget
       adequacy."
    )
  }
} else {
  cli::cli_alert_warning(
    "Licence FAILS -- certified_e and certified_m0 do not both hold, so this
     proxy does not support grid-exact sparsity for the flagship estimate.
     There is NO fallback to estimate_att_crossfit(), silent or explicit:
     estimate_att() remains the flagship (per config_estimator below) and
     estimate_att_crossfit() remains its companion, unconditionally. A failed
     licence is a signal to scrutinize the discrepancy reported below, not an
     instruction to report estimate_att_crossfit()'s number in its place."
  )
}

## ---- companion: estimate_att_crossfit() (no sparsity requirement) -----------
##
## Always run, whichever config_estimator names as flagship. See this file's
## header for why: a sparsity failure that is invisible in the flagship estimate
## shows up here as a discrepancy.

fit_cf <- doubletree::estimate_att_crossfit(
  X = X, A = A, Y = Y,
  K = config_crossfit_k,
  outcome_type = "continuous",
  seed = config_seed
)

cli::cli_h2("estimate_att_crossfit() -- companion, no sparsity requirement")
cli::cli_alert_info("theta = {signif(fit_cf$theta, 4)}  (SE {signif(fit_cf$sigma, 4)})")
cli::cli_alert_info("95% CI: [{signif(fit_cf$ci_95[1], 4)}, {signif(fit_cf$ci_95[2], 4)}]")

## ---- the discrepancy, reported explicitly -----------------------------------

delta <- fit_att$theta - fit_cf$theta
pooled_se <- sqrt(fit_att$sigma^2 + fit_cf$sigma^2)
cli::cli_h2("Discrepancy")
cli::cli_alert_info(
  "estimate_att() - estimate_att_crossfit() = {signif(delta, 4)}
   ({signif(delta / pooled_se, 2)} pooled SEs)"
)
cli::cli_alert_warning(
  "This is a DISAGREEMENT MEASURE, not a hypothesis test. The two estimates share
   the same data, so pooled_se overstates the standard error of their difference;
   and the estimators differ in more than the sparsity assumption (in-sample
   versus cross-fit nuisances, different regularization). A large value says
   'investigate', never 'sparsity is rejected'."
)

flagship <- switch(
  config_estimator,
  att = fit_att,
  crossfit = fit_cf,
  cli::cli_abort(
    "{.code config_estimator} must be {.val att} or {.val crossfit}, not
     {.val {config_estimator}}."
  )
)
cli::cli_alert_info(
  "Flagship per {.code config_estimator = \"{config_estimator}\"}: theta =
   {signif(flagship$theta, 4)}"
)

## ---- Persist the results ----------------------------------------------------
##
## Everything the manuscript's Empirical Application section reports is written
## here, with the run's provenance, so no number depends on a console session.
## The .tex tables and nuisance-tree figure are built from these files
## downstream; this script does not typeset anything.
##   output/att_fit_flagship.rds      checkpoint written right after estimate_att()
##   output/att_estimates.parquet     one row per estimator: theta, SE, CI, n
##   output/att_diagnostics.parquet   one row: sparsity proxies + licence flags
##   output/att_provenance.csv        config, versions, git state, data hash
##   output/att_fits.rds              both fit objects, plus sessionInfo()
##
## SERVER-ONLY. The .rds files hold patient-level quantities (score values,
## fitted propensities, outcome vectors) and the estimates files are only as
## safe as their small cells; nothing in output/ leaves the server except
## through smidata::smi_export(). The fitted tree structure is NOT persisted
## in a tidy form here: n_leaves_* in the diagnostics are all that is saved
## outside the RDS, so the nuisance-tree figure still needs its own step.

out_dir <- file.path(app_dir, "output")
fs::dir_create(out_dir)

licence_vacuous <- budget_check_vacuous   # set with the licence above

ci_of <- function(fit) {
  tibble::tibble(ci_lo = fit$ci_95[[1]], ci_hi = fit$ci_95[[2]])
}
att_estimates <- dplyr::bind_rows(
  dplyr::bind_cols(
    tibble::tibble(estimator = "estimate_att", theta = fit_att$theta,
                   sigma = fit_att$sigma),
    ci_of(fit_att)
  ),
  dplyr::bind_cols(
    tibble::tibble(estimator = "estimate_att_crossfit", theta = fit_cf$theta,
                   sigma = fit_cf$sigma),
    ci_of(fit_cf)
  )
) |>
  dplyr::mutate(
    flagship = (.data$estimator == "estimate_att") == (config_estimator == "att"),
    ## The saved sigma/CI carry these caveats into any table built from them.
    se_method = c("in-sample EIF plug-in; downward-biased, O((2*leaf_budget+1)/n)",
                  "cross-fit EIF"),
    ci_method = "Wald", ci_level = 0.95,
    n = nrow(X), n_treated = sum(A == 1L), n_control = sum(A == 0L)
  )

## tibble() silently drops NULL arguments, so a field a future doubletree
## stops returning would vanish from the output. Require every one to be a
## length-1 value instead.
diag_fields <- c(
  "certified_e", "certified_m0", "used_search_e", "used_search_m0",
  "n_leaves_e", "n_leaves_m0", "gap_e", "gap_m0", "any_truncated_e",
  "any_truncated_m0", "clip_rate_e",
  # The package-default values actually used (not set in _config.R), recorded
  # so a package upgrade that changes a default is visible in the output.
  "m_n", "lambda_n", "lambda_e", "lambda_m0", "clip_c"
)
bad_fields <- diag_fields[lengths(fit_att[diag_fields]) != 1L]
if (length(bad_fields) > 0L) {
  cli::cli_abort("estimate_att() did not return length-1 {.field {bad_fields}}.")
}
att_diagnostics <- tibble::as_tibble(fit_att[diag_fields]) |>
  dplyr::mutate(
    leaf_budget = config_leaf_budget, max_depth = config_max_depth,
    depth_restricted = TRUE,
    depth_sufficient = config_max_depth >= config_leaf_budget - 1L,
    flagship_licensed = flagship_licensed,
    budget_check_vacuous = licence_vacuous,
    delta = delta,
    pooled_se_upper_bound_not_a_test = pooled_se
  )

## git state of each code base the numbers depend on. Recorded, never silently
## skipped: where git is unavailable the value says so.
git_state <- function(dir) {
  sha <- tryCatch(
    suppressWarnings(system2("git", c("-C", dir, "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  if (length(sha) != 1L) return("git_unavailable")
  dirty <- length(suppressWarnings(
    system2("git", c("-C", dir, "status", "--porcelain"), stdout = TRUE,
            stderr = FALSE)
  )) > 0L
  paste0(sha, if (dirty) "+dirty" else "")
}
pkg_dir <- function(pkg) tryCatch(find.package(pkg), error = function(e) ".")
prov <- c(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  ingest_id = config_ingest_id, smi_run_id = Sys.getenv("SMI_RUN_ID", "none"),
  r_version = R.version.string,
  seed = config_seed, estimator_flagship = config_estimator,
  crossfit_k = config_crossfit_k, saa_max_month_gap = config_saa_max_month_gap,
  saa_achievement_threshold = config_saa_achievement_threshold,
  cost_setting_col = config_cost_setting_col,
  complete_case_months = config_complete_case_months,
  coverage_required_windows = paste(config_coverage_required_windows, collapse = ","),
  windows = paste(DOUBLETREE_WINDOWS$window, DOUBLETREE_WINDOWS$start_months,
                  DOUBLETREE_WINDOWS$end_months, sep = ":", collapse = ","),
  analytic_data = analytic_path,
  analytic_md5 = unname(tools::md5sum(analytic_path)),
  doubletree_version = as.character(utils::packageVersion("doubletree")),
  optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
  smidata_version = as.character(utils::packageVersion("smidata")),
  application_git = git_state(app_dir),
  doubletree_git = git_state(pkg_dir("doubletree")),
  optimaltrees_git = git_state(pkg_dir("optimaltrees")),
  smidata_git = git_state(pkg_dir("smidata"))
)
if (any(prov[grepl("_git$", names(prov))] == "git_unavailable")) {
  cli::cli_alert_warning(
    "git state could not be read for some code bases; recorded as
     {.val git_unavailable} in att_provenance.csv."
  )
}
att_provenance <- tibble::tibble(field = names(prov), value = unname(as.character(prov)))

arrow::write_parquet(att_estimates, file.path(out_dir, "att_estimates.parquet"))
arrow::write_parquet(att_diagnostics, file.path(out_dir, "att_diagnostics.parquet"))
readr::write_csv(att_provenance, file.path(out_dir, "att_provenance.csv"))
readr::write_rds(
  list(fit_att = fit_att, fit_cf = fit_cf, session = utils::sessionInfo()),
  file.path(out_dir, "att_fits.rds")
)
cli::cli_alert_success("Wrote the 07 result files to {.file {out_dir}}.")

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("07_estimate_att.R complete.")
}
