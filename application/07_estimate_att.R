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
## THIS SCRIPT ABORTS IMMEDIATELY TODAY, before loading anything. All four
## blocking decisions bear on it: msr_code (A is nondeterministic),
## cost_family_scope (Y may be double-counted, AND prior_cost's quartile
## boundaries with it), enrollment_source (the analytic sample is empty), and
## diagnostics (nothing specifies what licenses reporting the flagship estimate).
## Running an estimator on top of those would produce a number, and the number
## would be meaningless.
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
## Cost grows with the number of binary columns and the leaf budget; 15 columns at
## leaf_budget = 4 is modest, but the crossfit path fits 2K trees rather than 2 and
## its cv_regularization sweep multiplies that again. Budget accordingly.
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

source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))

library(doubletree)

set.seed(config_seed)

analytic_data <- arrow::read_parquet(
  file.path(app_dir, "output", "analytic_cohort.parquet")
)

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
  outcome_type = "continuous"
)

cli::cli_h2("estimate_att() -- flagship, requires grid-exact sparsity")
cli::cli_alert_info("theta = {signif(fit_att$theta, 4)}  (SE {signif(fit_att$sigma, 4)})")
cli::cli_alert_info("95% CI: [{signif(fit_att$ci_95[1], 4)}, {signif(fit_att$ci_95[2], 4)}]")

## Sparsity PROXIES, reported and never acted on automatically. Which of these
## licenses reporting this estimate is OPEN_DECISIONS$diagnostics -- one of the
## reasons the gate above blocks. Note what they can and cannot say: certified_*
## means the returned fit provably solves eq:select, i.e. it is the right tree
## for the budget; it does NOT mean the budget is large enough for the truth.
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
## Also worth stating in the write-up rather than discovering in review: this
## path's `sigma` is the in-sample EIF plug-in and is downward biased by a
## relative factor of order (2 * leaf_budget + 1) / n, so the nominal 95%
## interval can undercover at small n. estimate_att()'s own documentation says
## so; it is not implemented away because the fix would change the shared
## att_se() used by all three entry points.

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

flagship <- if (identical(config_estimator, "att")) fit_att else fit_cf
cli::cli_alert_info(
  "Flagship per {.code config_estimator = \"{config_estimator}\"}: theta =
   {signif(flagship$theta, 4)}"
)

## ---- Outputs for the manuscript's Empirical Application section -------------
##   application/output/att-estimate-table.tex     both estimators, side by side
##   application/output/sparsity-diagnostics.tex   certified_*, n_leaves_*, gap_*
##   application/output/prior-cost-cutpoints.tex   the realized quartile grid
##   application/output/nuisance-trees.pdf         the two fitted partitions --
##       the interpretability payoff, and the reason a 15-column grid was chosen
##       over dual-bounds' ~68.

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("07_estimate_att.R complete.")
}
