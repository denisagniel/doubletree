## ============================================================================
## application/helpers/discretize_prior_cost.R
##
## Pure function. No data reads, no side effects, no smidata dependency.
##
## prior_cost (continuous dollars, constructed over [INDEX_DT - 12mo, INDEX_DT)
## in 03_cost_windows.R) -> 3 binary dummies, Q1 as the implicit reference.
## Confirmed PI 2026-09-17 (smidata inst/analyses/doubletree__application.yml,
## covariates.grid_exact_sparsity: "discretized into QUARTILES (Q1-Q4, computed
## within the analytic sample at estimation time -- not fixed dollar cutoffs)").
## ============================================================================

#' Discretize `prior_cost` into quartile dummies for a grid-exact design matrix.
#'
#' @description
#' Assigns each patient to a within-sample quartile of `prior_cost` and returns
#' the quartile as three binary indicators (`prior_cost_q2`, `prior_cost_q3`,
#' `prior_cost_q4`), with Q1 the implicit reference category. The realized
#' cutpoints are attached as an attribute.
#'
#' @details
#' \strong{Why this must happen HERE, before the estimator, rather than being
#' left to the estimator.} Both `doubletree::estimate_att()` and
#' `doubletree::estimate_att_crossfit()` take `discretize_method` (default
#' `"quantiles"`) and `discretize_bins` (default `"adaptive"`) and forward them
#' to `optimaltrees::bisect_lambda_to_budget()` / `optimaltrees::fit_tree()`
#' (verified in `R/estimate_att.R` lines 384-385 and 402-403, and
#' `R/estimate_att_crossfit.R` lines 219-220). A continuous column handed to
#' either estimator un-preprocessed is therefore discretized *inside*
#' `optimaltrees`, adaptively, at cutpoints chosen by the fitting procedure and
#' never declared by the analyst.
#'
#' That destroys grid-exactness by construction. `estimate_att()`'s identifying
#' assumption (`theory.tex` `ass:sparsity`) is that a tree of at most
#' `leaf_budget` leaves represents both nuisances exactly *on the analyst's
#' pre-specified grid*, and `ass:finite` licenses a general covariate space only
#' when it is discretized "via a finite set of analyst-chosen cutpoints fixed
#' before seeing the data". A cutpoint chosen by the fitter, per fit, per fold,
#' is not that grid: the assumption would be stated over an object that does not
#' exist until after fitting, and would differ between the flagship and
#' cross-fit paths.
#'
#' The two estimators fail differently on the same mistake, and neither failure
#' is the one you want:
#' \itemize{
#'   \item `estimate_att()` REJECTS a non-binary column outright
#'     (`R/estimate_att.R`, the `non_binary` check: every column must satisfy
#'     `all(col %in% c(0, 1))`). Loud, and recoverable.
#'   \item `estimate_att_crossfit()` ACCEPTS it silently. Its only input
#'     validation is `check_att_data()`, which checks `NA`s, that `A` is binary,
#'     and `Y` against `outcome_type` -- it does not examine `X`'s coding at all.
#'     A continuous `prior_cost` therefore runs to completion and returns a
#'     number, computed on an undeclared grid.
#' }
#' Pre-discretizing removes the second failure mode, which is the dangerous one.
#'
#' \strong{Why the cutpoints are sample-dependent, and what follows.} The
#' confirmed design specifies quartiles "computed within the analytic sample",
#' not fixed dollar thresholds. So the grid is fixed before seeing the OUTCOME
#' but is a function of the analytic sample's own `prior_cost` distribution. Two
#' consequences, both belonging in the methods write-up rather than in a
#' comment: the realized cutpoints must be REPORTED (hence the `cutpoints`
#' attribute), and a sensitivity analysis at fixed dollar cutpoints is the
#' natural robustness check.
#'
#' \strong{Zero-inflation and the collapsed-boundary abort.} Pre-index cost is
#' heavily zero-inflated: a patient with no pre-index claims has
#' `prior_cost == 0` (a real zero, not a missing value). If more than 25% of the
#' sample is at zero, the 25th and 50th percentiles are both 0, "Q2" is empty,
#' and `prior_cost_q2` is an all-zero column. An all-zero covariate is not a
#' harmless no-op here: it is a declared grid coordinate with one realized
#' level, so the stated grid and the achieved grid differ. This function aborts
#' rather than returning it. The fix is a decision (fewer bins, a
#' zero-versus-positive indicator plus tertiles of the positive part, ...), and
#' decisions are not made silently by a helper.
#'
#' @param prior_cost Numeric vector of pre-index total cost in dollars, one
#'   element per patient. `NA` is not accepted: absence of a pre-index claim is
#'   a real zero (see `04_covariates.R`), so an `NA` here means an upstream join
#'   failed and must not be quietly binned.
#' @param probs Numeric vector of quantile probabilities defining the interior
#'   cutpoints. Default `c(0.25, 0.5, 0.75)` (quartiles), the confirmed design.
#'   Length `k - 1` produces `k - 1` dummies over `k` bins.
#'
#' @return A [tibble::tibble()] with `length(prior_cost)` rows and one binary
#'   integer column per interior cutpoint, named `prior_cost_q2`,
#'   `prior_cost_q3`, ... The realized cutpoints are attached as attribute
#'   `"cutpoints"` (a named numeric vector) and the probabilities used as
#'   `"probs"`.
discretize_prior_cost <- function(prior_cost, probs = c(0.25, 0.5, 0.75)) {
  if (!is.numeric(prior_cost)) {
    cli::cli_abort("{.arg prior_cost} must be numeric, not {.cls {class(prior_cost)}}.")
  }
  if (length(prior_cost) == 0L) {
    cli::cli_abort("{.arg prior_cost} is empty; there is nothing to discretize.")
  }
  if (anyNA(prior_cost)) {
    cli::cli_abort(c(
      "{.arg prior_cost} contains {sum(is.na(prior_cost))} {.val NA} value{?s}.",
      "i" = "A patient with no pre-index claim has {.val 0}, not {.val NA}
             (04_covariates.R). {.val NA} here means an upstream join failed.",
      "x" = "Refusing to bin it: an {.val NA} silently assigned to a quartile is a
             fabricated covariate value."
    ))
  }
  if (any(!is.finite(prior_cost))) {
    cli::cli_abort("{.arg prior_cost} contains non-finite values.")
  }
  if (any(prior_cost < 0)) {
    cli::cli_abort(c(
      "{.arg prior_cost} has {sum(prior_cost < 0)} negative value{?s}.",
      "i" = "Negative paid amounts are claim reversals. Whether to net them out
             or drop them is a decision, not a binning detail."
    ))
  }
  if (!is.numeric(probs) || length(probs) < 1L || anyNA(probs) ||
      any(probs <= 0) || any(probs >= 1)) {
    cli::cli_abort("{.arg probs} must be numeric values strictly inside (0, 1).")
  }
  if (is.unsorted(probs, strictly = TRUE)) {
    cli::cli_abort("{.arg probs} must be strictly increasing.")
  }

  cuts <- stats::quantile(prior_cost, probs = probs, names = FALSE, type = 7)

  ## A collapsed boundary means an empty bin, hence an all-zero dummy. Report
  ## every collapsed pair and the zero share, because zero-inflation is the
  ## overwhelmingly likely cause and the reader needs the number to decide what
  ## to do instead.
  collapsed <- which(diff(cuts) <= 0)
  if (length(collapsed) > 0L) {
    zero_share <- mean(prior_cost == 0)
    collapsed_at <- paste(round(cuts[collapsed], 2), collapse = ", ")
    cli::cli_abort(c(
      "Quartile boundaries collapsed at cutpoint value(s) {.val {collapsed_at}}:
       not strictly less than the next.",
      "i" = "Realized cutpoints at probs {.val {probs}}: {.val {round(cuts, 2)}}.",
      "i" = "{round(100 * zero_share, 1)}% of the sample has
             {.field prior_cost} == 0.",
      "x" = "An empty bin yields an all-zero dummy -- a declared grid coordinate
             with one realized level, so the stated and achieved grids differ.",
      "i" = "Fix by DECISION, not by coercion: fewer bins, or a
             zero-versus-positive indicator plus quantiles of the positive part."
    ))
  }

  ## Bin index: `left.open = TRUE` makes the bins (-Inf, c1], (c1, c2], ...,
  ## so a patient exactly AT the 25th percentile falls in Q1. That is the
  ## conventional quantile-group assignment, and it is what keeps the zero-cost
  ## patients together in Q1 when the 25th percentile is 0.
  bin <- 1L + findInterval(prior_cost, cuts, left.open = TRUE)

  dummy_names <- paste0("prior_cost_q", seq_along(cuts) + 1L)
  out <- tibble::as_tibble(stats::setNames(
    lapply(seq_along(cuts) + 1L, function(k) as.integer(bin == k)),
    dummy_names
  ))

  ## Strictly increasing cutpoints are necessary but NOT sufficient for every
  ## bin to be occupied: with ties at a cutpoint, `left.open = TRUE` can push a
  ## whole bin's worth of mass down into its predecessor. e.g. prior_cost =
  ## c(1, 1, 2, 3) gives cuts (1, 1.5, 2.5) -- strictly increasing -- yet the
  ## second bin is empty because both 1s land in the first. So the occupancy
  ## invariant is checked directly rather than inferred from the cutpoints.
  empty_bins <- setdiff(seq_along(cuts) + 1L, unique(bin))
  if (length(empty_bins) > 0L) {
    dead <- paste(paste0("prior_cost_q", empty_bins), collapse = ", ")
    cli::cli_abort(c(
      "Empty bin(s) {.val {empty_bins}}: dummy/dummies {.field {dead}} would be
       all-zero.",
      "i" = "Realized cutpoints at probs {.val {probs}}: {.val {round(cuts, 2)}};
             cutpoints are strictly increasing, so the cause is TIES AT a
             cutpoint rather than a collapsed boundary.",
      "x" = "An all-zero dummy is a declared grid coordinate with one realized
             level, so the stated and achieved grids differ.",
      "i" = "Same remedy as a collapsed boundary: fewer bins, or a
             zero-versus-positive indicator plus quantiles of the positive part."
    ))
  }

  attr(out, "cutpoints") <- stats::setNames(cuts, paste0("p", probs * 100))
  attr(out, "probs") <- probs
  out
}
