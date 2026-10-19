## ============================================================================
## application/helpers/winsorize_outcome.R
##
## Upper winsorization of the cost outcome, applied ONCE in 06 so every
## downstream estimator (flagship tree, cross-fit anchor, OLS comparison) sees
## the same Y. This changes the estimand to the ATT on winsorized cost (a
## trimmed-cost ATT); the raw cost is kept as `Y_raw`.
## ============================================================================

#' Cap a non-negative outcome at an upper quantile of its own sample.
#'
#' The cap is a fixed function of the analytic sample (it DEFINES the
#' outcome), not a nuisance fitted per fold, so it is computed once on the
#' pooled final sample.
#'
#' @param y Numeric vector, no NA.
#' @param q Upper quantile in (0, 1]. `1` means no winsorization.
#' @return `y` capped at `quantile(y, q)`, with attributes `cap` (the cap, or
#'   `Inf` when `q == 1`), `q`, and `n_capped` (rows strictly above the cap).
winsorize_upper <- function(y, q) {
  if (!is.numeric(y) || anyNA(y)) {
    cli::cli_abort("{.arg y} must be numeric with no NA.")
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0 || q > 1) {
    cli::cli_abort("{.arg q} must be a single number in (0, 1], got {.val {q}}.")
  }
  cap <- if (q >= 1) Inf else unname(stats::quantile(y, probs = q, type = 7))
  out <- pmin(y, cap)
  attr(out, "cap") <- cap
  attr(out, "q") <- q
  attr(out, "n_capped") <- sum(y > cap)
  out
}
