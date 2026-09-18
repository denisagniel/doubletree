## ============================================================================
## application/helpers/design_matrix.R
##
## Assemble and validate X: the 15-column, all-binary design matrix
## doubletree::estimate_att() requires.
##
## Pure functions. No data reads, no smidata dependency.
##
## WHY THIS FILE EXISTS AT ALL -- i.e. what the package does NOT check for you.
## Read `R/utils.R`'s check_att_data() before editing: it validates only that
## X/A/Y are non-NULL, that lengths agree, that there are no NAs anywhere in X,
## A or Y, that A is 0/1, and that Y matches `outcome_type`. It does NOT examine
## X's CODING. The two entry points then diverge:
##
##   estimate_att()          DOES reject non-binary X, itself, after
##                           check_att_data() -- see R/estimate_att.R's
##                           `non_binary` check, `all(col %in% c(0, 1))`.
##   estimate_att_crossfit()  does NOT. check_att_data() is its only input
##                           validation on X, so a 1/2-coded or continuous
##                           column runs to completion and returns a number,
##                           discretized internally by optimaltrees at
##                           cutpoints the analyst never declared.
##
## Neither entry point checks, on either path: the number of columns, their
## ORDER, or whether any column is constant. So of the four properties this
## file asserts, exactly one is redundant with estimate_att() (binary coding),
## it is redundant on ONE of the two paths, and the redundancy is checked HERE
## at assembly time in 06 -- before an artifact is written, rather than after
## a multi-hour cost aggregation has already been spent.
## ============================================================================

#' Assemble the 15-column binary design matrix.
#'
#' Joins the 12 confirmed already-binary `_YN` columns from
#' `larger_smi_covariates` to the 3 constructed `prior_cost` quartile dummies,
#' in the fixed order declared by `design_matrix_columns()`.
#'
#' Nothing is coerced. A `_YN` column arriving as `"Y"`/`"N"`, or as 1/2, is a
#' finding about the source table (and a Tier-1 check -- see
#' `90_checks_tier1.R`), not something to fix silently inside an assembler: a
#' `"Y"` -> 1 mapping that guesses which level means "yes" can invert a
#' covariate's sign without ever erroring.
#'
#' @param covariates A data frame with one row per patient, containing `ID` and
#'   the 12 columns named by `covariate_blocks()`. Extra columns are ignored.
#' @param prior_cost_dummies A data frame with the 3 columns named by
#'   `prior_cost_dummy_names()`, with rows ALIGNED to `covariates` (same order,
#'   same length) -- the output of `discretize_prior_cost()` applied to
#'   `covariates`' own patients. Aligned by position, not joined by key, because
#'   `discretize_prior_cost()` is a vectorized transform of one column and
#'   carries no key of its own; the caller is responsible for having built it
#'   from this same table (see `04_covariates.R`).
#'
#' @return A [tibble::tibble()] with `ID` followed by the 15 design columns in
#'   `design_matrix_columns()` order. Any `"cutpoints"` attribute on
#'   `prior_cost_dummies` is carried through, so the realized quartile
#'   boundaries stay attached to the matrix they define.
assemble_design_matrix <- function(covariates, prior_cost_dummies) {
  if (!"ID" %in% names(covariates)) {
    cli::cli_abort("{.arg covariates} must have an {.field ID} column.")
  }
  yn_cols <- unlist(covariate_blocks(), use.names = FALSE)
  missing_yn <- setdiff(yn_cols, names(covariates))
  if (length(missing_yn) > 0L) {
    cli::cli_abort(c(
      "{.arg covariates} is missing confirmed covariate(s)
       {.field {paste(missing_yn, collapse = ', ')}}.",
      "i" = "The 15-column selection is PI-confirmed (2026-09-17); a missing
             column is a contract failure, not something to substitute for."
    ))
  }

  dummy_cols <- prior_cost_dummy_names()
  missing_dummies <- setdiff(dummy_cols, names(prior_cost_dummies))
  if (length(missing_dummies) > 0L) {
    cli::cli_abort(
      "{.arg prior_cost_dummies} is missing
       {.field {paste(missing_dummies, collapse = ', ')}}; expected the output of
       {.fn discretize_prior_cost}."
    )
  }

  if (nrow(covariates) != nrow(prior_cost_dummies)) {
    cli::cli_abort(c(
      "Row counts differ: {.arg covariates} has {nrow(covariates)},
       {.arg prior_cost_dummies} has {nrow(prior_cost_dummies)}.",
      "i" = "These are aligned BY POSITION. Build the dummies from this same
             table's {.field prior_cost} column."
    ))
  }
  if (anyDuplicated(covariates$ID) > 0L) {
    cli::cli_abort("{.arg covariates} must be one row per patient.")
  }

  x <- dplyr::bind_cols(
    dplyr::select(covariates, "ID", dplyr::all_of(yn_cols)),
    dplyr::select(prior_cost_dummies, dplyr::all_of(dummy_cols))
  )
  out <- dplyr::select(x, "ID", dplyr::all_of(design_matrix_columns()))

  attr(out, "cutpoints") <- attr(prior_cost_dummies, "cutpoints")
  out
}

#' Assert that a design matrix meets every precondition of the ATT estimators.
#'
#' Four properties, none of which both entry points check (see this file's
#' header for exactly which does what):
#'
#' 1. **Exactly the 15 declared columns, in the declared order.** Order matters
#'    because the fitted trees' split indices are read against it.
#' 2. **Every column coded exactly `{0, 1}`.** Not "has two distinct values" --
#'    a column valued 1/2, or `"Y"`/`"N"`, has two distinct values and is still
#'    wrong. It passes `check_att_data()` unexamined, is accepted by
#'    `estimate_att_crossfit()`, and changes the grid the sparsity assumption is
#'    stated over.
#' 3. **No `NA` anywhere.** `check_att_data()` does catch this, at call time.
#'    Catching it here means catching it before 06 writes anything.
#' 4. **No constant column.** A one-level covariate is a declared grid
#'    coordinate that does not exist in the sample, so the stated grid and the
#'    achieved grid differ -- and neither estimator complains.
#'
#' @param X A data frame or tibble of covariates. An `ID` column, if present, is
#'   ignored: it is a key, not a covariate, and is not passed to the estimator.
#' @param expected Character vector of expected column names in expected order.
#'   Default `design_matrix_columns()`.
#' @return `invisible(TRUE)`, or aborts.
assert_binary_design_matrix <- function(X, expected = design_matrix_columns()) {
  if (!is.data.frame(X)) {
    cli::cli_abort("{.arg X} must be a data frame, not {.cls {class(X)}}.")
  }
  x <- X[, setdiff(names(X), "ID"), drop = FALSE]

  ## (1) exact column set, in exact order
  if (!identical(names(x), expected)) {
    extra <- setdiff(names(x), expected)
    missing <- setdiff(expected, names(x))
    cli::cli_abort(c(
      "{.arg X} does not match the declared {length(expected)}-column design.",
      "x" = "Got {ncol(x)} column(s): {.field {paste(names(x), collapse = ', ')}}",
      if (length(missing) > 0L) {
        c("x" = "Missing: {.field {paste(missing, collapse = ', ')}}")
      },
      if (length(extra) > 0L) {
        c("x" = "Unexpected: {.field {paste(extra, collapse = ', ')}}")
      },
      if (length(missing) == 0L && length(extra) == 0L) {
        c("x" = "Same columns, WRONG ORDER. Expected:
                 {.field {paste(expected, collapse = ', ')}}")
      },
      "i" = "Order is part of the specification: fitted split indices are read
             against it."
    ))
  }

  ## (3) NAs, checked before coding so the message names the real problem
  na_cols <- names(x)[vapply(x, anyNA, logical(1))]
  if (length(na_cols) > 0L) {
    cli::cli_abort(c(
      "{length(na_cols)} design column(s) contain {.val NA}:
       {.field {paste(na_cols, collapse = ', ')}}.",
      "i" = "Every one of the 12 {.field _YN} columns is a flag whose absence
             means something specific; decide what, do not impute here."
    ))
  }

  ## (2) coded exactly {0, 1}
  not_binary <- names(x)[!vapply(x, function(col) {
    (is.numeric(col) || is.logical(col)) && all(col %in% c(0, 1))
  }, logical(1))]
  if (length(not_binary) > 0L) {
    observed <- vapply(not_binary, function(nm) {
      vals <- sort(unique(x[[nm]]))
      paste0(nm, "={", paste(utils::head(vals, 4L), collapse = ","),
             if (length(vals) > 4L) ",..." else "", "}")
    }, character(1))
    cli::cli_abort(c(
      "{length(not_binary)} design column(s) not coded {.val 0}/{.val 1}:
       {.field {paste(observed, collapse = '; ')}}.",
      "x" = "{.fn check_att_data} does not examine X's coding, and
             {.fn estimate_att_crossfit} performs no further check -- a 1/2-coded
             column would run to completion and return a number.",
      "i" = "{.fn estimate_att} would reject it, but only after the analytic data
             was already assembled and written."
    ))
  }

  ## (4) constant columns
  constant <- names(x)[vapply(x, function(col) length(unique(col)) < 2L, logical(1))]
  if (length(constant) > 0L) {
    cli::cli_abort(c(
      "{length(constant)} design column(s) constant:
       {.field {paste(constant, collapse = ', ')}}.",
      "x" = "A one-level covariate is a declared grid coordinate absent from the
             sample, so the grid the sparsity assumption is stated over is not
             the grid that was fit.",
      "i" = "Drop it and restate the covariate set, or find out why it is
             constant. Neither estimator will complain."
    ))
  }

  invisible(TRUE)
}
