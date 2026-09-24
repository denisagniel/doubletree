## ============================================================================
## application/helpers/msr_filter.R
##
## Pure function. No data reads, no side effects, no smidata dependency (same
## convention as discretize_prior_cost.R / complete_case.R).
##
## Extracted from 02_population_and_eligibility.R's MSR-contamination guard
## (2026-09-24) so the filter-and-report logic is independently testable with
## hand-built fixtures, rather than only exercisable end-to-end against a
## Tier-0 fixture or the real server table.
## ============================================================================

#' Filter an `msr_aap`-shaped frame to a single literal `MSR` code, reporting
#' any exclusion rather than trusting it silently.
#'
#' @description
#' `smidata::smi_select_aap_row()`'s winner-selection groups and joins on `ID`
#' only, and never reads `MSR` except inside its tie-diagnostic message -- its
#' abort-on-tie protects against two rows tying at the exact same winning
#' `month_gap`, NOT against a non-AAP row simply being closer in time to the
#' target month than any true AAP row (no tie forms in that case, so nothing
#' aborts, and the wrong row's `MSR_NUM`/`MSR_DEN` silently flow into
#' `elig_aap`/`aap_achieved`). Filtering to the confirmed AAP code before
#' calling that function is what actually closes this gap; see
#' `OPEN_DECISIONS$msr_code`'s note (`application/_config.R`) for the full
#' correction and design review (`oracle`, 2026-09-24).
#'
#' @details
#' `trimws()` is applied before comparison: SAS character columns are
#' commonly blank-padded on read, and a bare `==` against a padded value fails
#' silently to zero rows, which then looks like a data problem rather than a
#' comparison problem.
#'
#' @param aap_raw Data frame/tibble with an `MSR` column (character, possibly
#'   blank-padded).
#' @param msr_code Single, non-missing character value: the literal code to
#'   keep (e.g. `"AAP"` -- `application/_config.R`'s `AAP_MSR_CODE`, NOT
#'   `OPEN_DECISIONS$msr_code$value`, which is a human-readable annotation of
#'   the resolution, not a filter value).
#' @param quiet If `TRUE`, suppress the `cli::cli_warn()` on exclusion (the
#'   caller is expected to report it in its own words instead -- e.g. as an
#'   informational, non-alarming message under a Tier-0 fixture, where
#'   exclusion is expected and not a finding). Default `FALSE`.
#'
#' @return The filtered data frame/tibble (same class and columns as
#'   `aap_raw`), rows where `trimws(MSR) != msr_code` removed. Attributes
#'   `"n_excluded"` (integer count) and `"excluded_msr_values"` (character
#'   vector of the other codes found, possibly empty) are attached
#'   unconditionally, so a caller can act on the exclusion without
#'   re-deriving it.
filter_to_msr_code <- function(aap_raw, msr_code, quiet = FALSE) {
  if (!"MSR" %in% names(aap_raw)) {
    cli::cli_abort("{.arg aap_raw} is missing column {.field MSR}.")
  }
  if (!is.character(msr_code) || length(msr_code) != 1L || is.na(msr_code)) {
    cli::cli_abort("{.arg msr_code} must be a single, non-missing character value.")
  }

  msr_clean <- trimws(as.character(aap_raw$MSR))
  distinct_msr <- unique(msr_clean)
  n_before <- nrow(aap_raw)
  out <- aap_raw[msr_clean == msr_code, , drop = FALSE]
  n_after <- nrow(out)
  n_excluded <- n_before - n_after
  excluded_values <- setdiff(distinct_msr, msr_code)

  if (!quiet && n_excluded > 0L) {
    cli::cli_warn(c(
      "!" = "{n_excluded} of {n_before} row{?s} had {.field MSR} !=
             {.val {msr_code}} (value{?s}: {.val {excluded_values}}) and
             {?was/were} EXCLUDED.",
      "i" = "If {.field MSR} is not constant, other column semantics assumed
             constant may be wrong too, not only these rows."
    ))
  }

  attr(out, "n_excluded") <- n_excluded
  attr(out, "excluded_msr_values") <- excluded_values
  out
}
