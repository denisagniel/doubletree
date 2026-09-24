## ============================================================================
## application/helpers/complete_case.R
##
## Pure functions. No data reads, no side effects, no smidata dependency (see
## `discretize_prior_cost.R`'s header for why that matters: this file is
## Tier-0-testable on hand-built fixtures alone, with no server and no
## `smidata` install).
##
## Implements OPEN_DECISIONS$enrollment_source's PI resolution
## (`application/_config.R`, confirmed 2026-09-24): a patient is complete-case
## iff `MEDICAID_FLAG == 1` in `larger_smi_medicaid_monthly_flag` for EVERY one
## of `config_complete_case_months` (24) consecutive calendar months starting
## at `INDEX_DT`'s own month, with NO filtering on `TYPE` and ZERO gap
## tolerance.
##
## `05_complete_case.R` is the orchestration script: it reads the real table
## via `smidata::smi_require()`/`smi_read()` and calls `compute_complete_case()`
## below. That split mirrors `discretize_prior_cost()` / `06_assemble_analytic_
## data.R` exactly, for the same reason: a function that reads no data can be
## tested with hand-built fixtures, at Tier 0, with the same rigor a server run
## would get -- and the correctness of the WINDOWING/MEMBERSHIP logic below is
## exactly the part that must not be trusted to "it ran on the server and the
## numbers looked plausible."
## ============================================================================

#' Absolute-month index for a `Date`, floored to the first of its own month.
#'
#' Reproduces `smidata::smi_month_index(yyyy, mm)`'s definition
#' (`year*12 + month`) rather than depending on the package, per this file's
#' header. Doing all window arithmetic in this integer space (instead of on
#' `Date`s with `lubridate::add_with_rollback()`) removes any dependence on
#' date-comparison boundary behaviour -- there is no rollback question when
#' adding an integer to an integer. Flooring to month-start is one line of
#' base R (`format()`'s `"%Y"`/`"%m"` already discard the day), so this file
#' adds no dependency beyond what `discretize_prior_cost.R` already uses.
#'
#' @param date A `Date` vector.
#' @return An integer vector, one element per input.
ym_month_index <- function(date) {
  d <- as.Date(date)
  as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
}

#' Parse `YEAR_MONTH` into a first-of-month `Date`, without guessing the format.
#'
#' @description
#' `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` format is not fixed by any
#' contract this pipeline has seen. smidata's own evidence-gathering script
#' (`~/RAND/tools/smidata/inst/server/07_enrollment_coverage.R`, Section 3)
#' discovers it at runtime rather than assuming it, because a wrong guess does
#' not error -- it silently mis-parses. This function does the same detection,
#' but checks EVERY value, not a sample: that script characterizes a
#' distribution, where a handful of mis-parsed outlier rows would be diluted
#' noise; here, every row can decide whether a specific patient is included in
#' the analytic sample, so a format check on a sample that misses an outlier
#' row is exactly the failure that would flip an inclusion decision.
#'
#' @param year_month A vector, one element per row. Accepted forms: already
#'   `Date`/`POSIXct` (floored to month-start), or character in `"YYYY-MM"` or
#'   `"YYYYMM"` form (whichever form is used, ALL values must match it, or the
#'   function aborts rather than mixing formats row-by-row).
#' @return A `Date` vector, one element per input element, the first day of
#'   the parsed month.
parse_year_month <- function(year_month) {
  if (length(year_month) == 0L) {
    return(as.Date(character()))
  }
  if (inherits(year_month, "Date") || inherits(year_month, "POSIXct")) {
    d <- as.Date(year_month)
    return(as.Date(format(d, "%Y-%m-01")))
  }
  if (!is.character(year_month)) {
    cli::cli_abort(
      "{.arg year_month} must be character, {.cls Date}, or {.cls POSIXct},
       not {.cls {class(year_month)}}."
    )
  }
  if (anyNA(year_month)) {
    cli::cli_abort(
      "{.arg year_month} has {sum(is.na(year_month))} {.val NA} value{?s}.
       {.field YEAR_MONTH} is a key column; a missing value here means an
       upstream read failed and must not be silently dropped or coerced."
    )
  }

  fmt_yyyymm      <- "^\\d{4}-\\d{2}$"   # e.g. "2019-06"
  fmt_yyyymm_flat <- "^\\d{6}$"          # e.g. "201906"

  is_yyyymm      <- grepl(fmt_yyyymm, year_month)
  is_yyyymm_flat <- grepl(fmt_yyyymm_flat, year_month)

  if (all(is_yyyymm)) {
    yyyy <- as.integer(substr(year_month, 1L, 4L))
    mm   <- as.integer(substr(year_month, 6L, 7L))
  } else if (all(is_yyyymm_flat)) {
    yyyy <- as.integer(substr(year_month, 1L, 4L))
    mm   <- as.integer(substr(year_month, 5L, 6L))
  } else {
    offenders <- utils::head(unique(year_month[!is_yyyymm & !is_yyyymm_flat]), 10L)
    cli::cli_abort(c(
      "{.arg year_month} does not match either anticipated format
       ({.val YYYY-MM} or {.val YYYYMM}) for every value, or matches both
       inconsistently across rows.",
      "x" = "Offending value{?s}: {.val {offenders}}",
      "i" = "Do not coerce this to {.val NA} and proceed -- the real format
             may differ from what
             {.path ~/RAND/tools/smidata/inst/server/07_enrollment_coverage.R}
             last found. Update this parser once the real format is known,
             the same way {.fn smidata::smi_parse_msr_yr} aborts on an
             MSR_YR mismatch rather than guessing."
    ))
  }

  if (any(mm < 1L | mm > 12L)) {
    bad <- unique(year_month[mm < 1L | mm > 12L])
    cli::cli_abort(
      "Parsed {.arg year_month} month component out of range 1-12 for
       value{?s}: {.val {utils::head(bad, 10L)}}."
    )
  }
  ## Plausible year band -- catches a transposed/garbled field a bare format
  ## regex would still pass (e.g. a 4-digit year that is not actually a year).
  if (any(yyyy < 1900L | yyyy > 2100L)) {
    bad <- unique(year_month[yyyy < 1900L | yyyy > 2100L])
    cli::cli_abort(
      "Parsed {.arg year_month} year component outside a plausible range
       (1900-2100) for value{?s}: {.val {utils::head(bad, 10L)}}."
    )
  }

  as.Date(sprintf("%04d-%02d-01", yyyy, mm))
}

#' Per-patient complete-case indicator from monthly Medicaid-enrollment flags.
#'
#' @description
#' A patient is complete-case iff `MEDICAID_FLAG == 1` for EVERY one of
#' `window_months` consecutive calendar months starting at `INDEX_DT`'s own
#' month, i.e. the month set
#' `{floor_month(INDEX_DT), floor_month(INDEX_DT) + 1, ..., + window_months - 1}`.
#' `TYPE` is NOT filtered on (PI, 2026-09-24: "note this as a question and
#' just proceed ignoring this for now"): a month counts as covered if ANY row
#' for that `(ID, calendar month)`, of ANY `TYPE`, has `MEDICAID_FLAG == 1` --
#' the most PERMISSIVE reading of "ignore TYPE". Gap tolerance is zero (PI:
#' "let's just not have a gap for now"): every required month is checked by
#' SET MEMBERSHIP, not by comparing a covered-month COUNT to `window_months` --
#' a count alone cannot distinguish genuine full coverage from
#' `window_months` months of coverage concentrated in the wrong calendar
#' months (see Details for why that distinction is real, not hypothetical).
#'
#' @details
#' **Why this is not the claims-presence proxy this pipeline's
#' `05_complete_case.R` header forbids.** That header's two failure modes
#' (false negatives from zero-utilization complete patients; false positives
#' from partially-enrolled patients with a claim inside a short window) are
#' both about inferring enrollment from claims PRESENCE -- a service being
#' billed. This function reads `MEDICAID_FLAG`, an ENROLLMENT indicator
#' recorded independently of claims activity: a patient enrolled the whole
#' window with zero utilization still has `MEDICAID_FLAG == 1` every month.
#' Neither directional failure mode in that header applies to this input.
#'
#' **Why set membership, not a count.** Consider computing, per patient, the
#' number of distinct covered calendar months that fall ANYWHERE inside
#' `[INDEX_DT, INDEX_DT + window_months months)` when `INDEX_DT` is not first
#' floored to its own month start (e.g. `INDEX_DT = 2019-06-15`): the
#' half-open interval then spans `2019-06-15` to `2021-06-15`, and the
#' calendar months with a month-start date strictly inside it are
#' `2019-07-01 .. 2021-06-01` -- 24 of them, so a naive count check
#' (`== 24`) PASSES, yet the index month itself (`2019-06`) was silently
#' excluded and a patient uncovered exactly at index can be certified
#' complete. This function floors `INDEX_DT` to its own month first (via
#' [ym_month_index()]) and only ever compares INTEGER month indices built
#' from that floored anchor, so this shift cannot occur; membership is
#' checked explicitly (via an anti-join against the exact required set)
#' rather than inferred from a count, so the safety property is visible in
#' the code rather than something a future reader has to re-derive.
#'
#' **Missing-patient semantics.** `monthly_flag` and `cohort` are independent
#' inputs. A cohort patient with zero rows in `monthly_flag` -- at ALL, not
#' just outside their window -- is not silently dropped, and is not left as
#' an implicit `NA` for a caller's `left_join()` to coerce via `isTRUE()`:
#' every cohort patient gets exactly one output row, with `complete_case =
#' FALSE` and `reason == "absent_from_flag_table"`. \strong{This treats
#' absence as "not enrolled," which assumes `larger_smi_medicaid_monthly_flag`
#' enumerates Medicaid-months for the whole relevant population -- an
#' assumption this function does NOT verify.} Compare the `n_absent` tally
#' against the full cohort size and reconcile it before trusting a real
#' number; if it is large, this may mean the table is a partial extract, and
#' absence should instead read `NA`, not `FALSE`.
#'
#' **The TYPE self-check.** The PI's "ignore TYPE" default is provisional on
#' a prerequisite (`application/90_checks_tier1.R` CHECK 4, whole-table,
#' server-only, NOT YET RUN as of this writing): `TYPE` must never cause
#' `MEDICAID_FLAG` to disagree for the same `(ID, calendar month)`, or
#' "ignore TYPE" silently picks whichever row happens to be read first/last
#' rather than being a no-op. CHECK 4 verifies this at WHOLE-TABLE scope;
#' this function cannot substitute for that. But every row this function
#' reads is, by construction, restricted to cohort patients and their
#' in-window months -- exactly the rows that actually drive THIS paper's
#' inclusion decisions -- so it computes the identical disagreement check AT
#' THAT RESTRICTED SCOPE for free, and \strong{aborts if it ever finds
#' disagreement} rather than silently trusting an unverified default. If it
#' aborts, "ignore TYPE" is not neutral for this cohort and needs a NEW PI
#' decision (which `TYPE` wins, or require agreement, or ANY-vs-ALL), not a
#' default silently applied over a real disagreement.
#'
#' **`INDEX_DT` is the `cohort` argument's, not this table's own copy.**
#' `monthly_flag`'s `INDEX_DT` column is read only to CROSS-CHECK against
#' `cohort`; the window itself is always anchored on `cohort$INDEX_DT`, to
#' match `06_assemble_analytic_data.R`'s own use of `analysis_cohort$INDEX_DT`
#' (see 06's Step 1). Two things abort rather than being resolved by
#' `min()`-and-warn (smidata's evidence script's own leniency, appropriate for
#' a descriptive distribution, not for a hard inclusion filter): (1) more
#' than one distinct `INDEX_DT` for the same `ID` within `monthly_flag`
#' itself -- may mean this table pools multiple cohorts with different index
#' conventions, in which case `ID` may not even be unique across them, a risk
#' `min()` would silently absorb; (2) `monthly_flag`'s `INDEX_DT` disagreeing
#' with `cohort`'s for the same `ID` -- either an `ID` collision or a stale
#' copy, neither resolved correctly by picking one silently.
#'
#' @param monthly_flag Tibble/data.frame from `larger_smi_medicaid_monthly_flag`,
#'   approximately one row per `(ID, YEAR_MONTH, TYPE)`. Required columns:
#'   `ID`, `INDEX_DT`, `TYPE`, `YEAR_MONTH`, `MEDICAID_FLAG`.
#' @param cohort Tibble/data.frame, one row per patient: `ID`, `INDEX_DT`. The
#'   window anchor.
#' @param window_months Single positive integer: consecutive months required,
#'   starting at `INDEX_DT`'s own calendar month. Callers pass
#'   `config_complete_case_months` (24, confirmed).
#'
#' @return A [tibble::tibble()] with exactly `nrow(cohort)` rows: `ID`,
#'   `complete_case` (logical), `reason` (character: `"covered"`,
#'   `"insufficient_coverage"`, or `"absent_from_flag_table"`). Attribute
#'   `"tally"`: a named list, `n_cohort`, `n_complete`,
#'   `n_insufficient_coverage`, `n_absent`, `n_type_disagreement_pairs`
#'   (always `0L` if the function returns at all -- see the TYPE self-check
#'   above; a non-zero count aborts instead of appearing here).
compute_complete_case <- function(monthly_flag, cohort, window_months) {
  ## ---- input validation ----------------------------------------------------
  req_mf <- c("ID", "INDEX_DT", "TYPE", "YEAR_MONTH", "MEDICAID_FLAG")
  missing_mf <- setdiff(req_mf, names(monthly_flag))
  if (length(missing_mf) > 0L) {
    cli::cli_abort(
      "{.arg monthly_flag} is missing column{?s}: {.field {missing_mf}}."
    )
  }
  req_c <- c("ID", "INDEX_DT")
  missing_c <- setdiff(req_c, names(cohort))
  if (length(missing_c) > 0L) {
    cli::cli_abort("{.arg cohort} is missing column{?s}: {.field {missing_c}}.")
  }
  if (!is.numeric(window_months) || length(window_months) != 1L ||
      is.na(window_months) || window_months != as.integer(window_months) ||
      window_months < 1L) {
    cli::cli_abort("{.arg window_months} must be a single positive integer.")
  }
  window_months <- as.integer(window_months)

  cohort <- tibble::as_tibble(cohort)
  monthly_flag <- tibble::as_tibble(monthly_flag)

  if (nrow(cohort) == 0L) {
    cli::cli_abort(
      "{.arg cohort} is empty; nothing to compute complete-case status for."
    )
  }
  if (anyDuplicated(cohort$ID) > 0L) {
    dupes <- unique(cohort$ID[duplicated(cohort$ID)])
    cli::cli_abort(c(
      "{.arg cohort} has duplicate {.field ID}s: expected one row per patient.",
      "x" = "Duplicated: {.val {utils::head(dupes, 10L)}}"
    ))
  }
  if (anyNA(cohort$INDEX_DT)) {
    cli::cli_abort(
      "{.arg cohort}${.field INDEX_DT} has {sum(is.na(cohort$INDEX_DT))}
       {.val NA} value{?s}; {.field INDEX_DT} is the window anchor and cannot
       be missing."
    )
  }

  ## A join-key class mismatch is the single most likely SILENT catastrophe
  ## here: a character-vs-integer or zero-padding mismatch between the two
  ## tables' ID columns would make every join fail to match, certifying the
  ## entire cohort "absent_from_flag_table" -- indistinguishable, on its
  ## face, from a genuinely unenrolled cohort, without this check.
  if (!identical(class(cohort$ID), class(monthly_flag$ID))) {
    cli::cli_abort(c(
      "{.arg cohort}${.field ID} is class {.cls {class(cohort$ID)}} but
       {.arg monthly_flag}${.field ID} is class
       {.cls {class(monthly_flag$ID)}}.",
      "x" = "A join-key class mismatch can silently certify every cohort
             patient as absent from the flag table rather than erroring --
             fix the mismatch; do not proceed past it."
    ))
  }

  ## ---- MEDICAID_FLAG validity ------------------------------------------------
  flag <- monthly_flag$MEDICAID_FLAG
  flag_num <- suppressWarnings(as.integer(flag))
  bad_flag <- !is.na(flag) & !(flag_num %in% c(0L, 1L))
  if (any(bad_flag)) {
    cli::cli_abort(c(
      "{.arg monthly_flag}${.field MEDICAID_FLAG} has {sum(bad_flag)}
       value(s) outside {{0, 1}} that are not {.val NA}.",
      "x" = "Offending value{?s}: {.val {utils::head(unique(flag[bad_flag]), 10L)}}"
    ))
  }
  n_na_flag <- sum(is.na(flag))
  if (n_na_flag > 0L) {
    cli::cli_warn(
      "{n_na_flag} {.field MEDICAID_FLAG} value{?s} {?is/are} {.val NA};
       treated as NOT covering that row's month -- an {.val NA} flag is not
       evidence of enrollment, so it is excluded rather than counted."
    )
  }
  monthly_flag$MEDICAID_FLAG <- flag_num

  ## ---- INDEX_DT: cohort is the anchor; cross-check the table's own copy ----
  mf_index <- monthly_flag |>
    dplyr::distinct(.data$ID, .data$INDEX_DT) |>
    dplyr::mutate(INDEX_DT = as.Date(.data$INDEX_DT))
  dupe_index <- mf_index |>
    dplyr::count(.data$ID, name = "n_index_dt") |>
    dplyr::filter(.data$n_index_dt > 1L)
  if (nrow(dupe_index) > 0L) {
    cli::cli_abort(c(
      "{.arg monthly_flag} has more than one distinct {.field INDEX_DT} for
       {nrow(dupe_index)} patient(s).",
      "i" = "smidata's inst/server/07_enrollment_coverage.R warns and takes
             {.fn min} here; that is too lenient for a hard inclusion
             filter, and may mean this table pools multiple cohorts with
             different index conventions -- in which case {.field ID} may
             not even be unique across them. See this function's Details.",
      "x" = "Affected ID{?s}: {.val {utils::head(dupe_index$ID, 10L)}}"
    ))
  }

  cohort_anchor <- cohort |>
    dplyr::transmute(ID = .data$ID, cohort_index_dt = as.Date(.data$INDEX_DT))
  mismatch <- mf_index |>
    dplyr::inner_join(cohort_anchor, by = "ID") |>
    dplyr::filter(.data$INDEX_DT != .data$cohort_index_dt)
  if (nrow(mismatch) > 0L) {
    cli::cli_abort(c(
      "{nrow(mismatch)} patient(s) have a different {.field INDEX_DT} in
       {.arg monthly_flag} than in {.arg cohort}.",
      "i" = "{.arg cohort}'s {.field INDEX_DT} is the window anchor (matches
             06_assemble_analytic_data.R's use of {.field analysis_cohort});
             a disagreement means either an {.field ID} collision or a
             stale copy in one of the two tables. Neither is resolved
             correctly by picking one silently.",
      "x" = "Affected ID{?s}: {.val {utils::head(mismatch$ID, 10L)}}"
    ))
  }

  ## ---- window arithmetic, entirely in integer month-index space -----------
  monthly_flag$ym_idx <- ym_month_index(parse_year_month(monthly_flag$YEAR_MONTH))

  cohort_idx <- cohort |>
    dplyr::transmute(
      ID = .data$ID,
      anchor_idx = ym_month_index(as.Date(.data$INDEX_DT))
    )

  restricted <- monthly_flag |>
    dplyr::inner_join(cohort_idx, by = "ID") |>
    dplyr::filter(
      .data$ym_idx >= .data$anchor_idx,
      .data$ym_idx <  .data$anchor_idx + window_months
    )

  ## ---- TYPE self-check, at cohort/window-restricted scope ------------------
  disagreement <- restricted |>
    dplyr::filter(!is.na(.data$MEDICAID_FLAG)) |>
    dplyr::distinct(.data$ID, .data$ym_idx, .data$MEDICAID_FLAG) |>
    dplyr::count(.data$ID, .data$ym_idx, name = "n_flag_values") |>
    dplyr::filter(.data$n_flag_values > 1L)
  n_type_disagreement_pairs <- nrow(disagreement)
  if (n_type_disagreement_pairs > 0L) {
    cli::cli_abort(c(
      "{n_type_disagreement_pairs} (ID, month) pair(s), within the cohort's
       own coverage window, have {.field TYPE} rows that DISAGREE on
       {.field MEDICAID_FLAG}.",
      "x" = "OPEN_DECISIONS$enrollment_source's 'ignore TYPE' default (PI,
             2026-09-24) is a no-op only if TYPE never disagrees on the
             flag. It does here, for at least this many (ID, month)
             pair(s) -- see application/90_checks_tier1.R CHECK 4.",
      "i" = "This needs a NEW PI decision (which TYPE wins, require
             agreement, or ANY-vs-ALL) -- not a default silently applied
             over a real disagreement.",
      "i" = "Affected ID{?s}: {.val {utils::head(unique(disagreement$ID), 10L)}}"
    ))
  }

  ## ---- coverage: ANY TYPE with MEDICAID_FLAG == 1 covers the month --------
  covered_idx <- restricted |>
    dplyr::filter(.data$MEDICAID_FLAG == 1L) |>
    dplyr::distinct(.data$ID, .data$ym_idx)

  ## Long-format required (ID, req_idx) pairs -- exactly window_months rows
  ## per cohort patient, built explicitly rather than via a count-based proxy
  ## (see this function's Details on why membership, not a count, is checked).
  required_long <- cohort_idx[rep(seq_len(nrow(cohort_idx)), each = window_months), ,
                               drop = FALSE]
  required_long$req_idx <- required_long$anchor_idx +
    rep(0:(window_months - 1L), times = nrow(cohort_idx))

  uncovered <- required_long |>
    dplyr::anti_join(
      dplyr::rename(covered_idx, req_idx = "ym_idx"),
      by = c("ID", "req_idx")
    ) |>
    dplyr::count(.data$ID, name = "n_uncovered")

  present_ids <- unique(monthly_flag$ID)

  result <- cohort |>
    dplyr::transmute(ID = .data$ID) |>
    dplyr::left_join(uncovered, by = "ID") |>
    dplyr::mutate(
      is_present   = .data$ID %in% present_ids,
      ## A patient with EVERY required month covered has no row in
      ## `uncovered` at all (dplyr::count() has nothing to count for
      ## them), which becomes NA after this left_join -- that NA means
      ## "the anti-join found nothing missing" (0 uncovered), NOT "fully
      ## uncovered". `is_present` (checked separately, from the raw
      ## table, before any windowing) is what correctly flags a patient
      ## who is absent from monthly_flag entirely.
      n_uncovered  = dplyr::if_else(is.na(.data$n_uncovered), 0L, .data$n_uncovered),
      complete_case = .data$is_present & .data$n_uncovered == 0L,
      reason = dplyr::case_when(
        !.data$is_present     ~ "absent_from_flag_table",
        .data$complete_case   ~ "covered",
        TRUE                  ~ "insufficient_coverage"
      )
    ) |>
    dplyr::select("ID", "complete_case", "reason")

  if (all(result$reason == "absent_from_flag_table")) {
    cli::cli_abort(c(
      "Every cohort patient came back {.val absent_from_flag_table}.",
      "i" = "This is far more consistent with an {.field ID} join-key
             mismatch (format, zero-padding, or a mismatched extract) than
             with a genuinely unenrolled cohort -- the CLASS check above
             passed, but matching classes does not guarantee matching
             VALUES.",
      "x" = "Refusing to return a result that certifies 100% of the cohort
             incomplete without a human looking at this first."
    ))
  }

  attr(result, "tally") <- list(
    n_cohort = nrow(cohort),
    n_complete = sum(result$reason == "covered"),
    n_insufficient_coverage = sum(result$reason == "insufficient_coverage"),
    n_absent = sum(result$reason == "absent_from_flag_table"),
    n_type_disagreement_pairs = n_type_disagreement_pairs
  )
  result
}
