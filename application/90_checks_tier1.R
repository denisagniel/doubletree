## ============================================================================
## 90_checks_tier1.R
##
## Tier-1 checks. Every check requires the real files on the secure server, so
## the comment blocks specify each check and the EXECUTABLE section at the
## bottom implements the ones that are this paper's own (CHECKS 3, 5 and 6).
## Checks 1, 2 (partly) and 4 reference dual-bounds' specification.
##
## Run on the server:   Rscript application/90_checks_tier1.R
## It aborts on a local machine (Tier 0 is run_tests.R). Tables go out through
## smidata::smi_export(), which suppresses declared count columns with n < 11.
## That protects each exported cell, NOT differences between cells or tables
## (smidata documents this as out of scope), so the exports here are built from
## DISJOINT bands and declared complements. Console output on the server is
## unsuppressed and stays on the server.
##
## "Tier 1" is smidata's middle tier: run on the server against a real
## subsample, PLUS the whole-data checks a subsample cannot answer (row counts,
## key overlap, totals).
##
## TWO OF THE CHECKS BELOW ARE dual-bounds', REFERENCED RATHER THAN COPIED.
## This paper inherits dual-bounds' outcome and exposure design verbatim, so it
## inherits the same two risks, on the same files, with the same resolving
## procedure. Re-deriving them here would produce two specifications of one check
## that can drift apart and be answered differently -- which is exactly the
## failure mode the shared OPEN_DECISIONS entries (with their `shared_with` and
## `registry_ref` fields) exist to prevent. Run them ONCE, on the server, and
## record the answer in smidata's registry where both papers read it.
## ============================================================================


## ---------------------------------------------------------------------------
## CHECK 1 -- corroborate OPEN_DECISIONS$cost_family_scope [SHARED, REFERENCED]
##
## STATUS: RESOLVED (PI, 2026-09-17: the families are disjoint patient
## subsets, summed; since 2026-10-15 only the two SERVICE families). This check is CORROBORATION, not a gate: nothing in 03-07
## waits on it, but a disagreement would reopen the decision.
##
## SPECIFICATION: ~/RAND/rprojects/smi/dual-bounds/analysis/real_data/
##                90_checks_tier1.R, CHECK 1.
## REGISTRY:      smidata::inst/analyses/dual-bounds__application.yml#cost_family_scope
##
## Do not re-derive it. That file specifies the whole procedure: pick one shared
## suffix year, read both aim3_* and new_* with (ID, SRV_DT, AMOUNT_PAID, SEQ_ID),
## then compare (a) row counts, (b) distinct claim keys AND their overlap, which
## is the decisive number, (c) total dollars including the de-duplicated union as
## a consistency check on (b), and (d) whether SEQ_ID is a real cross-file claim
## identifier or merely a within-file row sequence -- if the latter, (b) is
## meaningless and the fallback match key is specified there too.
##
## WHY THIS PAPER CARES *MORE* THAN dual-bounds DOES, and the one thing to add
## when running it: double counting here corrupts the COVARIATE GRID, not just
## the outcome scale. prior_cost's quartile cutpoints are computed within the
## analytic sample (06_assemble_analytic_data.R), so inflated dollars move the
## boundaries, which changes which patients are in which quartile, which changes
## the 2^15 grid that estimate_att()'s grid-exact sparsity assumption is stated
## over. An outcome inflated by a constant factor rescales an ATT; a covariate
## whose bins are drawn on inflated values is a different covariate.
##
## SO: when running dual-bounds' CHECK 1, additionally report the QUARTILE
## CUTPOINTS of pre-index total cost under each candidate scope (both service
## families summed, versus aim3_* only, versus the de-duplicated union) on the
## SAA-eligible cohort. If the three cutpoint vectors agree to within rounding,
## the resolution is corroborated for this paper. If they disagree, say so to
## the PI: the confirmed scope would then need to be reopened for both papers.
##
## SUPPRESSION. Counts and totals over millions of records -- no small cells.
## Cutpoints are quantiles of a large sample, likewise. Still route anything
## leaving the server through smidata::smi_export().
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 2 -- msr_code ambiguity, and the SAA gap distribution
##                                                        [SHARED, REFERENCED]
##
## SPECIFICATION: ~/RAND/rprojects/smi/dual-bounds/analysis/real_data/
##                90_checks_tier1.R, CHECK 2 (written for AAP; the exposure is
##                now msr_saa -- confirm dual-bounds has made the same switch
##                before treating this check as shared again).
## REGISTRY:      smidata::inst/analyses/dual-bounds__application.yml#msr_code
##
## SUPERSEDED IN PART (PI, 2026-10-15): there is no gap cutoff to choose. The
## msr_saa period must START EXACTLY at s (first of the month on or after
## INDEX_DT) or the patient is dropped (smi_select_measure_row(anchor =
## "ceiling_month")), so the "month_gap distribution" and "cohort size at
## candidate cutoffs" parts of that specification collapse to: how many patients
## are matched at all, and the eligible / arm sizes under that single rule. Ties
## cannot occur across periods (a period is matched or not).
##
## Its output resolves OPEN_DECISIONS$msr_code for BOTH papers simultaneously,
## because exposure A here is `status: inherited` from dual-bounds -- the same
## smidata::smi_select_measure_row() call on the same table with the same threshold.
##
## ONE ADDITION FOR THIS PAPER. Report the resulting A = 1 / A = 0 counts, not
## just total cohort size. estimate_att() fits its mu
## tree on CONTROLS ALONE and requires at least leaf_budget * m_n control units
## (R/estimate_att.R's n_control check); a cohort that is large but
## thin in one arm can fail the estimator where it would not fail marbounds.
##
## SECOND ADDITION: the MSR_YR CONVENTION. 02 reads MSR_YR as the START of the
## 12-month measurement period and requires one starting exactly at s; if MSR_YR is really the END, A would be measured over a period
## ending at that date, i.e. up to 12 months earlier and overlapping the
## covariate window, which is not a post-index exposure. The start
## convention is "not yet confirmed against msr_saa itself" (02). The executable
## check below reports the nearest-period month_gap distribution under both
## readings (offset 0 and offset 12) as EVIDENCE. It cannot settle the question if each
## patient has a measurement period every year, since both readings then give
## similar gaps; in that case confirm from the SAS label or data dictionary.
##
## SUPPRESSION. The per-gap patient counts WILL have small cells in the tail.
## Suppress n < 11 (CMS/HIPAA Safe Harbor) via smidata::smi_suppress() before
## export. Do not collapse the tail by hand, and do not widen buckets to clear 11
## without saying so.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 3 -- are the 12 confirmed _YN columns genuinely {0, 1}?
##                                                      [THIS PAPER'S OWN]
##
## QUESTION. The PI-confirmed covariate selection (helpers/covariate_blocks.R)
## rests on the claim that these 12 columns are ALREADY BINARY and need no
## construction -- which is what makes a 15-column grid-exact design feasible at
## all. The fingerprint confirms the NAMES and the `_YN` suffix. It does not
## confirm the CODING. A SAS `_YN` flag is, in the wild, any of: 0/1 numeric,
## 1/2 numeric, "Y"/"N" character, "1"/"0" character, or TRUE/FALSE.
##
## WHY THIS IS NOT A COSMETIC CHECK. helpers/design_matrix.R's
## assert_binary_design_matrix() will catch a mis-coded column -- but the package
## itself will not, on the path that matters: check_att_data() (R/utils.R) does
## not examine X's coding at all, and estimate_att_crossfit() adds no further
## check, so a 1/2-coded column runs to completion and returns a number,
## discretized internally by optimaltrees at cutpoints nobody declared.
## estimate_att() does reject it, but only after the analytic data has been built
## -- i.e. after the multi-hour cost aggregation in 03. Knowing the coding BEFORE
## that run is the difference between a five-minute fix and a repeated run.
##
## THE CHECK. Read larger_smi_covariates with ID plus the 12 columns
## (unlist(covariate_blocks())). For each column report:
##
##   (a) CLASS and the full set of DISTINCT VALUES, including NA as its own
##       level. With 12 flags on 381,018 rows this is a 12-row table and is
##       cheap. Do NOT summarise to "n distinct = 2": 1/2 has two distinct
##       values and is wrong.
##
##   (b) MISSINGNESS RATE per column. An NA in a _YN flag is not automatically
##       "no": for a Z-code social-determinant flag (HOUSING_YN, FOOD_YN,
##       POVERTY_YN) it plausibly means "not screened", which is a different
##       thing from "screened, negative" and is itself correlated with
##       utilization. assert_binary_design_matrix() rejects any NA rather than
##       imputing, so a non-trivial rate on any of the 12 is a NEW PI decision
##       (recode, or replace that covariate) and would need its own registry
##       entry.
##
##   (c) PREVALENCE (mean, once coding is confirmed 0/1) per column, on the
##       SAA-eligible cohort rather than the full table. Two reasons. A flag at
##       0.1% prevalence is very nearly a constant column: it will not be
##       rejected by the constant-column check, but it spends a grid coordinate
##       to distinguish a handful of patients, which is a poor use of a scarce
##       budget. And prevalence within the A = 0 arm specifically bounds what the
##       mu tree can split on, since that tree is fit on controls alone.
##
##   (d) EXPOSURE-LEAKAGE CHECK on VAR9E_ANTIPSY_YN. This column is
##       antipsychotic use and the exposure is antipsychotic ADHERENCE (SAA).
##       It is included as a PRE-index confounder, which is correct only if its
##       measurement window precedes s (the first of the month on or after
##       INDEX_DT, so preceding INDEX_DT suffices). Confirm against the table's own
##       PRE_START / INDEX_DT columns that it does (a dictionary/column lookup,
##       NOT automated below). The executable part cross-tabulates
##       VAR9E_ANTIPSY_YN against A: near-perfect agreement would indicate the
##       column is measured over a window overlapping the exposure window, in
##       which case it is not a confounder and must be dropped. CAUTION: prior
##       antipsychotic use is also EXPECTED to correlate with later adherence
##       (that is why it is a confounder), and if the flag is ~1 in this cohort
##       the cross-tab has no variation. Agreement is consistent with window
##       overlap but is neither necessary nor sufficient evidence of it; only
##       the dictionary / PRE_START lookup settles it.
##       (VAR5G_AP_YN, another AP-related flag, is NOT in this paper's selection
##       -- but if a replacement covariate is ever needed, it carries the same
##       risk and needs the same check.)
##
##   (e) GLU_OR_LIPID_YN TIMING (NOT automated). Its window is anchored on the
##       index antipsychotic fill date, not INDEX_DT (helpers/covariate_blocks.R
##       audit note). Confirm that window ends at or before s; if the fill date
##       can fall after s the flag can be post-treatment and must be replaced.
##       (The fill date can fall after INDEX_DT but before s; that is fine.)
##
##   (f) VAR9E_ANTIPSY_YN PREVALENCE IN THE SAA COHORT. SAA's denominator
##       requires antipsychotic dispensing, so this flag may be near 1 within the
##       eligible cohort. A (near-)constant column triggers the constant-column
##       abort in 06 or spends a grid coordinate on nothing; the executable
##       part reports prevalence overall and by arm.
##
## SUPPRESSION. (a) is a value-set listing, no cells. (b) and (c) are rates over
## a large cohort, but (c) crossed with A gives 2x2 tables per column and (d) is
## explicitly a cross-tabulation -- these CAN contain cells under 11 for a rare
## flag. Route through smidata::smi_suppress() before export.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## CHECK 4 -- resolve whether larger_smi_medicaid_monthly_flag$TYPE can be
## safely ignored (OPEN_DECISIONS$enrollment_source's IP/OP sub-question)
##   [SHARED, REFERENCED -- identical check to dual-bounds' CHECK 3]
##
## QUESTION. PI (2026-09-24) deferred the IP-vs-OP choice and directed:
## proceed WITHOUT filtering on TYPE for now -- every
## larger_smi_medicaid_monthly_flag row counts as enrollment evidence
## (toward complete-case status) regardless of TYPE. That is only safe if
## TYPE does not fragment a single patient's own enrollment signal in a way
## that changes the answer: if a patient's IP and OP rows for the same
## (ID, YEAR_MONTH) never DISAGREE on MEDICAID_FLAG, ignoring TYPE is a
## no-op; if they can disagree, "ignore TYPE" silently picks whichever row
## happens to be read first/last, which is not a decision anyone made.
##
## THE CHECK. Identical to dual-bounds' CHECK 3 -- run once, applies to both
## papers (same table, same question). On the server, group
## larger_smi_medicaid_monthly_flag by (ID, YEAR_MONTH) and compute:
##
##   (a) HOW MANY PATIENTS carry BOTH TYPE values at all (any row with
##       TYPE=='IP' AND any row with TYPE=='OP' for the same ID, anywhere in
##       their history) vs. patients with only one TYPE ever.
##
##   (b) WITHIN an (ID, YEAR_MONTH) cell that has both TYPE=='IP' and
##       TYPE=='OP' rows, do MEDICAID_FLAG values ever DISAGREE (one row 1,
##       the other 0)? Tabulate agree/disagree counts. Disagreement means
##       "ignore TYPE" is not neutral, and the deferred IP/OP choice will
##       change results, not just documentation.
##
##   (c) If (b) finds disagreement, report its size (patient-months
##       affected) so the PI can see whether the deferred choice is a
##       rounding error or a real modeling decision before
##       OPEN_DECISIONS$enrollment_source's interim "ignore TYPE" default is
##       relied on for a real number.
##
## SUPPRESSION. (a) is a per-category patient count and (c) is a
## patient-month count -- suppress any cell with n < 11 (CMS/HIPAA Safe
## Harbor) via smidata::smi_suppress() before export.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 5 -- do the 18 service cost files cover what elig_coverage assumes?
##                                                      [THIS PAPER'S OWN]
##
## QUESTION. 02's elig_coverage filter (via smidata::smi_window_coverage())
## treats file suffix N as calendar year 2000 + N, Jan 1 to Dec 31, for both service
## families, and calls a patient-window "full" if it lies inside the pooled span
## 2016-01-01 .. 2025-01-01. Nothing has verified that. If suffixes are fiscal
## years, if new_* starts later than aim3_*, or if late-2024 claims are not yet
## complete, patients marked "full" can have truncated windows -- and truncation
## looks like LOW COST, not missing cost.
##
## THE CHECK (executable below). Stream all 18 files once and tabulate claims by
## month. Report, per file, the first and last claim month, the number of claims
## with a missing date, and whether both fall in 20yy; per family, the final
## year's last-three-months claim volume against that year's earlier months
## (a run-out heuristic, not a test); and the Managed Care Invoice claims
## matched on config_cost_setting_col, as corroboration that the exclusion in 03
## is doing something (a total of zero would mean it matches nothing).
##
## SUPPRESSION. Counts over millions of claims; the exported tables are
## per-file and per-family aggregates. smi_export() suppresses n < 11.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 6 -- exposure convention and exact-start sizes (msr_saa)  [THIS PAPER'S OWN]
## CHECK 3 (a)-(c), (d) cross-tab only, (f) as noted above   [THIS PAPER'S OWN]
##
## Specified in CHECK 2 and CHECK 3. Implemented below.
## ---------------------------------------------------------------------------


## ============================================================================
## EXECUTABLE SECTION -- server only
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "msr_filter.R"))
source(file.path(app_dir, "helpers", "period_start.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "disclosure.R"))

## Which of the executable checks to run, in order. The two that read only the
## covariate and msr_saa tables come first; CHECK 5 streams all 18 service cost files
## (about the cost of 03's file pass), so a late failure there costs nothing.
config_tier1_checks <- c("saa_convention", "yn_columns", "file_coverage")

tier1_export <- function(x, name, n_cols, derived_cols = NULL) {
  smidata::smi_export(
    x, name = name, kind = "table", n_cols = n_cols, derived_cols = derived_cols
  )
}

tier1_show <- function(x) cli::cli_verbatim(utils::capture.output(print(x, n = Inf)))

## The two tables CHECKS 3 and 6 share, read once. Mirrors 02's guards (zero
## rows, NA INDEX_DT, one row per patient, zero SAA rows) and, unlike 02's
## quiet-on-local read, always reports an MSR exclusion. It does NOT apply
## 02's elig_coverage filter, so "SAA-eligible" below means eligible before the
## coverage filter, and nothing here applies the complete-case restriction.
tier1_data <- new.env()
tier1_load <- function() {
  if (!is.null(tier1_data$base)) return(invisible(tier1_data))
  yn_cols <- unlist(covariate_blocks(), use.names = FALSE)
  base <- smi_read_pinned(SMI_KEYS$covariates, columns = c(SMI_COLS$covariates, yn_cols))
  if (nrow(base) == 0L) cli::cli_abort("{.val {SMI_KEYS$covariates}} returned zero rows.")
  if (anyDuplicated(base$ID) > 0L) {
    cli::cli_abort("{.val {SMI_KEYS$covariates}} is not one row per patient.")
  }
  if (anyNA(base$INDEX_DT)) cli::cli_abort("{.field INDEX_DT} has NA.")
  saa <- filter_to_msr_code(smi_read_pinned(SMI_KEYS$saa, columns = SMI_COLS$saa), SAA_MSR_CODE)
  if (nrow(saa) == 0L) {
    cli::cli_abort("No {.val {SAA_MSR_CODE}} rows in {.val {SMI_KEYS$saa}}.")
  }
  tier1_data$base <- base
  tier1_data$saa <- saa
  tier1_data$yn_cols <- yn_cols
  invisible(tier1_data)
}

## The selection 02 uses (PI, 2026-10-15): the msr_saa period starting EXACTLY at
## s = the first of the month on or after INDEX_DT; unmatched otherwise.
tier1_select_saa <- function() {
  d <- tier1_load()
  smidata::smi_select_measure_row(
    d$saa,
    index_dt = dplyr::transmute(d$base, ID = .data$ID, INDEX_DT = as.Date(.data$INDEX_DT)),
    achievement_threshold = config_saa_achievement_threshold,
    anchor = config_period_anchor
  )
}

## Nearest period to a target month by ABSOLUTE gap, before or after INDEX_DT.
## Not what 02 uses; only for the MSR_YR-convention evidence in CHECK 6.
tier1_select_saa_any_period <- function(target_offset_months, max_month_gap) {
  d <- tier1_load()
  smidata::smi_select_measure_row(
    d$saa,
    index_dt = dplyr::transmute(d$base, ID = .data$ID, INDEX_DT = as.Date(.data$INDEX_DT)),
    target_offset_months = target_offset_months,
    achievement_threshold = config_saa_achievement_threshold,
    max_month_gap = max_month_gap
  )
}

## ---- CHECK 6 ----------------------------------------------------------------

tier1_check_saa_convention <- function() {
  cli::cli_h2("CHECK 6 -- msr_saa exposure convention and exact-start sizes")
  d <- tier1_load()

  ## (i) EVIDENCE on the MSR_YR convention. If periods are anchored to each
  ## patient's index, the distribution of (period start - INDEX_DT month) has its
  ## mass at 0 under START and at -12 under END. If periods are calendar-anchored
  ## and recur every year, it cannot tell the two apart.
  ## msr_saa rows for IDs outside the cohort are expected and ignored, exactly as
  ## in 02; they have no INDEX_DT to compare against.
  idx <- dplyr::select(d$base, "ID", "INDEX_DT")
  in_cohort <- d$saa$ID %in% idx$ID
  cli::cli_alert_info(
    "{sum(!in_cohort)} msr_saa row{?s} are for IDs outside the cohort (excluded, as in 02)."
  )
  saa_in <- d$saa[in_cohort, , drop = FALSE]
  if (nrow(saa_in) == 0L) cli::cli_abort("No msr_saa rows for any cohort ID.")
  parsed <- smidata::smi_parse_msr_yr(as.character(saa_in$MSR_YR))
  joined <- dplyr::inner_join(dplyr::select(saa_in, "ID"), idx, by = "ID")
  rel <- smidata::smi_month_index(parsed$yyyy, parsed$mm) -
    smidata::smi_month_index(lubridate::year(joined$INDEX_DT), lubridate::month(joined$INDEX_DT))
  rel_bands <- tibble::tibble(
    band = cut(rel, c(-Inf, -13, -12, -2, -1, 0, 1, 11, 12, Inf),
               labels = c("<= -13", "-12", "-11..-2", "-1", "0", "1", "2..11", "12", ">= 13")),
  ) |>
    dplyr::count(.data$band, name = "n_rows") |>
    dplyr::mutate(band = as.character(.data$band))
  cli::cli_alert_info(
    "All msr_saa rows by (period start - INDEX_DT month). Mass at 0 supports the
     START reading, mass at -12 the END reading. If rows spread across all
     bands (a period every year, calendar-anchored) this cannot tell them
     apart, and the convention must be confirmed from the SAS label or the
     dictionary."
  )
  tier1_show(rel_bands)

  by_offset <- lapply(c(0L, 12L), function(off) {
    sel <- tier1_select_saa_any_period(off, Inf)
    tibble::tibble(
      target_offset_months = as.character(off), n_matched = sum(sel$matched),
      median_gap = stats::median(sel$month_gap, na.rm = TRUE),
      p90_gap = unname(stats::quantile(sel$month_gap, 0.9, na.rm = TRUE)),
      share_abs_gap_le_1 = mean(sel$month_gap <= 1L, na.rm = TRUE)
    )
  }) |> dplyr::bind_rows()
  cli::cli_alert_info(
    "|month_gap| of the nearest period under each reading of MSR_YR (offset 0 =
     START, 12 = END). Informative only if the readings differ materially."
  )
  tier1_show(by_offset)

  ## (ii) Cohort and arm sizes under the exact-start rule 02 uses. There is no
  ## gap cutoff to tune any more: a patient either has a period starting exactly
  ## at s or is dropped, so this is one DECLARED row, not a table of bands.
  chosen <- tier1_select_saa()
  n_with_rows <- sum(d$base$ID %in% d$saa$ID)
  by_band <- tibble::tibble(
    n_with_saa_rows = n_with_rows,
    n_matched = sum(chosen$matched),
    n_elig = sum(chosen$elig %in% TRUE),
    n_treated = sum(chosen$elig %in% TRUE & chosen$achieved %in% TRUE),
    n_control = sum(chosen$elig %in% TRUE & chosen$achieved %in% FALSE),
    n_achieved_na = sum(chosen$elig %in% TRUE & is.na(chosen$achieved))
  )
  cli::cli_alert_info(
    "SAA-matched (a period starting exactly at s), eligible (MSR_DEN >= 1) and
     arm sizes, BEFORE the coverage filter and the complete-case restriction.
     Everyone else is dropped. estimate_att() needs at least leaf_budget * m_n
     = {config_leaf_budget} * m_n controls AFTER both. A low matched share means
     msr_saa's periods do not start on the first of the month on or after
     INDEX_DT for most patients."
  )
  tier1_show(by_band)

  ## Exported form. with_rows >= matched >= elig is a chain, and treated + control
  ## + na = elig is a partition; smi_export() masks cells, not their differences.
  by_band_out <- tier1_disclosable_sizes(by_band)
  tier1_export(by_band_out, "doubletree_tier1_exact_start_sizes",
               n_cols = c("n_with_saa_rows", "n_matched", "n_elig", "n_treated",
                          "n_control", "n_achieved_na"))
  tier1_export(rel_bands, "doubletree_tier1_msr_yr_offsets", n_cols = "n_rows")
  tier1_export(by_offset, "doubletree_tier1_gap_by_offset",
               n_cols = "n_matched",
               derived_cols = c("median_gap", "p90_gap", "share_abs_gap_le_1"))
  invisible(list(by_offset = by_offset, by_band = by_band, rel_bands = rel_bands))
}

## ---- CHECK 3 (a)-(c), (d) cross-tab, (f) -------------------------------------

tier1_check_yn_columns <- function() {
  cli::cli_h2("CHECK 3 -- the 12 _YN columns")
  d <- tier1_load()
  yn_cols <- d$yn_cols
  cov <- d$base
  missing_cols <- setdiff(yn_cols, names(cov))
  if (length(missing_cols) > 0L) {
    cli::cli_abort("larger_smi_covariates is missing {.field {missing_cols}}.")
  }

  coding <- tibble::tibble(
    column = yn_cols,
    class = vapply(cov[yn_cols], function(x) class(x)[[1]], character(1)),
    n_missing = vapply(cov[yn_cols], function(x) sum(is.na(x)), integer(1))
  ) |> dplyr::mutate(na_rate = .data$n_missing / nrow(cov))
  value_sets <- vapply(cov[yn_cols], function(x) {
    paste(sort(unique(as.character(x)), na.last = TRUE), collapse = "|")
  }, character(1))
  cli::cli_alert_info(
    "(a)/(b) Coding and missingness over the full table. Anything other than
     0|1 is a mis-coded flag (1|2 has two distinct values and is wrong); a
     non-trivial missing rate on a Z-code flag is a new PI decision. The value
     sets are console-only (a rare stray code would be a small cell)."
  )
  tier1_show(dplyr::mutate(coding, values = unname(value_sets)))

  cohort <- tier1_select_saa() |>
    dplyr::filter(.data$elig %in% TRUE) |>
    dplyr::select("ID", A = "achieved") |>
    dplyr::inner_join(cov[c("ID", yn_cols)], by = "ID")
  binary <- yn_cols[vapply(yn_cols, function(col) {
    all(stats::na.omit(as.character(cohort[[col]])) %in% c("0", "1"))
  }, logical(1))]
  if (length(binary) < length(yn_cols)) {
    cli::cli_alert_warning(
      "{.field {setdiff(yn_cols, binary)}} {?is/are} not coded 0/1; prevalence
       skipped for {?it/them}."
    )
  }
  ## Counts of the 3-way split (yes / no / NA) by arm, unmasked (console only).
  count_cell <- function(v, arm, value) {
    in_arm <- if (is.na(arm)) rep(TRUE, length(v)) else cohort$A %in% arm
    sum(in_arm & if (is.na(value)) is.na(v) else (!is.na(v) & v == value))
  }
  prevalence <- lapply(binary, function(col) {
    v <- as.integer(as.character(cohort[[col]]))
    tibble::tibble(
      column = col,
      n_yes_control = count_cell(v, FALSE, 1L), n_no_control = count_cell(v, FALSE, 0L),
      n_na_control = count_cell(v, FALSE, NA_integer_),
      n_yes_treated = count_cell(v, TRUE, 1L), n_no_treated = count_cell(v, TRUE, 0L),
      n_na_treated = count_cell(v, TRUE, NA_integer_)
    )
  }) |> dplyr::bind_rows() |>
    dplyr::mutate(
      prevalence_control = .data$n_yes_control / (.data$n_yes_control + .data$n_no_control),
      prevalence_treated = .data$n_yes_treated / (.data$n_yes_treated + .data$n_no_treated)
    )
  ## Exported form: disjoint cells only (yes / no / NA within each arm). The
  ## overall counts are not exported (they are the sums of the two arms), and the
  ## arm sizes are published in doubletree_tier1_exact_start_sizes, so a small
  ## cell is masked together with a complementary cell of its (column, arm)
  ## partition, and the prevalence derived from a masked cell is masked.
  prevalence_out <- tier1_disclosable_prevalence(prevalence)
  cli::cli_alert_info(
    "(c)/(f) Prevalence in the SAA-eligible cohort ({nrow(cohort)} patients,
     before the coverage filter), by arm (console only: unmasked counts). Near 0 or near 1 spends a
     grid coordinate on (almost) nothing; the mu tree is fit on controls alone,
     so prevalence_control is the one that bounds it. VAR9E_ANTIPSY_YN near 1
     is expected to be the problem case."
  )
  tier1_show(prevalence)

  ## Console-only: the cells duplicate the prevalence table's yes/no counts.
  cli::cli_alert_info(
    "(d) VAR9E_ANTIPSY_YN x A (console only). Agreement is consistent with window
     overlap but also with strong confounding; see the CHECK 3(d) caution."
  )
  tier1_show(dplyr::count(
    cohort, A = as.character(.data$A),
    VAR9E_ANTIPSY_YN = as.character(.data$VAR9E_ANTIPSY_YN)
  ))

  tier1_export(coding, "doubletree_tier1_yn_coding",
               n_cols = "n_missing", derived_cols = "na_rate")
  count_cols <- grep("^n_", names(prevalence_out), value = TRUE)
  tier1_export(prevalence_out, "doubletree_tier1_yn_prevalence", n_cols = count_cols,
               derived_cols = c("prevalence_control", "prevalence_treated"))
  invisible(list(coding = coding, prevalence = prevalence))
}

## ---- CHECK 5 ----------------------------------------------------------------

tier1_check_file_coverage <- function(cost_sources = COST_SOURCES) {
  cli::cli_h2("CHECK 5 -- cost-file coverage")
  monthly <- smidata::smi_stream_cost_files(
    cost_sources,
    read_fn = function(key) smi_read_pinned(key, columns = SMI_COLS$cost_claims),
    process_fn = function(file_data, dataset_key) {
      d <- as.Date(file_data$SRV_DT)
      is_mci <- smidata::smi_classify_is_mci(file_data[[config_cost_setting_col]])
      tibble::tibble(
        key = lubridate::year(d) * 100L + lubridate::month(d),   # integer YYYYMM
        is_mci = is_mci, amount = file_data$AMOUNT_PAID
      ) |>
        dplyr::summarise(
          n_claims = dplyr::n(), n_mci = sum(.data$is_mci),
          mci_dollars = sum(.data$amount[.data$is_mci], na.rm = TRUE),
          .by = "key"
        ) |>
        dplyr::mutate(dataset_key = dataset_key, .before = 1)
    }
  ) |>
    dplyr::left_join(
      dplyr::select(cost_sources, "dataset_key", "family", "suffix_year"),
      by = "dataset_key"
    ) |>
    dplyr::mutate(
      expected_year = 2000L + .data$suffix_year,
      in_year = !is.na(.data$key) & .data$key %/% 100L == .data$expected_year
    )

  spans <- monthly |>
    dplyr::summarise(
      n_claims = sum(.data$n_claims),
      n_missing_date = sum(.data$n_claims[is.na(.data$key)]),
      n_outside_year = sum(.data$n_claims[!is.na(.data$key) & !.data$in_year]),
      months_with_claims = as.character(sum(.data$in_year)),   # of the 12 expected
      n_mci = sum(.data$n_mci), mci_dollars = sum(.data$mci_dollars),
      .by = c("dataset_key", "family", "suffix_year")
    ) |>
    dplyr::arrange(.data$family, .data$suffix_year)
  tier1_show(spans)

  bad_files <- spans$n_outside_year > 0L | spans$months_with_claims != "12"
  if (any(bad_files)) {
    cli::cli_alert_danger(
      "{sum(bad_files)} of {nrow(spans)} files have claims outside calendar year
       20yy (suffix yy) or fewer than 12 months with claims. elig_coverage's
       'full' is NOT trustworthy for those until this is understood: the suffix
       year may not be the service year, or the file has gaps."
    )
  } else {
    cli::cli_alert_success(
      "All {nrow(spans)} files hold claims in all 12 months of their own 20yy only."
    )
  }
  if (any(spans$n_missing_date > 0L)) {
    cli::cli_alert_warning(
      "{sum(spans$n_missing_date)} claims have a missing SRV_DT; they fall in no
       window and are silently absent from every cost sum."
    )
  }

  last_year <- max(cost_sources$suffix_year)
  runout <- monthly |>
    dplyr::filter(.data$in_year, .data$suffix_year == last_year) |>
    dplyr::arrange(.data$family, .data$key) |>
    dplyr::summarise(
      median_earlier = stats::median(.data$n_claims[seq_len(max(dplyr::n() - 3L, 1L))]),
      mean_last3 = mean(utils::tail(.data$n_claims, 3L)),
      .by = "family"
    ) |>
    dplyr::mutate(ratio = .data$mean_last3 / .data$median_earlier)
  cli::cli_alert_info(
    "Final-year ({2000L + last_year}) run-out (console only): mean claims in the
     last three months over the median of the earlier months, by family. A ratio
     well below 1 suggests incomplete run-out, so late-index patients' outcomes
     would be understated while marked 'full'."
  )
  tier1_show(runout)

  cli::cli_alert_info(
    "Managed Care Invoice claims matched on {.field {config_cost_setting_col}}:
     {sum(spans$n_mci)} claims, ${signif(sum(spans$mci_dollars), 4)}."
  )
  if (sum(spans$n_mci) == 0L) {
    cli::cli_alert_danger(
      "Zero Managed Care Invoice claims matched: the exclusion in 03 matches
       nothing on this column."
    )
  }

  tier1_export(spans, "doubletree_tier1_file_coverage",
               n_cols = c("n_claims", "n_missing_date", "n_outside_year", "n_mci"),
               derived_cols = "mci_dollars")
  invisible(list(spans = spans, runout = runout))
}

if (!interactive() && sys.nframe() == 0L) {
  if (!smidata::smi_is_server()) {
    cli::cli_abort(c(
      "The Tier-1 checks need the real files and cannot run in a local
       environment.",
      "i" = "Tier 0 is {.file application/run_tests.R}."
    ))
  }
  runners <- list(
    saa_convention = tier1_check_saa_convention,
    yn_columns = tier1_check_yn_columns,
    file_coverage = tier1_check_file_coverage
  )
  unknown <- setdiff(config_tier1_checks, names(runners))
  if (length(unknown) > 0L) {
    cli::cli_abort("Unknown entries in config_tier1_checks: {.val {unknown}}.")
  }
  ran <- lapply(runners[config_tier1_checks], function(f) f())
  cli::cli_alert_success(
    "Tier-1 checks run: {.val {names(ran)}}. Exports are in the outbox; list
     them with smidata::smi_outbox_manifest()."
  )
}
