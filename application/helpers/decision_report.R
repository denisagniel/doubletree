## ============================================================================
## application/helpers/decision_report.R
##
## Everything 08_decision_report.R needs that is not a smidata builder:
##   * report_record()        how 02/03/05/06 save the counts the report reads
##   * summarise_design()     the 15-column binary design: prevalence, prior_cost
##                            bins, floor counts, realized grid (raw counts, local
##                            intermediate only)
##   * design_patterns()      distinct covariate patterns (pure, tested by hand)
##   * report_levers()        the levers in force, read FROM the config objects
##   * build_decision_report() intermediates -> smidata::smi_report()
##   * toy_report_inputs()    hand-calculated toy intermediates (off-server run
##                            and tests/test-decision-report.R)
##
## DISCLOSURE. The saved intermediates (output/report_inputs/*.rds) hold RAW,
## unsuppressed counts and some patient-level rows; they live under
## application/output/ (gitignored, never exported) exactly like the parquet
## outputs 06 already writes. Every cell of every table returned by
## build_decision_report() is suppressed (cells < 11), through the smidata
## builders or the masks below. Cell masking does NOT stop a masked count being
## recovered as the difference of two published ones, so (helpers/disclosure.R)
## a small attrition step is merged into the next step, a count that is a small
## gap from its neighbour is hidden, and a table that decomposes a hidden count
## is withheld. Stand-alone counts: report_mask(); partitions:
## smidata::smi_suppress(); a 0/1 split of a published total:
## mask_binary_split(). Cross-table coverage is a guard (tests/
## test-decision-report.R), not a proof: read the rendered report.
##
## Needs smidata >= 0.5.0 with R/decision_report.R (smi_report_*). Loaded after
## _config.R (levers are read from its objects).
## ============================================================================

## config_report_min_cell and the gap helpers live in helpers/disclosure.R.

#' Directory the pipeline stages save report intermediates to.
report_inputs_dir <- function(app_dir) file.path(app_dir, "output", "report_inputs")

#' Save one stage's report intermediate.
#'
#' @param name File stem, e.g. `"02_population"`.
#' @param value A list. A `.meta` element (environment, ingest id, time) is added
#'   so the report can refuse to mix stages from different runs or label a
#'   Tier-0 fixture run as such.
#' @param app_dir The application directory.
#' @return The path, invisibly.
report_record <- function(name, value, app_dir) {
  dir <- report_inputs_dir(app_dir)
  fs::dir_create(dir)
  value$.meta <- list(
    env = if (requireNamespace("smidata", quietly = TRUE)) smidata::smi_env() else "unknown",
    ingest_id = if (exists("effective_ingest_id")) effective_ingest_id() else NA_character_,
    recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  path <- file.path(dir, paste0(name, ".rds"))
  saveRDS(value, path)
  invisible(path)
}

#' Mask stand-alone counts: 1..(min_cell - 1) become NA, zero stays zero.
report_mask <- function(x, min_cell = config_report_min_cell) {
  x[!is.na(x) & x > 0 & x < min_cell] <- NA
  x
}

#' Mask a 0/1 split of a published total.
#'
#' `n_yes + n_no` is the published sample size, so a small `n_yes` or a small
#' `n_no` discloses the other; both are masked, with the rates derived from them.
#' Zero is not masked (it is a degeneracy flag, not a small cell).
mask_binary_split <- function(tb, min_cell = config_report_min_cell) {
  small <- function(v) !is.na(v) & v > 0 & v < min_cell
  hit <- small(tb$n_yes) | small(tb$n_no)
  for (cl in intersect(c("n_yes", "n_no", "prevalence", "minority_share"), names(tb))) {
    tb[[cl]][hit] <- NA
  }
  tb
}

#' Distinct covariate patterns of a binary design matrix.
#'
#' @param X Data frame of 0/1 columns.
#' @return A list: `n_patterns` (distinct rows) and `n_in_small_patterns`
#'   (patients whose pattern is shared by fewer than `min_cell` patients).
design_patterns <- function(X, min_cell = config_report_min_cell) {
  key <- do.call(paste, c(unname(as.list(as.data.frame(X))), sep = "|"))
  size <- table(key)
  list(
    n_patterns = length(size),
    n_in_small_patterns = sum(size[size < min_cell])
  )
}

#' Raw counts describing the final analytic design (unsuppressed).
#'
#' @param analytic_data 06's analytic data: `A`, `prior_cost`, and the 15 design
#'   columns.
#' @param n_floored Patients whose negative `prior_cost` was floored before
#'   binning (06's `n_floored`).
#' @param floor_value The floor (the OPEN `prior_cost_negative_floor`).
#' @param cutpoints The realized quartile cutpoints (named numeric).
#' @return A list of plain data frames and scalars; **raw, never export**.
summarise_design <- function(analytic_data, n_floored, floor_value, cutpoints) {
  cols <- design_matrix_columns()
  X <- as.data.frame(analytic_data[, cols, drop = FALSE])
  n <- nrow(X)
  blocks <- covariate_blocks()
  block_of <- stats::setNames(
    rep(names(blocks), lengths(blocks)), unlist(blocks, use.names = FALSE)
  )
  n_yes <- vapply(X, function(v) sum(v == 1L), numeric(1))
  prevalence <- tibble::tibble(
    covariate = cols,
    block = unname(ifelse(cols %in% names(block_of), block_of[cols], "prior_cost quartile")),
    n_yes = as.integer(n_yes),
    n_no = as.integer(n - n_yes)
  )
  q_cols <- prior_cost_dummy_names()
  n_q <- c(n - sum(n_yes[q_cols]), n_yes[q_cols])
  bins <- tibble::tibble(
    bin = c("Q1 (reference)", "Q2", "Q3", "Q4"), n = as.integer(n_q)
  )
  pat <- design_patterns(X)
  list(
    n = n,
    prevalence = prevalence,
    bins = bins,
    cutpoints = cutpoints,
    floor = tibble::tibble(
      quantity = c("n_final", "n_negative_prior_cost_floored",
                   "n_prior_cost_at_floor_after_flooring"),
      n = as.integer(c(
        n, n_floored,
        sum(pmax(analytic_data$prior_cost, floor_value) == floor_value)
      ))
    ),
    floor_value = floor_value,
    grid = tibble::tibble(
      quantity = c("grid_atoms (2^15)", "distinct_patterns_realized",
                   "patients_in_patterns_shared_by_lt11"),
      n = c(2^length(cols), pat$n_patterns, pat$n_in_small_patterns)
    ),
    arms = tibble::tibble(
      arm = c("A=0 (control)", "A=1 (treated)"),
      n = as.integer(c(sum(analytic_data$A == 0L), sum(analytic_data$A == 1L)))
    )
  )
}

## ---- levers: read from the config objects, never retyped -------------------

#' The levers in force, built from `_config.R` so the report cannot drift.
#'
#' @param meta The `.meta` of the intermediates (environment, ingest id).
#' @return A named list for `smidata::smi_report(levers =)`.
report_levers <- function(meta = list(env = NA_character_, ingest_id = NA_character_)) {
  es <- confirmed_value("enrollment_source")
  off <- function(a, k) if (k == 0L) a else sprintf("%s%+dmo", a, as.integer(k))
  w <- DOUBLETREE_WINDOWS
  anchor_label <- ifelse(w$anchor == "s", "s", "INDEX_DT")
  win <- c(
    sprintf("%s [%s, %s)", w$window,
            mapply(off, anchor_label, w$start_months),
            mapply(off, anchor_label, w$end_months)),
    "exposure: msr_saa period, [s, s+12mo) (identified, not summed)",
    "covariates (12 _YN flags): [INDEX_DT-12mo, INDEX_DT) (INDEX_DT-anchored, as prior_cost)",
    sprintf("complete case: [s, s+%dmo)", config_complete_case_months)
  )
  open_status <- vapply(OPEN_DECISIONS, function(d) d$status, character(1))
  list(
    data_environment = sprintf("%s (ingest %s)", meta$env, meta$ingest_id),
    period_rule = sprintf(
      paste("s = first of the month on or after INDEX_DT; the %s period must START",
            "exactly at s (anchor = \"%s\", no nearest period, no gap tolerance);",
            "no such row, or MSR_DEN < 1 -> patient DROPPED"),
      SAA_MSR_CODE, config_period_anchor
    ),
    measure_code = sprintf("%s (msr_saa$MSR), rows filtered to this literal before selection",
                           SAA_MSR_CODE),
    achievement_threshold = sprintf("MSR_NUM / MSR_DEN >= %s", config_saa_achievement_threshold),
    windows = paste(trimws(win), collapse = "; "),
    negative_amount_policy = sprintf(
      "%s (negative AMOUNT_PAID lines kept; window cost = gross positive + negative)",
      config_negative_amount_policy
    ),
    cost_families = sprintf(
      "service only: %s (%d files, suffix years %d-%d); pharmacy excluded",
      paste(unique(COST_SOURCES$family), collapse = ", "), nrow(COST_SOURCES),
      min(config_cost_years), max(config_cost_years)
    ),
    mci_exclusion = sprintf("lines whose %s matches \"%s\" excluded from every sum",
                            config_cost_setting_col, config_exclude_setting_regex),
    prices = "nominal dollars (no inflation adjustment)",
    coverage_eligibility = sprintf(
      "patients whose %s window(s) extend past the cost files' date coverage are excluded in 02",
      paste(config_coverage_required_windows, collapse = " and ")
    ),
    enrollment_rule = sprintf(
      "complete case over [s, s+%dmo): type_policy = %s, TYPE filter = %s, gap_months = %d, gap_mode = %s",
      config_complete_case_months, es$type_policy, es$type_filter,
      as.integer(es$gap_months), es$gap_mode
    ),
    prior_cost_floor = sprintf(
      "net prior_cost below 0 set to %s before the quartile cut (the raw dollars are not floored); decision prior_cost_negative_floor",
      format(OPEN_DECISIONS$prior_cost_negative_floor$value)
    ),
    prior_cost_discretization = paste(
      "within-sample quartiles (p25, p50, p75) of the FINAL complete-case sample;",
      "bins (-Inf, c1], (c1, c2], (c2, c3], (c3, Inf); Q1 reference; 3 dummies"
    ),
    covariate_set = paste(design_matrix_columns(), collapse = ", "),
    winsorizing = "not applicable: no winsorizing or trimming of Y or prior_cost",
    provider_attribution = "not applicable: patient-level ATT, no provider attribution",
    open_decisions = paste0(names(open_status), "=", open_status, collapse = "; ")
  )
}

## ---- the report ------------------------------------------------------------

.mask_partition <- function(tb, n_col, derived = NULL) {
  smidata::smi_suppress(tb, n_cols = n_col, derived_cols = derived,
                        threshold = config_report_min_cell)
}

#' Build the review report from the stages' saved intermediates.
#'
#' @param inputs List with `p02`, `p03`, `p06` and optionally `p05` (the
#'   contents of the four `.rds` files).
#' @return A `smidata::smi_report()`.
build_decision_report <- function(inputs, title = "doubletree application: decision report") {
  p02 <- inputs$p02; p03 <- inputs$p03; p05 <- inputs$p05; p06 <- inputs$p06
  mc <- config_report_min_cell

  ## Stage-consistency guards: a report that mixes two runs is worse than none.
  n_cohort_02 <- utils::tail(p02$flow$n, 1L)
  if (!identical(as.integer(p06$n_pre), as.integer(n_cohort_02))) {
    cli::cli_abort(c(
      "The 02 and 06 intermediates disagree on the eligible cohort size.",
      "x" = "02 recorded {n_cohort_02}; 06 recorded {p06$n_pre}.",
      "i" = "They come from different runs; re-run the pipeline end to end."
    ))
  }
  envs <- unique(c(p02$.meta$env, p03$.meta$env, p06$.meta$env, p05$.meta$env))
  if (length(envs) != 1L) {
    cli::cli_abort("Intermediates come from different environments: {.val {envs}}.")
  }

  ## 2. attrition: 02's steps, then 06's. Asserted-zero steps are listed so the
  ## reviewer sees that nothing else was dropped.
  flow <- data.frame(
    step = c(p02$flow$step,
             sprintf("complete case: fully enrolled over [s, s+%dmo)",
                     config_complete_case_months),
             "no missing Y, prior_cost, or covariate (asserted; the run aborts otherwise)"),
    n = c(p02$flow$n, p06$n_complete_case, p06$n_final),
    reason = c(p02$flow$reason,
               "not fully enrolled (insufficient coverage or absent from the monthly flag table)",
               "none dropped"),
    stringsAsFactors = FALSE
  )
  ## A small step (1..10 patients) is not a cell the reader can be denied while
  ## both neighbouring counts are shown, so such steps are merged into the next
  ## published step instead of masked.
  flow_pub <- coarsen_flow(flow, mc)

  ## Exposure-period counts. patients / matched / unmatched / den_below_1 are
  ## published ONCE, in the flow; here only eligible and achieved remain, and
  ## every table that decomposes `eligible` is withheld when the flow hid it.
  period_full <- smidata::smi_report_measure_period(p02$selection, p02$index_dt)
  n_eligible <- sum(p02$selection$elig)
  n_achieved <- sum(p02$selection$achieved, na.rm = TRUE)
  eligible_hidden <- !(n_eligible %in% flow_pub$n)
  achieved_hidden <- !eligible_hidden && n_eligible - n_achieved > 0L &&
    n_eligible - n_achieved < mc
  withheld <- function(what) tibble::tibble(
    note = sprintf("withheld: %s would let a masked count (< %d) be recovered by subtraction.",
                   what, mc)
  )
  period_summary <- tibble::tibble(
    quantity = c("eligible", "achieved"),
    n = c(if (eligible_hidden) NA_integer_ else as.integer(n_eligible),
          if (achieved_hidden) NA_integer_ else as.integer(n_achieved))
  )
  period_summary$n <- report_mask(period_summary$n, mc)
  decomposed <- eligible_hidden || achieved_hidden
  period <- list(
    summary = period_summary,
    by_period_start = if (decomposed) {
      withheld("the per-period eligible/achieved counts")
    } else period_full$by_period_start,
    index_day = if (eligible_hidden) {
      withheld("the per-day-of-month eligible counts")
    } else period_full$index_day
  )

  ## The two windows overlap (not a partition), so primary suppression only.
  coverage_by_window <- p02$coverage_counts
  coverage_by_window$n_matched <- report_mask(as.integer(coverage_by_window$n_matched))
  coverage_by_window$n_not_full <- report_mask(as.integer(coverage_by_window$n_not_full))

  msr_filter <- tibble::tibble(
    quantity = c("msr_saa_rows_not_SAA_excluded", "msr_saa_rows_for_ids_outside_cohort"),
    n = report_mask(as.integer(c(p02$n_msr_rows_excluded, p02$n_saa_rows_out_of_cohort)))
  )

  ## 4. costs
  final_ids <- p06$final_ids
  cost_final <- p03$cost[p03$cost$id %in% final_ids, , drop = FALSE]
  cost_known <- p03$known
  ## cohort - known (patients with no cost) and cohort - final are drops too.
  unknown_gap <- cost_known$n_cohort - cost_known$n_known
  cost_known$n_known[unknown_gap > 0L & unknown_gap < mc] <- NA
  cost_known$n_cohort <- report_mask(as.integer(cost_known$n_cohort))
  cost_known$n_known <- report_mask(as.integer(cost_known$n_known))
  n_cohort_cost <- max(p03$known$n_cohort)
  n_excluded_cost <- n_cohort_cost - length(final_ids)
  costs_final_sample <- if (n_excluded_cost > 0L && n_excluded_cost < mc) {
    withheld("the final-sample cost table next to the cohort cost table")
  } else {
    smidata::smi_report_costs(cost_final, min_cell = mc)
  }

  ## 5. enrollment
  if (is.null(p05)) {
    enrollment <- tibble::tibble(
      note = "No 05 intermediate: complete-case status was not computed in this run (Tier-0 fixture limitation or 05 not run)."
    )
    enrollment_by_arm <- enrollment
  } else {
    enrollment <- smidata::smi_report_enrollment(p05$coverage, p05$window_months)
    ## One arm-level partition only. The by-lag-by-arm cells are NOT published:
    ## they are cells of the same cross-tab, so a masked one is the difference
    ## between this table and the lag table. n is not published either (it is
    ## the sum of the columns, and the rate derived from them is masked together).
    enrollment_by_arm <- p05$by_arm |>
      dplyr::mutate(
        n_complete = as.integer(.data$n_complete),
        n_not_complete = as.integer(.data$n - .data$n_complete - .data$n_unresolved),
        n_unresolved = as.integer(.data$n_unresolved)
      ) |>
      dplyr::select("A", "n_complete", "n_not_complete", "n_unresolved", "rate") |>
      .mask_partition(c("n_complete", "n_not_complete", "n_unresolved"), "rate")
  }

  ## 6. project-specific: prior_cost floor, bins, cut points, design degeneracy
  d <- p06$design
  ## at_floor - floored = patients whose raw prior_cost is exactly the floor: a
  ## small value there is recoverable from the two counts, so at_floor is hidden.
  floor_n <- d$floor$n
  zero_gap <- floor_n[3] - floor_n[2]
  if (zero_gap > 0L && zero_gap < mc) floor_n[3] <- NA
  prior_cost_floor <- tibble::tibble(quantity = d$floor$quantity, n = report_mask(floor_n))
  prior_cost_bins <- d$bins |>
    dplyr::mutate(pct = round(100 * .data$n / sum(.data$n), 2)) |>
    .mask_partition("n", "pct")
  cuts <- d$cutpoints
  prior_cost_cutpoints <- tibble::tibble(
    percentile = names(cuts),
    cutpoint_usd = if (d$n >= mc) round(unname(cuts), 2) else rep(NA_real_, length(cuts))
  )
  ## The three quartile dummies are the Q2-Q4 counts of prior_cost_bins; they are
  ## published there once (masked jointly with Q1), not again here.
  design_prevalence <- d$prevalence |>
    dplyr::filter(.data$block != "prior_cost quartile") |>
    dplyr::mutate(
      prevalence = round(.data$n_yes / (.data$n_yes + .data$n_no), 4),
      minority_share = pmin(.data$prevalence, 1 - .data$prevalence)
    ) |>
    mask_binary_split()
  design_grid <- tibble::tibble(quantity = d$grid$quantity, n = report_mask(d$grid$n))
  final_arms <- .mask_partition(d$arms, "n")

  smidata::smi_report(
    title,
    levers = report_levers(p06$.meta),
    flow = smidata::smi_report_flow(flow_pub, min_cell = mc),
    measure_period = period,
    coverage_by_window = coverage_by_window,
    msr_filter = msr_filter,
    costs_cohort = smidata::smi_report_costs(p03$cost, min_cell = mc),
    costs_final_sample = costs_final_sample,
    cost_categories = smidata::smi_report_cost_categories(p03$cost, min_cell = mc),
    cost_known = cost_known,
    enrollment = enrollment,
    enrollment_by_arm = enrollment_by_arm,
    prior_cost_floor = prior_cost_floor,
    prior_cost_bins = prior_cost_bins,
    prior_cost_cutpoints = prior_cost_cutpoints,
    design_prevalence = design_prevalence,
    design_grid = design_grid,
    final_arms = final_arms
  )
}

## ---- toy intermediates (hand-calculated; see tests/test-decision-report.R) ---

#' Toy report intermediates with hand-calculable values.
#'
#' Same shape 02/03/05/06 record. Not data, not a stand-in for data: it exists so
#' 08 runs end to end off-server and so the test can pin counts by arithmetic.
#' The sections are independent of each other except where the test says so.
#'
#' @return `list(p02, p03, p05, p06)`.
toy_report_inputs <- function() {
  meta <- list(env = "toy", ingest_id = NA_character_, recorded = NA_character_)

  ## ---- 02: 200 patients p001..p200 -------------------------------------
  ## p001-p100 INDEX_DT 2020-03-01 (a 1st: s = INDEX_DT), p101-p200 2020-03-15
  ## (s = 2020-04-01). SAA rows: p001-p180. Period starting exactly at s:
  ## p001-p160. MSR_DEN >= 1: p001-p150. Achieved: p001-p060 and p101-p130.
  ids <- sprintf("p%03d", 1:200)
  index_dt <- tibble::tibble(
    ID = ids,
    INDEX_DT = as.Date(c(rep("2020-03-01", 100), rep("2020-03-15", 100)))
  )
  matched <- seq_along(ids) <= 160L
  elig <- seq_along(ids) <= 150L
  achieved_idx <- c(1:60, 101:130)
  selection <- tibble::tibble(
    ID = ids,
    matched = matched,
    elig = elig,
    MSR_DEN = dplyr::if_else(elig, 1, 0),
    MSR_YR_selected = dplyr::case_when(
      !matched ~ NA_character_,
      seq_along(ids) <= 100L ~ "RY2020-03",
      TRUE ~ "RY2020-04"
    ),
    achieved = dplyr::if_else(elig, seq_along(ids) %in% achieved_idx, NA)
  )
  ## Coverage: of the 150 eligible, 14 (p137-p150) fall outside the cost files.
  flow02 <- data.frame(
    step = c("base cohort (larger_smi_covariates)",
             "has >= 1 msr_saa row for the SAA measure",
             "a period starts exactly at s",
             "MSR_DEN >= 1 (eligible for SAA)",
             "cost-file coverage of the outcome and prior_cost windows"),
    n = c(200L, 180L, 160L, 150L, 136L),
    reason = c("", "no SAA rows", "rows exist but none starts exactly at s",
               "MSR_DEN < 1", "a gating window extends past the cost files"),
    stringsAsFactors = FALSE
  )
  p02 <- list(
    flow = flow02, selection = selection, index_dt = index_dt,
    coverage_counts = tibble::tibble(
      window = c("outcome", "prior_cost"), gating = c(TRUE, TRUE),
      n_matched = c(160L, 160L), n_not_full = c(20L, 14L)
    ),
    n_msr_rows_excluded = 7L, n_saa_rows_out_of_cohort = 40L,
    .meta = meta
  )

  ## ---- 03: 136 analysis-cohort patients p001..p136 -----------------------
  ## T0 = INDEX_DT = 2019-06-01: outcome [2020-06-01, 2021-06-01), prior_cost
  ## [2018-06-01, 2019-06-01). Outcome file (2020-07-01): +200 Inpatient and +100
  ## Outpatient each; p001-p012 also -50 Outpatient; p001-p012 also a 999 Managed
  ## Care Invoice line. Prior file (2018-09-01): +100 Inpatient each; p001-p019
  ## also -150 Inpatient.
  cids <- sprintf("p%03d", 1:136)
  sel <- function(k) cids[seq_len(k)]
  claims20 <- dplyr::bind_rows(
    tibble::tibble(ID = cids, SRV_DT = as.Date("2020-07-01"),
                   AMOUNT_PAID = 200, SS_DESC = "Inpatient"),
    tibble::tibble(ID = cids, SRV_DT = as.Date("2020-07-01"),
                   AMOUNT_PAID = 100, SS_DESC = "Outpatient"),
    tibble::tibble(ID = sel(12), SRV_DT = as.Date("2020-07-01"),
                   AMOUNT_PAID = -50, SS_DESC = "Outpatient"),
    tibble::tibble(ID = sel(12), SRV_DT = as.Date("2020-07-01"),
                   AMOUNT_PAID = 999, SS_DESC = "Managed Care Invoice")
  )
  claims18 <- dplyr::bind_rows(
    tibble::tibble(ID = cids, SRV_DT = as.Date("2018-09-01"),
                   AMOUNT_PAID = 100, SS_DESC = "Inpatient"),
    tibble::tibble(ID = sel(19), SRV_DT = as.Date("2018-09-01"),
                   AMOUNT_PAID = -150, SS_DESC = "Inpatient")
  )
  src <- tibble::tibble(
    dataset_key = c("toy_svc_18", "toy_svc_20"), family = "aim3_svc_cost",
    suffix_year = c(18L, 20L), amount_col = "AMOUNT_PAID", date_col = "SRV_DT",
    setting_col = "SS_DESC"
  )
  cohort <- tibble::tibble(ID = cids, INDEX_DT = as.Date("2019-06-01"),
                           T0 = as.Date("2019-06-01"))
  cost <- suppressMessages(suppressWarnings(smidata::smi_sum_cost_in_windows(
    doubletree_patient_windows(cohort), src,
    negative_amount_policy = config_negative_amount_policy,
    exclude_setting_regex = config_exclude_setting_regex,
    read_fn = function(dataset_key) {
      if (dataset_key == "toy_svc_18") claims18 else claims20
    },
    progress = FALSE
  )))
  p03 <- list(
    cost = cost,
    known = tibble::tibble(window = c("outcome", "prior_cost"),
                           n_cohort = c(136L, 136L), n_known = c(136L, 136L)),
    .meta = meta
  )

  ## ---- 05: complete case over the 136 ------------------------------------
  ## 100 covered, 30 insufficient, 6 absent.
  cov <- data.frame(n_covered = c(rep(24L, 100), rep(20L, 30), rep(0L, 6)))
  attr(cov, "tally") <- list(n_cohort = 136L, n_covered = 100L, n_insufficient = 30L,
                             n_absent = 6L, n_uncovered_at_index = 0L)
  lag_arm <- tibble::tribble(
    ~t0_lag_months, ~A, ~reason, ~n,
    0L, 1L, "covered", 40L,
    0L, 1L, "insufficient_coverage", 10L,
    0L, 0L, "covered", 35L,
    0L, 0L, "insufficient_coverage", 9L,
    0L, 0L, "absent_from_flag_table", 6L,
    1L, 1L, "covered", 15L,
    1L, 1L, "insufficient_coverage", 11L,
    1L, 0L, "covered", 10L
  )
  by_arm <- tibble::tibble(
    A = c(0L, 1L), n = c(60L, 76L), n_complete = c(45L, 55L), n_unresolved = c(0L, 0L)
  ) |> dplyr::mutate(rate = .data$n_complete / .data$n)
  p05 <- list(coverage = cov, window_months = 24L, by_lag_arm = lag_arm,
              by_arm = by_arm, .meta = meta)

  ## ---- 06: final sample of 100 (p001-p100) --------------------------------
  ## Design column j (j = 1..12) is 1 for row i (i = 0..99) iff i %% (j + 1) == 0.
  yn <- unlist(covariate_blocks(), use.names = FALSE)
  i <- 0:99
  X <- as.data.frame(lapply(seq_along(yn), function(j) as.integer(i %% (j + 1L) == 0L)))
  names(X) <- yn
  ## prior_cost = (row - 20) * 100 for row 1..100: rows 1-19 negative (floored),
  ## row 20 is 0, rows 21-100 are 100..8000.
  raw_pc <- (seq_len(100) - 20) * 100
  floored <- floor_negative_prior_cost(raw_pc, 0)
  dummies <- discretize_prior_cost(floored)
  analytic <- dplyr::bind_cols(
    tibble::tibble(ID = sprintf("p%03d", 1:100),
                   A = as.integer(seq_len(100) <= 55L), prior_cost = raw_pc),
    X, dummies
  )
  p06 <- list(
    n_pre = 136L, n_complete_case = 100L, n_final = 100L,
    final_ids = analytic$ID,
    design = summarise_design(analytic, n_floored = attr(floored, "n_floored"),
                              floor_value = 0, cutpoints = attr(dummies, "cutpoints")),
    .meta = meta
  )
  list(p02 = p02, p03 = p03, p05 = p05, p06 = p06)
}
