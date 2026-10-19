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
## DISCLOSURE. The saved intermediates (output/private/report_inputs/*.rds) hold RAW,
## unsuppressed counts and some patient-level rows; they live under
## application/output/ (gitignored, never exported) exactly like the parquet
## outputs 06 already writes. build_decision_report() masks every cell below 11
## (smidata builders; helpers/disclosure.R for this paper's own tables), DECLARES
## the arithmetic between published counts as smidata identities and gaps
## (smi_identity(), smi_gap(): a count published in several tables, a total beside
## its breakdown, a cohort count beside a subset), and hands them to
## smidata::smi_report(), which masks further cells until no declared identity
## returns a masked value and then AUDITS the result. The audit guarantees that
## no masked value is determined by the published cells THROUGH THE DECLARED
## IDENTITIES. It does not see undeclared identities, does not reason about
## dollar dominance, and does not bound a masked cell inside 1..10. Read the
## rendered report's header line, and review before releasing.
##
## Needs smidata >= 0.5.0 with R/decision_report.R (smi_report_*). Loaded after
## _config.R (levers are read from its objects).
## ============================================================================

## config_report_min_cell and the masking helpers live in helpers/disclosure.R.

#' Directory the pipeline stages save report intermediates to.
report_inputs_dir <- function(app_dir) file.path(app_dir, "output", "private", "report_inputs")

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
    winsorizing = if (OPEN_DECISIONS$outcome_winsor_quantile$value >= 1) {
      "none: Y is raw cost (outcome_winsor_quantile = 1); prior_cost is not winsorized or trimmed"
    } else {
      sprintf(
        "Y capped above its p%s quantile over the final complete-case sample (pooled across arms); estimand is the ATT on winsorized cost; prior_cost is not winsorized or trimmed; decision outcome_winsor_quantile",
        format(100 * OPEN_DECISIONS$outcome_winsor_quantile$value)
      )
    },
    provider_attribution = "not applicable: patient-level ATT, no provider attribution",
    open_decisions = paste0(names(open_status), "=", open_status, collapse = "; ")
  )
}

## ---- the report ------------------------------------------------------------

#' Build the review report from the stages' saved intermediates.
#'
#' Tables are masked first (cells < 11), the arithmetic between them is declared
#' (see "identities" below), and `smidata::smi_report()` masks further cells until
#' no declared identity returns a masked value, then audits.
#'
#' @param inputs List with `p02`, `p03`, `p06` and optionally `p05` (the
#'   contents of the four `.rds` files).
#' @return An audited `smidata::smi_report()`.
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
  ## The winsorization lever is read from the config; make sure 06 used it.
  if (!is.null(p06$winsor) &&
      !isTRUE(all.equal(p06$winsor$q, OPEN_DECISIONS$outcome_winsor_quantile$value))) {
    cli::cli_abort(c(
      "06 winsorized Y at q = {p06$winsor$q} but the config now says
       {OPEN_DECISIONS$outcome_winsor_quantile$value}.",
      "i" = "Re-run 06 (and 07) so the report describes the level actually used."
    ))
  }

  envs <- unique(c(p02$.meta$env, p03$.meta$env, p06$.meta$env, p05$.meta$env))
  if (length(envs) != 1L) {
    cli::cli_abort("Intermediates come from different environments: {.val {envs}}.")
  }

  ## ---- identities: collected as the tables are built -------------------------
  ## Every identity is declared from the TRUE (unmasked) counts, and only when the
  ## recorded counts actually satisfy it: a stage mismatch is reported, never
  ## declared away. `member()` is one side of a cluster of cells that publish the
  ## same quantity; `cluster()` links every member to the first single cell.
  ids <- list()
  to_blank <- list()
  member <- function(cells, val) list(cells = cells, val = as.numeric(val))
  ## `value_hidden`: the flow does not publish this count (a small step was merged
  ## away), so no other table may publish it either, or the hidden drop is the
  ## difference of the two. Its single-cell copies are blanked after the tables
  ## are built instead of linked by identities.
  cluster <- function(label, members, optional = FALSE, value_hidden = FALSE) {
    members <- members[!vapply(members, is.null, NA)]
    if (value_hidden) {
      for (m in members) {
        if (inherits(m$cells, "smi_cell") && !m$cells$all) {
          to_blank[[length(to_blank) + 1L]] <<- m$cells
        }
      }
      return(invisible())
    }
    single <- which(vapply(members, function(m) inherits(m$cells, "smi_cell"), NA))
    if (length(single) == 0L || length(members) < 2L) return(invisible())
    ref <- members[[single[[1L]]]]
    for (m in members[-single[[1L]]]) {
      if (!isTRUE(all.equal(ref$val, m$val))) {
        if (!optional) {
          cli::cli_alert_warning(
            "Identity {.val {label}} not declared: recorded counts disagree ({ref$val} vs {m$val}); the stages may come from different runs."
          )
        }
        next
      }
      ids[[length(ids) + 1L]] <<- smidata::smi_identity(ref$cells, m$cells, label = label)
    }
    invisible()
  }
  gap <- function(big, small, label) {
    ids[[length(ids) + 1L]] <<- smidata::smi_gap(big, small, label = label)
  }
  cell <- smidata::smi_cell
  all_rows <- smidata::smi_cells_all

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
  ## Eligible (MSR_DEN >= 1) = achieved + not achieved is published below, so the
  ## flow must not publish any count within 1..10 of it either: the difference is
  ## a small cell (for example matched - eligible = patients with MSR_DEN < 1).
  sel <- p02$selection
  ach <- sel$elig & sel$achieved %in% TRUE
  n_eligible <- sum(sel$elig); n_matched <- sum(sel$matched)
  n_achieved <- sum(ach); n_not_achieved <- sum(sel$elig & !ach)
  ## A small step (1..10 patients) is merged into the next published step
  ## (smi_coarsen_flow), and smi_report_flow() declares the chain as identities.
  flow_pub <- smidata::smi_coarsen_flow(flow, mc)
  repeat {
    near <- which(flow_pub$n != n_eligible & abs(flow_pub$n - n_eligible) < mc)
    near <- setdiff(near, 1L)          # the starting population cannot be merged away
    if (length(near) == 0L) break
    r <- near[[1L]]
    to <- if (r < nrow(flow_pub)) r + 1L else r - 1L
    for (col in c("step", "reason")) {
      parts <- if (to > r) c(flow_pub[[col]][r], flow_pub[[col]][to]) else c(flow_pub[[col]][to], flow_pub[[col]][r])
      flow_pub[[col]][to] <- paste(parts[!is.na(parts) & nzchar(parts)], collapse = "; ")
    }
    flow_pub <- smidata::smi_coarsen_flow(flow_pub[-r, , drop = FALSE], mc)
  }
  flow_tbl <- smidata::smi_report_flow(flow_pub, min_cell = mc)
  flow_members <- function(value) {
    lapply(which(flow_pub$n == value), function(k) {
      member(cell("flow", "n", step_no = k), value)
    })
  }
  hidden_by_flow <- function(value) value %in% flow$n && !(value %in% flow_pub$n)

  ## Exposure-period counts. smidata's summary also republishes patients, not
  ## matched and MSR_DEN < 1, which the flow carries already (as steps and
  ## drops); publishing them twice is how a hidden step comes back, so only the
  ## split of the eligible patients is kept. Its by-month and by-INDEX_DT-day
  ## tables are smidata's own.
  period <- smidata::smi_report_measure_period(sel, p02$index_dt, min_cell = mc)
  period$summary <- mask_counts_primary(
    tibble::tibble(
      quantity = c("eligible", "eligible_not_achieved", "eligible_achieved"),
      n = as.integer(c(n_eligible, n_not_achieved, n_achieved))
    ), "n", min_cell = mc
  )
  attr(period, "smi_identities") <- Filter(
    function(i) !(i$label %in% c("measure_period partition", "patients by INDEX_DT day")),
    attr(period, "smi_identities")
  )
  mp <- function(k) cell("measure_period.summary", "n", quantity = k)
  ids[[length(ids) + 1L]] <- smidata::smi_identity(
    mp("eligible"), list(mp("eligible_not_achieved"), mp("eligible_achieved")),
    label = "eligible = achieved + not achieved"
  )
  cluster("eligible: flow vs measure_period",
          c(flow_members(n_eligible), list(member(mp("eligible"), n_eligible))))
  cluster("patients: flow vs INDEX_DT day", list(
    if (nrow(flow_pub) > 0L) member(cell("flow", "n", step_no = 1L), flow_pub$n[1]),
    member(all_rows("measure_period.index_day", "n_patients"), nrow(sel))
  ))
  for (k in seq_len(nrow(flow_pub))) {
    if (flow_pub$n[k] > n_eligible) {
      gap(cell("flow", "n", step_no = k), mp("eligible"), paste0("flow step ", k, " - eligible"))
    } else if (flow_pub$n[k] < n_eligible) {
      gap(mp("eligible"), cell("flow", "n", step_no = k), paste0("eligible - flow step ", k))
    }
  }

  ## The two windows overlap (not a partition), so primary suppression. The
  ## SAA-matched count they are measured against (`n_matched`) is the flow's
  ## exact-start step and is published there once, not repeated here.
  cov_counts <- p02$coverage_counts
  coverage_by_window <- mask_counts_primary(
    dplyr::transmute(cov_counts, window = .data$window, gating = .data$gating,
                     n_not_full = as.integer(.data$n_not_full)),
    "n_not_full", min_cell = mc
  )

  msr_filter <- mask_counts_primary(
    tibble::tibble(
      quantity = c("msr_saa_rows_not_SAA_excluded", "msr_saa_rows_for_ids_outside_cohort"),
      n = as.integer(c(p02$n_msr_rows_excluded, p02$n_saa_rows_out_of_cohort))
    ), "n", min_cell = mc
  )

  ## 4. costs. cost_final is the final sample's cost cells.
  final_ids <- p06$final_ids
  cost_final <- p03$cost[p03$cost$id %in% final_ids, , drop = FALSE]
  cost_n <- function(cost) {
    dplyr::summarise(cost, n_cells = dplyr::n(), n_patients = dplyr::n_distinct(.data$id),
                     .by = "window")
  }
  cn <- cost_n(p03$cost); fn <- cost_n(cost_final)
  costs_cohort <- smidata::smi_report_costs(p03$cost, min_cell = mc)
  costs_final_sample <- smidata::smi_report_costs(cost_final, min_cell = mc)
  cost_categories <- smidata::smi_report_cost_categories(p03$cost, min_cell = mc)
  cost_ids <- smidata::smi_report_cost_identities(
    p03$cost, costs = "costs_cohort", categories = "cost_categories"
  )
  ## smidata declares the "excluded lines" identity for every window, but a window
  ## with no Managed Care Invoice rows has no excluded category row, and a cell-set
  ## over no rows aborts. Drop only those identities (their count is 0 on both
  ## sides, so nothing is lost).
  cost_ids <- Filter(function(idn) {
    if (!inherits(idn, "smi_identity")) return(TRUE)
    for (r in idn$rhs) {
      if (r$all && identical(r$table, "cost_categories")) {
        ct <- cost_categories
        keep <- rep(TRUE, nrow(ct))
        for (k in names(r$key)) keep <- keep & !is.na(ct[[k]]) & ct[[k]] == r$key[[k]]
        if (!any(keep)) return(FALSE)
      }
    }
    TRUE
  }, cost_ids)
  ## Cohort-minus-final is the number of patients (and claim lines, ...) the
  ## complete-case restriction removed: a count in its own right, so it must not
  ## be the difference of two published ones.
  gap_cols <- c("n_cells", "n_patients", "n_claims", "n_negative_lines",
                "n_cells_negative_total", "n_cells_zero_total", "n_na_amount_lines",
                "n_excluded_mci_lines")
  for (w in intersect(costs_cohort$window, costs_final_sample$window)) {
    for (cl in gap_cols) {
      gap(cell("costs_cohort", cl, window = w), cell("costs_final_sample", cl, window = w),
          paste0("cohort - final sample: ", cl, ", ", w))
    }
  }

  cost_known <- mask_counts_primary(
    dplyr::mutate(p03$known, n_cohort = as.integer(.data$n_cohort),
                  n_known = as.integer(.data$n_known)),
    c("n_cohort", "n_known"), min_cell = mc
  )
  for (i in seq_len(nrow(p03$known))) {
    w <- p03$known$window[i]
    gap(cell("cost_known", "n_cohort", window = w), cell("cost_known", "n_known", window = w),
        paste0("cohort - known cost: ", w))
    ## n_known is the cohort count again when nobody lacks a cost: it is then
    ## hidden with it if the flow does not publish that count.
    nk <- p03$known$n_known[i]
    if (hidden_by_flow(nk)) {
      cluster(paste0("known cost count is hidden by the flow, ", w), value_hidden = TRUE,
              list(member(cell("cost_known", "n_known", window = w), nk)))
    } else if (w %in% cn$window) {
      cluster(paste0("patients with a known cost, ", w), optional = TRUE, list(
        member(cell("costs_cohort", "n_patients", window = w), cn$n_patients[cn$window == w]),
        member(cell("cost_known", "n_known", window = w), p03$known$n_known[i])
      ))
    }
  }

  ## 5. enrollment
  tally <- attr(p05$coverage, "tally")
  if (is.null(p05)) {
    enrollment <- tibble::tibble(
      note = "No 05 intermediate: complete-case status was not computed in this run (Tier-0 fixture limitation or 05 not run)."
    )
    enrollment_by_arm <- enrollment
  } else {
    enrollment <- smidata::smi_report_enrollment(p05$coverage, p05$window_months,
                                                 min_cell = mc)
    ## smidata's own masks pool a complement across the column (the smallest other
    ## positive cell), which here is a count that is also published in the flow,
    ## the cost tables and the arm table. Primary suppression only; the declared
    ## identities choose the complement (and the flow decides which copies of the
    ## cohort count exist at all).
    enrollment$summary <- mask_counts_primary(
      tibble::tibble(quantity = names(tally), n = as.integer(unlist(tally))),
      "n", min_cell = mc
    )
    ## Split by reason, as smidata's builder does: an absent patient sits in the
    ## n_covered = window_months bucket, so undivided buckets would publish the
    ## absent count as bucket - covered. Identities are per reason (declared by
    ## smi_report_enrollment()); primary suppression only, for the same reason as
    ## the summary above.
    enrollment$by_months_covered <- mask_counts_primary(
      tibble::as_tibble(p05$coverage[c("reason", "n_covered")]) |>
        dplyr::summarise(n_patients = dplyr::n(), .by = c("reason", "n_covered")) |>
        dplyr::arrange(.data$reason, .data$n_covered) |>
        dplyr::mutate(window_months = p05$window_months,
                      pct = round(100 * .data$n_patients / sum(.data$n_patients), 2)),
      "n_patients", "pct", label_cols = c("n_covered", "window_months"), min_cell = mc
    )
    ## One arm-level partition only. The by-lag-by-arm cells are NOT published:
    ## they are cells of the same cross-tab as the arm table.
    enrollment_by_arm <- p05$by_arm |>
      dplyr::mutate(
        n_complete = as.integer(.data$n_complete),
        n_not_complete = as.integer(.data$n - .data$n_complete - .data$n_unresolved),
        n_unresolved = as.integer(.data$n_unresolved)
      ) |>
      dplyr::select("A", "n_complete", "n_not_complete", "n_unresolved", "rate") |>
      mask_counts_partition(c("n_complete", "n_not_complete", "n_unresolved"), "rate",
                            label_cols = "A", min_cell = mc)
  }
  by_arm <- function(col) all_rows("enrollment_by_arm", col)
  arm_total <- if (!is.null(p05)) {
    member(list(by_arm("n_complete"), by_arm("n_not_complete"), by_arm("n_unresolved")),
           sum(p05$by_arm$n))
  }
  arm_complete <- if (!is.null(p05)) member(by_arm("n_complete"), sum(p05$by_arm$n_complete))
  es <- function(k) cell("enrollment.summary", "n", quantity = k)

  ## The analysis cohort (after the coverage filter) is published in the flow, the
  ## cost tables, the coverage tally and the arm table.
  cluster("analysis cohort size", value_hidden = hidden_by_flow(n_cohort_02), c(
    flow_members(n_cohort_02),
    lapply(seq_len(nrow(p03$known)), function(i) {
      member(cell("cost_known", "n_cohort", window = p03$known$window[i]), p03$known$n_cohort[i])
    }),
    unlist(lapply(seq_len(nrow(cn)), function(i) list(
      member(cell("costs_cohort", "n_patients", window = cn$window[i]), cn$n_patients[i]),
      member(cell("costs_cohort", "n_cells", window = cn$window[i]), cn$n_cells[i])
    )), recursive = FALSE),
    if (!is.null(p05)) list(member(es("n_cohort"), tally$n_cohort)),
    list(arm_total)
  ))
  ## Complete-case count (covered): flow, tally, arm table, months-covered cells.
  cluster("complete-case (covered) count", value_hidden = hidden_by_flow(p06$n_complete_case), c(
    flow_members(p06$n_complete_case),
    if (!is.null(p05)) list(member(es("n_covered"), tally$n_covered)),
    list(arm_complete)
  ))
  ## smi_report_enrollment() already declares covered = sum of the covered-reason
  ## months-covered cells, so no second link to that table is needed here.

  ## 6. project-specific: prior_cost floor, bins, cut points, design degeneracy
  d <- p06$design
  prior_cost_floor <- mask_counts_primary(
    tibble::tibble(quantity = d$floor$quantity, n = as.integer(d$floor$n)), "n",
    min_cell = mc
  )
  pcf <- function(k) cell("prior_cost_floor", "n", quantity = k)
  ## Patients exactly at the floor before flooring, and those not floored, are
  ## differences of published counts.
  gap(pcf("n_prior_cost_at_floor_after_flooring"), pcf("n_negative_prior_cost_floored"),
      "at the floor - floored")
  gap(pcf("n_final"), pcf("n_prior_cost_at_floor_after_flooring"), "final - at the floor")
  gap(pcf("n_final"), pcf("n_negative_prior_cost_floored"), "final - floored")

  prior_cost_bins <- d$bins |>
    dplyr::mutate(pct = round(100 * .data$n / sum(.data$n), 2)) |>
    mask_counts_partition("n", "pct", min_cell = mc)
  cuts <- d$cutpoints
  prior_cost_cutpoints <- tag_report_table(
    tibble::tibble(
      percentile = names(cuts),
      cutpoint_usd = if (d$n >= mc) round(unname(cuts), 2) else rep(NA_real_, length(cuts))
    ),
    n_cols = character(), derived_cols = "cutpoint_usd"
  )
  ## The three quartile dummies are the Q2-Q4 counts of prior_cost_bins; they are
  ## published there once (masked jointly with Q1), not again here.
  design_prevalence <- d$prevalence |>
    dplyr::filter(.data$block != "prior_cost quartile") |>
    dplyr::mutate(
      prevalence = round(.data$n_yes / (.data$n_yes + .data$n_no), 4),
      minority_share = pmin(.data$prevalence, 1 - .data$prevalence)
    ) |>
    mask_counts_primary(c("n_yes", "n_no"), c("prevalence", "minority_share"),
                        min_cell = mc)
  ## n_yes + n_no is the final sample size: a small split would return its other
  ## side, so each row is declared against n_final.
  for (cv in d$prevalence$covariate[d$prevalence$block != "prior_cost quartile"]) {
    row <- d$prevalence[d$prevalence$covariate == cv, ]
    cluster(paste0("n_yes + n_no = final sample, ", cv), list(
      member(pcf("n_final"), d$n),
      member(list(cell("design_prevalence", "n_yes", covariate = cv),
                  cell("design_prevalence", "n_no", covariate = cv)),
             row$n_yes + row$n_no)
    ))
  }
  design_grid <- mask_counts_primary(
    tibble::tibble(quantity = d$grid$quantity, n = d$grid$n), "n", min_cell = mc
  )
  final_arms <- mask_counts_partition(d$arms, "n", label_cols = character(), min_cell = mc)

  cluster("final sample size", value_hidden = hidden_by_flow(p06$n_final), c(
    flow_members(p06$n_final),
    list(member(pcf("n_final"), d$n)),
    unlist(lapply(seq_len(nrow(fn)), function(i) list(
      member(cell("costs_final_sample", "n_patients", window = fn$window[i]), fn$n_patients[i]),
      member(cell("costs_final_sample", "n_cells", window = fn$window[i]), fn$n_cells[i])
    )), recursive = FALSE),
    list(member(all_rows("final_arms", "n"), sum(d$arms$n)),
         member(all_rows("prior_cost_bins", "n"), sum(d$bins$n)))
  ))
  if (!is.null(p05)) {
    ## No drop after the complete case: the arm sizes are the complete-case counts.
    arm_val <- c(sum(p05$by_arm$n_complete[p05$by_arm$A == 0L]),
                 sum(p05$by_arm$n_complete[p05$by_arm$A == 1L]))
    for (i in 1:2) {
      cluster(paste0("final arm size = complete-case count, A=", i - 1L), optional = TRUE, list(
        member(cell("enrollment_by_arm", "n_complete", A = i - 1L), arm_val[i]),
        member(cell("final_arms", "n", arm = d$arms$arm[i]), d$arms$n[i])
      ))
    }
  }

  tabs <- list(
    flow = flow_tbl,
    measure_period = period,
    coverage_by_window = coverage_by_window,
    msr_filter = msr_filter,
    costs_cohort = costs_cohort,
    costs_final_sample = costs_final_sample,
    cost_categories = cost_categories,
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
  for (cl in to_blank) tabs <- blank_cell(tabs, cl)

  do.call(smidata::smi_report, c(
    list(title, levers = report_levers(p06$.meta)), tabs,
    list(identities = c(ids, cost_ids), min_cell = mc)
  ))
}

#' Blank one published cell (and the derived columns of its row).
#'
#' Used when the flow does not publish a count: every other copy of it is removed
#' rather than linked, so the hidden drop is not the difference of two cells.
#' `cl` is a `smidata::smi_cell()`; its table may be nested (`"enrollment.summary"`).
blank_cell <- function(tabs, cl) {
  path <- strsplit(cl$table, ".", fixed = TRUE)[[1L]]
  tb <- tryCatch(tabs[[path]], error = function(e) NULL)
  if (!is.data.frame(tb) || !cl$col %in% names(tb)) return(tabs)
  keep <- rep(TRUE, nrow(tb))
  for (k in names(cl$key)) {
    keep <- keep & !is.na(tb[[k]]) & as.character(tb[[k]]) == as.character(cl$key[[k]])
  }
  tb[[cl$col]][keep] <- NA
  from <- attr(tb, "smi_derived_from")
  for (d in attr(tb, "smi_derived_cols")) {
    src <- from[[d]]
    if (d %in% names(tb) && (is.null(src) || cl$col %in% src)) tb[[d]][keep] <- NA
  }
  tabs[[path]] <- tb
  tabs
}

## ---- toy intermediates (hand-calculated; see tests/test-decision-report.R) ---

#' Toy `coverage` as smi_enrollment_coverage() shapes it (`n_covered`, `reason`).
#'
#' Covered patients have every month (`window_months`), insufficient ones
#' `insufficient_months`, and absent ones `window_months` too (they have no
#' uncovered months to count), with a `tally` attribute.
toy_coverage <- function(n_covered, n_insufficient, n_absent, window_months = 24L,
                         insufficient_months = 20L) {
  cov <- data.frame(
    n_covered = c(rep(window_months, n_covered), rep(insufficient_months, n_insufficient),
                  rep(window_months, n_absent)),
    reason = c(rep("covered", n_covered), rep("insufficient_coverage", n_insufficient),
               rep("absent_from_flag_table", n_absent))
  )
  attr(cov, "tally") <- list(
    n_cohort = n_covered + n_insufficient + n_absent, n_covered = n_covered,
    n_insufficient = n_insufficient, n_absent = n_absent, n_uncovered_at_index = 0L
  )
  cov
}

#' Toy report intermediates with hand-calculable values.
#'
#' Same shape 02/03/05/06 record. Not data, not a stand-in for data: it exists so
#' 08 runs end to end off-server and so the test can pin counts by arithmetic.
#' The sections are independent of each other except where the test says so.
#'
#' @param extra_outcome_claims Optional extra rows (`ID`, `SRV_DT`, `AMOUNT_PAID`,
#'   `SS_DESC`) appended to the outcome-window claims, so a test can add a service
#'   category of a chosen size.
#' @return `list(p02, p03, p05, p06)`.
toy_report_inputs <- function(extra_outcome_claims = NULL) {
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
                   AMOUNT_PAID = 999, SS_DESC = "Managed Care Invoice"),
    extra_outcome_claims
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
  ## An absent patient has n_covered = window_months (smi_enrollment_coverage()).
  cov <- toy_coverage(100L, 30L, 6L)
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
