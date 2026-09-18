## ============================================================================
## 03_cost_windows.R
##
## Build the 36-file cost_sources declaration, the per-patient windows, and the
## coverage cross-check. Two of the three steps run for real on toy patients;
## the third needs the files.
##
## The OUTCOME window and the PRIOR_COST covariate window are two slices of ONE
## pass over the cost files. That is not an optimization -- it is what makes them
## the same measurement by construction, with the same MCI exclusion and the same
## total_cost definition, differing only in which months are summed.
##
## Design (inherited from dual-bounds, smidata
## inst/analyses/doubletree__application.yml sections 5-6, `status: inherited`):
##   outcome.primary               = total_cost
##   outcome.ascertainment_window  = [INDEX_DT + 12mo, INDEX_DT + 24mo)
##   covariates.measurement_window = [INDEX_DT - 12mo, INDEX_DT)   -> prior_cost
##   outcome.exclusions            = Managed Care Invoice records
##
## THIS PAPER'S ONE ADDITION: prior_cost does not enter X as dollars. It is
## discretized into quartile dummies in 04_covariates.R, because
## doubletree::estimate_att() requires every covariate to be individually binary.
## The dollars are still computed here -- discretization needs them.
##
## Run time when unblocked: the cost-claims family is 171.5 GB across 36 files
## -- 91.8% of every byte in the secure tree (smidata fingerprint 2026-09-15).
## smidata::smi_sum_cost_in_windows() streams one file at a time for exactly that
## reason. Budget HOURS, not minutes, and expect the per-file cost to scale with
## file size (aim3_svc_cost_20 alone is 23,117,933 rows).
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("03 -- cost windows")

## ---- cost_sources: 4 families x 9 years = 36 files -------------------------
##
## Column names are verbatim from the real fingerprint (all 36 files share
## ID / SRV_DT / AMOUNT_PAID / SEQ_ID / SS_DESC; the pharm files additionally
## carry NDC, which this pipeline does not use).
##
## WHETHER ALL FOUR FAMILIES SHOULD BE SUMMED IS UNRESOLVED, and it is the SAME
## unresolved question dual-bounds has -- same files, same consequence. See
## OPEN_DECISIONS$cost_family_scope (blocking, shared) and its registry_ref. The
## resolving check is specified once, in dual-bounds'
## analysis/real_data/90_checks_tier1.R CHECK 1; it is referenced, not re-derived,
## from this project's 90_checks_tier1.R.

cost_family_scope <- open_value("cost_family_scope")

cost_sources <- tidyr::expand_grid(
  family = unlist(SMI_KEYS$cost_claims, use.names = FALSE),
  suffix_year = config_cost_years
) |>
  dplyr::mutate(
    dataset_key = paste0(.data$family, "_", .data$suffix_year),
    amount_col = "AMOUNT_PAID",
    date_col = "SRV_DT",
    setting_col = config_cost_setting_col,
    claim_key_col = "SEQ_ID"
  ) |>
  dplyr::select("dataset_key", "family", "suffix_year",
                "amount_col", "date_col", "setting_col", "claim_key_col")

if (nrow(cost_sources) != 36L) {
  cli::cli_abort("Expected 36 cost files (4 families x 9 years); built {nrow(cost_sources)}.")
}
cli::cli_alert_info(
  "cost_sources: {nrow(cost_sources)} files across
   {dplyr::n_distinct(cost_sources$family)} families, suffix years
   {min(cost_sources$suffix_year)}-{max(cost_sources$suffix_year)}."
)
if (identical(cost_family_scope, "all_four")) {
  cli::cli_alert_warning(
    "Summing ALL FOUR families on an UNCONFIRMED default. If aim3_* and new_*
     overlap, both Y AND prior_cost are inflated -- and prior_cost's inflation is
     worse than a scale error, because the QUARTILE BOUNDARIES are computed from
     the inflated values, so the covariate's grid changes too."
  )
}

## ---- windows ---------------------------------------------------------------
##
## The window table itself lives in _config.R as DOUBLETREE_WINDOWS, not here:
## smidata::smi_patient_windows() deliberately has no default `windows`
## argument, so the offsets are stated exactly once for this paper and passed
## explicitly at every call site. Half-open, start inclusive, end exclusive;
## prior_cost's end IS INDEX_DT, so the index day itself belongs to neither
## window and no claim is counted twice.
##
## The VALUES equal dual-bounds' DUAL_BOUNDS_WINDOWS, because the windows are
## inherited. The CONSTANT is still local -- see _config.R's header for why an
## inherited design decision must not become a cross-project code dependency.

## ---- INTENDED PIPELINE (needs real data) -----------------------------------
##
## patients comes from 02_population_and_eligibility.R's analysis_cohort,
## renamed to the (id, index_dt) contract smi_patient_windows() expects:
##
## patients <- analysis_cohort |>
##   dplyr::transmute(id = .data$ID, index_dt = as.Date(.data$INDEX_DT))
##
## windows <- smidata::smi_patient_windows(patients, windows = DOUBLETREE_WINDOWS)
##
## coverage <- smidata::smi_window_coverage(windows, cost_sources)
## if (any(coverage$coverage != "full")) {
##   cli::cli_warn(
##     "{sum(coverage$coverage != 'full')} patient-windows extend beyond the
##      2016-2024 file coverage; their cost sums are TRUNCATED, which looks like
##      low cost. Decide explicitly whether to exclude them."
##   )
## }
##
## ## The expensive call. Streams 36 files; see the run-time note above.
## cost_by_window <- smidata::smi_sum_cost_in_windows(
##   patient_windows = windows,
##   cost_sources = cost_sources,
##   exclude_setting_regex = "Managed Care Invoice",  # -> smi_classify_is_mci()
##   detail = "total",
##   read_fn = smidata::smi_read,
##   id_col = "ID"
## )
##
## ## Y (outcome) and prior_cost are two slices of the same object.
## outcome_cost <- dplyr::filter(cost_by_window, .data$window == "outcome")
## prior_cost   <- dplyr::filter(cost_by_window, .data$window == "prior_cost")
##
## ## TRUNCATION IS A COVARIATE PROBLEM HERE, NOT ONLY AN OUTCOME PROBLEM.
## ## A patient whose PRE-index window starts before 2016 has a truncated
## ## prior_cost that is indistinguishable from a genuinely low one -- and because
## ## the quartile cutpoints in 04 are computed WITHIN the analytic sample, a
## ## cluster of truncated-low values shifts the boundaries for everyone. Join
## ## `coverage` in and decide explicitly; do not let truncation pass as data.

## ---- smoke-check: the pure functions, on two hand-written patients ---------
##
## Confirms THIS script's own cost_sources declaration and _config.R's
## DOUBLETREE_WINDOWS compose with smidata's promoted functions. Needs no data,
## but DOES need smidata installed -- the window and coverage logic lives there
## (>= 0.2.0), so the check is skipped, loudly, when it is absent rather than
## reimplemented here.

if (!interactive() && sys.nframe() == 0L) {
  if (!config_has_smidata) {
    cli::cli_alert_warning(
      "Smoke-check SKIPPED: {.pkg smidata} (>= 0.2.0) is not installed, and
       {.fn smi_patient_windows} / {.fn smi_window_coverage} live there."
    )
    cli::cli_alert_info("03_cost_windows.R complete (smoke-check skipped).")
  } else {
    cli::cli_h2("Smoke-check on {config_smoke_n_patients} hand-written patients")

    smoke_patients <- tibble::tibble(
      id = paste0("smoke_", seq_len(config_smoke_n_patients)),
      ## One late-indexed patient deliberately runs past the 2024 file coverage,
      ## so the coverage check has something to catch.
      index_dt = as.Date(c("2018-06-15", "2023-06-15"))[seq_len(config_smoke_n_patients)]
    )

    smoke_windows <- smidata::smi_patient_windows(
      smoke_patients, windows = DOUBLETREE_WINDOWS
    )
    print(smoke_windows)

    smoke_coverage <- smidata::smi_window_coverage(smoke_windows, cost_sources)
    print(dplyr::select(smoke_coverage, "id", "window", "window_start", "window_end", "coverage"))

    n_incomplete <- sum(smoke_coverage$coverage != "full")
    cli::cli_alert_info(
      "Smoke-check: {n_incomplete} of {nrow(smoke_coverage)} patient-windows are not
       fully covered by the 2016-2024 files."
    )
    cli::cli_alert_info(
      "smi_sum_cost_in_windows() NOT called -- it needs a real read_fn and the 36
       files. See the commented block above."
    )
    cli::cli_alert_info("03_cost_windows.R complete.")
  }
}
