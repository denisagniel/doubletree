## ============================================================================
## 03_cost_windows.R
##
## Build the 36-file cost_sources declaration, the per-patient windows, the
## coverage cross-check, and the per-patient cost sums. REAL as of 2026-09-24:
## the previously-commented "intended pipeline" is now live code, cascade-
## guarded on 02's analysis_cohort (require_stage()) -- runs end to end
## against a Tier-0 fixture locally, or the real 36 files on the server.
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
## Run time: ~10 seconds locally against 36 Tier-0 fixtures (verified by
## `oracle`, 2026-09-24, before this was wired live). On the server: the
## cost-claims family is 171.5 GB across 36 files -- 91.8% of every byte in
## the secure tree (smidata fingerprint 2026-09-15).
## smidata::smi_sum_cost_in_windows() streams one file at a time for exactly
## that reason. Budget HOURS, not minutes, and expect the per-file cost to
## scale with file size (aim3_svc_cost_20 alone is 23,117,933 rows). smidata's
## own docs flag this as "real streaming code, not a stub, but exercised only
## through in-memory fake readers... treat the first server run as a Tier-1
## check, not a re-run" -- true of the Tier-0 fixture run here too.
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
## WHETHER ALL FOUR FAMILIES SHOULD BE SUMMED was the SAME unresolved question
## dual-bounds had -- same files, same consequence -- and is now CONFIRMED
## (PI, 2026-09-17; see OPEN_DECISIONS$cost_family_scope): the four families
## are disjoint patient subsets (aim3_* has 12mo post-index enrollment,
## new_* does not), so summing all four does not double-count. STALE GATE
## FIXED 2026-09-24 (same bug class as msr_code's, found the same day): this
## script called `open_value("cost_family_scope")`, which always warns
## "Unconfirmed value" regardless of actual status -- misleading for a
## decision that has, in fact, been confirmed since 2026-09-17. Not yet
## independently verified against real claim-level data; the
## aim3_svc_cost_20-vs-new_svc_cost_20 comparison in dual-bounds'
## 90_checks_tier1.R CHECK 1 remains worth running as corroboration.

cost_family_scope <- confirmed_value("cost_family_scope")

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
  cli::cli_alert_info(
    "Summing all four families -- PI-confirmed (2026-09-17) they are
     disjoint patient subsets, not overlapping extracts. Not yet
     independently verified against real claim-level data; see
     OPEN_DECISIONS$cost_family_scope's note for the corroborating check
     still worth running."
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

## ---- windows/coverage/cost sums: REAL as of 2026-09-24 ----------------------
##
## Design-reviewed by `oracle` before writing (session
## ses_f2b12eccdffel35nxc8t7UppvR): actually RAN this against a Tier-0 fixture
## first (10.4s, nonzero sums) to find failure modes before trusting the code,
## not after. Two real gaps that finding surfaced, fixed here:
##
##   (a) `read_fn = smidata::smi_read` with no column restriction reads every
##       column of every file (20-21 cols incl. NDC) -- invisible on a 200-row
##       fixture, real cost on aim3_svc_cost_20's 23,117,933 rows on the
##       server. Restricted via a wrapping closure below.
##   (b) MCI exclusion (n_excluded_mci) is UNVERIFIABLE at Tier 0: the census
##       carries no real SS_DESC levels, so the fixture's SS_DESC values never
##       match "Managed Care Invoice" and n_excluded_mci is always 0 --
##       passing cleanly regardless of whether the real regex/column is
##       right. Checked only when NOT local, below; 90_checks_tier1.R remains
##       the real corroboration either way.
##
## TRUNCATION (file-date coverage) is resolved into `coverage`, joined onto
## BOTH `outcome_cost` and `prior_cost` as an explicit column, rather than
## left for 06 to rediscover: a patient absent from `smi_sum_cost_in_windows()`'s
## output because they have coverage == "full" and genuinely zero claims is
## Y = 0 (or prior_cost = 0) -- a real, informative value. A patient absent
## because their window falls outside the 2016-2024 file coverage has an
## UNKNOWABLE cost, and coalescing that to 0 would silently look like the
## former. `06` must branch on the carried `coverage` column, not blanket
## `coalesce(..., 0)` -- verified there is a real difference: at Tier-0
## fixture scale, 17/200 cohort patients were absent from `outcome_cost`
## despite full file coverage (genuine zero-utilization), and separately
## 26/200 outcome-windows were only PARTIALLY covered by file years (the
## truncation case) -- two different populations, both real at fixture scale.

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; skipping the real windows/coverage/cost
     pipeline entirely (same convention as 01/02)."
  )
} else {
  require_stage("analysis_cohort", "02_population_and_eligibility.R")
  is_local <- identical(smidata::smi_env(), "local")

  patients <- analysis_cohort |>
    dplyr::transmute(id = .data$ID, index_dt = as.Date(.data$INDEX_DT))

  windows <- smidata::smi_patient_windows(patients, windows = DOUBLETREE_WINDOWS)

  coverage <- smidata::smi_window_coverage(windows, cost_sources) |>
    dplyr::select("id", "window", "coverage")
  n_incoverage <- sum(coverage$coverage != "full")
  if (n_incoverage > 0L) {
    ## cli_alert_warning(), not cli_warn(): cli_warn() raises an R warning
    ## CONDITION, which Rscript defers and buckets ("There were 50 or more
    ## warnings...") -- this notice got silently buried under
    ## smi_fixture()'s own ~50 "no numeric summary" warnings the first time
    ## this ran, found while testing 2026-09-24. cli_alert_warning() prints
    ## immediately, matching every other notice in this file.
    cli::cli_alert_warning(
      "{n_incoverage} of {nrow(coverage)} patient-windows extend beyond the
       2016-2024 file coverage; their cost sums would be TRUNCATED, which
       looks like low cost, not missing cost. Carried forward as an explicit
       {.field coverage} column on {.field outcome_cost}/{.field prior_cost}
       below -- 06 must branch on it, not coalesce blindly."
    )
  }

  ## Column-restricted read_fn: only what smi_sum_cost_in_windows() actually
  ## uses (id_col, date_col, amount_col, setting_col), not every column of
  ## every file. cost_sources' columns are the confirmed contract (all 36
  ## files share ID/SRV_DT/AMOUNT_PAID/SEQ_ID/SS_DESC; SEQ_ID is read too,
  ## even though this function does not use it, because it is this table's
  ## claim-key column per SMI_COLS$cost_claims and dropping it here would
  ## silently diverge this read from every other read of the same files).
  cost_read_fn <- function(key) {
    smi_read_pinned(key, columns = SMI_COLS$cost_claims)
  }

  ## The expensive call on the server (budget hours, per this file's header);
  ## ~10s locally against 36 Tier-0 fixtures (verified by oracle before this
  ## was written).
  cost_by_window <- smidata::smi_sum_cost_in_windows(
    patient_windows = windows,
    cost_sources = cost_sources,
    exclude_setting_regex = "Managed Care Invoice",  # -> smi_classify_is_mci()
    detail = "total",
    read_fn = cost_read_fn,
    id_col = "ID"
  )

  if (!is_local) {
    ## Tier-0 fixtures cannot exercise this (see (b) above) -- checked only
    ## against real data, where a total of zero MCI exclusions across 36
    ## files and ~24-month windows is implausible enough to be a
    ## misconfigured config_cost_setting_col or regex, not a real finding.
    n_excluded_mci_total <- sum(cost_by_window$n_excluded_mci, na.rm = TRUE)
    if (n_excluded_mci_total == 0L) {
      cli::cli_alert_warning(
        "Zero Managed Care Invoice records excluded across all cost windows.
         Plausible, but check {.field config_cost_setting_col}
         ({.val {config_cost_setting_col}}) and the exclude_setting_regex
         before trusting this -- see 90_checks_tier1.R's calibrating check."
      )
    }
  }

  ## Y (outcome) and prior_cost are two slices of the same object. Renamed to
  ## ID (matching every other table's join key) and LEFT-JOINED onto the full
  ## per-patient `coverage` rows (not the other way around), so every cohort
  ## patient gets a row here -- absent-with-full-coverage becomes a real 0;
  ## absent-with-partial/none coverage stays NA, on purpose (see the
  ## TRUNCATION note above -- 06 relies on this column, not a re-derivation).
  build_window_cost <- function(window_name) {
    coverage |>
      dplyr::filter(.data$window == window_name) |>
      dplyr::left_join(
        dplyr::filter(cost_by_window, .data$window == window_name),
        by = c("id", "window")
      ) |>
      dplyr::mutate(
        ## coverage != "full" means NA regardless of whether cost_by_window
        ## happened to return a non-NA sum: a "partial" patient CAN still
        ## have a claim inside the covered portion, and that partial sum is
        ## NOT the intended full-window total -- passing it through as if it
        ## were complete would silently misrepresent a truncated observation
        ## as a real one. Only "full" coverage gets its absence-is-zero
        ## treatment; anything else is unknowable, full stop.
        cost = dplyr::if_else(
          .data$coverage == "full", dplyr::coalesce(.data$cost, 0), NA_real_
        )
      ) |>
      dplyr::transmute(ID = .data$id, cost = .data$cost, coverage = .data$coverage)
  }
  outcome_cost <- build_window_cost("outcome")
  prior_cost   <- build_window_cost("prior_cost")

  cli::cli_alert_info(
    "outcome_cost: {sum(!is.na(outcome_cost$cost))} of {nrow(outcome_cost)}
     patients have a known Y (rest: coverage != 'full', UNKNOWABLE)."
  )
  cli::cli_alert_info(
    "prior_cost: {sum(!is.na(prior_cost$cost))} of {nrow(prior_cost)}
     patients have a known prior_cost."
  )
  if (is_local) {
    cli::cli_alert_warning(
      "The counts above are Tier-0 FIXTURE counts -- proof of execution, not
       a real cost distribution."
    )
  }
}

## ---- smoke-check: the pure functions, on two hand-written patients ---------
##
## Confirms THIS script's own cost_sources declaration and _config.R's
## DOUBLETREE_WINDOWS compose with smidata's promoted functions. Needs no data,
## but DOES need smidata installed -- the window and coverage logic lives there
## (>= 0.2.0), so the check is skipped, loudly, when it is absent rather than
## reimplemented here. Separate from the real pipeline above (which now uses
## the real analysis_cohort, not two hand-written patients) -- kept because it
## deliberately constructs a LATE-indexed patient to demonstrate the coverage
## check catches something even when the real cohort's own coverage happens
## not to need showing.

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
      "smi_sum_cost_in_windows() NOT called on these two hand-written patients
       (it already ran above, for real, on analysis_cohort -- this
       smoke-check is only demonstrating smi_patient_windows()/
       smi_window_coverage() composing correctly, not re-running the
       expensive call)."
    )
    cli::cli_alert_info("03_cost_windows.R complete.")
  }
}
