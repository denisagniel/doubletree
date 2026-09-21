#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# enrollment_coverage.R -- evidence for the `enrollment_source` open decision
# (SHARED by dual-bounds and doubletree; see application/_config.R's
# OPEN_DECISIONS$enrollment_source and smidata's inst/analyses/*.yml).
# Characterizes, from larger_smi_medicaid_monthly_flag, how much of each
# patient's [INDEX_DT, INDEX_DT+24mo) window is covered by MEDICAID_FLAG == 1,
# and what the TYPE column actually partitions.
#
# WHERE THIS SHOULD LIVE: this is package-level infrastructure, not a
# doubletree-specific pipeline step -- it belongs alongside smidata's own
# inst/server/05_cohort_overlap.R (same role: server-side EVIDENCE for a
# cross-project open decision, not a ruling). It was authored here, in
# doubletree's application/server/, because the authoring environment could
# not write to ~/RAND/tools/smidata directly (outside this session's sandboxed
# workspace). PROMOTE this file into smidata's inst/server/ as the next
# available number (07_enrollment_coverage.R, after 06_ingest.R) the next time
# you have write access there, and delete it from here once promoted --
# resolving enrollment_source here should resolve it for both dual-bounds and
# doubletree at once, per the registry's own registry_ref convention, and that
# only works cleanly if the evidence script lives in the shared package, not
# duplicated per consumer.
#
# THIS IS EVIDENCE, NOT A RULING. It does not pick a threshold for "enough"
# months of coverage to count as complete-case -- that is a PI call. It exists
# so the PI can make that call by looking at the real distribution instead of
# guessing at one.
#
# COST WARNING: larger_smi_medicaid_monthly_flag is 5.8 GB, 24M rows. Per
# smidata's SERVER_SESSION_RUNBOOK.md, opening a .sas7bdat file costs hours
# regardless of how much is subsequently read. This is almost certainly your
# ONE hours-scale, file-touching operation for the session -- if you also
# planned to run 01_census.R (structure/distribution capture for
# larger_smi_covariates / msr_aap / the cost-claims family), that is a
# SEPARATE hours-scale commitment and the two probably do not both fit in one
# sitting. Pick one.
#
# Structured in sections so you can run it semi-interactively (paste into a
# Citrix R console) and stop after Section 2 if the format checks in Section 2
# don't match what Section 3 expects -- Section 3 aborts cleanly rather than
# silently mis-parse if that happens, but you shouldn't need to find that out
# the hard way after paying the file-open cost twice.
#
# Nothing here calls smi_export() unconditionally. Everything printed to the
# console stays on the server. Section 5 is an OPT-IN final step for anyone
# who wants a suppressed summary copied into the SharePoint outbox.
#
# Safe to paste straight into the R console in a Citrix window, or to run with
#   Rscript enrollment_coverage.R
#
# Requires smidata installed on the server (remotes::install_github(
# "denisagniel/smidata@v0.1.0") or later, per SERVER_SESSION_RUNBOOK.md step 1).
#
# CONFIGURATION: plain script variables below (edit before running), not
# environment variables -- smidata's project convention (PI, 2026-09-17).
# ---------------------------------------------------------------------------

library(smidata)
library(dplyr)
library(lubridate)

## ---- 0 gates --------------------------------------------------------------

dataset_key   <- "larger_smi_medicaid_monthly_flag"
window_months <- 24L          # the complete-case window this evidence targets
min_cell      <- 11L          # smidata's own small-cell suppression threshold
data_dir      <- smi_data_dir()

cli::cli_h1("enrollment_coverage")
cli::cli_bullets(c(
  "*" = "Environment: {.strong {smi_env()}}",
  "*" = "Data directory: {.path {data_dir}}",
  "*" = "Dataset: {.val {dataset_key}}",
  "*" = "Window: INDEX_DT to INDEX_DT + {window_months} months",
  "!" = "This will read the FULL raw file if no ingest cache exists yet.
         Expect this to take a while -- see the header comment."
))

if (!smi_is_server()) {
  cli::cli_abort(c(
    "This does not look like the secure server.",
    "i" = 'Set the data directory and re-run:
           Sys.setenv(SMI_DATA_DIR = "H:/CiTER-ER/ETS/Rand_OMH_SMI_Project/DATA")'
  ))
}

## Suppress any count below min_cell before it is printed or exported.
## Applied by hand here (not a smi_export() call) because everything in this
## script stays on the server -- smi_export()'s own suppression only fires on
## the export path itself, and this script's tables are printed to console,
## not exported, unless Section 5 is explicitly run.
suppress_small <- function(df, count_col) {
  df[[count_col]][df[[count_col]] < min_cell & df[[count_col]] > 0L] <- NA_integer_
  df
}

## ---- 1 read ----------------------------------------------------------------

started <- Sys.time()
cli::cli_h2("1. Reading {.val {dataset_key}}")
d <- smi_read(dataset_key)
cli::cli_alert_success(
  "Read {nrow(d)} rows, {ncol(d)} columns in
   {round(difftime(Sys.time(), started, units = 'mins'), 1)} minutes."
)

expected_cols <- c("ID", "INDEX_DT", "TYPE", "COHORT", "YEAR_MONTH", "MEDICAID_FLAG")
missing_cols <- setdiff(expected_cols, names(d))
if (length(missing_cols) > 0L) {
  cli::cli_abort(c(
    "{.val {dataset_key}} is missing expected column{?s}: {.field {missing_cols}}.",
    "i" = "Confirmed schema was ID, INDEX_DT, TYPE, COHORT, YEAR_MONTH,
           MEDICAID_FLAG per the 2026-09-15 fingerprint. The source may have
           changed -- do not proceed on a guess."
  ))
}

## ---- 2 structural diagnostics (cheap; do these first) ----------------------

cli::cli_h2("2. Structural diagnostics")

cli::cli_h3("2a. TYPE -- what does it partition?")
type_counts <- d |>
  count(TYPE, name = "n") |>
  suppress_small("n") |>
  arrange(desc(n))
print(type_counts, n = Inf)

cli::cli_h3("2b. COHORT -- does this table span more than the target cohort?")
cohort_counts <- d |>
  count(COHORT, name = "n") |>
  suppress_small("n") |>
  arrange(desc(n))
print(cohort_counts, n = Inf)

cli::cli_h3("2c. Does every (ID) have exactly one row per (YEAR_MONTH, TYPE)?")
dupe_check <- d |>
  count(ID, YEAR_MONTH, TYPE, name = "n_rows") |>
  count(n_rows, name = "n_id_month_type_combos")
print(dupe_check, n = Inf)
cli::cli_alert_info(
  "If n_rows is ever > 1 above, (ID, YEAR_MONTH, TYPE) is not a unique key and
   Section 4's per-patient coverage counts below would silently double-count
   without a fix -- stop and look before proceeding."
)

cli::cli_h3("2d. Is INDEX_DT constant within patient (as expected)?")
index_dt_check <- d |>
  group_by(ID) |>
  summarise(n_distinct_index_dt = n_distinct(INDEX_DT), .groups = "drop") |>
  count(n_distinct_index_dt, name = "n_patients")
print(index_dt_check, n = Inf)
if (any(index_dt_check$n_distinct_index_dt > 1L)) {
  cli::cli_alert_danger(
    "Some patients have more than one distinct INDEX_DT within this table.
     Section 4 below uses the MINIMUM per patient -- confirm that is the right
     choice before trusting the coverage numbers; it may instead mean this
     table pools multiple cohorts/papers with different index-date
     conventions (cross-check against Section 2b's COHORT breakdown)."
  )
}

cli::cli_h3("2e. YEAR_MONTH -- observed format (unconfirmed; do not assume)")
sample_values <- unique(d$YEAR_MONTH)
cli::cli_alert_info("{length(sample_values)} distinct YEAR_MONTH value(s).")
cli::cli_bullets(c("*" = "First few, sorted: {.val {utils::head(sort(sample_values), 8L)}}"))
cli::cli_bullets(c("*" = "Last few, sorted: {.val {utils::tail(sort(sample_values), 8L)}}"))

fmt_yyyymm      <- "^\\d{4}-\\d{2}$"   # e.g. "2019-06"
fmt_yyyymm_flat <- "^\\d{6}$"          # e.g. "201906"
fmt_ry          <- "^RY\\d{4}-\\d{2}$" # msr_aap's MSR_YR convention, for contrast

n_match <- c(
  "YYYY-MM"  = sum(grepl(fmt_yyyymm, sample_values)),
  "YYYYMM"   = sum(grepl(fmt_yyyymm_flat, sample_values)),
  "RY<yyyy>-<mm>" = sum(grepl(fmt_ry, sample_values))
)
cli::cli_bullets(c("*" = "Pattern match counts (out of {length(sample_values)} distinct values): "))
print(n_match)

## STOP HERE if you are running this semi-interactively and none of the three
## patterns above accounts for all distinct values. Section 3 aborts if the
## format is not confirmed, rather than guessing -- but confirming visually
## here first saves re-paying Section 1's read cost after a wasted trip
## through Section 3's abort.

## ---- 3 parse YEAR_MONTH (aborts, does not guess, on a format mismatch) ----

cli::cli_h2("3. Parsing YEAR_MONTH")

if (all(grepl(fmt_yyyymm, sample_values))) {
  cli::cli_alert_success("Format confirmed: YYYY-MM.")
  d <- d |> mutate(
    ym_year  = as.integer(substr(YEAR_MONTH, 1L, 4L)),
    ym_month = as.integer(substr(YEAR_MONTH, 6L, 7L))
  )
} else if (all(grepl(fmt_yyyymm_flat, sample_values))) {
  cli::cli_alert_success("Format confirmed: YYYYMM.")
  d <- d |> mutate(
    ym_year  = as.integer(substr(YEAR_MONTH, 1L, 4L)),
    ym_month = as.integer(substr(YEAR_MONTH, 5L, 6L))
  )
} else {
  offenders <- utils::head(sample_values[!grepl(fmt_yyyymm, sample_values) &
                                          !grepl(fmt_yyyymm_flat, sample_values)], 10L)
  cli::cli_abort(c(
    "YEAR_MONTH does not match either anticipated format (YYYY-MM or YYYYMM).",
    "x" = "Offending value{?s}: {.val {offenders}}",
    "i" = "Do not coerce this to NA and proceed -- update this script's parser
           once you know the real format, the same way smi_parse_msr_yr()
           aborts on MSR_YR mismatches rather than guessing."
  ))
}

if (any(d$ym_month < 1L | d$ym_month > 12L, na.rm = TRUE)) {
  cli::cli_abort("Parsed YEAR_MONTH month component out of range 1-12 -- the
                  format guess above was wrong despite matching the regex.")
}

## ---- 4 per-patient coverage over [INDEX_DT, INDEX_DT + window_months) -----

cli::cli_h2("4. Per-patient enrollment coverage")

patients <- d |>
  group_by(ID) |>
  summarise(index_dt = min(as.Date(INDEX_DT)), .groups = "drop") |>
  mutate(
    window_start = as.Date(index_dt),
    window_end   = lubridate::add_with_rollback(
      as.Date(index_dt), lubridate::period(months = window_months)
    )
  )

d2 <- d |>
  mutate(ym_date = as.Date(sprintf("%04d-%02d-01", ym_year, ym_month))) |>
  inner_join(patients, by = "ID") |>
  filter(ym_date >= window_start, ym_date < window_end)

## months_present: any row at all for this (ID, calendar month) in the window,
## regardless of TYPE or flag value.
## months_covered: same, but MEDICAID_FLAG == 1.
## The two are reported separately on purpose -- a month absent from this
## table is NOT the same claim as a month present with MEDICAID_FLAG == 0, and
## conflating them would silently answer the enrollment_source question
## incorrectly.
coverage <- d2 |>
  group_by(ID) |>
  summarise(
    months_present = n_distinct(ym_date),
    months_covered = n_distinct(ym_date[MEDICAID_FLAG == 1]),
    .groups = "drop"
  ) |>
  mutate(pct_covered = months_covered / window_months)

cli::cli_h3("4a. Distribution of months_covered (out of {window_months})")
coverage_dist <- coverage |>
  count(months_covered, name = "n_patients") |>
  suppress_small("n_patients") |>
  arrange(months_covered)
print(coverage_dist, n = Inf)

cli::cli_h3("4b. Quantiles of pct_covered (safe -- computed over all patients, not a small cell)")
print(quantile(coverage$pct_covered, probs = c(0, .05, .1, .25, .5, .75, .9, .95, 1), na.rm = TRUE))

cli::cli_h3("4c. months_present vs months_covered -- how much is 'absent' vs 'present but flag=0'?")
presence_vs_flag <- coverage |>
  mutate(months_absent_or_unflagged = window_months - months_covered) |>
  summarise(
    n_patients = n(),
    mean_months_present = mean(months_present),
    mean_months_covered = mean(months_covered)
  )
print(presence_vs_flag)

cli::cli_h3("4d. Coverage broken out by TYPE (does TYPE matter for this question?)")
coverage_by_type <- d2 |>
  group_by(ID, TYPE) |>
  summarise(months_covered_this_type = n_distinct(ym_date[MEDICAID_FLAG == 1]), .groups = "drop") |>
  group_by(TYPE) |>
  summarise(
    n_patients = n(),
    mean_months_covered = mean(months_covered_this_type),
    .groups = "drop"
  ) |>
  filter(n_patients >= min_cell) |>
  arrange(desc(n_patients))
print(coverage_by_type, n = Inf)

## ---- 5 OPTIONAL: export a suppressed summary off-server -------------------
## Not run automatically. Uncomment and run this section only if you want a
## copy of the (already-suppressed) summary tables landed in the SharePoint
## outbox for the PI conversation / session notes.

# smi_export(coverage_dist,    name = "enrollment_coverage_distribution", kind = "table")
# smi_export(coverage_by_type, name = "enrollment_coverage_by_type",      kind = "table")

## ---- 6 verdict --------------------------------------------------------------

elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)
cli::cli_h1("DONE")
cli::cli_bullets(c(
  "v" = "Elapsed: {elapsed} minute{?s}",
  "i" = "This is EVIDENCE for the enrollment_source open decision, not a
         resolution. The PI still needs to pick a months_covered (or
         pct_covered) threshold from Section 4a/4b's real distribution, and
         decide what TYPE means for Section 4d's breakdown.",
  "i" = "Record the answer in application/_config.R's
         OPEN_DECISIONS$enrollment_source$value AND dual-bounds' equivalent
         entry -- this question is shared; resolving it here resolves both,
         per smidata's own registry_ref convention.",
  "*" = "diagnostics (the OTHER blocking open decision for doubletree) is not
         addressed by this script at all -- it needs a live PI conversation,
         not data."
))
