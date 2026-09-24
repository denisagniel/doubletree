## ============================================================================
## application/run_pipeline.R
##
## Convenience only, not plumbing: sources 01-07 in order, in one session, so
## the whole chain runs with a single command. Every numbered script is
## independently runnable standalone (each cascades to its own upstream
## dependencies via `require_stage()` in `_config.R` -- see e.g. 06's own
## `require_stage()` calls) -- this file exists only to save typing seven
## commands, not because any script actually needs it.
##
##   Rscript application/run_pipeline.R
##
## Stops wherever the real chain stops. Locally (Tier 0), that is inside
## 05_complete_case.R's own YEAR_MONTH limitation (larger_smi_medicaid_
## monthly_flag has no real-data census levels for that column, so the
## fixture cannot exercise compute_complete_case() -- see that file's
## header): 06 and 07 will correctly abort just past that point, which is an
## honest Tier-0 boundary, not a bug in this script. On the server, the full
## chain runs end to end.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."

for (script in c(
  "01_declare_requirements.R",
  "02_population_and_eligibility.R",
  "03_cost_windows.R",
  "04_covariates.R",
  "05_complete_case.R",
  "06_assemble_analytic_data.R",
  "07_estimate_att.R"
)) {
  cli::cli_rule("Sourcing {.file {script}}")
  source(file.path(app_dir, script))
}
