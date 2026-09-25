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
##
## RUN TRACKING. smi_run_start()/smi_run_finish() bracket the WHOLE pipeline
## and are called from HERE ONLY -- never inside an individual stage script
## (01-07): a stage is independently runnable on its own (require_stage()'s
## whole point), so a run record that started fresh every time a downstream
## stage got sourced standalone would not describe "one pipeline execution"
## at all. Idempotent on the START side: if a run is ALREADY active (this
## file sourced twice in one session, or an outer caller already started
## one), a second smi_run_start() is skipped rather than repeated, because it
## would silently overwrite the run-id bookkeeping smi_run_finish() needs to
## restore SMI_RUN_ID correctly afterward. Guarded on config_has_smidata,
## matching every numbered script's own convention (_config.R's header):
## this repo has no data and, on some machines, no smidata, and this file
## must stay sourceable end to end without it, the same way 01-07 already
## degrade gracefully on their own. `ingest_id = config_ingest_id` is exactly
## the pin `smi_read_pinned()` reads at every call site inside 02-05, so the
## run record and the reads it produced are pinned to the same ingest.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))

pipeline_scripts <- c(
  "01_declare_requirements.R",
  "02_population_and_eligibility.R",
  "03_cost_windows.R",
  "04_covariates.R",
  "05_complete_case.R",
  "06_assemble_analytic_data.R",
  "07_estimate_att.R"
)

## Named, not inline in a loop, so it can be called from exactly one place
## below regardless of which branch (tracked or untracked) is taken.
source_pipeline_stages <- function() {
  for (script in pipeline_scripts) {
    cli::cli_rule("Sourcing {.file {script}}")
    source(file.path(app_dir, script))
  }
}

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed; sourcing stages with no run tracking
     (same convention as every numbered script's own smidata guard)."
  )
  source_pipeline_stages()
} else {
  if (is.null(smidata::smi_run_info())) {
    smidata::smi_run_start(
      script    = "doubletree/application/run_pipeline.R",
      tier      = if (identical(smidata::smi_env(), "local")) "local" else "smoke",
      ingest_id = config_ingest_id
    )
  } else {
    cli::cli_alert_info(
      "A run is already active ({.val {smidata::smi_run_info()$run_id}});
       not starting a second one."
    )
  }

  ## A failed stage must still close out the run record as "error" rather
  ## than leaving run.json stuck at status "running" forever -- see
  ## smi_run_finish()'s own docs (status: "ok" | "error" | "cancelled"). The
  ## original condition is re-signalled after finishing, never swallowed: no
  ## bare tryCatch(error = function(e) NULL) here -- a real pipeline error
  ## must still reach the caller/Rscript exit code, not just this run's log.
  tryCatch({
    source_pipeline_stages()
    smidata::smi_run_finish("ok")
  }, error = function(e) {
    smidata::smi_run_finish("error")
    stop(e)
  })
}
