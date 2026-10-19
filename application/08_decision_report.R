## ============================================================================
## 08_decision_report.R
##
## REVIEW REPORT: every lever in force and what it did, in tables the PI can
## check and quote in the paper. Builds from the intermediates 02, 03, 05 and 06
## SAVE as they run (application/output/private/report_inputs/*.rds); it recomputes
## nothing heuristically and reads no data.
##
##   Rscript application/08_decision_report.R
##
## Levers come from _config.R's own objects (report_levers()), so the report
## cannot drift from the code. Tables come from smidata's decision-report
## builders (smi_report_flow / _measure_period / _costs / _cost_categories /
## _enrollment) plus this paper's own (prior_cost floor and bins, design
## degeneracy), all suppressed at n < 11. smidata::smi_report() then resolves and
## AUDITS the identities declared between them (helpers/decision_report.R); the
## report's header line states what was audited. The audit guarantees that no
## masked value is determined by the published cells through those declared
## identities, and nothing more (undeclared identities, dollar dominance and
## narrowing inside 1..10 are not covered).
##
## OUTPUT. The Markdown goes to application/output/export/decision_report/ (local,
## gitignored). The tables leave the secure server ONLY through
## smidata::smi_report_export() -> smi_export(), and only when
## smidata::smi_is_server() and the intermediates were recorded on the server.
## Off-server this script still runs end to end: with no intermediates it uses
## helpers/decision_report.R::toy_report_inputs() (hand-calculated toy counts),
## says so in the report title, writes to output/private/decision_report_toy/, and
## exports nothing.
##
## Run after 06 (it does not need 07).
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "discretize_prior_cost.R"))
source(file.path(app_dir, "helpers", "disclosure.R"))
source(file.path(app_dir, "helpers", "decision_report.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("08 -- decision report")

if (!config_has_smidata ||
    !"smi_report" %in% getNamespaceExports("smidata")) {
  cli::cli_alert_info(
    "{.pkg smidata} with {.fn smi_report} (decision-report builders) is not
     installed; skipping (same convention as 01-07)."
  )
} else {
  inputs_dir <- report_inputs_dir(app_dir)
  stage_files <- c(p02 = "02_population", p03 = "03_costs",
                   p05 = "05_enrollment", p06 = "06_final")
  required <- c("p02", "p03", "p06")
  have <- file.exists(file.path(inputs_dir, paste0(stage_files, ".rds")))
  names(have) <- names(stage_files)
  use_toy <- !any(have)
  ## Some intermediates but not all is a half-finished or mixed run, never "no
  ## pipeline output": do not discard what exists and swap in toy numbers.
  if (!use_toy && !all(have[required])) {
    cli::cli_abort(c(
      "Partial report intermediates in {.path {inputs_dir}}.",
      "x" = "Have {.file {unname(stage_files[have])}}; missing {.file {unname(stage_files[required][!have[required]])}}.",
      "i" = "Re-run the pipeline end to end (or delete {.path {inputs_dir}} to build the toy report)."
    ))
  }
  if (!use_toy && have[["p06"]] && !have[["p05"]]) {
    cli::cli_abort(c(
      "A 06 intermediate exists without 05's: 06 requires 05, so these are from different runs.",
      "i" = "Re-run the pipeline end to end."
    ))
  }

  if (use_toy) {
    if (identical(smidata::smi_env(), "server")) {
      cli::cli_abort(c(
        "Missing report intermediate{?s}: {.file {unname(stage_files[required][!have[required]])}}.",
        "i" = "02, 03 and 06 save them as they run; run the pipeline on the
               server first. A toy report is never produced on the server."
      ))
    }
    cli::cli_alert_warning(
      "No pipeline intermediates in {.path {inputs_dir}}: building the report
       from the TOY inputs (hand-calculated numbers, NOT results)."
    )
    inputs <- toy_report_inputs()
  } else {
    inputs <- lapply(stage_files, function(f) {
      path <- file.path(inputs_dir, paste0(f, ".rds"))
      if (file.exists(path)) readRDS(path) else NULL
    })
    names(inputs) <- names(stage_files)
  }

  report <- build_decision_report(
    inputs,
    title = if (use_toy) {
      "doubletree application: decision report (TOY INPUTS, not results)"
    } else if (!identical(inputs$p06$.meta$env, "server")) {
      sprintf("doubletree application: decision report (%s run, not real data)",
              inputs$p06$.meta$env)
    } else {
      "doubletree application: decision report"
    }
  )

  ## The audited report is shareable (export/); a toy report is not a result and
  ## stays in private/.
  md_path <- output_file(if (use_toy) "private" else "export",
                         if (use_toy) "decision_report_toy" else "decision_report",
                         "decision_report.md", app_dir = app_dir)
  writeLines(smidata::smi_report_md(report), md_path)
  cli::cli_alert_success(
    "Wrote {.file {md_path}} ({length(report$tables)} table{?s}; cells < 11 masked; audited against {length(report$identities)} declared identit{?y/ies}, {report$resolved} further cell{?s} masked by the repair step)."
  )

  ## Export: server only, real intermediates only, recorded on this environment.
  recorded_on_server <- !use_toy && identical(inputs$p06$.meta$env, "server")
  if (smidata::smi_is_server() && recorded_on_server) {
    manifest <- smidata::smi_report_export(report, prefix = "doubletree_application")
    cli::cli_alert_success("Exported {nrow(manifest)} table{?s} through smi_export().")
  } else {
    cli::cli_alert_info(
      "Not exporting: {if (!smidata::smi_is_server()) 'not on the server' else 'intermediates were not recorded on the server'}."
    )
  }
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("08_decision_report.R complete.")
}
