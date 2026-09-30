## ============================================================================
## 91_leaf_budget_feasibility.R
##
## STANDALONE. Not part of the 01-07 cascade, and deliberately NOT gated by
## assert_no_blocking_open(): this measures WALL-CLOCK FEASIBILITY of
## estimate_att() at config_leaf_budget, on the real analytic data -- it does
## NOT produce, print as final, or write a reportable theta. It exists so
## that "the gate passes" never gets conflated with "the estimator can even
## finish at this budget," which is a separate, purely computational question
## this script answers instead.
##
## WHY THIS EXISTS. Benchmarked on the 400-row Tier-0 toy fixture (helper-toy.R,
## same 15-column shape as the real design matrix), estimate_att() runtime grew
## steeply, not linearly, with leaf_budget: 4.8s at leaf_budget=4, 94.5s at
## leaf_budget=6, still running past 5+ minutes at leaf_budget=15 (killed).
## That is on 400 TOY rows -- optimaltrees' GOSDT-style search cost is driven
## by the tree-structure space (leaf budget x covariate count), not
## meaningfully eased by n. Running the real ~250k-patient cohort at
## config_leaf_budget straight into 07's full estimate_att() +
## estimate_att_crossfit(K=5) pass, with no prior timing signal, risks
## discovering an infeasible budget only AFTER 03's multi-hour, 171.5GB
## cost-streaming pass has already been paid for -- 07 is the LAST stage.
## This script is the five-minute check that avoids that.
##
## SCOPE, DELIBERATELY NARROW. Only estimate_att() (the flagship), not
## estimate_att_crossfit() (K=5 trees plus its own cv_regularization sweep --
## strictly more expensive; if the flagship alone does not finish in budget,
## the crossfit companion will not either, and there is no need to pay for
## both signals). Only ONE leaf_budget value: config_leaf_budget, read from
## _config.R exactly as 07 would. This is a feasibility probe, not a tuning
## sweep -- re-run with a different config_leaf_budget to probe another value.
##
## THE RELIABLE TIME BOUND IS THE SHELL, NOT THIS SCRIPT'S OWN setTimeLimit().
## optimaltrees' solver is compiled C++ (Rcpp). A long-running C++ loop that
## does not call Rcpp::checkUserInterrupt() internally will NOT be
## interrupted by R's own setTimeLimit()/elapsed-time mechanism until control
## returns to R -- i.e. setTimeLimit() may report the overrun only AFTER the
## call finally finishes, not stop it early. setTimeLimit() below is
## best-effort defense in depth, not the real guarantee. THE REAL GUARANTEE IS
## AN OS-LEVEL KILL, invoked from the shell:
##
##   timeout 900 Rscript application/91_leaf_budget_feasibility.R
##
## (GNU coreutils' `timeout`, standard on the server's Linux; `gtimeout` via
## `brew install coreutils` if ever run on macOS.) `timeout` SIGTERMs the
## whole process tree, including any C++ call in progress, and this script's
## own exit code is irrelevant at that point -- `timeout` itself exits 124.
## Pick the outer budget to match how much server time you are willing to
## spend on a probe before concluding "not feasible at this leaf_budget."
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))

cli::cli_h1("91 -- leaf_budget feasibility probe (NOT a reportable estimate)")
cli::cli_alert_warning(
  "This script deliberately does NOT check assert_no_blocking_open(). It
   measures whether estimate_att() can even FINISH at
   {.code config_leaf_budget = {config_leaf_budget}}. This probe applies no
   diagnostics licence and reports no CI -- 07_estimate_att.R is the only
   place a reportable ATT is produced."
)

if (!config_has_smidata) {
  cli::cli_abort(
    "{.pkg smidata} is not installed. This probe needs the real analytic
     data written by 06_assemble_analytic_data.R on the server -- there is
     nothing to time locally."
  )
}

## ---- locate the analytic data: REPORTABLE first, PROVISIONAL as a warned
## fallback -----------------------------------------------------------------
##
## Order reversed 2026-09-30 when `diagnostics` resolved: 06 now writes the
## un-suffixed name once nothing blocks, and any
## analytic_cohort_PROVISIONAL.parquet still on disk predates that
## resolution. Preferring it -- as this script originally did, back when
## `diagnostics` was still open and only the PROVISIONAL file could exist --
## would silently time a stale cohort, the exact silent fallback this
## directory's conventions forbid. Falling back to it here is still allowed
## (this is a feasibility probe, not a reportable run) but is now WARNED,
## not silent.
output_dir <- file.path(app_dir, "output")
candidate_paths <- file.path(
  output_dir, c("analytic_cohort.parquet", "analytic_cohort_PROVISIONAL.parquet")
)
data_path <- candidate_paths[fs::file_exists(candidate_paths)][1]
if (is.na(data_path) || length(data_path) == 0L) {
  cli::cli_abort(c(
    "No analytic data found at {.file {candidate_paths}}.",
    "i" = "Run {.file 01_declare_requirements.R} through
           {.file 06_assemble_analytic_data.R} (or
           {.file run_pipeline.R}) first -- this script only TIMES
           estimate_att() on their output, it does not build it."
  ))
}
if (grepl("_PROVISIONAL", data_path, fixed = TRUE)) {
  cli::cli_alert_warning(
    "Timing the PROVISIONAL cohort -- no reportable
     {.file analytic_cohort.parquet} exists yet. Acceptable for a feasibility
     probe; re-run {.file 06_assemble_analytic_data.R} to time the reportable
     sample instead."
  )
}
cli::cli_alert_info("Reading {.file {data_path}}.")

analytic_data <- arrow::read_parquet(data_path)
x_cols <- design_matrix_columns()
X <- as.data.frame(analytic_data[, x_cols, drop = FALSE])
A <- as.integer(analytic_data$A)
Y <- as.numeric(analytic_data$Y)

## Same re-assert 07 does at the estimator boundary -- a parquet round-trip
## is where a 0/1 integer column can come back as something else.
assert_binary_design_matrix(X)

n_control <- sum(A == 0L)
cli::cli_alert_info(
  "n = {nrow(X)}; A = 1: {sum(A == 1L)}; A = 0 (control): {n_control}.
   estimate_att()'s own control-count floor at m_n = 1L needs n_control >=
   leaf_budget = {config_leaf_budget}; got {n_control} -- {if (n_control >= config_leaf_budget) 'satisfied' else 'WOULD ABORT on this alone, before any search cost is paid'}."
)

## ---- the timed call ---------------------------------------------------------
##
## setTimeLimit() is best-effort (see header) -- it is still set, because it
## is free and DOES work for the R-level bookkeeping around the C++ call
## even when it cannot interrupt the C++ call itself. transient = TRUE so a
## limit set here does not leak into any later call in the same session.
probe_elapsed_budget_secs <- 3600  ## 1 hour; edit before running if a
                                   ## different in-R signal is wanted --
                                   ## the OUTER shell `timeout` (see header)
                                   ## is what actually bounds wall-clock cost.
setTimeLimit(elapsed = probe_elapsed_budget_secs, transient = TRUE)

cli::cli_alert_info(
  "Starting estimate_att() at leaf_budget = {config_leaf_budget}. In-R
   elapsed guard: {probe_elapsed_budget_secs}s (best-effort only -- see this
   file's header on why an outer shell {.code timeout} is the real bound)."
)

t0 <- Sys.time()
fit <- tryCatch(
  {
    doubletree::estimate_att(
      X = X, A = A, Y = Y,
      leaf_budget = config_leaf_budget,
      outcome_type = "continuous"
    )
  },
  error = function(e) {
    cli::cli_alert_danger("estimate_att() errored: {conditionMessage(e)}")
    NULL
  }
)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
setTimeLimit(elapsed = Inf, transient = FALSE)  ## clear the guard explicitly

cli::cli_h2("Result")
if (is.null(fit)) {
  cli::cli_alert_danger(
    "estimate_att() did not return a fit (errored after {round(elapsed, 1)}s
     -- see the message above). NOT a timing result; see the error."
  )
} else {
  cli::cli_alert_success(
    "estimate_att() finished in {round(elapsed, 1)}s at leaf_budget =
     {config_leaf_budget} on n = {nrow(X)} real patients."
  )
  cli::cli_alert_info(
    "Realized leaves: e = {fit$n_leaves_e}, m0 = {fit$n_leaves_m0} (budget
     {fit$leaf_budget}); certified: e = {fit$certified_e}, m0 =
     {fit$certified_m0}."
  )
  cli::cli_alert_warning(
    "INFORMATIONAL ONLY -- not a reportable estimate. This script applies no
     diagnostics licence and computes no CI; 07_estimate_att.R is the only
     place a reportable ATT is produced. theta from this run:
     {signif(fit$theta, 4)}."
  )
  cli::cli_alert_info(
    "estimate_att_crossfit(K = {config_crossfit_k}) in 07 fits
     {config_crossfit_k} trees per nuisance (vs. this probe's 2 total) plus
     its own cv_regularization sweep -- budget accordingly relative to this
     number; it will not be faster."
  )
}

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("91_leaf_budget_feasibility.R complete.")
}
