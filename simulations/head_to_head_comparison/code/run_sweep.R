# ============================================================
# run_sweep.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md  (§6)
#
# Run:  Rscript simulations/head_to_head_comparison/code/run_sweep.R [regimes]
#       (from the doubletree package root; default: all regimes)
#       e.g.  Rscript .../run_sweep.R R1,R3
#
# THE FULL-SCALE RUN, at spec §6's nsim targets over the FULL n grid. Sized by
# code/run_pilot.R, which must be run FIRST -- its section 4 projects the wall
# time per regime and this script refuses to guess at it.
#
# Sequential, one cell at a time, checkpointing after every cell. A kill mid-sweep
# costs at most the cell in flight. Already-complete checkpoints are SKIPPED, so
# re-running after an interruption resumes rather than restarts.
#
# ENVIRONMENT GATES (declared here, in one block, with defaults, and echoed):
#   H2H_REGIMES   comma-separated regime ids            (default: all)
#   H2H_NSIM      override nsim for EVERY regime        (default: spec §6 targets)
#   H2H_FORCE     "1" to recompute existing checkpoints (default: "0")
#
# THREAD PINNING. BLAS/OpenMP threads are pinned to 1 below. Without this, a
# multithreaded BLAS makes the per-method wall-clock comparison of spec §5
# meaningless AND makes the run's memory footprint unpredictable
# (simulations/docs/MEMORY_SAFE_SIMULATIONS.md). ranger's own threads are pinned
# inside att_forest(); common.R asserts that.
# ============================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           RCPP_PARALLEL_NUM_THREADS = "1")

source(file.path("simulations", "head_to_head_comparison", "code", "common.R"))
source(file.path("simulations", "head_to_head_comparison", "code", "analyze.R"))

suppressPackageStartupMessages(library(knitr))

# ---- 0 gates -------------------------------------------------------------

.args <- commandArgs(trailingOnly = TRUE)
SWEEP_REGIMES <- local({
  from_arg <- if (length(.args) >= 1L && nzchar(.args[[1L]])) .args[[1L]] else ""
  raw <- if (nzchar(from_arg)) from_arg else Sys.getenv("H2H_REGIMES", "")
  if (nzchar(raw)) trimws(strsplit(raw, ",", fixed = TRUE)[[1L]]) else REGIME_IDS
})
NSIM_OVERRIDE <- local({
  v <- Sys.getenv("H2H_NSIM", "")
  if (nzchar(v)) as.integer(v) else NA_integer_
})
FORCE <- Sys.getenv("H2H_FORCE", "0") == "1"

unknown <- setdiff(SWEEP_REGIMES, REGIME_IDS)
if (length(unknown)) {
  cli::cli_abort("Unknown regime(s) {.val {unknown}}; expected {.val {REGIME_IDS}}.")
}

cli::cli_h1("Full sweep")
cli::cli_inform(c(
  "*" = "regimes = {paste(SWEEP_REGIMES, collapse = ', ')}",
  "*" = "nsim = {if (is.na(NSIM_OVERRIDE)) 'spec §6 targets per regime' else NSIM_OVERRIDE}",
  "*" = "force recompute = {FORCE}",
  "*" = "BLAS/OpenMP threads pinned to 1; worker_limit = {WORKER_LIMIT}; parallel_cv = {PARALLEL_CV}"
))

#' nsim for a regime (spec §6, unless overridden)
nsim_for <- function(regime) {
  if (!is.na(NSIM_OVERRIDE)) return(NSIM_OVERRIDE)
  NSIM_FULL[[regime]]
}

# ---- 1 the gate that must pass before any replication --------------------

cli::cli_h1("DGP gate")
rem <- verify_main_effects_remainder(dgp_spec_shared_interaction())
cli::cli_alert_success(
  "DGP-A remainder = {signif(rem$bilinear_remainder_P, 6)} (nonzero, matches spec §3)."
)

# ---- 2 the full grid ----------------------------------------------------

# design_grid() returns the PILOT n values; the sweep uses N_GRID_FULL.
grid <- do.call(rbind, lapply(SWEEP_REGIMES, function(r) {
  spec <- REGIMES[[r]]
  g <- expand.grid(n = N_GRID_FULL[[r]], dgp = spec$dgps,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g <- g[order(match(g$dgp, spec$dgps), g$n), , drop = FALSE]
  data.frame(regime = r, dgp = g$dgp, n = g$n, reps = nsim_for(r),
             stringsAsFactors = FALSE)
}))
rownames(grid) <- NULL
grid$out <- mapply(cell_path, grid$regime, grid$dgp, grid$n, grid$reps)
grid$exists <- file.exists(grid$out)

cli::cli_h2("{nrow(grid)} cell(s); {sum(grid$exists)} already on disk")
print(knitr::kable(grid[, c("regime", "dgp", "n", "reps", "exists")],
                   format = "simple", row.names = FALSE))

todo <- if (FORCE) grid else grid[!grid$exists, , drop = FALSE]
if (!nrow(todo)) {
  cli::cli_alert_success("Nothing to do; every checkpoint exists. Set H2H_FORCE=1 to recompute.")
} else {
  t_start <- Sys.time()
  for (i in seq_len(nrow(todo))) {
    g <- todo[i, ]
    cli::cli_inform("--- cell {i}/{nrow(todo)} ---")
    run_cell(g$regime, g$dgp, g$n, reps = g$reps, out = g$out,
             on_error = "record", heartbeat = 50L)
    gc(full = TRUE, verbose = FALSE)
  }
  cli::cli_alert_success(
    "Sweep complete: {nrow(todo)} cell(s) in {signif(as.numeric(difftime(Sys.time(), t_start, units = 'hours')), 4)} h."
  )
}

# ---- 3 summarise everything on disk -------------------------------------

STAMP <- format(Sys.time(), "%Y%m%d-%H%M%S")
res <- load_results(SWEEP_REGIMES)
smry <- summarise_cell(res)
utils::write.csv(smry, file.path(DIR_TABLES, sprintf("sweep_summary_%s.csv", STAMP)),
                 row.names = FALSE)

ratio_tabs <- list()
for (r in SWEEP_REGIMES) {
  cli::cli_h2("{r} -- {REGIMES[[r]]$label}")
  print_regime_summary(smry, r)
  rt <- print_ratio_summary(res, r)
  if (!is.null(rt)) ratio_tabs[[r]] <- rt
}
if (length(ratio_tabs)) {
  utils::write.csv(do.call(rbind, ratio_tabs),
                   file.path(DIR_TABLES, sprintf("sweep_rmse_ratios_%s.csv", STAMP)),
                   row.names = FALSE)
}
cli::cli_alert_success("Wrote sweep tables to {.path {DIR_TABLES}} (stamp {STAMP}).")
cli::cli_h1("Read spec §0 before quoting any number above. This is not a leaderboard.")
