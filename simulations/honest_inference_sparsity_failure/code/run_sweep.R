# ============================================================
# run_sweep.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# Drives the full grid by launching ONE Rscript PROCESS PER CELL (see run_cell.R's
# header for why process-per-cell). Resumable: a cell whose checkpoint already
# exists is skipped, so a kill mid-sweep costs at most the cell in flight.
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/code/run_sweep.R
#       (from the doubletree package root)
#
# Gates (all defaulted; a bare invocation runs the real study):
#   HIS_REPS        replications per cell for n < max(N_GRID)          [300]
#   HIS_REPS_MAX_N  replications per cell at n = max(N_GRID)           [150]
#   HIS_CELLS       comma-separated "dgp:n" cells to run           [all 15]
#   HIS_FORCE       "1" to re-run cells that already have a checkpoint  [0]
#
# The two-valued R is the spec's own sanctioned exception, measured by run_pilot.R
# and documented at REPS_DEFAULT in common.R: R = 300 at n in {500, 2000},
# R = 150 at n = 8000 (which alone is 65% of the wall time, because the anchor
# runs nested CV). It is stated here, in the metadata of every checkpoint, and in
# analyze.R's output -- never applied silently.
#
# Thread pinning: BLAS / OpenMP / Accelerate thread counts are forced to 1 in
# every child process. Not cosmetic -- it is one leg of the "no parallelism"
# memory rule; the others are worker_limit = 1L on both estimator calls and
# parallel = FALSE on the cross-fit call (both in common.R). A multi-threaded BLAS
# would reintroduce exactly the concurrency this study's memory posture avoids.
# ============================================================

source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                 "common.R"))

CELLS_GATE <- Sys.getenv("HIS_CELLS", "")
FORCE <- identical(Sys.getenv("HIS_FORCE", "0"), "1")
if (is.na(REPS_DEFAULT) || REPS_DEFAULT < 1L ||
    is.na(REPS_AT_MAX_N) || REPS_AT_MAX_N < 1L) {
  cli::cli_abort("HIS_REPS and HIS_REPS_MAX_N must be positive integers.")
}

grid <- design_grid()
if (nzchar(CELLS_GATE)) {
  wanted <- trimws(strsplit(CELLS_GATE, ",", fixed = TRUE)[[1]])
  grid <- grid[paste(grid$dgp, grid$n, sep = ":") %in% wanted, , drop = FALSE]
  if (nrow(grid) == 0L) {
    cli::cli_abort("HIS_CELLS={CELLS_GATE} matched no cell of the design grid.")
  }
}
grid$out <- vapply(seq_len(nrow(grid)),
                   function(i) cell_path(grid$dgp[[i]], grid$n[[i]], grid$reps[[i]]),
                   character(1))
grid$exists <- file.exists(grid$out)

cli::cli_h1("honest_inference_sparsity_failure sweep (S2)")
cli::cli_inform(c(
  "*" = "gates: HIS_REPS={REPS_DEFAULT} HIS_REPS_MAX_N={REPS_AT_MAX_N} HIS_CELLS={if (nzchar(CELLS_GATE)) CELLS_GATE else '<all>'} HIS_FORCE={FORCE}",
  "*" = "{nrow(grid)} cells, {sum(grid$reps)} replications = {2 * sum(grid$reps)} estimator calls (estimate_att + estimate_att_crossfit per rep)",
  "*" = "R = {REPS_DEFAULT} at n < {max(N_GRID)}, R = {REPS_AT_MAX_N} at n = {max(N_GRID)} (spec-sanctioned reduction; see common.R)",
  "*" = "{sum(grid$exists)} cell(s) already checkpointed{if (FORCE) ' (will be re-run: HIS_FORCE=1)' else ' (will be skipped)'}"
))

# Smallest cells first so a configuration problem surfaces cheaply.
grid <- grid[order(grid$n, match(grid$dgp, DGP_IDS)), , drop = FALSE]

rscript <- file.path(R.home("bin"), "Rscript")
child_env <- c(
  "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
  "VECLIB_MAXIMUM_THREADS=1", "NUMEXPR_NUM_THREADS=1"
)

t_all <- Sys.time()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  if (g$exists && !FORCE) {
    cli::cli_alert_info("skip {g$dgp} n={g$n} r={g$reps} (checkpoint present)")
    next
  }
  cli::cli_alert("[{i}/{nrow(grid)}] {g$dgp} eps={g$eps} n={g$n} r={g$reps} -- elapsed {signif(as.numeric(difftime(Sys.time(), t_all, units = 'mins')), 3)} min")
  st <- system2(
    rscript,
    args = c("--vanilla",
             file.path("simulations", "honest_inference_sparsity_failure",
                       "code", "run_cell.R"),
             g$dgp, g$n, g$reps),
    env = child_env
  )
  if (st != 0L) {
    # Do NOT press on: a nonzero exit means that cell has no checkpoint, and
    # continuing would produce a sweep that looks complete but is not.
    cli::cli_abort(c(
      "Cell {g$dgp} n={g$n} exited with status {st}.",
      i = "Completed cells are already checkpointed; fix the cause and re-run (finished cells will be skipped)."
    ))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t_all, units = "mins"))
missing <- grid$out[!file.exists(grid$out)]
if (length(missing)) {
  cli::cli_abort(c("{length(missing)} cell checkpoint(s) missing after the sweep.",
                   stats::setNames(basename(missing), rep("*", length(missing)))))
}
cli::cli_alert_success(
  "sweep complete: {nrow(grid)} cells in {signif(elapsed, 3)} min ({signif(elapsed / 60, 3)} h). Next: Rscript .../code/analyze.R"
)
