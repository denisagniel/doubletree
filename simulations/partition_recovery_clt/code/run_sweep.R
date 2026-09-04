# ============================================================
# run_sweep.R
# Study: partition_recovery_clt  (doubletree, S1)
#
# Drives the full grid by launching ONE Rscript PROCESS PER CELL (see
# run_cell.R's header for why process-per-cell). Resumable: a cell whose
# checkpoint already exists is skipped, so a kill mid-sweep costs at most the
# cell in flight.
#
# Run:  Rscript simulations/partition_recovery_clt/code/run_sweep.R
#       (from the doubletree package root)
#
# Gates (all defaulted; a bare invocation runs the real study):
#   PRC_REPS   replications per cell                            [1000, spec §5]
#   PRC_CELLS  comma-separated "dgp:n" cells to run             [all 18]
#   PRC_FORCE  "1" to re-run cells that already have a checkpoint  [0]
#
# Thread pinning: BLAS / OpenMP / Accelerate thread counts are forced to 1 in
# every child process. Not cosmetic -- it is half of the "no parallelism" memory
# rule, the other half being worker_limit = 1L inside common.R. A multi-threaded
# BLAS would reintroduce exactly the concurrency this study's memory posture
# exists to avoid.
# ============================================================

source(file.path("simulations", "partition_recovery_clt", "code", "common.R"))

N_REPS <- as.integer(Sys.getenv("PRC_REPS", "1000"))
CELLS_GATE <- Sys.getenv("PRC_CELLS", "")
FORCE <- identical(Sys.getenv("PRC_FORCE", "0"), "1")
if (is.na(N_REPS) || N_REPS < 1L) {
  cli::cli_abort("PRC_REPS must be a positive integer, got {.val {Sys.getenv('PRC_REPS')}}.")
}

grid <- design_grid()
if (nzchar(CELLS_GATE)) {
  wanted <- trimws(strsplit(CELLS_GATE, ",", fixed = TRUE)[[1]])
  grid <- grid[paste(grid$dgp, grid$n, sep = ":") %in% wanted, , drop = FALSE]
  if (nrow(grid) == 0L) {
    cli::cli_abort("PRC_CELLS={CELLS_GATE} matched no cell of the design grid.")
  }
}
grid$out <- vapply(seq_len(nrow(grid)),
                   function(i) cell_path(grid$dgp[[i]], grid$n[[i]], N_REPS),
                   character(1))
grid$exists <- file.exists(grid$out)

cli::cli_h1("partition_recovery_clt sweep")
cli::cli_inform(c(
  "*" = "gates: PRC_REPS={N_REPS} PRC_CELLS={if (nzchar(CELLS_GATE)) CELLS_GATE else '<all>'} PRC_FORCE={FORCE}",
  "*" = "{nrow(grid)} cells x {N_REPS} reps = {nrow(grid) * N_REPS} estimate_att() calls",
  "*" = "{sum(grid$exists)} cell(s) already checkpointed{if (FORCE) ' (will be re-run: PRC_FORCE=1)' else ' (will be skipped)'}"
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
    cli::cli_alert_info("skip {g$dgp} n={g$n} (checkpoint present)")
    next
  }
  cli::cli_alert("[{i}/{nrow(grid)}] {g$regime} {g$dgp} n={g$n} r={N_REPS}")
  st <- system2(
    rscript,
    args = c("--vanilla",
             file.path("simulations", "partition_recovery_clt", "code", "run_cell.R"),
             g$dgp, g$n, N_REPS),
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
  "sweep complete: {nrow(grid)} cells in {signif(elapsed, 3)} min. Next: Rscript .../code/analyze.R"
)
