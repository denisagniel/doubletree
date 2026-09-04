# ============================================================
# run_cell.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# Runs ONE (dgp, n) cell and writes ONE checkpoint. Deliberately one cell per
# PROCESS, invoked from the shell:
#
#     Rscript simulations/honest_inference_sparsity_failure/code/run_cell.R <dgp> <n> <reps>
#
# (from the doubletree package root). Defaults: eps0 500 5 -- a bare invocation is
# the micro-test, not the biggest job.
#
# WHY ONE PROCESS PER CELL rather than a loop over cells in one session, matching
# partition_recovery_clt/code/run_cell.R and threshold_superconsistency's: this
# codebase has previously exhausted host memory during tree simulations (see
# ../docs/MEMORY_SAFE_SIMULATIONS.md and MEMORY.md's [LEARN:rashomon-memory]).
# Process-per-cell means (a) nothing accumulates across cells -- the OS reclaims
# everything at exit, (b) a crash costs at most one cell, and (c) a staged
# rollout is enforced by how the script is invoked rather than by discipline
# inside a loop.
#
# MEMORY RULES OBSERVED (see common.R for where each is set):
#   * NO Rashomon enumeration. estimate_att() -> bisect_lambda_to_budget() ->
#     fit_tree(); estimate_att_crossfit() -> fit_tree_with_cv() ->
#     cv_regularization_adaptive() -> fit_tree(). fit_rashomon() is never called
#     on either path. Rashomon set enumeration was the documented cause of every
#     previous memory-exhaustion incident, not ordinary single-tree fitting.
#   * NO parallelism. worker_limit = 1L on BOTH estimator calls, and
#     parallel = FALSE on the cross-fit call so cv_regularization_adaptive()'s
#     furrr/future fold-parallelism (ON by default) is switched off explicitly
#     rather than relying on no future plan happening to be set. BLAS/OpenMP
#     thread counts are pinned to 1 by run_sweep.R in every child process.
#   * gc(full = TRUE) on a heartbeat, and a checkpoint at the end of the cell.
#
# Sourced-and-callable: run_cell() in common.R is a plain function with working
# defaults, so this file can be sourced in a live session and stepped through.
# The invocation below fires only for a genuine top-level Rscript run --
# `!interactive()` alone is NOT enough, since it is also FALSE when another
# script sources this one under Rscript; sys.nframe() == 0L is the discriminator.
# ============================================================

source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                 "common.R"))

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  dgp_arg  <- if (length(args) >= 1L) args[[1L]] else "eps0"
  n_arg    <- if (length(args) >= 2L) as.integer(args[[2L]]) else 500L
  reps_arg <- if (length(args) >= 3L) as.integer(args[[3L]]) else 5L
  if (!dgp_arg %in% DGP_IDS) {
    cli::cli_abort("Unknown dgp {.val {dgp_arg}}; expected one of {.val {DGP_IDS}}.")
  }
  if (is.na(n_arg) || is.na(reps_arg)) {
    cli::cli_abort("Usage: Rscript run_cell.R <dgp> <n> <reps>")
  }
  run_cell(dgp_id = dgp_arg, n = n_arg, reps = reps_arg,
           on_error = "record", heartbeat = 25L)
}
