#!/bin/bash
# =============================================================================
# slurm/run_simulations.slurm -- SLURM array job for head_to_head_comparison (S5)
# =============================================================================
# NOT run directly. launch_all_simulations.sh / launch_subset.sh submit it once
# per REGIME, supplying --array / --time / --mem / --partition / --output on the
# sbatch command line (sizing comes from slurm/units.R, so the submitter and the
# worker cannot disagree about it) and exporting:
#
#   PKG_ROOT     absolute path to the doubletree package root (has DESCRIPTION)
#   SCRATCH_DIR  run-specific scratch dir; this regime's shards land here
#   REGIME       regime id (R1..R6); selects the per-regime unit table
#   TARGET_SECS  target wall cost per task -- MUST be the value the launcher used
#                to compute --array, or the two would build different maps
#
# One array task = one contiguous block of replications of one (regime, dgp, n)
# cell. slurm/units.R holds the map and the reasoning behind the block sizing.
# =============================================================================
#SBATCH --job-name=h2h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# --array, --time, --mem, --partition, --output, --error come from the launcher.
# cpus-per-task=1 is not a default to be tuned: every arm runs single-threaded by
# construction (worker_limit = 1L, ranger num.threads = 1L, BLAS pinned in
# run_single_replication.R) because spec §5's per-method wall-clock comparison is
# only meaningful if no arm is quietly given more cores than another.

set -euo pipefail

# --- Module environment -------------------------------------------------------
# Same versions this project's other cluster studies use today (copied from
# six-approach-arbitration/slurm/{array.slurm,submit_per_method.sh}, 2026-07).
# NOTE: functional_consistency/slurm/README_O2.md still documents gcc/9.2.0 and
# R/4.2.1 -- that is a stale April-2026 pairing and is NOT what this study uses.
#
# Module R carries no working devtools/pkgload dev-load, which is exactly why
# run_single_replication.R sets H2H_USE_INSTALLED=1 and the packages must be
# R CMD INSTALLed first (README_O2.md, prerequisites).
module load gcc/14.2.0 2>/dev/null || module load gcc || true
module load R/4.4.2   2>/dev/null || module load R
# No `|| true` on the R load: without R this task cannot run at all, and failing
# here names the cause, whereas continuing would surface it as "Rscript: command
# not found" in 134 logs at once.

: "${PKG_ROOT:?PKG_ROOT not exported}"
: "${SCRATCH_DIR:?SCRATCH_DIR not exported}"
: "${REGIME:?REGIME not exported}"
TARGET_SECS="${TARGET_SECS:-1800}"

mkdir -p "${SCRATCH_DIR}"

echo "[$(date '+%F %T')] host=$(hostname) regime=${REGIME} array_task=${SLURM_ARRAY_TASK_ID} target_secs=${TARGET_SECS}"

# --- SIGTERM trap -------------------------------------------------------------
# On a wall-time TIMEOUT SLURM sends SIGTERM and then, after a grace period,
# SIGKILL. run_single_replication.R flushes a partial every FLUSH_EVERY
# replications, so a timeout loses at most that many replications. Report how many
# partials/shards this scratch dir holds so the over-run is visibly recoverable:
# resubmitting the SAME run-id resumes each task from its last flush. The
# "received SIGTERM" text is the sentinel check_progress.sh greps for.
on_term() {
  local n_shard n_part
  n_shard=$(find "${SCRATCH_DIR}" -maxdepth 1 -name 'batch_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
  n_part=$(find "${SCRATCH_DIR}" -maxdepth 1 -name 'part_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
  echo "[$(date '+%F %T')] regime=${REGIME} task=${SLURM_ARRAY_TASK_ID} received SIGTERM (likely TIMEOUT); ${n_shard} shard(s) and ${n_part} partial(s) on disk are safe. Resubmit the same run-id to resume." >&2
  exit 143
}
trap on_term TERM

# --- Run this task's block ----------------------------------------------------
# `--no-init-file` is LOAD-BEARING, not tidiness. The doubletree repository ships a
# TRACKED .Rprofile that runs devtools::load_all("../optimaltrees") at startup, and
# R reads ./.Rprofile as the user profile whenever R_PROFILE_USER is unset. Without
# this flag a task launched anywhere near the package root would dev-load
# optimaltrees from source, common.R's isNamespaceLoaded() guard would find it
# already present, and the whole regime would silently run on UNINSTALLED source
# code while reporting pkg_load_path = "installed". run_single_replication.R's
# assert_installed_load() is the backstop that turns that into a hard failure; this
# flag is what stops it happening. .Renviron is deliberately still read, so
# R_LIBS_USER set per-machine keeps working.
#
# Backgrounded and `wait`ed so the trap can fire while Rscript is still running;
# a foreground child would swallow the signal until it returned.
Rscript --no-init-file "${PKG_ROOT}/simulations/head_to_head_comparison/slurm/run_single_replication.R" \
  --regime      "${REGIME}" \
  --task-id     "${SLURM_ARRAY_TASK_ID}" \
  --pkg-root    "${PKG_ROOT}" \
  --scratch-dir "${SCRATCH_DIR}" \
  --target-secs "${TARGET_SECS}" &
wait $!
status=$?

echo "[$(date '+%F %T')] regime=${REGIME} task=${SLURM_ARRAY_TASK_ID} finished with status ${status}"
exit ${status}
