#!/bin/bash
# =============================================================================
# array.slurm -- SLURM array job for the "six-approach-arbitration" study
# =============================================================================
# NOT run directly. submit.sh sets the sizing (#SBATCH --array, --time, --mem)
# on the sbatch command line and exports the variables below. This script maps
# its array index to a GLOBAL task id, prepares scratch, and runs one task.
#
# Exported by submit.sh (via --export):
#   STUDY_DIR    absolute path to the study dir (config/, R/, slurm/)
#   SCRATCH_DIR  run-specific scratch dir (task files + logs live here)
#   REPS_PER_JOB units per task
#   ARRAY_OFFSET added to SLURM_ARRAY_TASK_ID to get the global task id
#                (arrays are chunked at <=1000; wave k covers a different range)
# =============================================================================
#SBATCH --job-name=six-approach-arbitration
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# --array, --time, --mem, --output, --error are supplied by submit.sh.

set -euo pipefail

# --- Module environment (module R has no devtools; use library()) ------------
module load gcc/14.2.0 2>/dev/null || module load gcc || true
module load R/4.4.2   2>/dev/null || module load R

# --- Global task id -----------------------------------------------------------
GLOBAL_TASK_ID=$(( SLURM_ARRAY_TASK_ID + ARRAY_OFFSET ))

# --- Per-method block offsets -------------------------------------------------
# submit_per_method.sh exports UNIT_OFFSET (global start of this method's unit
# block) and MAX_UNIT (its end); run_replication.R maps GLOBAL_TASK_ID to units as
# UNIT_OFFSET + (task-1)*reps + 1 .. min(UNIT_OFFSET + task*reps, MAX_UNIT). These
# MUST be forwarded as CLI flags -- without them run_replication.R defaults to
# unit_offset=0, so EVERY method's array recomputes the FIRST block (`full`) and
# mislabels its output (the 2026-07-10 contamination bug). The stock global
# submit.sh does not set them; default to 0 / (unset) so it still works.
UNIT_OFFSET="${UNIT_OFFSET:-0}"
OFFSET_ARGS=(--unit-offset "${UNIT_OFFSET}")
if [[ -n "${MAX_UNIT:-}" ]]; then
  OFFSET_ARGS+=(--max-unit "${MAX_UNIT}")
fi

echo "[$(date '+%F %T')] host=$(hostname) array_idx=${SLURM_ARRAY_TASK_ID} array_offset=${ARRAY_OFFSET} global_task=${GLOBAL_TASK_ID} unit_offset=${UNIT_OFFSET} max_unit=${MAX_UNIT:-<none>}"

mkdir -p "${SCRATCH_DIR}"

# --- SIGTERM trap (M3): on a wall-time TIMEOUT SLURM sends SIGTERM, then (after a
# grace period) SIGKILL. run_replication.R flushes each unit to partials/ the
# moment it finishes, so a timeout loses at most the one in-flight unit. Report
# how many per-unit partials this task has flushed so the over-run is visibly
# recoverable (resubmitting the same RUN_ID skips every completed unit). The
# "received SIGTERM" line is a sentinel monitor.sh keys on to flag timeouts.
on_term() {
  local pdir="${SCRATCH_DIR}/partials"
  local n=0
  [[ -d "${pdir}" ]] && n=$(find "${pdir}" -maxdepth 1 -name 'unit_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  echo "[$(date '+%F %T')] global_task=${GLOBAL_TASK_ID} received SIGTERM (likely TIMEOUT); ${n} unit partial(s) flushed and safe. Resubmit the same run-id to resume." >&2
  exit 143
}
trap on_term TERM

# --- Run one task's block of units (checkpointing per unit) ------------------
# Run in the BACKGROUND and `wait` so the SIGTERM trap can fire while Rscript is
# still running (a foreground child would swallow the signal until it returns).
Rscript "${STUDY_DIR}/slurm/run_replication.R" \
  --task-id      "${GLOBAL_TASK_ID}" \
  --reps-per-job "${REPS_PER_JOB}" \
  --study-dir    "${STUDY_DIR}" \
  --scratch-dir  "${SCRATCH_DIR}" \
  "${OFFSET_ARGS[@]}" &
wait $!

echo "[$(date '+%F %T')] global_task=${GLOBAL_TASK_ID} finished with status $?"
