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

# --- Run one task's block of units -------------------------------------------
Rscript "${STUDY_DIR}/slurm/run_replication.R" \
  --task-id      "${GLOBAL_TASK_ID}" \
  --reps-per-job "${REPS_PER_JOB}" \
  --study-dir    "${STUDY_DIR}" \
  --scratch-dir  "${SCRATCH_DIR}" \
  "${OFFSET_ARGS[@]}"

echo "[$(date '+%F %T')] global_task=${GLOBAL_TASK_ID} finished with status $?"
