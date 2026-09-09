#!/bin/bash
# =============================================================================
# slurm/check_progress.sh -- monitor one run of the head_to_head_comparison sweep
# =============================================================================
# Run ON O2, from anywhere:
#
#   bash slurm/check_progress.sh                      # newest run under scratch
#   RUN_ID=20260909-131500_abc1234 bash slurm/check_progress.sh
#   WATCH=1 bash slurm/check_progress.sh              # re-print every 60 s
#
# Reports, per regime: queued/running task counts from squeue, finished shards on
# disk against the number slurm/units.R says the array should produce, partials in
# flight, and any task that was killed at the wall limit.
#
# TIMEOUTS ARE SURFACED EXPLICITLY. A task killed at --time leaves its finished
# replications on disk as a partial and its log carries the "received SIGTERM"
# sentinel run_simulations.slurm's trap prints. Those tasks are recoverable by
# resubmitting the same run id -- but only if someone notices, and SLURM's own
# accounting states them as plain failures.
# =============================================================================

set -euo pipefail

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="head_to_head_comparison"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

SCRATCH_ROOT="${SCRATCH_ROOT:-/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}}"
H2H_TARGET_SECS="${H2H_TARGET_SECS:-1800}"
WATCH="${WATCH:-0}"
export H2H_TARGET_SECS

# Default to the newest run directory rather than making the caller paste a run id
# out of the manifest; the id is echoed below either way so the report is never
# ambiguous about which run it describes.
if [[ -z "${RUN_ID:-}" ]]; then
  if [[ ! -d "${SCRATCH_ROOT}" ]]; then
    echo "No scratch root at ${SCRATCH_ROOT}; nothing has been launched." >&2
    exit 1
  fi
  RUN_ID="$(ls -1t "${SCRATCH_ROOT}" | head -n 1)"
  [[ -n "${RUN_ID}" ]] || { echo "No run directories under ${SCRATCH_ROOT}." >&2; exit 1; }
fi
SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
LOG_DIR="${SCRATCH_DIR}/logs"

[[ -d "${SCRATCH_DIR}" ]] || { echo "No such run: ${SCRATCH_DIR}" >&2; exit 1; }

report_once() {
  echo "=============================================================="
  echo " ${STUDY_NAME}   run id: ${RUN_ID}"
  echo " scratch: ${SCRATCH_DIR}"
  echo " time   : $(date '+%F %T')"
  echo "=============================================================="

  # --- SLURM's view ----------------------------------------------------------
  echo
  echo "--- squeue ---"
  if command -v squeue >/dev/null 2>&1; then
    squeue -u "${USER}" -o "%.12i %.28j %.10T %.10M %.10l %.6D %R" || true
    local n_run n_pend
    n_run=$(squeue -u "${USER}" -h -t RUNNING -n "${STUDY_NAME}-R1,${STUDY_NAME}-R2,${STUDY_NAME}-R3,${STUDY_NAME}-R4,${STUDY_NAME}-R5,${STUDY_NAME}-R6" 2>/dev/null | wc -l | tr -d ' ' || true)
    n_pend=$(squeue -u "${USER}" -h -t PENDING -n "${STUDY_NAME}-R1,${STUDY_NAME}-R2,${STUDY_NAME}-R3,${STUDY_NAME}-R4,${STUDY_NAME}-R5,${STUDY_NAME}-R6" 2>/dev/null | wc -l | tr -d ' ' || true)
    echo "  running: ${n_run}   pending: ${n_pend}"
  else
    echo "  (squeue not on PATH -- not an O2 login node?)"
  fi

  # --- shards on disk vs the array size ------------------------------------
  # Expected counts come from the same table the workers use, so "83/134" cannot
  # drift from what was actually submitted.
  echo
  echo "--- shards on disk vs slurm/units.R (target ${H2H_TARGET_SECS} s/task) ---"
  local expected_file
  expected_file="$(mktemp)"
  Rscript "${SLURM_DIR}/print_sizing.R" --format env --pkg-root "${PKG_ROOT}" \
    > "${expected_file}" 2>/dev/null || { echo "  (could not compute expected task counts)"; }

  local total_done=0 total_expected=0
  printf "  %-7s %8s %8s %9s %9s\n" regime shards expected pct partials
  while read -r REGIME N_TASKS _WALLTIME _MEM _PART; do
    [[ -z "${REGIME}" ]] && continue
    local d done part pct
    d="${SCRATCH_DIR}/${REGIME}"
    done=0; part=0
    if [[ -d "${d}" ]]; then
      done=$(find "${d}" -maxdepth 1 -name 'batch_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
      part=$(find "${d}" -maxdepth 1 -name 'part_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
    fi
    pct=$(( N_TASKS > 0 ? 100 * done / N_TASKS : 0 ))
    printf "  %-7s %8s %8s %8s%% %9s\n" "${REGIME}" "${done}" "${N_TASKS}" "${pct}" "${part}"
    total_done=$(( total_done + done ))
    total_expected=$(( total_expected + N_TASKS ))
  done < "${expected_file}"
  rm -f "${expected_file}"
  if (( total_expected > 0 )); then
    printf "  %-7s %8s %8s %8s%%\n" TOTAL "${total_done}" "${total_expected}" \
      "$(( 100 * total_done / total_expected ))"
  fi

  # --- wall-limit kills -----------------------------------------------------
  # `|| true` on every counting pipeline below is load-bearing, not decoration:
  # under `set -o pipefail` a grep/find that matches nothing returns 1, the
  # pipeline inherits it, and `set -e` would abort this report mid-print. "No
  # timeouts yet" is the normal case, so without it the monitor would only survive
  # runs that had already gone wrong.
  echo
  echo "--- wall-limit kills (SIGTERM sentinel in the logs) ---"
  if [[ -d "${LOG_DIR}" ]]; then
    local n_term
    n_term=$(grep -rl "received SIGTERM" "${LOG_DIR}" 2>/dev/null | wc -l | tr -d ' ' || true)
    if (( n_term > 0 )); then
      echo "  ${n_term} task(s) hit their wall limit. Their finished replications are"
      echo "  safe as partials; resubmit the SAME run id to resume:"
      echo "    RUN_ID=${RUN_ID} bash slurm/launch_subset.sh <regime...>"
      echo "  If a whole regime times out repeatedly, raise SAFETY_FACTOR or lower"
      echo "  H2H_TARGET_SECS in slurm/units.R -- the cost model under-states the"
      echo "  crossfit-heavy regimes (R2, R5) by construction."
      grep -rl "received SIGTERM" "${LOG_DIR}" 2>/dev/null | head -n 10 | sed 's/^/    /' || true
    else
      echo "  none"
    fi

    echo
    echo "--- R-level errors in the logs ---"
    local n_err
    n_err=$(grep -rl "^Error" "${LOG_DIR}" 2>/dev/null | wc -l | tr -d ' ' || true)
    if (( n_err > 0 )); then
      echo "  ${n_err} log file(s) contain an R error:"
      grep -rh "^Error" "${LOG_DIR}" 2>/dev/null | sort | uniq -c | sort -rn | head -n 10 | sed 's/^/    /' || true
    else
      echo "  none"
    fi
  else
    echo "  (no log dir yet: ${LOG_DIR})"
  fi

  echo
  echo "--- next step ---"
  if (( total_expected > 0 && total_done >= total_expected )); then
    echo "  All shards present. Combine:"
    echo "    Rscript slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} \\"
    echo "      --pkg-root ${PKG_ROOT} --run-id ${RUN_ID}"
  else
    echo "  Still running. Re-check, or:"
    echo "    Rscript slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} --pkg-root ${PKG_ROOT} --dry-run"
    echo "  (a dry run lists exactly which cells are short and which task ids to resubmit)"
  fi
  echo
}

if [[ "${WATCH}" == "1" ]]; then
  while true; do
    clear
    report_once
    sleep 60
  done
else
  report_once
fi
