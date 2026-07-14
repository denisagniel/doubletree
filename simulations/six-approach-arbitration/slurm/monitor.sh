#!/bin/bash
# =============================================================================
# monitor.sh -- progress, failed-task detection, and log discovery
# =============================================================================
# Run ON O2 from the study directory. With no arguments it reports on the most
# recent run (via ./logs/latest). Pass a run-id to inspect an older run.
#
# Solves "hard to find recent logs": it resolves the current run's scratch log
# dir, shows the newest log files, counts completed vs expected tasks, and tails
# the logs of any FAILED array tasks so you can debug fast.
#
# Usage:
#   bash slurm/monitor.sh                 # latest run
#   bash slurm/monitor.sh <run-id>        # specific run
#   bash slurm/monitor.sh --tail-failures # also print tails of failed task logs
# =============================================================================

set -euo pipefail

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="six-approach-arbitration"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
SCRATCH_ROOT="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}"

TAIL_FAILURES=0
RUN_ID=""
for arg in "$@"; do
  case "${arg}" in
    --tail-failures) TAIL_FAILURES=1 ;;
    *) RUN_ID="${arg}" ;;
  esac
done

# --- Resolve run dir ----------------------------------------------------------
if [[ -n "${RUN_ID}" ]]; then
  SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
  LOG_DIR="${SCRATCH_DIR}/logs"
elif [[ -L "${STUDY_DIR}/logs/latest" ]]; then
  LOG_DIR="$(readlink -f "${STUDY_DIR}/logs/latest")"
  SCRATCH_DIR="$(dirname "${LOG_DIR}")"
  RUN_ID="$(basename "${SCRATCH_DIR}")"
else
  echo "No run-id given and ./logs/latest missing. Recent runs in scratch:" >&2
  ls -1t "${SCRATCH_ROOT}" 2>/dev/null | head -10 >&2 || echo "  (none)" >&2
  exit 1
fi

if [[ ! -d "${SCRATCH_DIR}" ]]; then
  echo "ERROR: scratch dir not found: ${SCRATCH_DIR}" >&2
  exit 1
fi

echo "=============================================================="
echo " run-id : ${RUN_ID}"
echo " scratch: ${SCRATCH_DIR}"
echo " logs   : ${LOG_DIR}"
echo "=============================================================="

# --- Progress: completed task files vs expected -------------------------------
# Recurse: the per-method launcher writes task_*.rds into per-method SUBDIRS
# (SCRATCH_DIR/<method>/), so a -maxdepth 1 count would always be 0. No -maxdepth.
DONE=$(find "${SCRATCH_DIR}" -name 'task_*.rds' 2>/dev/null | wc -l | tr -d ' ')
# EXPECTED = sum of TOTAL_TASKS across per-method sizing files if present (per-method
# run); else the single global sizing.env (stock single-array run).
EXPECTED="?"
per_method_sizes=("${STUDY_DIR}"/config/sizing_*.env)
if [[ -e "${per_method_sizes[0]}" ]]; then
  EXPECTED=$(grep -h '^TOTAL_TASKS=' "${STUDY_DIR}"/config/sizing_*.env | cut -d= -f2 \
             | paste -sd+ - | bc 2>/dev/null || echo "?")
elif [[ -f "${STUDY_DIR}/config/sizing.env" ]]; then
  EXPECTED=$(grep '^TOTAL_TASKS=' "${STUDY_DIR}/config/sizing.env" | cut -d= -f2)
fi
echo "Completed task files: ${DONE} / ${EXPECTED}"

# --- Queue state for this user's jobs ----------------------------------------
# Match by name PREFIX, not exact: the per-method launcher names jobs
# "${STUDY_NAME}-<method>" (e.g. six-approach-arbitration-full), so an exact
# --name=${STUDY_NAME} filter matches nothing. Filter our own rows by prefix.
echo
echo "Queue (squeue) for ${HMS_ID}, jobs matching ${STUDY_NAME}*:"
SQ_FMT="%.18i %.9P %.28j %.8T %.10M %.6D %R"
squeue -u "${HMS_ID}" --format="${SQ_FMT}" 2>/dev/null | awk -v s="${STUDY_NAME}" \
  'NR==1 || index($3, s)==1' || echo "  (squeue unavailable)"
# Count running/pending among rows whose job name starts with STUDY_NAME.
RUNNING=$(squeue -u "${HMS_ID}" -h -t RUNNING --format="%j" 2>/dev/null | grep -c "^${STUDY_NAME}" || true)
PENDING=$(squeue -u "${HMS_ID}" -h -t PENDING --format="%j" 2>/dev/null | grep -c "^${STUDY_NAME}" || true)
echo "Running: ${RUNNING}   Pending: ${PENDING}"

# --- Most recent log files (log discovery) -----------------------------------
echo
echo "Newest log files:"
# find + sort (S4), never `ls *.out` glob: at ~100k logs the glob overflows the
# argv limit / errors, silently truncating log discovery. `-printf '%T@ %p'`
# sorts numerically by mtime with no dependence on locale-formatted `ls` output.
# Matches BOTH log layouts: task_%A_%a.{out,err} (submit.sh) and
# <method>_%A_%a.{out,err} (submit_per_method.sh).
find "${LOG_DIR}" -maxdepth 1 \( -name '*.out' -o -name '*.err' \) -printf '%T@ %p\n' 2>/dev/null \
  | sort -rn | head -8 | while read -r ts f; do
      printf "  %s  %s\n" "$(date -d "@${ts%.*}" '+%F %T' 2>/dev/null || echo '?')" "${f}"
    done || echo "  (no logs yet)"

# --- Failed-task detection (A3: sentinel/exit-code, not grep-the-.err) --------
# array.slurm writes "finished with status N" to the .out on completion. The
# project's tasks `library(doubletree)`/`library(optimaltrees)`, which emit S7/
# jsonlite import WARNINGS into .err on every HEALTHY task -- so grepping .err for
# error|cannot|killed|oom yields constant FALSE POSITIVES. Instead trust the
# finish sentinel and its exit status:
#   FAILED     = definite failure (non-zero exit status in the finish sentinel)
#   INCOMPLETE = no "finished with status 0" sentinel AND no result file yet
#                -> still running, or killed/OOM'd (no chance to print the sentinel)
# Scan BOTH .out layouts (task_*.out and <method>_*.out).
echo
echo "Scanning for failed/incomplete tasks (exit-status sentinels)..."
FAILED=()
INCOMPLETE=()
while IFS= read -r -d '' outf; do
  if grep -q 'received SIGTERM' "${outf}" 2>/dev/null; then
    FAILED+=("${outf}")                                    # wall-time TIMEOUT (array.slurm trap)
  elif grep -qE 'finished with status [1-9]' "${outf}" 2>/dev/null; then
    FAILED+=("${outf}")                                    # non-zero exit
  elif ! grep -q 'finished with status 0' "${outf}" 2>/dev/null; then
    INCOMPLETE+=("${outf}")                                # no sentinel: running OR killed
  fi
done < <(find "${LOG_DIR}" -maxdepth 1 -name '*.out' -print0 2>/dev/null)

if (( ${#INCOMPLETE[@]} > 0 )); then
  echo "  ${#INCOMPLETE[@]} task log(s) have no finish sentinel yet (still running, or killed/OOM if the job is gone from squeue)."
fi

if (( ${#FAILED[@]} == 0 )); then
  echo "  No failed tasks detected via exit sentinels."
else
  echo "  ${#FAILED[@]} task log(s) failed or timed out:"
  for f in "${FAILED[@]}"; do echo "    ${f}"; done
  if (( TAIL_FAILURES == 1 )); then
    echo
    echo "---- tails of failed task logs ----"
    for f in "${FAILED[@]}"; do
      echo ">>> ${f}"
      tail -n 15 "${f}"
      echo
    done
  else
    echo "  (re-run with --tail-failures to see the tails)"
  fi
fi

echo
if [[ "${DONE}" == "${EXPECTED}" ]]; then
  echo "All tasks complete. Combine with:"
  echo "  Rscript slurm/combine.R --run-id ${RUN_ID} --scratch-dir ${SCRATCH_DIR} --study-dir ${STUDY_DIR}"
fi
