#!/bin/bash
# =============================================================================
# slurm/check_progress.sh -- monitor one run of the honest_inference_sparsity_failure sweep
# =============================================================================
# Run ON O2, from anywhere:
#
#   bash slurm/check_progress.sh                      # newest run under scratch
#   RUN_ID=20260909-131500_abc1234 bash slurm/check_progress.sh
#   WATCH=1 bash slurm/check_progress.sh              # re-print every 60 s
#
# Reports, per wave: queued/running task counts from squeue, finished shards on disk
# against the number slurm/units.R says the array should produce, partials in
# flight, any task that was killed at the wall limit, and R-level errors. Then the
# spec-grid view: which of the 15 cells are on disk in the study's results/ dir.
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
STUDY_NAME="honest_inference_sparsity_failure"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

SCRATCH_ROOT="${SCRATCH_ROOT:-/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}}"
HIS_TARGET_SECS="${HIS_TARGET_SECS:-900}"
HIS_REPS="${HIS_REPS:-300}"
HIS_REPS_MAX_N="${HIS_REPS_MAX_N:-150}"
WATCH="${WATCH:-0}"
export HIS_TARGET_SECS HIS_REPS HIS_REPS_MAX_N

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
  echo " ${STUDY_NAME} (S2)   run id: ${RUN_ID}"
  echo " scratch: ${SCRATCH_DIR}"
  echo " gates  : HIS_TARGET_SECS=${HIS_TARGET_SECS} HIS_REPS=${HIS_REPS} HIS_REPS_MAX_N=${HIS_REPS_MAX_N}"
  echo " time   : $(date '+%F %T')"
  echo "=============================================================="

  # --- SLURM's view ----------------------------------------------------------
  echo
  echo "--- squeue ---"
  if command -v squeue >/dev/null 2>&1; then
    squeue -u "${USER}" -o "%.12i %.34j %.10T %.10M %.10l %.6D %R" || true
    local n_run n_pend
    n_run=$(squeue -u "${USER}" -h -t RUNNING -n "${STUDY_NAME}-n500,${STUDY_NAME}-n2000,${STUDY_NAME}-n8000" 2>/dev/null | wc -l | tr -d ' ' || true)
    n_pend=$(squeue -u "${USER}" -h -t PENDING -n "${STUDY_NAME}-n500,${STUDY_NAME}-n2000,${STUDY_NAME}-n8000" 2>/dev/null | wc -l | tr -d ' ' || true)
    echo "  running: ${n_run}   pending: ${n_pend}"
  else
    echo "  (squeue not on PATH -- not an O2 login node?)"
  fi

  # --- shards on disk vs the array size ------------------------------------
  # Expected counts come from the same table the workers use, so "13/20" cannot
  # drift from what was actually submitted. Only the waves this run actually has a
  # scratch subdir for are reported.
  echo
  echo "--- shards on disk vs slurm/units.R (target ${HIS_TARGET_SECS} s/task) ---"
  local expected_file present_waves
  expected_file="$(mktemp)"
  present_waves="$(find "${SCRATCH_DIR}" -maxdepth 1 -type d -name 'n[0-9]*' -exec basename {} \; 2>/dev/null | sort | paste -sd, - || true)"
  if [[ -z "${present_waves}" ]]; then
    echo "  (no wave subdirectories under ${SCRATCH_DIR} yet)"
    rm -f "${expected_file}"
  else
    Rscript --no-init-file "${SLURM_DIR}/print_sizing.R" --format env \
      --waves "${present_waves}" --pkg-root "${PKG_ROOT}" \
      > "${expected_file}" 2>/dev/null || { echo "  (could not compute expected task counts)"; }

    local total_done=0 total_expected=0
    printf "  %-7s %8s %8s %9s %9s\n" wave shards expected pct partials
    while read -r WAVE N_TASKS _WALLTIME _MEM _PART; do
      [[ -z "${WAVE}" ]] && continue
      local d done part pct
      d="${SCRATCH_DIR}/${WAVE}"
      done=0; part=0
      if [[ -d "${d}" ]]; then
        done=$(find "${d}" -maxdepth 1 -name 'batch_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
        part=$(find "${d}" -maxdepth 1 -name 'part_*.rds' 2>/dev/null | wc -l | tr -d ' ' || true)
      fi
      pct=$(( N_TASKS > 0 ? 100 * done / N_TASKS : 0 ))
      printf "  %-7s %8s %8s %8s%% %9s\n" "${WAVE}" "${done}" "${N_TASKS}" "${pct}" "${part}"
      total_done=$(( total_done + done ))
      total_expected=$(( total_expected + N_TASKS ))
    done < "${expected_file}"
    rm -f "${expected_file}"
    if (( total_expected > 0 )); then
      printf "  %-7s %8s %8s %8s%%\n" TOTAL "${total_done}" "${total_expected}" \
        "$(( 100 * total_done / total_expected ))"
    fi
  fi

  # --- wall-limit kills -----------------------------------------------------
  # `|| true` on every counting pipeline below is load-bearing, not decoration:
  # under `set -o pipefail` a grep/find that matches nothing returns 1, the pipeline
  # inherits it, and `set -e` would abort this report mid-print. "No timeouts yet" is
  # the normal case, so without it the monitor would only survive runs that had
  # already gone wrong.
  echo
  echo "--- wall-limit kills (SIGTERM sentinel in the logs) ---"
  if [[ -d "${LOG_DIR}" ]]; then
    local n_term
    n_term=$(grep -rl "received SIGTERM" "${LOG_DIR}" 2>/dev/null | wc -l | tr -d ' ' || true)
    if (( n_term > 0 )); then
      echo "  ${n_term} task(s) hit their wall limit. Their finished replications are"
      echo "  safe as partials; resubmit the SAME run id to resume:"
      echo "    RUN_ID=${RUN_ID} bash slurm/launch_subset.sh <wave...>"
      echo "  If a whole wave times out repeatedly, raise SAFETY_FACTOR or lower"
      echo "  HIS_TARGET_SECS in slurm/units.R. Note the cost anchors there were"
      echo "  measured on the DEV BOX, not on an O2 node: a slower core is the most"
      echo "  likely cause of a systematic over-run, and re-measuring one anchor is"
      echo "  the fix, not raising the safety factor forever."
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
      echo
      echo "  If the message mentions DEV-LOADED packages, the tracked .Rprofile at the"
      echo "  package root ran: that is assert_installed_load() doing its job. The task"
      echo "  must be launched with 'Rscript --no-init-file' (run_simulations.slurm does)."
    else
      echo "  none"
    fi

    echo
    echo "--- dev-load guard trips ---"
    local n_dev
    n_dev=$(grep -rl "DEV-LOADED" "${LOG_DIR}" 2>/dev/null | wc -l | tr -d ' ' || true)
    echo "  ${n_dev} task(s) refused to run on dev-loaded packages (0 is correct)"
  else
    echo "  (no log dir yet: ${LOG_DIR})"
  fi

  # --- the spec grid, in the study's own results/ dir ------------------------
  # The deployment's actual goal. Shards in scratch are means; these 15 files are
  # the end, and analyze.R aborts unless every one of them exists.
  echo
  echo "--- spec grid in results/ (what analyze.R requires) ---"
  Rscript --no-init-file -e '
    setwd(commandArgs(trailingOnly = TRUE)[[1]]); Sys.setenv(HIS_USE_INSTALLED = "1")
    suppressMessages(source(file.path("simulations","honest_inference_sparsity_failure","code","common.R")))
    g <- design_grid()
    g$path <- vapply(seq_len(nrow(g)),
                     function(i) cell_path(g$dgp[[i]], g$n[[i]], g$reps[[i]]), character(1))
    g$on_disk <- file.exists(g$path)
    for (nn in sort(unique(g$n))) {
      s <- g[g$n == nn, ]
      cat(sprintf("  n = %-5d %d/%d cells (R = %d)%s\n", nn, sum(s$on_disk), nrow(s),
                  s$reps[[1]],
                  if (all(s$on_disk)) "" else paste0("   missing: ",
                    paste(s$dgp[!s$on_disk], collapse = ", "))))
    }
    cat(sprintf("  TOTAL     %d/%d\n", sum(g$on_disk), nrow(g)))
  ' "${PKG_ROOT}" 2>/dev/null || echo "  (could not read the results dir)"

  echo
  echo "--- next step ---"
  if (( ${total_expected:-0} > 0 && ${total_done:-0} >= ${total_expected:-1} )); then
    echo "  All shards present. Combine:"
    echo "    Rscript --no-init-file slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} \\"
    echo "      --pkg-root ${PKG_ROOT} --run-id ${RUN_ID}"
  else
    echo "  Still running. Re-check, or:"
    echo "    Rscript --no-init-file slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} --pkg-root ${PKG_ROOT} --dry-run"
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
