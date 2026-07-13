#!/bin/bash
# =============================================================================
# submit_per_method.sh -- submit ONE SLURM array PER METHOD, each sized for that
# method's cost (produced by profile_per_method.R).
# =============================================================================
# WHY: the seven methods differ ~50x in cost; a single global sizing would time
# out the slow M-split tasks. Each method occupies a contiguous unit block; we
# submit one array per method with method-appropriate --time/--mem, into a
# per-method scratch SUBDIR so task_NNNNNN.rds files never collide across methods.
#
# Run ON O2 from the study dir, AFTER: Rscript slurm/profile_per_method.R
#   bash slurm/submit_per_method.sh
#
# Mirrors submit.sh's run-id, scratch/home discipline, code-hash guard, chunking
# (<=1000/array) and wave throttling (<=10000 queued), but per method.
# =============================================================================

set -euo pipefail

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="six-approach-arbitration"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
CONFIG_DIR="${STUDY_DIR}/config"

METHODS=(doubletree dt_averaged single_tree)

# Require every method's sizing file up front (fail loudly, no partial submit).
for m in "${METHODS[@]}"; do
  if [[ ! -f "${CONFIG_DIR}/sizing_${m}.env" ]]; then
    echo "ERROR: ${CONFIG_DIR}/sizing_${m}.env not found. Run profile_per_method.R first." >&2
    exit 1
  fi
done

GIT_SHA="$(git -C "${STUDY_DIR}" rev-parse --short HEAD 2>/dev/null || echo nogit)"
# RUN_ID is normally minted fresh. To RESUME/backfill an existing run (re-run only
# the tasks whose scratch task_*.rds is missing -- e.g. after timeouts), export
# RUN_ID=<existing-run-id> before calling; array.slurm skips tasks whose file exists.
RUN_ID="${RUN_ID:-$(date '+%Y%m%d-%H%M%S')_${GIT_SHA}}"
SCRATCH_ROOT="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}"
SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${SCRATCH_DIR}" "${LOG_DIR}"

# Code-hash guard (byte-identical to combine.R's code_hash()).
Rscript -e "sd<-'${STUDY_DIR}'; \
  files<-c('config/grid.R','R/dgp.R','R/estimators.R','R/run_one.R'); \
  h<-0; for(f in files){bytes<-readBin(file.path(sd,f),'raw',n=file.size(file.path(sd,f))); \
  for(b in as.integer(bytes)){h<-(h*257+b)%%2147483647}}; \
  cat(sprintf('%.0f', h))" > "${SCRATCH_DIR}/GRID_HASH"

mkdir -p "${STUDY_DIR}/logs"
ln -sfn "${LOG_DIR}" "${STUDY_DIR}/logs/latest"

echo "=============================================================="
echo " study  : ${STUDY_NAME}   run-id: ${RUN_ID}"
echo " scratch: ${SCRATCH_DIR}"
echo "=============================================================="

declare -a ALL_JOB_IDS=()

for m in "${METHODS[@]}"; do
  # Fresh env per method.
  unset TOTAL_TASKS REPS_PER_JOB MAX_ARRAY_SIZE MAX_CONCURRENT_JOBS \
        CONCURRENCY_CAP WALLTIME MEM_GB UNIT_OFFSET N_UNITS_METHOD METHOD
  # shellcheck disable=SC1090
  source "${CONFIG_DIR}/sizing_${m}.env"
  : "${TOTAL_TASKS:?}" "${REPS_PER_JOB:?}" "${MAX_ARRAY_SIZE:?}" \
    "${MAX_CONCURRENT_JOBS:?}" "${CONCURRENCY_CAP:?}" "${WALLTIME:?}" "${MEM_GB:?}" \
    "${UNIT_OFFSET:?}" "${N_UNITS_METHOD:?}"

  MAX_UNIT=$(( UNIT_OFFSET + N_UNITS_METHOD ))   # last global unit for this method
  METHOD_SCRATCH="${SCRATCH_DIR}/${m}"           # per-method subdir: no task-file collision
  mkdir -p "${METHOD_SCRATCH}"

  echo "--- ${m}: ${TOTAL_TASKS} tasks x ${REPS_PER_JOB} units  --time ${WALLTIME} --mem ${MEM_GB}G ---"

  ARRAYS_PER_WAVE=$(( MAX_CONCURRENT_JOBS / MAX_ARRAY_SIZE )); (( ARRAYS_PER_WAVE < 1 )) && ARRAYS_PER_WAVE=1
  offset=0; arrays_in_wave=0; prev_wave_last=""
  while (( offset < TOTAL_TASKS )); do
    remaining=$(( TOTAL_TASKS - offset ))
    chunk=$(( remaining < MAX_ARRAY_SIZE ? remaining : MAX_ARRAY_SIZE ))
    cap=$(( CONCURRENCY_CAP < chunk ? CONCURRENCY_CAP : chunk ))
    dep_args=()
    if (( arrays_in_wave == 0 )) && [[ -n "${prev_wave_last}" ]]; then
      dep_args=(--dependency=afterany:"${prev_wave_last}")
    fi
    jobid=$(sbatch --parsable \
      --job-name="${STUDY_NAME}-${m}" \
      --array=1-"${chunk}"%"${cap}" \
      --time="${WALLTIME}" --mem="${MEM_GB}G" \
      --output="${LOG_DIR}/${m}_%A_%a.out" \
      --error="${LOG_DIR}/${m}_%A_%a.err" \
      "${dep_args[@]}" \
      --export=ALL,STUDY_DIR="${STUDY_DIR}",SCRATCH_DIR="${METHOD_SCRATCH}",REPS_PER_JOB="${REPS_PER_JOB}",ARRAY_OFFSET="${offset}",UNIT_OFFSET="${UNIT_OFFSET}",MAX_UNIT="${MAX_UNIT}" \
      "${SLURM_DIR}/array.slurm")
    ALL_JOB_IDS+=("${jobid}")
    echo "    array ${jobid}: tasks $((offset+1))-$((offset+chunk)) (%${cap})"
    offset=$(( offset + chunk )); arrays_in_wave=$(( arrays_in_wave + 1 ))
    if (( arrays_in_wave >= ARRAYS_PER_WAVE )); then
      prev_wave_last="${jobid}"; arrays_in_wave=0
    fi
  done
done

# MANIFEST
{
  echo "# Run Manifest -- ${STUDY_NAME} (per-method)"
  echo
  echo "- run-id: \`${RUN_ID}\`   git SHA: \`${GIT_SHA}\`"
  echo "- grid hash: \`$(cat "${SCRATCH_DIR}/GRID_HASH")\`"
  echo "- submitted: $(date '+%F %T %Z')"
  echo "- scratch dir: \`${SCRATCH_DIR}\`  (per-method subdirs)"
  echo "- SLURM job ids: ${ALL_JOB_IDS[*]}"
  echo
  echo "## Combine (after all arrays finish)"
  echo '```bash'
  echo "Rscript slurm/combine.R --run-id ${RUN_ID} \\"
  echo "  --scratch-dir ${SCRATCH_DIR} --study-dir ${STUDY_DIR}"
  echo '```'
} > "${STUDY_DIR}/MANIFEST.md"

echo "Submitted ${#ALL_JOB_IDS[@]} array job(s) across ${#METHODS[@]} methods."
echo "Monitor: bash slurm/monitor.sh   Manifest: ${STUDY_DIR}/MANIFEST.md"
