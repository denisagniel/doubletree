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

METHODS=(full crossfit doubletree dt_averaged msplit msplit_averaged single_tree)

# Require every method's sizing file up front (fail loudly, no partial submit).
for m in "${METHODS[@]}"; do
  if [[ ! -f "${CONFIG_DIR}/sizing_${m}.env" ]]; then
    echo "ERROR: ${CONFIG_DIR}/sizing_${m}.env not found. Run profile_per_method.R first." >&2
    exit 1
  fi
done

# =============================================================================
# PREFLIGHT (S3): fail loudly BEFORE launching thousands of tasks.
# =============================================================================
preflight_fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 1; }

# (a) Required staged input files. This study's DGP is self-contained (generate_data
# simulates from parameters), so no external inputs are listed. EDIT THIS if the
# tasks ever read a staged .rds/.csv from the cluster filesystem.
REQUIRED_INPUTS=(
  # "${STUDY_DIR}/inputs/real_data.rds"
)
# ${arr[@]+"${arr[@]}"} expands safely even when the array is empty under set -u.
for f in ${REQUIRED_INPUTS[@]+"${REQUIRED_INPUTS[@]}"}; do
  [[ -e "${f}" ]] || preflight_fail "required input not found: ${f}"
done

# (b) Installed-package freshness. The agent cannot `git push` from the dev box, so
# a forgotten R CMD INSTALL means the cluster silently runs STALE package code. The
# tasks require doubletree (the estimator under study) and optimaltrees; verify both
# are installed on this node before submitting.
module load gcc/14.2.0 2>/dev/null || module load gcc || true
module load R/4.4.2   2>/dev/null || module load R
Rscript -e '
  pkgs <- c("doubletree", "optimaltrees")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(sprintf("project package(s) NOT installed on this node: %s -- R CMD INSTALL them first (agent cannot git push).",
                 paste(miss, collapse = ", ")))
  for (p in pkgs) cat(sprintf("preflight: %s %s OK\n", p, utils::packageVersion(p)))
' || preflight_fail "required package(s) missing (see message above)."
echo "Preflight OK."

GIT_SHA="$(git -C "${STUDY_DIR}" rev-parse --short HEAD 2>/dev/null || echo nogit)"
# RUN_ID is normally minted fresh. To RESUME/backfill an existing run (re-run only
# the tasks whose scratch task_*.rds is missing -- e.g. after timeouts), export
# RUN_ID=<existing-run-id> before calling; array.slurm skips tasks whose file exists.
RUN_ID="${RUN_ID:-$(date '+%Y%m%d-%H%M%S')_${GIT_SHA}}"
# SCRATCH_ROOT defaults to O2 scratch; overridable (e.g. for a local pipeline smoke
# test with a fake sbatch, where /n/scratch is not writable). Production leaves it unset.
SCRATCH_ROOT="${SCRATCH_ROOT:-/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}}"
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
        CONCURRENCY_CAP WALLTIME MEM_GB UNIT_OFFSET N_UNITS_METHOD METHOD PARTITION
  # shellcheck disable=SC1090
  source "${CONFIG_DIR}/sizing_${m}.env"
  : "${TOTAL_TASKS:?}" "${REPS_PER_JOB:?}" "${MAX_ARRAY_SIZE:?}" \
    "${MAX_CONCURRENT_JOBS:?}" "${CONCURRENCY_CAP:?}" "${WALLTIME:?}" "${MEM_GB:?}" \
    "${UNIT_OFFSET:?}" "${N_UNITS_METHOD:?}"
  # PARTITION is optional: older sizing envs (pre target-tasks) omit it. Default to
  # `short` so those still submit unchanged; the target-tasks path always writes it.
  PARTITION="${PARTITION:-short}"

  MAX_UNIT=$(( UNIT_OFFSET + N_UNITS_METHOD ))   # last global unit for this method
  METHOD_SCRATCH="${SCRATCH_DIR}/${m}"           # per-method subdir: no task-file collision
  mkdir -p "${METHOD_SCRATCH}"

  echo "--- ${m}: ${TOTAL_TASKS} tasks x ${REPS_PER_JOB} units  --time ${WALLTIME} --mem ${MEM_GB}G -p ${PARTITION} ---"

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
      --partition="${PARTITION}" \
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
