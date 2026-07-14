#!/bin/bash
# =============================================================================
# submit.sh -- submit the "six-approach-arbitration" study as chunked, throttled SLURM arrays
# =============================================================================
# Run ON O2 from the study directory after profiling. It:
#   1. Reads config/sizing.env (produced by profile_timing.R).
#   2. Mints a run-id: <timestamp>_<gitSHA>  (isolates this run's scratch/logs).
#   3. Creates a run-specific scratch dir and records the grid hash there.
#   4. Splits TOTAL_TASKS into arrays of <= MAX_ARRAY_SIZE (1000).
#   5. Throttles so no more than MAX_CONCURRENT_JOBS (10000) are queued at once,
#      submitting in WAVES with dependencies when the total would exceed it.
#   6. Writes MANIFEST.md and points logs/latest at this run's log dir.
#
# Usage:  bash slurm/submit.sh
# =============================================================================

set -euo pipefail

# --- Identity (filled by the skill) ------------------------------------------
HMS_ID="dma12"                         # e.g. dma12
PROJECT_NAME="global-scholars"
STUDY_NAME="six-approach-arbitration"

# Resolve the study dir as the parent of this script's slurm/ dir.
SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"

# --- Load sizing --------------------------------------------------------------
SIZING_ENV="${STUDY_DIR}/config/sizing.env"
if [[ ! -f "${SIZING_ENV}" ]]; then
  echo "ERROR: ${SIZING_ENV} not found. Run profile_timing.R locally first." >&2
  exit 1
fi
# shellcheck disable=SC1090
source "${SIZING_ENV}"
: "${TOTAL_TASKS:?}" "${REPS_PER_JOB:?}" "${MAX_ARRAY_SIZE:?}" \
  "${MAX_CONCURRENT_JOBS:?}" "${CONCURRENCY_CAP:?}" "${WALLTIME:?}" "${MEM_GB:?}"

# =============================================================================
# PREFLIGHT (S3): fail loudly BEFORE launching thousands of tasks.
# =============================================================================
preflight_fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 1; }

# (a) Required staged input files. This study's DGP is self-contained, so none are
# listed. EDIT THIS if tasks ever read a staged .rds/.csv from the cluster FS.
REQUIRED_INPUTS=(
  # "${STUDY_DIR}/inputs/real_data.rds"
)
# ${arr[@]+"${arr[@]}"} expands safely even when the array is empty under set -u.
for f in ${REQUIRED_INPUTS[@]+"${REQUIRED_INPUTS[@]}"}; do
  [[ -e "${f}" ]] || preflight_fail "required input not found: ${f}"
done

# (b) Installed-package freshness. The agent cannot `git push` from the dev box, so
# a forgotten R CMD INSTALL means the cluster silently runs STALE package code. The
# tasks require doubletree and optimaltrees; verify both are installed on this node.
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

# --- Run identity -------------------------------------------------------------
GIT_SHA="$(git -C "${STUDY_DIR}" rev-parse --short HEAD 2>/dev/null || echo nogit)"
RUN_ID="$(date '+%Y%m%d-%H%M%S')_${GIT_SHA}"

# Scratch: intermediate task files + logs (NOT backed up). Home holds only the
# final combined result (written later by combine.R).
SCRATCH_ROOT="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}"
SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${SCRATCH_DIR}" "${LOG_DIR}"

# Record a hash of the study SOURCE FILES so combine.R can detect stale code.
# Must produce byte-identical output to code_hash() in combine.R (same files, order).
Rscript -e "sd<-'${STUDY_DIR}'; \
  files<-c('config/grid.R','R/dgp.R','R/estimators.R','R/run_one.R'); \
  h<-0; for(f in files){bytes<-readBin(file.path(sd,f),'raw',n=file.size(file.path(sd,f))); \
  for(b in as.integer(bytes)){h<-(h*257+b)%%2147483647}}; \
  cat(sprintf('%.0f', h))" \
  > "${SCRATCH_DIR}/GRID_HASH"

# Point logs/latest at this run (convenient log discovery from the study dir).
mkdir -p "${STUDY_DIR}/logs"
ln -sfn "${LOG_DIR}" "${STUDY_DIR}/logs/latest"

echo "=============================================================="
echo " study        : ${STUDY_NAME}   (project ${PROJECT_NAME})"
echo " run-id       : ${RUN_ID}"
echo " total tasks  : ${TOTAL_TASKS}  (${REPS_PER_JOB} units/task)"
echo " scratch      : ${SCRATCH_DIR}"
echo " logs         : ${LOG_DIR}   (also ./logs/latest)"
echo " --time ${WALLTIME}  --mem ${MEM_GB}G  concurrency %${CONCURRENCY_CAP}"
echo "=============================================================="

# --- Chunk into arrays of <= MAX_ARRAY_SIZE, submit in throttled waves --------
# Global task ids run 1..TOTAL_TASKS. Each array covers a contiguous block; the
# array's local index 1..CHUNK maps to global id via ARRAY_OFFSET.
# To respect MAX_CONCURRENT_JOBS, arrays beyond the first wave are chained with
# --dependency so the queued count never exceeds the cap.

# How many arrays fit in one wave without exceeding the concurrent-jobs cap.
ARRAYS_PER_WAVE=$(( MAX_CONCURRENT_JOBS / MAX_ARRAY_SIZE ))
if (( ARRAYS_PER_WAVE < 1 )); then ARRAYS_PER_WAVE=1; fi

submitted=0
offset=0
wave_index=0
arrays_in_wave=0
prev_wave_last_jobid=""
this_wave_last_jobid=""
declare -a JOB_IDS=()

while (( offset < TOTAL_TASKS )); do
  remaining=$(( TOTAL_TASKS - offset ))
  chunk=$(( remaining < MAX_ARRAY_SIZE ? remaining : MAX_ARRAY_SIZE ))

  # Concurrency within the array: never exceed the chunk size.
  cap=$(( CONCURRENCY_CAP < chunk ? CONCURRENCY_CAP : chunk ))

  dep_args=()
  if (( arrays_in_wave == 0 )) && [[ -n "${prev_wave_last_jobid}" ]]; then
    # First array of a new wave waits for the previous wave to finish.
    dep_args=(--dependency=afterany:"${prev_wave_last_jobid}")
  fi

  out_pat="${LOG_DIR}/task_%A_%a.out"
  err_pat="${LOG_DIR}/task_%A_%a.err"

  jobid=$(sbatch --parsable \
    --array=1-"${chunk}"%"${cap}" \
    --time="${WALLTIME}" \
    --mem="${MEM_GB}G" \
    --output="${out_pat}" \
    --error="${err_pat}" \
    "${dep_args[@]}" \
    --export=ALL,STUDY_DIR="${STUDY_DIR}",SCRATCH_DIR="${SCRATCH_DIR}",REPS_PER_JOB="${REPS_PER_JOB}",ARRAY_OFFSET="${offset}" \
    "${SLURM_DIR}/array.slurm")

  JOB_IDS+=("${jobid}")
  this_wave_last_jobid="${jobid}"
  echo "  submitted array job ${jobid}: global tasks $((offset+1))-$((offset+chunk)) (%${cap})"

  offset=$(( offset + chunk ))
  submitted=$(( submitted + chunk ))
  arrays_in_wave=$(( arrays_in_wave + 1 ))

  if (( arrays_in_wave >= ARRAYS_PER_WAVE )); then
    prev_wave_last_jobid="${this_wave_last_jobid}"
    wave_index=$(( wave_index + 1 ))
    arrays_in_wave=0
    echo "  --- wave ${wave_index} full (${ARRAYS_PER_WAVE} arrays); next wave waits on ${prev_wave_last_jobid} ---"
  fi
done

echo "Submitted ${submitted} tasks across ${#JOB_IDS[@]} array job(s)."

# --- Write MANIFEST.md --------------------------------------------------------
{
  echo "# Run Manifest -- ${STUDY_NAME}"
  echo
  echo "- run-id: \`${RUN_ID}\`"
  echo "- git SHA: \`${GIT_SHA}\`"
  echo "- grid hash: \`$(cat "${SCRATCH_DIR}/GRID_HASH")\`"
  echo "- submitted: $(date '+%F %T %Z')"
  echo "- total tasks: ${TOTAL_TASKS} (${REPS_PER_JOB} units/task, ${TOTAL_UNITS:-?} units)"
  echo "- --time: ${WALLTIME}   --mem: ${MEM_GB}G   concurrency: %${CONCURRENCY_CAP}"
  echo "- scratch dir: \`${SCRATCH_DIR}\`"
  echo "- log dir: \`${LOG_DIR}\` (also \`./logs/latest\`)"
  echo "- SLURM job ids: ${JOB_IDS[*]}"
  echo
  echo "## Next steps"
  echo '```bash'
  echo "bash slurm/monitor.sh                 # watch progress / find failed tasks"
  echo "Rscript slurm/combine.R --run-id ${RUN_ID} \\"
  echo "  --scratch-dir ${SCRATCH_DIR} --study-dir ${STUDY_DIR}"
  echo "bash slurm/clean.sh                   # remove superseded runs when done"
  echo '```'
} > "${STUDY_DIR}/MANIFEST.md"

echo "Wrote ${STUDY_DIR}/MANIFEST.md"
echo "Monitor with: bash slurm/monitor.sh"
