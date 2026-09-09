#!/bin/bash
# =============================================================================
# slurm/launch_subset.sh -- submit ONE SLURM array PER WAVE, for a chosen subset
# =============================================================================
# THE WORKHORSE LAUNCHER. launch_all_simulations.sh is a one-line wrapper that
# calls this with the waves that are actually missing; there is only one copy of
# the submit logic.
#
# Run ON O2, from anywhere:
#   bash slurm/launch_subset.sh n2000            # the cheaper missing wave first
#   bash slurm/launch_subset.sh n8000            # the expensive one on its own
#   bash slurm/launch_subset.sh                  # n2000 n8000 (same as launch_all)
#   bash slurm/launch_subset.sh n500             # backfill/re-run an n = 500 cell
#
# Environment gates (declared here, in one block, with defaults, and echoed):
#   RUN_ID            resume/backfill an existing run instead of minting a new id
#   SCRATCH_ROOT      override the O2 scratch root (used by the local pipeline test)
#   HIS_TARGET_SECS   target wall cost per array task, seconds      (default 900)
#   HIS_REPS          replications per cell for n < 8000            (default 300)
#   HIS_REPS_MAX_N    replications per cell at n = 8000             (default 150)
#   CONCURRENCY_CAP   max simultaneously-running tasks per array    (default 200)
#   DRY_RUN           "1" prints the sbatch lines without submitting
#
# HIS_REPS / HIS_REPS_MAX_N are exported rather than merely documented: slurm/units.R
# builds its map from code/common.R's reps_for_n(), and combine_results.R writes to
# cell_<dgp>_n<n>_r<reps>.rds. If the launcher and the combine step disagreed about R
# the shards would assemble into a cell whose FILENAME lied about its replication
# count, and code/analyze.R globs by that exact name.
#
# ONE ARRAY PER WAVE, not one array for the sweep. n changes the per-replication cost
# 3x (3.83 s at n = 2000 against 11.25 s at n = 8000, measured 2026-09-09) and the
# replication count from 300 to 150, so a single --time/--mem sizing would either
# time out the n = 8000 tasks or reserve an hour for n = 2000 tasks that finish in
# ten minutes. Per-wave arrays also make the task-id -> block map offset-free: see
# slurm/units.R.
# =============================================================================

set -euo pipefail

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="honest_inference_sparsity_failure"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

# --- 0 gates, in one block, echoed -------------------------------------------
HIS_TARGET_SECS="${HIS_TARGET_SECS:-900}"
HIS_REPS="${HIS_REPS:-300}"
HIS_REPS_MAX_N="${HIS_REPS_MAX_N:-150}"
CONCURRENCY_CAP="${CONCURRENCY_CAP:-200}"
DRY_RUN="${DRY_RUN:-0}"
export HIS_TARGET_SECS HIS_REPS HIS_REPS_MAX_N

WAVES=("$@")
if [[ ${#WAVES[@]} -eq 0 ]]; then
  # The two waves the spec asks for that are not already on disk. n500 is a valid
  # argument but never a default, so a backfill is possible and never accidental.
  WAVES=(n2000 n8000)
fi

echo "=============================================================="
echo " study   : ${STUDY_NAME}  (S2)"
echo " spec    : quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md"
echo " pkg root: ${PKG_ROOT}"
echo " waves   : ${WAVES[*]}"
echo " gates   : HIS_TARGET_SECS=${HIS_TARGET_SECS} HIS_REPS=${HIS_REPS} HIS_REPS_MAX_N=${HIS_REPS_MAX_N}"
echo "           CONCURRENCY_CAP=${CONCURRENCY_CAP} DRY_RUN=${DRY_RUN}"
echo "=============================================================="

# --- 1 preflight: fail loudly BEFORE submitting anything ---------------------
preflight_fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 1; }

[[ -f "${PKG_ROOT}/DESCRIPTION" ]] || preflight_fail "no DESCRIPTION at ${PKG_ROOT}"

# Guarded on `command -v module` so DRY_RUN=1 can be exercised on a dev box, where
# no module system exists. On O2 the loads run normally, and a genuinely missing R
# module still fails loudly at the Rscript preflight below rather than being
# swallowed here.
if command -v module >/dev/null 2>&1; then
  module load gcc/14.2.0 2>/dev/null || module load gcc || true
  module load R/4.4.2   2>/dev/null || module load R   || true
fi

# Installed-package freshness. A forgotten `R CMD INSTALL` after a git pull means
# every task silently runs STALE estimator code and the whole wave has to be thrown
# away; this is the cheapest possible place to catch it. optimaltrees >= 0.4.1 is
# checked explicitly because this study was BLOCKED on exactly that requirement from
# 2026-09-04 to 2026-09-09 (doubletree's DESCRIPTION carries the constraint, but a
# manually-installed sibling can satisfy `requireNamespace` at 0.4.0).
#
# --no-init-file for the same reason run_simulations.slurm uses it: the tracked
# .Rprofile at the package root would dev-load optimaltrees from source and this
# check would then report the SOURCE version, not the installed one.
Rscript --no-init-file -e '
  need <- c(doubletree = "0.0.0.9000", optimaltrees = "0.4.1")
  for (p in names(need)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("package %s is NOT installed on this node -- R CMD INSTALL it first (README_O2.md, prerequisites).", p))
    }
    v <- utils::packageVersion(p)
    if (v < need[[p]]) {
      stop(sprintf("%s %s is installed but this study needs >= %s.", p, v, need[[p]]))
    }
    cat(sprintf("preflight: %s %s OK (%s)\n", p, v, dirname(find.package(p))))
  }
' || preflight_fail "required package(s) missing or too old (see message above)."

# The DGP gate. Cheap (the enumeration at leaf_budget = 2 is 6 partitions) and it is
# the one construction the whole study rests on: the blind-spot variant must have a
# bilinear inner product of essentially zero -- so prop:spectest's documented blind
# spot is actually being exercised -- while the severe variant must have a nonzero
# predicted bias, or there is nothing for the anchor interval to guard against. A
# broken DGP here would produce 10 perfectly well-formed cells that measure nothing.
export PKG_ROOT
Rscript --no-init-file -e '
  setwd(Sys.getenv("PKG_ROOT")); Sys.setenv(HIS_USE_INSTALLED = "1")
  suppressMessages(source(file.path("simulations","honest_inference_sparsity_failure","code","common.R")))
  sev  <- make_dgp("eps_sev")
  blind <- make_dgp("orth_blind")
  cat(sprintf("preflight: eps_sev  predicted bias = %.6g (must be nonzero)\n", sev$pseudo$bias))
  cat(sprintf("preflight: orth_blind <g,h>_nu    = %.3e, predicted bias = %.3e (both must be ~0)\n",
              blind$ip_gh$ip, blind$pseudo$bias))
  if (abs(sev$pseudo$bias) < 1e-6)
    stop("eps_sev has no predicted bias; the eps dial is broken.")
  if (abs(blind$ip_gh$ip) > 1e-8)
    stop("orth_blind is not orthogonal under the (1-e_0)-weighted inner product; prop:spectest blind spot not exercised.")
' || preflight_fail "DGP gate failed."
echo "Preflight OK."

# --- 2 run id, scratch, logs --------------------------------------------------
GIT_SHA="$(git -C "${PKG_ROOT}" rev-parse --short HEAD 2>/dev/null || echo nogit)"
# RUN_ID is normally minted fresh. Export an existing RUN_ID to RESUME: every task
# whose shard already exists exits 0 immediately, and every task with a partial
# resumes from its last flush.
RUN_ID="${RUN_ID:-$(date '+%Y%m%d-%H%M%S')_${GIT_SHA}}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}}"
SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
LOG_DIR="${SCRATCH_DIR}/logs"
mkdir -p "${SCRATCH_DIR}" "${LOG_DIR}"

mkdir -p "${STUDY_DIR}/logs"
ln -sfn "${LOG_DIR}" "${STUDY_DIR}/logs/latest"

echo " run id  : ${RUN_ID}   git SHA: ${GIT_SHA}"
echo " scratch : ${SCRATCH_DIR}"
echo " logs    : ${LOG_DIR}  (symlinked as ${STUDY_DIR}/logs/latest)"

# --- 3 sizing, from the SAME table the workers use ---------------------------
# print_sizing.R --format env emits: WAVE N_TASKS WALLTIME MEM_GB PARTITION
SIZING_FILE="${SCRATCH_DIR}/sizing.txt"
Rscript --no-init-file "${SLURM_DIR}/print_sizing.R" --format env \
  --waves "$(IFS=,; echo "${WAVES[*]}")" --pkg-root "${PKG_ROOT}" \
  > "${SIZING_FILE}" || preflight_fail "could not compute sizing."
echo
echo "Sizing (from slurm/units.R, target ${HIS_TARGET_SECS} s/task):"
cat "${SIZING_FILE}"
echo

# --- 4 submit one array per wave ---------------------------------------------
declare -a ALL_JOB_IDS=()
declare -a SUBMIT_LINES=()

while read -r WAVE N_TASKS WALLTIME MEM_GB PARTITION; do
  [[ -z "${WAVE}" ]] && continue
  CAP=$(( CONCURRENCY_CAP < N_TASKS ? CONCURRENCY_CAP : N_TASKS ))
  WAVE_SCRATCH="${SCRATCH_DIR}/${WAVE}"
  mkdir -p "${WAVE_SCRATCH}"

  echo "--- ${WAVE}: ${N_TASKS} task(s)  --time ${WALLTIME} --mem ${MEM_GB}G -p ${PARTITION} (%${CAP}) ---"

  if [[ "${DRY_RUN}" == "1" ]]; then
    SUBMIT_LINES+=("${WAVE}: 1-${N_TASKS}%${CAP} ${WALLTIME} ${MEM_GB}G ${PARTITION}")
    continue
  fi

  jobid=$(sbatch --parsable \
    --job-name="${STUDY_NAME}-${WAVE}" \
    --array=1-"${N_TASKS}"%"${CAP}" \
    --partition="${PARTITION}" \
    --time="${WALLTIME}" \
    --mem="${MEM_GB}G" \
    --output="${LOG_DIR}/${WAVE}_%A_%a.out" \
    --error="${LOG_DIR}/${WAVE}_%A_%a.err" \
    --export=ALL,PKG_ROOT="${PKG_ROOT}",SCRATCH_DIR="${WAVE_SCRATCH}",WAVE="${WAVE}",TARGET_SECS="${HIS_TARGET_SECS}",HIS_REPS="${HIS_REPS}",HIS_REPS_MAX_N="${HIS_REPS_MAX_N}" \
    "${SLURM_DIR}/run_simulations.slurm")
  ALL_JOB_IDS+=("${jobid}")
  echo "    submitted array ${jobid}"
done < "${SIZING_FILE}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo
  echo "DRY RUN -- nothing submitted. Would have submitted:"
  printf '  %s\n' "${SUBMIT_LINES[@]}"
  exit 0
fi

# --- 5 manifest ---------------------------------------------------------------
# One file that says exactly which run produced which numbers, and holds the combine
# command with the run's own paths already substituted -- so recovering a run three
# weeks later is copy-paste, not archaeology.
{
  echo "# Run Manifest -- ${STUDY_NAME} (S2)"
  echo
  echo "- run id: \`${RUN_ID}\`   git SHA: \`${GIT_SHA}\`"
  echo "- submitted: $(date '+%F %T %Z')"
  echo "- waves: ${WAVES[*]}"
  echo "- target seconds per task: ${HIS_TARGET_SECS}"
  echo "- reps per cell: HIS_REPS=${HIS_REPS}, HIS_REPS_MAX_N=${HIS_REPS_MAX_N}"
  echo "- scratch dir: \`${SCRATCH_DIR}\` (one subdir per wave)"
  echo "- log dir: \`${LOG_DIR}\`"
  echo "- SLURM array job ids: ${ALL_JOB_IDS[*]}"
  echo
  echo '## Sizing'
  echo
  echo '```'
  cat "${SIZING_FILE}"
  echo '```'
  echo
  echo '## Monitor'
  echo
  echo '```bash'
  echo "RUN_ID=${RUN_ID} bash slurm/check_progress.sh"
  echo '```'
  echo
  echo '## Combine (after every array finishes)'
  echo
  echo '```bash'
  echo "Rscript --no-init-file slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} \\"
  echo "  --pkg-root ${PKG_ROOT} --run-id ${RUN_ID}"
  echo '```'
  echo
  echo '## Then'
  echo
  echo '```bash'
  echo "Rscript --no-init-file simulations/${STUDY_NAME}/code/analyze.R"
  echo '```'
} > "${STUDY_DIR}/MANIFEST.md"

echo
echo "Submitted ${#ALL_JOB_IDS[@]} array job(s) across ${#WAVES[@]} wave(s)."
echo "Manifest: ${STUDY_DIR}/MANIFEST.md"
echo "Monitor : RUN_ID=${RUN_ID} bash slurm/check_progress.sh"
