#!/bin/bash
# =============================================================================
# slurm/launch_subset.sh -- submit ONE SLURM array PER REGIME, for a chosen subset
# =============================================================================
# THE WORKHORSE LAUNCHER. launch_all_simulations.sh is a one-line wrapper that
# calls this with every regime; there is only one copy of the submit logic.
#
# Run ON O2, from anywhere:
#   bash slurm/launch_subset.sh R1 R4 R6        # the three cheap regimes first
#   bash slurm/launch_subset.sh R2              # the expensive one on its own
#   bash slurm/launch_subset.sh                 # all six (same as launch_all)
#
# Environment gates (declared here, in one block, with defaults, and echoed):
#   RUN_ID            resume/backfill an existing run instead of minting a new id
#   SCRATCH_ROOT      override the O2 scratch root (used by the local pipeline test)
#   H2H_TARGET_SECS   target wall cost per array task, seconds     (default 1800)
#   CONCURRENCY_CAP   max simultaneously-running tasks per array   (default 200)
#   DRY_RUN           "1" prints the sbatch lines without submitting
#
# ONE ARRAY PER REGIME, not one array for the sweep. The regimes differ ~350x in
# per-replication cost (R6 at 0.019 s/rep vs R5 at 6.53 s/rep, both measured in the
# 2026-09-09 pilot), so a single --time/--mem sizing would either time out R2/R5 or
# reserve 1.7 h and 8 GB for R6 tasks that finish in 2 seconds. Per-regime arrays
# also make the task-id -> block map offset-free: see slurm/units.R.
# =============================================================================

set -euo pipefail

# --- Login-node environment bootstrap (added 2026-09-23) ----------------------
# This launcher runs on a LOGIN node, and its preflight calls Rscript. `module` is a
# bash function Lmod `export -f`s into the environment, so it is INHERITED: present for
# a human at an interactive prompt, ABSENT under `ssh host 'bash -s'`. Without this,
# agent-driven submission dies at the preflight with "Rscript: command not found" while
# the identical command works when typed by hand. Bootstrap before strict mode --
# /etc/profile.d/* scripts reference unset vars and return non-zero.
set +eu
if ! command -v module >/dev/null 2>&1; then
  for profile_script in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /etc/profile; do
    [ -r "${profile_script}" ] && . "${profile_script}" && break
  done
fi
set -euo pipefail
if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load gcc/14.2.0 2>/dev/null || module load gcc || true
  module load R/4.4.2 2>/dev/null || module load R || true
fi
command -v Rscript >/dev/null 2>&1 || {
  echo "ERROR: Rscript not on PATH after module bootstrap; cannot run preflight." >&2
  exit 1
}
export R_LIBS_USER="${R_LIBS_USER:-${HOME}/R/x86_64-pc-linux-gnu-library/4.4}"
echo " R       : $(command -v Rscript)  R_LIBS_USER=${R_LIBS_USER}"

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="head_to_head_comparison"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

# --- 0 gates, in one block, echoed -------------------------------------------
H2H_TARGET_SECS="${H2H_TARGET_SECS:-1800}"
CONCURRENCY_CAP="${CONCURRENCY_CAP:-200}"
DRY_RUN="${DRY_RUN:-0}"
export H2H_TARGET_SECS

REGIMES=("$@")
if [[ ${#REGIMES[@]} -eq 0 ]]; then
  REGIMES=(R1 R2 R3 R4 R5 R6)
fi

echo "=============================================================="
echo " study   : ${STUDY_NAME}"
echo " pkg root: ${PKG_ROOT}"
echo " regimes : ${REGIMES[*]}"
echo " gates   : H2H_TARGET_SECS=${H2H_TARGET_SECS} CONCURRENCY_CAP=${CONCURRENCY_CAP} DRY_RUN=${DRY_RUN}"
echo "=============================================================="

# --- 1 preflight: fail loudly BEFORE submitting anything ---------------------
preflight_fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 1; }

[[ -f "${PKG_ROOT}/DESCRIPTION" ]] || preflight_fail "no DESCRIPTION at ${PKG_ROOT}"

# Same module pairing as six-approach-arbitration (2026-07). Module R has no
# working pkgload dev-load, which is the whole reason the packages must be
# INSTALLED -- see README_O2.md.
#
# Guarded on `command -v module` so DRY_RUN=1 can be exercised on a dev box, where
# no module system exists. On O2 the loads run normally, and a genuinely missing R
# module still fails loudly at the Rscript preflight two lines down rather than
# being swallowed here.
if command -v module >/dev/null 2>&1; then
  module load gcc/14.2.0 2>/dev/null || module load gcc || true
  module load R/4.4.2   2>/dev/null || module load R   || true
fi

# Installed-package freshness. A forgotten `R CMD INSTALL` after a git pull means
# every task silently runs STALE estimator code and the whole sweep has to be
# thrown away; this is the cheapest possible place to catch it.
Rscript -e '
  pkgs <- c("doubletree", "optimaltrees")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(sprintf("package(s) NOT installed on this node: %s -- R CMD INSTALL them first (README_O2.md, prerequisites).",
                 paste(miss, collapse = ", ")))
  for (p in pkgs) cat(sprintf("preflight: %s %s OK\n", p, utils::packageVersion(p)))
' || preflight_fail "required package(s) missing (see message above)."

# The DGP gate run_pilot.R/run_sweep.R both run before any replication: DGP-A's
# bilinear remainder must be nonzero, or R1 is measuring nothing. Cheap (seconds)
# and it fails the whole submit rather than 48 arrays' worth of tasks.
export PKG_ROOT
Rscript -e '
  setwd(Sys.getenv("PKG_ROOT")); Sys.setenv(H2H_USE_INSTALLED = "1")
  suppressMessages(source(file.path("simulations","head_to_head_comparison","code","common.R")))
  rem <- verify_main_effects_remainder(dgp_spec_shared_interaction())
  cat(sprintf("preflight: DGP-A bilinear remainder = %.6g (must be nonzero)\n",
              rem$bilinear_remainder_P))
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
# print_sizing.R --format env emits: REGIME N_TASKS WALLTIME MEM_GB PARTITION
SIZING_FILE="${SCRATCH_DIR}/sizing.txt"
Rscript "${SLURM_DIR}/print_sizing.R" --format env \
  --regimes "$(IFS=,; echo "${REGIMES[*]}")" --pkg-root "${PKG_ROOT}" \
  > "${SIZING_FILE}" || preflight_fail "could not compute sizing."
echo
echo "Sizing (from slurm/units.R, target ${H2H_TARGET_SECS} s/task):"
cat "${SIZING_FILE}"
echo

# --- 4 submit one array per regime -------------------------------------------
declare -a ALL_JOB_IDS=()
declare -a SUBMIT_LINES=()

while read -r REGIME N_TASKS WALLTIME MEM_GB PARTITION; do
  [[ -z "${REGIME}" ]] && continue
  CAP=$(( CONCURRENCY_CAP < N_TASKS ? CONCURRENCY_CAP : N_TASKS ))
  REGIME_SCRATCH="${SCRATCH_DIR}/${REGIME}"
  mkdir -p "${REGIME_SCRATCH}"

  echo "--- ${REGIME}: ${N_TASKS} task(s)  --time ${WALLTIME} --mem ${MEM_GB}G -p ${PARTITION} (%${CAP}) ---"

  if [[ "${DRY_RUN}" == "1" ]]; then
    SUBMIT_LINES+=("${REGIME}: 1-${N_TASKS}%${CAP} ${WALLTIME} ${MEM_GB}G ${PARTITION}")
    continue
  fi

  # Shell-exported, then a plain --export=ALL. The combined --export=ALL,KEY=value
  # form is CANCELLED BY ROOT on O2 within seconds with NO output written at all --
  # it destroyed run 20260918-130255_88ba1bd. Isolated 2026-09-23; see
  # O2_SSH_GOTCHAS.md section 12. Do not collapse this back into the flag.
  export PKG_ROOT="${PKG_ROOT}"
  export SCRATCH_DIR="${REGIME_SCRATCH}"
  export REGIME="${REGIME}"
  export TARGET_SECS="${H2H_TARGET_SECS}"
  jobid=$(sbatch --parsable \
    --job-name="${STUDY_NAME}-${REGIME}" \
    --array=1-"${N_TASKS}"%"${CAP}" \
    --partition="${PARTITION}" \
    --time="${WALLTIME}" \
    --mem="${MEM_GB}G" \
    --output="${LOG_DIR}/${REGIME}_%A_%a.out" \
    --error="${LOG_DIR}/${REGIME}_%A_%a.err" \
    --export=ALL \
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
# One file that says exactly which run produced which numbers, and holds the
# combine command with the run's own paths already substituted -- so recovering a
# run three weeks later is copy-paste, not archaeology.
{
  echo "# Run Manifest -- ${STUDY_NAME}"
  echo
  echo "- run id: \`${RUN_ID}\`   git SHA: \`${GIT_SHA}\`"
  echo "- submitted: $(date '+%F %T %Z')"
  echo "- regimes: ${REGIMES[*]}"
  echo "- target seconds per task: ${H2H_TARGET_SECS}"
  echo "- scratch dir: \`${SCRATCH_DIR}\` (one subdir per regime)"
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
  echo "Rscript slurm/combine_results.R --scratch-dir ${SCRATCH_DIR} \\"
  echo "  --pkg-root ${PKG_ROOT} --run-id ${RUN_ID}"
  echo '```'
} > "${STUDY_DIR}/MANIFEST.md"

echo
echo "Submitted ${#ALL_JOB_IDS[@]} array job(s) across ${#REGIMES[@]} regime(s)."
echo "Manifest: ${STUDY_DIR}/MANIFEST.md"
echo "Monitor : RUN_ID=${RUN_ID} bash slurm/check_progress.sh"
