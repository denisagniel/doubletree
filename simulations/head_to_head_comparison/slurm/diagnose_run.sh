#!/bin/bash
# =============================================================================
# slurm/diagnose_run.sh -- one-shot post-mortem for the 2026-09-18 S5 submission
# =============================================================================
#   ON O2 (VPN OFF), from the study dir:
#       cd ~/doubletree/simulations/head_to_head_comparison
#       bash slurm/diagnose_run.sh
#
# Answers, in ONE session (= one Duo approval, per O2_SSH_GOTCHAS.md §1):
#   (a) why run 20260918-130255_88ba1bd died within ~1s of starting
#   (b) whether compute nodes have outbound egress, or whether SLURM mail is the
#       only viable completion channel
#   (c) that run's final state and whether ANY results survived on scratch
#
# Leading hypothesis under test (O2_SSH_GOTCHAS.md §4, §11): `module` is a bash
# FUNCTION that Lmod `export -f`s into the environment, so `sbatch --export=ALL`
# carries it only from an INTERACTIVE submitting shell. This run was submitted from
# `ssh host 'bash -s'`, where it is absent -- so `module load R/4.4.2` failed and
# `set -euo pipefail` (run_simulations.slurm:29) aborted every task. The probeALL
# vs probeNONE pair below is the discriminating experiment.
#
# Everything runs with `set +e`: this script diagnoses failures, so it must not
# abort on one. It deliberately does NOT source any _env.sh -- whether that
# bootstrap works is part of what is being tested.
# =============================================================================
set +e

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="head_to_head_comparison"
RUN_ID="20260918-130255_88ba1bd"
ARRAYS="53840731,53840732,53840733,53840734,53840735,53840736"
MAIL_USER="dagniel@rand.org"
REGIMES="R1 R2 R3 R4 R5 R6"

SCRATCH_ROOT="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}"
SCRATCH_DIR="${SCRATCH_ROOT}/${RUN_ID}"
OUT_DIR="${HOME}/o2diag_${RUN_ID}"
mkdir -p "${OUT_DIR}"

{
echo "=== 0 IDENTITY (this login shell) ==="
date; hostname; echo "user=${USER} shell=${SHELL}"
# `type -t` distinguishes an inherited shell FUNCTION from a binary or nothing;
# `command -v` would report a function and a binary identically.
echo "module_is=$(type -t module 2>/dev/null || echo UNDEFINED)"
echo "bash_func_count=$(env | grep -c '^BASH_FUNC_')"
echo "LOADEDMODULES=${LOADEDMODULES:-<unset>}"
echo "R_LIBS_USER=${R_LIBS_USER:-<unset>}"

echo
echo "=== 1 SACCT: what the six arrays actually did ==="
# ExitCode is the highest-information field and needs no log reading:
#   FAILED 1:0  = script exited non-zero (set -e on a failed `module load`)
#   FAILED 127  = command not found
#   CANCELLED   = killed externally  <- the "O2 noticed the SSHs" theory
#   OUT_OF_MEMORY / NODE_FAIL / TIMEOUT = none of the above
sacct -j "${ARRAYS}" \
  -o JobID%20,JobName%10,State%18,ExitCode,DerivedExitCode,Reason%26,Elapsed,Start \
  2>&1 | head -50
echo "-- same window, in case the job ids aged out of the default range --"
sacct -S 2026-09-18 -E 2026-09-19 -u "${USER}" --name=h2h \
  -o JobID%20,State%18,ExitCode,Reason%26,Elapsed 2>&1 | head -20

echo
echo "=== 2 SCRATCH: did anything survive? ==="
echo "scratch=${SCRATCH_DIR}"
ls -ld "${SCRATCH_DIR}" 2>&1
# This study writes batch_*.rds shards and part_*.rds partials into a PER-REGIME
# subdir (launch_subset.sh:134) -- NOT the scaffold's task_*.rds / partials/unit_*.rds.
# find, never a glob: a full sweep can hold ~100k files.
for R in ${REGIMES}; do
  d="${SCRATCH_DIR}/${R}"
  printf '  %-3s shards=%-6s partials=%-6s dir=%s\n' "${R}" \
    "$(find "${d}" -maxdepth 1 -name 'batch_*.rds' 2>/dev/null | wc -l | tr -d ' ')" \
    "$(find "${d}" -maxdepth 1 -name 'part_*.rds'  2>/dev/null | wc -l | tr -d ' ')" \
    "$([ -d "${d}" ] && echo present || echo MISSING)"
done
echo "logs dir: ${SCRATCH_DIR}/logs"
echo "  .out files=$(find "${SCRATCH_DIR}/logs" -name '*.out' 2>/dev/null | wc -l | tr -d ' ')"
echo "  .err files=$(find "${SCRATCH_DIR}/logs" -name '*.err' 2>/dev/null | wc -l | tr -d ' ')"
echo "  non-empty .err=$(find "${SCRATCH_DIR}/logs" -name '*.err' -size +0 2>/dev/null | wc -l | tr -d ' ')"
echo "-- first NON-EMPTY .err (THE payload: an empty one means death before any output) --"
ERRF=$(find "${SCRATCH_DIR}/logs" -name '*.err' -size +0 2>/dev/null | head -1)
echo "err_file=${ERRF:-<none non-empty>}"
[ -n "${ERRF}" ] && sed -n '1,25p' "${ERRF}"
echo "-- first .out --"
OUTF=$(find "${SCRATCH_DIR}/logs" -name '*.out' 2>/dev/null | head -1)
echo "out_file=${OUTF:-<none>}"
[ -n "${OUTF}" ] && sed -n '1,15p' "${OUTF}"

echo
echo "=== 3 PROBE SCRIPT ==="
cat > "${OUT_DIR}/probe.sh" <<'PROBE'
#!/bin/bash
set +eu
echo "PROBE host=$(hostname) xmode=${XMODE:-unset}"
echo "P_module_before=$(type -t module 2>/dev/null || echo UNDEFINED)"
echo "P_bash_func_count=$(env | grep -c '^BASH_FUNC_')"
echo "P_LOADEDMODULES_in=${LOADEDMODULES:-<unset>}"
echo "P_R_LIBS_USER_in=${R_LIBS_USER:-<unset>}"
if ! command -v module >/dev/null 2>&1; then
  for f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /etc/profile; do
    [ -r "$f" ] && . "$f" && echo "P_sourced=$f" && break
  done
fi
echo "P_module_after=$(type -t module 2>/dev/null || echo UNDEFINED)"
module purge >/dev/null 2>&1
module load gcc/14.2.0 >/dev/null 2>&1; echo "P_gcc_rc=$?"
module load R/4.4.2    >/dev/null 2>&1; echo "P_R_rc=$?"
echo "P_Rscript=$(command -v Rscript || echo NONE)"
# --no-init-file matches how this study's real tasks run: the doubletree repo ships
# a TRACKED .Rprofile that devtools::load_all()s optimaltrees, which would otherwise
# shadow the installed packages under test.
Rscript --no-init-file -e '
  cat("P_Rver=", R.version.string, "\n", sep="")
  cat("P_libpaths=", paste(.libPaths(), collapse=":"), "\n", sep="")
  cat("P_optparse=", requireNamespace("optparse", quietly=TRUE), "\n", sep="")
  cat("P_doubletree=", requireNamespace("doubletree", quietly=TRUE), "\n", sep="")
  cat("P_optimaltrees=", requireNamespace("optimaltrees", quietly=TRUE), "\n", sep="")
  cat("P_blas=", tryCatch(La_library(), error=function(e) "unknown"), "\n", sep="")' 2>&1 | head -8
echo "P_scratch_write=$( { [ -n "${SCRDIR}" ] && touch "${SCRDIR}/.wtest" && rm -f "${SCRDIR}/.wtest" && echo OK; } 2>&1 || echo FAIL )"
# Egress: decides whether a compute node could ever push results outward, or whether
# SLURM mail is the only completion channel available to an unattended pipeline.
echo "P_curl=$(curl -sS -m 10 -o /dev/null -w '%{http_code}' https://example.com 2>&1 || echo CURLFAIL)"
echo "P_dns_transfer=$(getent hosts transfer.rc.hms.harvard.edu >/dev/null 2>&1 && echo OK || echo FAIL)"
echo "P_rsync=$(rsync -e 'ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=8' \
  -q /etc/hostname transfer.rc.hms.harvard.edu:/tmp/o2egress.$$ 2>&1 | head -1 || true)"
echo "PROBE_DONE"
PROBE
chmod +x "${OUT_DIR}/probe.sh"
echo "wrote ${OUT_DIR}/probe.sh"

echo
echo "=== 4 SUBMIT THE DISCRIMINATING PAIR (differ ONLY in --export) ==="
# Shell-exported, never via --export=ALL,K=V -- that form is CANCELLED BY ROOT on O2
# and killed these very probes on the first run (O2_SSH_GOTCHAS.md §12).
export XMODE=ALL
export SCRDIR="${OUT_DIR}"

J_ALL=$(sbatch --parsable -p short -t 0:05:00 --mem=2G -J probeALL \
  --export=ALL \
  --mail-type=END --mail-user="${MAIL_USER}" \
  -o "${OUT_DIR}/probeALL.out" -e "${OUT_DIR}/probeALL.err" \
  "${OUT_DIR}/probe.sh" 2>&1)
echo "J_ALL=${J_ALL}   (mimics a human-submitted job: module inherited)"
# --export=NONE may be rejected by some SLURM builds; nothing below depends on it.
J_NONE=$(sbatch --parsable -p short -t 0:05:00 --mem=2G -J probeNONE \
  --export=NONE \
  -o "${OUT_DIR}/probeNONE.out" -e "${OUT_DIR}/probeNONE.err" \
  "${OUT_DIR}/probe.sh" 2>&1)
echo "J_NONE=${J_NONE}  (mimics the agent-submitted job: module NOT inherited)"

echo
echo "=== 5 POLL (up to 150s) ==="
for _ in $(seq 1 30); do
  grep -q PROBE_DONE "${OUT_DIR}/probeALL.out" 2>/dev/null && break
  sleep 5
done
squeue -u "${USER}" -o '%.14i %.10P %.12T %.24R' 2>&1 | head
for f in probeALL.out probeALL.err probeNONE.out probeNONE.err; do
  echo "--- ${f} ---"; cat "${OUT_DIR}/${f}" 2>&1 | head -26
done

echo
echo "=== 6 SACCT FOR THE PROBES ==="
sacct -j "${J_ALL},${J_NONE}" -o JobID%16,JobName%12,State%18,ExitCode,Reason%24 2>&1 | head
echo
echo "=== END DIAG ==="
} 2>&1 | tee "${OUT_DIR}/DIAG.txt"

cat <<EOF

================== PASTE FROM HERE ==================
$(sed -n '/=== 0 IDENTITY/,/=== END DIAG/p' "${OUT_DIR}/DIAG.txt")
=================== TO HERE =========================

Saved to ${OUT_DIR}/DIAG.txt -- re-readable without spending another Duo approval.

HOW TO READ IT
  P_module_before: 'function' in probeALL but 'UNDEFINED' in probeNONE
      -> hypothesis CONFIRMED; the step-1 bootstrap fix is the right fix.
  Both 'UNDEFINED'
      -> module is never inherited here; the fix still works, the mechanism differs.
  sacct Reason on the six arrays is anything other than 'None'
      -> read that instead; the mechanism above is not the cause.
  State=CANCELLED on the arrays
      -> something DID kill them externally; the SSH theory is back in play.
  P_optparse / P_doubletree / P_optimaltrees = FALSE
      -> a second, independent instant-death path (library resolution).
  P_curl=200  -> compute nodes have egress.  CURLFAIL -> mail is the only channel.
  A "SLURM Job_id=... Ended" email to ${MAIL_USER} settles the mail question
  independently of curl.
EOF
