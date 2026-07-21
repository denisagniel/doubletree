#!/bin/bash
# =============================================================================
# smoke_autosubmit.sh -- build an ISOLATED tiny "smoke" study to test the
# profile -> gate -> AUTO_SUBMIT -> array-dispatch -> combine chain ON O2, in
# minutes, with dummy code -- BEFORE committing the real multi-day run.
# =============================================================================
# WHY: the AUTO_SUBMIT path (profile.slurm chain-launching submit_per_method.sh
# from inside a running job) depends on O2 permitting sbatch-from-within-a-job --
# a cluster policy we cannot verify by reading code. This harness exercises the
# REAL slurm scripts (symlinked, not copied) against a 7-method, ~42-unit dummy
# study whose run_one() just sleeps ~0.5s and returns a valid row, so the whole
# chain runs fast and cannot hit the Rashomon blow-up.
#
# ISOLATION: the smoke study lives in its own dir with its own config/R, and the
# submit uses a SMOKE-specific SCRATCH_ROOT, so it never touches the real study's
# config/, results/, or /n/scratch/.../six-approach-arbitration/ run-ids.
#
# USAGE (ON O2, from the real study dir):
#   bash dev-scripts/smoke_autosubmit.sh            # generate the smoke study
#   # then follow the printed commands (cd into the smoke dir and sbatch).
# Override the smoke dir:  SMOKE_DIR=/path bash dev-scripts/smoke_autosubmit.sh
# =============================================================================

set -euo pipefail

SLURM_DEV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DEV_DIR}")"
REAL_SLURM="${STUDY_DIR}/slurm"

# Smoke study dir (isolated). Default: sibling of the study dir so it is obvious
# and easy to delete; overridable.
SMOKE_DIR="${SMOKE_DIR:-${STUDY_DIR}/../six-approach-arbitration-SMOKE}"
SMOKE_DIR="$(mkdir -p "${SMOKE_DIR}" && cd "${SMOKE_DIR}" && pwd)"

# Sanity: the real scripts we symlink must exist.
for f in profile.slurm profile_per_method.R submit_per_method.sh array.slurm \
         run_replication.R combine.R monitor.sh; do
  [[ -f "${REAL_SLURM}/${f}" ]] || { echo "ERROR: ${REAL_SLURM}/${f} not found." >&2; exit 1; }
done

echo "Building smoke study in: ${SMOKE_DIR}"
mkdir -p "${SMOKE_DIR}/config" "${SMOKE_DIR}/R" "${SMOKE_DIR}/slurm" "${SMOKE_DIR}/logs" "${SMOKE_DIR}/results"

# --- Symlink the REAL slurm infra (test the actual code, not copies) ---------
for f in profile.slurm profile_per_method.R submit_per_method.sh array.slurm \
         run_replication.R combine.R monitor.sh; do
  ln -sfn "${REAL_SLURM}/${f}" "${SMOKE_DIR}/slurm/${f}"
done

# --- Tiny grid.R: same 7 method names, 1 n x 1 dgp, few reps => ~42 units -----
# Keeps the METHODS contract submit_per_method.sh hardcodes; empties INFEASIBLE.
cat > "${SMOKE_DIR}/config/grid.R" <<'RGRID'
# SMOKE grid -- tiny, for infra testing only. Mirrors the real grid.R API.
STUDY_NAME   <- "six-approach-arbitration-SMOKE"
PROJECT_NAME <- "global-scholars"
BASE_SEED    <- 20240101L
TOTAL_REPS   <- 6L                       # 6 reps x 7 methods = 42 units

GRID <- expand.grid(
  n      = 200L,
  dgp    = "simple",
  method = c("full", "crossfit", "doubletree", "dt_averaged",
             "msplit", "msplit_averaged", "single_tree"),
  escalate = FALSE,
  stringsAsFactors = FALSE
)
.method_levels <- c("full", "crossfit", "doubletree", "dt_averaged",
                    "msplit", "msplit_averaged", "single_tree")
GRID <- GRID[order(factor(GRID$method, levels = .method_levels), GRID$escalate), , drop = FALSE]
rownames(GRID) <- NULL
GRID$config_id <- seq_len(nrow(GRID))

.unit_seed <- function(config_id, rep_id, total_reps = TOTAL_REPS, base_seed = BASE_SEED) {
  MOD <- 2147483647
  offset <- (as.double(config_id) - 1) * total_reps + rep_id
  as.integer((base_seed + offset) %% MOD)
}
unit_table <- function(grid = GRID, total_reps = TOTAL_REPS, base_seed = BASE_SEED) {
  reps <- seq_len(total_reps)
  ut <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    data.frame(config_id = grid$config_id[i], rep_id = reps,
               grid[i, setdiff(names(grid), "config_id"), drop = FALSE],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  ut$unit <- seq_len(nrow(ut))
  ut$seed <- .unit_seed(ut$config_id, ut$rep_id, total_reps, base_seed)
  ut[, c("unit", "config_id", "rep_id", "seed",
         setdiff(names(ut), c("unit", "config_id", "rep_id", "seed")))]
}
n_units <- function() nrow(GRID) * TOTAL_REPS
INFEASIBLE_CELLS <- data.frame(method = character(), n = integer(),
                               dgp = character(), stringsAsFactors = FALSE)
is_feasible <- function(df) rep(TRUE, nrow(df))
RGRID

# --- Dummy dgp/estimators/run_one: fast, valid schema, no heavy deps ----------
cat > "${SMOKE_DIR}/R/dgp.R" <<'RDGP'
# SMOKE dgp -- trivial.
generate_data <- function(config) list(n = config$n)
true_value    <- function(config) 0.0
RDGP

cat > "${SMOKE_DIR}/R/estimators.R" <<'REST'
# SMOKE estimator -- returns a deterministic-ish value + full est schema.
estimate <- function(data, config) {
  Sys.sleep(0.5)                          # simulate a little work (measurable med)
  est <- stats::rnorm(1, 0, 0.1)          # seeded upstream in run_one()
  se  <- 0.1
  z   <- stats::qnorm(0.975)
  list(estimate = est, std_error = se, ci_lower = est - z*se, ci_upper = est + z*se,
       converged = TRUE, theta_crossfit = est, se_crossfit = se,
       delta = 0, delta_over_se = 0, intersection_nonempty = TRUE,
       rashomon_c_e = NA_real_, rashomon_c_m0 = NA_real_)
}
REST

cat > "${SMOKE_DIR}/R/run_one.R" <<'RONE'
# SMOKE run_one -- mirrors the real one-row schema (see R/run_one.R contract).
run_one <- function(unit_row) {
  set.seed(unit_row$seed)
  config <- unit_row
  data <- generate_data(config); est <- estimate(data, config); truth <- true_value(config)
  z <- stats::qnorm(0.975)
  data.frame(
    unit = unit_row$unit, config_id = unit_row$config_id, rep_id = unit_row$rep_id,
    n = config$n, dgp = config$dgp, method = config$method,
    escalate = isTRUE(config$escalate),
    estimate = est$estimate, std_error = est$std_error,
    ci_lower = est$ci_lower, ci_upper = est$ci_upper, truth = truth,
    error = est$estimate - truth,
    covered = as.integer(truth >= est$ci_lower & truth <= est$ci_upper),
    converged = as.integer(est$converged),
    theta_crossfit = est$theta_crossfit, se_crossfit = est$se_crossfit,
    delta = est$delta, delta_over_se = est$delta_over_se,
    covered_crossfit = 1L,
    intersection_nonempty = as.integer(est$intersection_nonempty),
    rashomon_c_e = est$rashomon_c_e, rashomon_c_m0 = est$rashomon_c_m0,
    stringsAsFactors = FALSE)
}
RONE

# --- Smoke-specific scratch root (isolate from the real study's scratch) ------
HMS_ID="dma12"
SMOKE_SCRATCH="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME:-global-scholars}/SMOKE-six-approach-arbitration"

cat <<EOF

=============================================================================
 Smoke study ready:  ${SMOKE_DIR}
 7 methods x 6 reps = 42 dummy units (~0.5s each). Symlinks the REAL slurm code.
=============================================================================

RUN THE FULL AUTO_SUBMIT CHAIN (from the SMOKE dir, so SLURM_SUBMIT_DIR is set):

  cd "${SMOKE_DIR}"
  sbatch --partition=short --mem=4G --time=0-01:00 \\
    --export=ALL,SCRATCH_ROOT=${SMOKE_SCRATCH},TARGET_TASKS=14,AUTO_SUBMIT=1,N_UNITS_PROBE=2 \\
    slurm/profile.slurm

  # The --partition/--mem/--time flags shrink the PROFILING job (override its
  # #SBATCH medium/48G/24h). TARGET_TASKS=14 -> ~2 tasks/method. AUTO_SUBMIT=1
  # chains submit_per_method.sh IF all 7 sizing envs come out valid.

WATCH IT:
  squeue -u ${HMS_ID}                       # profile job, then 7 method arrays
  tail -f logs/profile_*.out                # gate decision + "Submitted N array job(s)"
  bash slurm/monitor.sh                      # completed task files vs expected

VERIFY (success = all of):
  1. logs/profile_*.out shows "AUTO_SUBMIT: all 7 sizing envs valid; launching..."
     and submit_per_method.sh's per-method "array <jobid>" lines.
     -> confirms sbatch-from-within-a-job is PERMITTED on O2 (the key unknown).
  2. config/sizing_*.env exist for all 7 methods, each with PARTITION=.
  3. find ${SMOKE_SCRATCH} -name 'task_*.rds' | wc -l   == 14 (all tasks wrote).
  4. Combine:
       Rscript slurm/combine.R --run-id <run-id> \\
         --scratch-dir ${SMOKE_SCRATCH}/<run-id> --study-dir "${SMOKE_DIR}"
     -> writes results/<run-id>.rds with 42 unique unit rows, no errors.

TEST THE GATE'S REFUSAL (optional, proves it won't launch a partial study):
  rm config/sizing_single_tree.env
  sbatch --partition=short --mem=4G --time=0-01:00 \\
    --export=ALL,SCRATCH_ROOT=${SMOKE_SCRATCH},AUTO_SUBMIT=1 slurm/profile.slurm
  # (re-profiles + resubmits; but to test refusal in isolation, inspect the log:
  #  it should say "NOT submitting" only if a method's env is missing/invalid.)

CLEANUP when done:
  rm -rf "${SMOKE_DIR}" "${SMOKE_SCRATCH}"
=============================================================================
EOF
