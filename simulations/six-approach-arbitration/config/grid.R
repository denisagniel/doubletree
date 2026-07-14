# =============================================================================
# grid.R -- single source of truth for the "six-approach-arbitration" simulation study
# =============================================================================
# Sourced by run_replication.R, profile_timing.R, combine.R, and (via Rscript)
# submit.sh. Defines the parameter grid, the number of replications, and study
# identity. Nothing here should have side effects beyond assigning objects.
#
# Work model:
#   - GRID has one row per configuration (a "cell" of the design).
#   - Each configuration is run TOTAL_REPS times.
#   - A "unit" is one (config, rep) pair. Total units = nrow(GRID) * TOTAL_REPS.
#   - Units are enumerated deterministically (see unit_table()) so that the
#     same unit index always maps to the same (config, rep) and the same seed.
# =============================================================================

STUDY_NAME   <- "six-approach-arbitration"
PROJECT_NAME <- "global-scholars"

# Base seed for the whole study. Per-(config, rep) seeds are derived from this
# via .unit_seed() so every replication is independent, reproducible, parallel-safe
# -- and INDEPENDENT of unit ordering (see .unit_seed / S6 note below).
BASE_SEED <- 20240101L

# Replications per configuration. Total Monte Carlo work = nrow(GRID) * TOTAL_REPS.
TOTAL_REPS <- 1000L

# -----------------------------------------------------------------------------
# Parameter grid. EDIT THIS: one row per configuration you want to run.
# Use named columns with documented meanings -- no magic numbers.
# Include at least one STRESS regime (Constitution Section 9): a setting where
# the method under study is expected to struggle.
# -----------------------------------------------------------------------------
# Main grid: seven estimators at the FIXED theory Rashomon tolerance (escalate = FALSE).
#   full            : full-sample single tree, in-sample (biased baseline)
#   crossfit        : K separate trees, out-of-sample (valid; no single tree)
#   doubletree      : Rashomon-intersection struct, cross-fit leaves (valid twin)
#   dt_averaged     : intersection struct, averaged leaves, single tree + honest CI
#   msplit          : modal struct, cross-fit predictions (valid twin; SE caveat)
#   msplit_averaged : modal struct, averaged leaves, single tree + honest CI
#   single_tree     : Alt A -- intersection struct, all-n leaves, single tree
GRID_MAIN <- expand.grid(
  n      = c(500L, 1000L, 2000L),        # sample size
  # DGP regime (see dgp.R). "continuous" is the STRESS regime (Constitution S9):
  # x4^2 nonlinearity requires many thresholds; tree discretization struggles.
  dgp    = c("simple", "moderate", "complex", "continuous"),
  method = c("full", "crossfit", "doubletree", "dt_averaged",
             "msplit", "msplit_averaged", "single_tree"),
  escalate = FALSE,                      # fixed theory epsilon_n = log(n)/n
  stringsAsFactors = FALSE
)

# Escalation coverage sweep (escalate = TRUE): only the Rashomon-intersection methods
# (whose intersection can be empty at the theory tolerance). Widening epsilon_n = c*
# log(n)/n TRADES the fixed-tolerance validity guarantee for a non-empty intersection;
# run_one records the realized rashomon_c_e/m0 so data-analysis can plot coverage vs c
# and answer empirically whether large c degrades the CLT (o(n^{-1/2}) question). The
# "continuous" DGP is omitted here: at the fixed tolerance it already exceeds 64 GB
# (INFEASIBLE_CELLS below), and escalation only enlarges the Rashomon set further.
GRID_ESCALATE <- expand.grid(
  n      = c(500L, 1000L, 2000L),
  dgp    = c("simple", "moderate", "complex"),
  method = c("doubletree", "dt_averaged", "single_tree"),
  escalate = TRUE,
  stringsAsFactors = FALSE
)

GRID <- rbind(GRID_MAIN, GRID_ESCALATE)

# METHOD-MAJOR ORDER (required by the per-method submission path). profile_per_method.R
# and submit_per_method.sh assume each method occupies a CONTIGUOUS block of global
# units (UNIT_OFFSET = min(unit[method]); N_UNITS_METHOD = sum(method)). The escalation
# sweep adds a second block of rows for doubletree/dt_averaged/single_tree, so without
# this sort each of those methods would span two disjoint unit ranges and the per-method
# array would run the wrong units. Order by method (grid appearance order), then keep
# escalate=FALSE before TRUE and the expand.grid n/dgp order within each method so
# offsets stay deterministic. Stable sort preserves within-group order.
.method_levels <- c("full", "crossfit", "doubletree", "dt_averaged",
                    "msplit", "msplit_averaged", "single_tree")
GRID <- GRID[order(factor(GRID$method, levels = .method_levels), GRID$escalate), ,
             drop = FALSE]
rownames(GRID) <- NULL

# Stable configuration id (1..nrow(GRID)); do not reorder GRID after a run
# without cleaning stale scratch, or unit->config mapping will change.
GRID$config_id <- seq_len(nrow(GRID))

# -----------------------------------------------------------------------------
# .unit_seed() -- deterministic seed for a (config, rep) pair (S6).
# Depends ONLY on (config_id, rep_id), never on the `unit` index or unit ordering,
# so re-ordering / re-stratifying the grid never changes a replication's seed. The
# old `base_seed + unit` scheme was ordering-DEPENDENT: adding/reordering configs
# shifted every downstream unit's seed, silently re-randomizing untouched cells.
# Arithmetic is done in DOUBLES then reduced mod 2^31-1 into the valid 32-bit range
# set.seed() requires -- base-R integer arithmetic overflows to NA past 2^31 (see
# MEMORY base-r-hashing-gotcha), so we never form an integer beyond that. Distinct
# (config, rep) pairs get distinct seeds as long as rep_id <= total_reps.
# -----------------------------------------------------------------------------
.unit_seed <- function(config_id, rep_id, total_reps = TOTAL_REPS, base_seed = BASE_SEED) {
  MOD    <- 2147483647                                   # 2^31 - 1 (Mersenne prime)
  offset <- (as.double(config_id) - 1) * total_reps + rep_id
  as.integer((base_seed + offset) %% MOD)
}

# -----------------------------------------------------------------------------
# unit_table() -- deterministic enumeration of all (config, rep) work units.
# Returns a data frame with columns: unit, config_id, rep_id, seed, plus the
# grid columns. unit runs 1..(nrow(GRID) * TOTAL_REPS).
# -----------------------------------------------------------------------------
unit_table <- function(grid = GRID, total_reps = TOTAL_REPS, base_seed = BASE_SEED) {
  reps <- seq_len(total_reps)
  # rep varies fastest within a config, so a config's reps are contiguous.
  ut <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    data.frame(
      config_id = grid$config_id[i],
      rep_id    = reps,
      grid[i, setdiff(names(grid), "config_id"), drop = FALSE],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
  ut$unit <- seq_len(nrow(ut))
  # Ordering-invariant per-(config, rep) seed (S6). NOTE: this changes RNG streams
  # relative to the old base_seed+unit scheme -- intended and accepted.
  ut$seed <- .unit_seed(ut$config_id, ut$rep_id, total_reps, base_seed)
  ut[, c("unit", "config_id", "rep_id", "seed",
         setdiff(names(ut), c("unit", "config_id", "rep_id", "seed")))]
}

# Convenience: total number of work units.
n_units <- function() nrow(GRID) * TOTAL_REPS

# -----------------------------------------------------------------------------
# INFEASIBLE_CELLS -- (method, n, dgp) cells excluded because the Rashomon set
# explodes beyond feasible memory (>64 GB) at the stress DGP. Determined by
# slurm/probe_tier2_boundary.R (2026-07-11): at dgp="continuous" the many-threshold
# discretization makes the per-fold Rashomon set past 64 GB even at the FIXED theory
# tolerance epsilon_n = log(n)/n. (The probe predates the escalation fix and ran with
# escalation inert, so this boundary is the fixed-tolerance cost, not escalation.
# Escalation cells, when added, are expected to be at least this memory-intensive.)
# This is a documented STRESS-REGIME LIMITATION of the Rashomon-intersection methods,
# not a bug (Constitution S9: report the boundary).
# NOTE: this does NOT change unit numbering -- unit_table() is unchanged, so unit
# ids/seeds/offsets are preserved. is_feasible() only tells the runner which units
# to SKIP; combine.R records them as not-run (distinct from converged=FALSE).
INFEASIBLE_CELLS <- rbind(
  data.frame(method = "doubletree",  n = 2000L, dgp = "continuous"),
  data.frame(method = "dt_averaged", n = 2000L, dgp = "continuous"),
  data.frame(method = "single_tree", n = 1000L, dgp = "continuous"),
  data.frame(method = "single_tree", n = 2000L, dgp = "continuous"),
  stringsAsFactors = FALSE
)

# is_feasible(config_row) -> TRUE unless (method,n,dgp) is in INFEASIBLE_CELLS.
# Vectorized over a data frame of config rows (the run block).
is_feasible <- function(df) {
  key  <- paste(df$method, df$n, df$dgp, sep = "\r")
  bad  <- paste(INFEASIBLE_CELLS$method, INFEASIBLE_CELLS$n,
                INFEASIBLE_CELLS$dgp, sep = "\r")
  !(key %in% bad)
}
