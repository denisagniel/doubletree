# =============================================================================
# run_one.R -- run a single work unit (one config x one replication)
# =============================================================================
# The heart of the study. run_one() is pure with respect to its seed: given the
# same unit row it always returns the same one-row data frame. It is called both
# locally (profile_timing.R) and on the cluster (run_replication.R), so it must
# not depend on any cluster-specific state.
#
# Contract:
#   input : unit_row -- one row of unit_table() (has config_id, rep_id, seed,
#                       and the grid columns)
#   output: a ONE-ROW data.frame (never a list) so results rbind cleanly.
#           MUST include an `estimate` column: slurm/profile_timing.R treats a
#           unit whose `estimate` is all-NA as a FAILED replication and refuses
#           to size the job array if every probe unit fails. Keep this column
#           (or generalize that check) if you change the output schema.
# =============================================================================

# Assumes generate_data(), estimate(), true_value() are already sourced
# (dgp.R, estimators.R) and grid.R for column names.

run_one <- function(unit_row) {
  # Deterministic, unit-specific seed -> reproducible & independent across units.
  set.seed(unit_row$seed)

  config <- unit_row  # unit_row carries all grid columns plus ids/seed

  data <- generate_data(config)
  est  <- estimate(data, config)
  truth <- true_value(config)

  # Derived performance quantities computed per-unit; aggregation happens later.
  # covered: goal-(i) coverage of the reported estimator's CI (single tree for
  # method="single_tree"). covered_crossfit: goal-(ii) coverage of the valid
  # cross-fit twin's CI (NA where the approach has no twin). Both recorded so the
  # analysis can compare single-tree vs twin coverage directly.
  z <- stats::qnorm(0.975)
  covered <- as.integer(truth >= est$ci_lower & truth <= est$ci_upper)

  cf_lower <- est$theta_crossfit - z * est$se_crossfit
  cf_upper <- est$theta_crossfit + z * est$se_crossfit
  covered_crossfit <- as.integer(truth >= cf_lower & truth <= cf_upper)

  # One-row data frame. Keep identity columns first for easy grouping.
  data.frame(
    unit      = unit_row$unit,
    config_id = unit_row$config_id,
    rep_id    = unit_row$rep_id,
    n         = config$n,
    dgp       = config$dgp,
    method    = config$method,
    estimate  = est$estimate,
    std_error = est$std_error,
    ci_lower  = est$ci_lower,
    ci_upper  = est$ci_upper,
    truth     = truth,
    error     = est$estimate - truth,
    covered   = covered,
    converged = as.integer(est$converged),
    # Goal-(ii) fidelity diagnostic + twin coverage (NA where not applicable).
    theta_crossfit        = est$theta_crossfit,
    se_crossfit           = est$se_crossfit,
    delta                 = est$delta,
    delta_over_se         = est$delta_over_se,
    covered_crossfit      = covered_crossfit,
    intersection_nonempty = as.integer(est$intersection_nonempty),
    # Rashomon tolerance multiplier selected by escalation per nuisance
    # (epsilon_n = c*log(n)/n; c=1 = theory value, larger = intersection had to be
    # widened). NA for non-Rashomon methods. Enables the inference-vs-c analysis.
    rashomon_c_e          = est$rashomon_c_e,
    rashomon_c_m0         = est$rashomon_c_m0,
    stringsAsFactors = FALSE
  )
}
