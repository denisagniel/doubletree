# ============================================================
# run_pilot.R
# Study: single_tree_corollaries  (doubletree)
# Spec:  quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md  (§6)
#
# Run:  Rscript simulations/single_tree_corollaries/code/run_pilot.R [reps] [n] [regimes]
#       (from the doubletree package root; defaults reps = 50, n = 1000, A,B,C,E)
#       e.g.  Rscript .../run_pilot.R 50 1000 A,C
#
# Spec §6's LOCAL PILOT, and nothing beyond it. The cluster stage (nsim = 1000,
# n in {1000, 4000}, 17 cells) is a SEPARATE decision and is deliberately not
# runnable from this file.
#
# Does four things and then stops:
#   1. Re-runs code/dgps.R's population gates, so no replication can run against
#      DGP constants that no longer reproduce spec §4's pre-computed numbers.
#   2. Runs `reps` replications of every requested regime at `n`.
#   3. Prints, per regime, spec §5's metrics with Monte Carlo SEs, plus Regime E's
#      partition-stratified variances.
#   4. Prints the gate-by-gate reading against spec §4's table.
#
# WHAT THE PILOT CAN AND CANNOT SETTLE (spec §6, stated before the run):
#   CAN  -- (i) the cor:single-tree identity to 1e-10; (ii) fitted-arm partition
#           recovery >= 95%; (iii) whether empirical variances track the analytic
#           V_tau/V ratios in DIRECTION and rough magnitude.
#   CANNOT -- resolve a coverage or variance-ratio gate to spec §4's precision. At
#           R = 50 a variance ratio carries ~20% relative MC SE (against the 6.3%
#           spec §4 pre-registers at R = 1000) and a coverage estimate carries
#           ~3-7 percentage points. Regime A's [0.88, 1.14] band is registered for
#           R = 1000; a pilot excursion outside it is NOT a violation.
# ============================================================

source(file.path("simulations", "single_tree_corollaries", "code", "common.R"))
source(file.path("simulations", "single_tree_corollaries", "code", "analyze.R"))

suppressPackageStartupMessages(library(knitr))

# ---- 0 arguments, with working defaults ---------------------------------

.args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(.args) >= 1L && !is.na(suppressWarnings(as.integer(.args[[1L]])))) {
  as.integer(.args[[1L]])
} else {
  PILOT_NSIM
}
N_PILOT <- if (length(.args) >= 2L && !is.na(suppressWarnings(as.integer(.args[[2L]])))) {
  as.integer(.args[[2L]])
} else {
  PILOT_N
}
PILOT_REGIMES <- if (length(.args) >= 3L && nzchar(.args[[3L]])) {
  trimws(strsplit(.args[[3L]], ",", fixed = TRUE)[[1L]])
} else {
  ST_REGIME_IDS
}
unknown <- setdiff(PILOT_REGIMES, ST_REGIME_IDS)
if (length(unknown)) {
  cli::cli_abort(c(
    "Unknown regime(s) {.val {unknown}}; expected {.val {ST_REGIME_IDS}}.",
    i = "Regime D has no DGP: spec §3 folds rem:single-tree-se into A-C as the naive-vs-corrected SE metric."
  ))
}

STAMP <- format(Sys.time(), "%Y%m%d-%H%M%S")

cli::cli_h1("Pilot: {REPS} reps, n = {N_PILOT}, regimes {paste(PILOT_REGIMES, collapse = ', ')}")
cli::cli_inform(c(
  "!" = "PILOT SCALE ONLY. Variance-ratio relative MC SE ~ {signif(100 * st_mcse_var_rel(REPS), 3)}% (spec §4 pre-registers 6.3% at nsim = 1000); coverage MC SE ~ {signif(st_mcse_prop(0.95, REPS), 3)}-{signif(st_mcse_prop(0.7, REPS), 3)}.",
  "i" = "Regime A's [0.88, 1.14] band is registered for nsim = 1000. A pilot excursion is not a violation."
))

# ---- 1 the population gate ----------------------------------------------

cli::cli_h1("1. Population gates (spec §4), recomputed from code/dgps.R")
gates_pop <- verify_population_gates()
print(knitr::kable(gates_pop, digits = 10, format = "simple", row.names = FALSE))
cli::cli_alert_success("All {nrow(gates_pop)} population gate(s) reproduce spec §4.")

pop_all <- do.call(rbind, lapply(ST_REGIME_IDS, function(r) {
  st_population_table(make_st_dgp(ST_REGIMES[[r]]$dgp))
}))
cli::cli_h2("Analytic population quantities, all regimes")
print(knitr::kable(
  pop_all[, c("dgp", "partition", "pi_pop", "theta0", "Tr", "C_untied", "C_tied",
              "V_num", "V_tau_num", "ratio_Vtau_V", "V_full", "V_tau_full",
              "naive_over_correct", "implied_coverage_naive")],
  digits = 6, format = "simple", row.names = FALSE))
utils::write.csv(pop_all,
                 file.path(DIR_TABLES, sprintf("population_%s.csv", STAMP)),
                 row.names = FALSE)

# ---- 2 run every regime -------------------------------------------------

cli::cli_h1("2. Running {length(PILOT_REGIMES)} regime(s)")
t_start <- Sys.time()
all_res <- lapply(PILOT_REGIMES, function(r) {
  st_run_cell(r, n = N_PILOT, reps = REPS, on_error = "record")
})
res <- do.call(rbind, all_res)
rownames(res) <- NULL
elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

saveRDS(res, file.path(DIR_RESULTS, sprintf("pilot_%s.rds", STAMP)))

n_failed <- length(unique(paste(res$regime, res$rep)[!is.na(res$error_message)]))
if (n_failed > 0L) {
  cli::cli_warn("{n_failed} replication(s) failed; see the per-cell messages above.")
} else {
  cli::cli_alert_success("All {length(PILOT_REGIMES)} regime(s) x {REPS} reps executed with zero failures.")
}

# ---- 3 metrics ----------------------------------------------------------

cli::cli_h1("3. Metrics (spec §5), with Monte Carlo SEs")
smry <- st_summarise(res)
utils::write.csv(smry, file.path(DIR_TABLES, sprintf("pilot_summary_%s.csv", STAMP)),
                 row.names = FALSE)

for (r in PILOT_REGIMES) {
  cli::cli_h2("Regime {r} -- {ST_REGIMES[[r]]$claims}")
  st_print_regime(smry, r)
}

if ("E" %in% PILOT_REGIMES) {
  cli::cli_h2("Regime E -- variance STRATIFIED by realised partition (spec §3, §5)")
  eE <- st_summarise_regime_E(res)
  if (is.null(eE)) {
    cli::cli_alert_warning("No Regime E fitted-arm rows.")
  } else {
    print(knitr::kable(eE, digits = 6, format = "simple", row.names = FALSE))
    utils::write.csv(eE, file.path(DIR_TABLES, sprintf("pilot_regimeE_%s.csv", STAMP)),
                     row.names = FALSE)
    cli::cli_inform(c(
      "i" = "The finding lem:single-tree-linear's caveat asks for is that the stratum-specific variances DIFFER. Read `emp_var_scaled` against each stratum's own `V_tau_full`, and note the analytic contrast between the two ties is only {signif(100 * abs(diff(range(pop_all$V_tau_full[pop_all$dgp == 'E']))) / min(pop_all$V_tau_full[pop_all$dgp == 'E']), 3)}% at these constants -- far below what R = {REPS} can resolve."
    ))
  }
}

# ---- 4 the gates --------------------------------------------------------

cli::cli_h1("4. Gate-by-gate reading against spec §4")
st_print_gates(smry)
gates_run <- do.call(rbind, lapply(
  intersect(c("oracle_tie", "fitted_tie"), unique(smry$arm)),
  function(a) st_gate_table(smry, arm = a)
))
utils::write.csv(gates_run, file.path(DIR_TABLES, sprintf("pilot_gates_%s.csv", STAMP)),
                 row.names = FALSE)

cli::cli_inform(c(
  "*" = "Pilot wall time: {signif(elapsed_min, 4)} min for {length(PILOT_REGIMES)} regime(s) x {REPS} reps at n = {N_PILOT}.",
  "*" = "Spec §6's cluster stage (nsim = {CLUSTER_NSIM}, n in {.val {CLUSTER_N}}, 17 cells) is a SEPARATE decision and is not run here.",
  "i" = "Sequential and single-threaded (worker_limit = {WORKER_LIMIT}). Parallelise over CELLS, never within a fit."
))

cli::cli_h1("Pilot complete. Read the MC SEs before quoting any gate above.")
