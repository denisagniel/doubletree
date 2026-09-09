# ============================================================
# run_pilot.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md
#
# Run:  Rscript simulations/head_to_head_comparison/code/run_pilot.R [reps] [regimes]
#       (from the doubletree package root; default reps = 25, all regimes)
#       e.g.  Rscript .../run_pilot.R 25 R1,R3
#
# Prototype BEFORE committing to the full grid, per this project's r-code
# conventions ("Prototype on a small subset or fewer iterations first; estimate
# run time") and spec §6 (nsim ~1000/2000 per cell is a cluster/batch job, not
# an interactive one).
#
# Does FOUR things and then stops. It deliberately does not roll into the full
# sweep: the sizing decision, and the directional sanity check, are the output.
#
#   1. Re-runs the DGP verification gate (code/verify_dgps.R's section 1) so the
#      pilot cannot run on a DGP whose bilinear remainder has silently cancelled.
#   2. Runs `reps` replications in every cell of every requested regime.
#   3. Prints, per regime, the metrics of spec §5 WITH Monte Carlo SEs, plus the
#      paired RMSE-ratio table where the regime has a reference arm.
#   4. Extrapolates to spec §6's nsim targets over the FULL n grid and prints a
#      projected wall time per regime, so the real sweep can be sized.
#
# READ SPEC §0 BEFORE QUOTING ANY NUMBER BELOW. This is not a leaderboard: under
# standard regularity a correctly-specified GLM, a correctly-specified black-box
# AIPW/DML estimator and doubletree under exact sparsity all attain the SAME
# efficiency bound, and no arm can win asymptotically. See common.R's header for
# the three claims (C1/C2/C3) this study is actually scoped to check, and for the
# regimes where a competitor is EXPECTED to do as well or better.
# ============================================================

source(file.path("simulations", "head_to_head_comparison", "code", "common.R"))
source(file.path("simulations", "head_to_head_comparison", "code", "analyze.R"))

suppressPackageStartupMessages(library(knitr))

# ---- 0 arguments, with working defaults ---------------------------------

.args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(.args) >= 1L && !is.na(suppressWarnings(as.integer(.args[[1L]])))) {
  as.integer(.args[[1L]])
} else {
  25L
}
PILOT_REGIMES <- if (length(.args) >= 2L && nzchar(.args[[2L]])) {
  trimws(strsplit(.args[[2L]], ",", fixed = TRUE)[[1L]])
} else {
  REGIME_IDS
}
unknown <- setdiff(PILOT_REGIMES, REGIME_IDS)
if (length(unknown)) {
  cli::cli_abort("Unknown regime(s) {.val {unknown}}; expected {.val {REGIME_IDS}}.")
}

STAMP <- format(Sys.time(), "%Y%m%d-%H%M%S")

cli::cli_h1("Pilot: {REPS} reps/cell, regimes {paste(PILOT_REGIMES, collapse = ', ')}")
cli::cli_inform(c(
  "!" = "PILOT SCALE ONLY. Spec §6's targets are nsim = {NSIM_FULL[['R1']]} (coverage regimes) and {NSIM_FULL[['R2']]} (RMSE-ratio regimes); at {REPS} reps a coverage MC SE is ~{signif(sqrt(0.95 * 0.05 / REPS), 3)}, so a 0.95-vs-0.90 distinction is NOT resolvable here. Directional sanity only.",
  "i" = "Every reported number carries its own MC SE (spec §6). Read them."
))

# ---- 1 the DGP gate ------------------------------------------------------

cli::cli_h1("1. DGP gate: DGP-A's bilinear remainder must be nonzero")
rem <- verify_main_effects_remainder(dgp_spec_shared_interaction())
cli::cli_alert_success(
  "remainder = {signif(rem$bilinear_remainder_P, 6)}, implied bias {signif(rem$implied_bias_P, 5)} ({signif(100 * rem$relative_bias_P, 3)}% of theta_0)."
)

# ---- 2 run every cell ----------------------------------------------------

grid <- design_grid(PILOT_REGIMES)
cli::cli_h1("2. Design grid: {nrow(grid)} cell(s)")
print(knitr::kable(grid, format = "simple", row.names = FALSE))

t_start <- Sys.time()
all_res <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  all_res[[i]] <- run_cell(g$regime, g$dgp, g$n, reps = REPS, on_error = "record")
}
res <- do.call(rbind, all_res)
rownames(res) <- NULL
elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

saveRDS(res, file.path(DIR_RESULTS, sprintf("pilot_%s.rds", STAMP)))

n_failed <- length(unique(paste(res$regime, res$dgp, res$n, res$rep)[
  !is.na(res$error_message)]))
if (n_failed > 0L) {
  cli::cli_warn("{n_failed} replication(s) failed across the pilot; see the per-cell messages above.")
} else {
  cli::cli_alert_success("All {nrow(grid)} cell(s) x {REPS} reps executed with zero failures.")
}

# ---- 3 metrics, per regime, with MC SEs ---------------------------------

cli::cli_h1("3. Metrics (spec §5), with Monte Carlo SEs")

smry <- summarise_cell(res)
utils::write.csv(smry, file.path(DIR_TABLES, sprintf("pilot_summary_%s.csv", STAMP)),
                 row.names = FALSE)

ratio_tabs <- list()
for (r in PILOT_REGIMES) {
  cli::cli_h2("{r} -- {REGIMES[[r]]$label}")
  cli::cli_inform("evidence for: {REGIMES[[r]]$evidence}")
  print_regime_summary(smry, r)
  rt <- print_ratio_summary(res, r)
  if (!is.null(rt)) ratio_tabs[[r]] <- rt
}
if (length(ratio_tabs)) {
  utils::write.csv(do.call(rbind, ratio_tabs),
                   file.path(DIR_TABLES, sprintf("pilot_rmse_ratios_%s.csv", STAMP)),
                   row.names = FALSE)
}

# ---- 4 sizing the real sweep --------------------------------------------

cli::cli_h1("4. Projected wall time for spec §6's targets over the FULL n grid")

# Per-replication cost scales roughly linearly in n for the tree fits and worse
# than linearly for the nested-CV crossfit arm, so the projection uses the MEASURED
# per-rep cost at the pilot's largest n in each regime and scales linearly in n
# from there. That UNDER-states the crossfit-heavy regimes; stated rather than
# silently optimistic.
proj <- lapply(PILOT_REGIMES, function(r) {
  sub <- res[res$regime == r, , drop = FALSE]
  per_rep <- tapply(sub$secs, list(sub$n, sub$dgp), function(z) sum(z, na.rm = TRUE))
  per_rep <- per_rep / REPS
  n_pilot_max <- max(as.integer(rownames(per_rep)))
  secs_at_max <- mean(per_rep[as.character(n_pilot_max), ], na.rm = TRUE)
  n_full <- N_GRID_FULL[[r]]
  n_dgps <- length(REGIMES[[r]]$dgps)
  target <- NSIM_FULL[[r]]
  secs_total <- sum(secs_at_max * (n_full / n_pilot_max)) * n_dgps * target
  data.frame(
    regime = r, nsim_target = target,
    n_cells_full = length(n_full) * n_dgps,
    secs_per_rep_at_n = secs_at_max, at_n = n_pilot_max,
    projected_hours = secs_total / 3600,
    stringsAsFactors = FALSE
  )
})
proj <- do.call(rbind, proj)
print(knitr::kable(proj, digits = 4, format = "simple", row.names = FALSE))
utils::write.csv(proj, file.path(DIR_TABLES, sprintf("pilot_projection_%s.csv", STAMP)),
                 row.names = FALSE)

cli::cli_inform(c(
  "*" = "Pilot wall time: {signif(elapsed_min, 4)} min for {nrow(grid)} cells x {REPS} reps.",
  "*" = "Projected FULL sweep (spec §6 nsim over the full n grid, sequential): {signif(sum(proj$projected_hours), 4)} h.",
  "i" = "Sequential and single-threaded throughout (worker_limit = 1, ranger num.threads = 1, parallel_cv = FALSE). Parallelising over CELLS -- not within a fit -- is the safe axis; see simulations/docs/MEMORY_SAFE_SIMULATIONS.md."
))

cli::cli_h1("Pilot complete. Read spec §0 before quoting any number above.")
