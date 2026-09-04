# ============================================================
# run_pilot.R
# Study: partition_recovery_clt  (doubletree, S1)
#
# Prototype BEFORE committing to the full grid, per spec §5 ("confirm actual
# per-fit cost empirically on a handful of fits before committing to the full
# 1000-replication grid") and this project's own r-code conventions ("Prototype
# on a small subset or fewer iterations first; estimate run time").
#
# Run:  Rscript simulations/partition_recovery_clt/code/run_pilot.R [reps]
#       (from the doubletree package root; default reps = 5)
#
# Does FOUR things and then stops. It deliberately does not roll into the full
# sweep: the sizing decision is the output.
#
#   1. Verifies, for EVERY DGP, that the population spec agrees with the sample
#      draw and that each nuisance depends on exactly the coordinates it
#      declares -- so the enumeration is anchored to the DGP the data comes from.
#   2. Prints each DGP's exact population structure -- |T_Lbar|, |S_e|, |S_mu|,
#      Delta_j in both currencies -- next to lambda_n = log(n)/n at every n, so
#      the "penalty excludes the necessary split" cells are identified BEFORE any
#      replication runs rather than diagnosed afterwards.
#   3. Runs `reps` replications in every cell of the real design grid, measuring
#      per-fit wall time.
#   4. Extrapolates to the spec's target of 1000 reps/cell and prints the
#      projected wall time, with an explicit verdict against a 2-3 hour budget.
# ============================================================

source(file.path("simulations", "partition_recovery_clt", "code", "common.R"))

suppressPackageStartupMessages(library(knitr))

REPS <- local({
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1L && !is.na(suppressWarnings(as.integer(a[[1L]])))) {
    as.integer(a[[1L]])
  } else {
    5L
  }
})
TARGET_REPS <- 1000L   # spec §5

set.seed(SEED_MASTER)

# ---- 1 anchor every DGP's enumeration to its own sample draw -------------

cli::cli_h1("1. Population spec vs sample draw, and declared coordinate roles")
dgps <- lapply(DGP_IDS, make_dgp)   # build_dgp() runs check_variable_roles()
names(dgps) <- DGP_IDS
for (id in DGP_IDS) {
  cmp <- check_spec_vs_draw(dgps[[id]], n = 20000L)
  cli::cli_alert_success(
    "{id}: spec matches its draw on all {nrow(cmp)} cells (tol 1e-12); roles e={paste(dgps[[id]]$e_vars, collapse='+')}, mu={paste(dgps[[id]]$mu_vars, collapse='+')} verified."
  )
}

# ---- 2 exact population structure, before any replication ---------------

cli::cli_h1("2. Exact population structure per DGP (closed form, no simulation)")
pop <- do.call(rbind, lapply(dgps, dgp_population_summary))
rownames(pop) <- NULL
print(knitr::kable(
  pop[, c("dgp", "nuisance", "theta0", "p", "card_T", "n_sufficient",
          "min_sufficient_leaves", "n_distinct_values", "delta_sq", "delta_kl",
          "e0_min", "e0_max")],
  digits = 6, format = "simple"
))

cli::cli_h2("lambda_n = log(n)/n vs Delta_j(kl): is the PENALISED optimum in S_j?")
# The comparison that separates two failure mechanisms the spec insists must not
# be conflated: lambda_n >= Delta_j means the penalised population optimum is
# outside S_j (structural failure, no noise required); lambda_n < Delta_j means
# any observed failure is sampling noise, which is what lem:selection is about.
lam <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  do.call(rbind, lapply(N_GRID, function(n) {
    data.frame(
      dgp = id, regime = regime_of(id, n), n = n,
      lambda_n = log(n) / n,
      delta_e_kl = d$class_e$delta_kl,
      delta_mu_kl = d$class_mu$delta_kl,
      penalised_opt_in_S_e = log(n) / n < d$class_e$delta_kl,
      penalised_opt_in_S_mu = log(n) / n < d$class_mu$delta_kl,
      stringsAsFactors = FALSE
    )
  }))
}))
print(knitr::kable(lam, digits = 5, format = "simple"))

# ---- 3 pilot replications across the whole grid --------------------------

grid <- design_grid()
cli::cli_h1("3. Pilot: {REPS} reps in each of {nrow(grid)} cells")
t_all <- Sys.time()
pilot <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  d <- dgps[[g$dgp]]
  # on_error = "stop": in the pilot a failure must surface at its origin.
  rows <- lapply(seq_len(REPS), function(r) run_one_rep(d, g$n, r, on_error = "stop"))
  pilot[[i]] <- do.call(rbind, rows)
  gc(full = TRUE, verbose = FALSE)
  cur <- pilot[[i]]
  cli::cli_inform(paste0(
    "  ", g$regime, " ", g$dgp, " n=", g$n,
    "  rec_e=", signif(mean(cur$recovered_e), 3),
    " rec_mu=", signif(mean(cur$recovered_mu), 3),
    " cov=", signif(mean(cur$covered), 3),
    " leaves_e=", signif(mean(cur$n_leaves_e), 3),
    " clipfrac=", signif(mean(cur$clip_frac), 3),
    "  ", signif(mean(cur$secs), 3), " s/fit"
  ))
}
pilot <- do.call(rbind, pilot)
elapsed_min <- as.numeric(difftime(Sys.time(), t_all, units = "mins"))

# Integrity: the fitted partition must lie in the class the margin was computed
# over, or the bound comparison in analyze.R is meaningless.
if (!all(pilot$in_class_e) || !all(pilot$in_class_mu)) {
  cli::cli_abort(c(
    "Some fitted partitions fall OUTSIDE the enumerated T_Lbar.",
    x = "e: {sum(!pilot$in_class_e)} of {nrow(pilot)}; mu: {sum(!pilot$in_class_mu)} of {nrow(pilot)}.",
    i = "Delta_j was computed over T_Lbar, so the selection-bound comparison would be invalid."
  ))
}
cli::cli_alert_success("All {nrow(pilot)} fitted partitions lie inside the enumerated T_Lbar.")

wm <- unique(stats::na.omit(pilot$warning_messages))
if (length(wm)) {
  cli::cli_inform(c("i" = "{length(wm)} distinct muffled warning text(s) across the pilot:",
                    stats::setNames(substr(wm, 1, 150), rep("*", length(wm)))))
}

# ---- 4 timing projection ------------------------------------------------

cli::cli_h1("4. Timing projection to the spec's {TARGET_REPS} reps/cell")
timing <- do.call(rbind, lapply(
  split(pilot, list(pilot$dgp, pilot$n), drop = TRUE),
  function(z) data.frame(
    dgp = z$dgp[[1]], regime = z$regime[[1]], n = z$n[[1]], reps = nrow(z),
    secs_per_fit = mean(z$secs),
    projected_min = mean(z$secs) * TARGET_REPS / 60,
    stringsAsFactors = FALSE
  )
))
timing <- timing[order(match(timing$dgp, DGP_IDS), timing$n), ]
rownames(timing) <- NULL
print(knitr::kable(timing, digits = 4, format = "simple"))

total_h <- sum(timing$projected_min) / 60
cli::cli_inform(c(
  "*" = "pilot wall time: {signif(elapsed_min, 3)} min for {nrow(pilot)} fits",
  "*" = "mean {signif(mean(pilot$secs), 4)} s/fit overall; slowest cell {signif(max(timing$secs_per_fit), 4)} s/fit",
  "!" = "PROJECTED FULL GRID at {TARGET_REPS} reps/cell: {signif(total_h, 3)} hours ({nrow(grid) * TARGET_REPS} fits)"
))
if (total_h > 3) {
  cli::cli_alert_danger(
    "Projected {signif(total_h, 3)} h exceeds the 2-3 h budget. STOP and report before sweeping."
  )
} else {
  cli::cli_alert_success(
    "Projected {signif(total_h, 3)} h is within the 2-3 h budget; the full sweep can proceed."
  )
}

out <- file.path(DIR_RESULTS, sprintf("pilot_r%d_%s.rds", REPS,
                                      format(Sys.time(), "%Y%m%d-%H%M%S")))
saveRDS(list(results = pilot, population = pop, lambda_vs_delta = lam,
             timing = timing, target_reps = TARGET_REPS,
             elapsed_min = elapsed_min,
             meta = cell_metadata(dgps[["f1"]], NA_integer_, REPS, NA_real_,
                                  "pilot", pilot)),
        out)
cli::cli_alert_success("wrote {.file {out}}")
