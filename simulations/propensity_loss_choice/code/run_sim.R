# ============================================================
# Propensity-tree loss ablation -- FULL RUN
#
# Stem:    2026-08-20_propensity-loss-choice
# Purpose: 2 DGPs x n in {500, 2000} x R reps x 2 arms of `propensity_loss`,
#          paired on the data (identical dataset per replicate, identical
#          leaf_budget / lambda_n / m_n / outcome tree).
# Inputs:  code/common.R
# Outputs: results/results_<runid>.rds
#
# Run from the doubletree package root:
#   Rscript simulations/propensity_loss_choice/code/run_sim.R
#
# Gates (all defaulted; sourcing this file runs the real study):
#   PLC_REPS   number of replicates per (dgp, n) cell   [300]
#   PLC_CELLS  comma-separated "dgp:n" cells to run     [all four]
# ============================================================

# ---- 0 setup ----------------------------------------------------------------

source(file.path("simulations", "propensity_loss_choice", "code", "common.R"))

set.seed(SEED_MASTER)

N_REPS      <- as.integer(Sys.getenv("PLC_REPS", "300"))
CELLS_GATE  <- Sys.getenv("PLC_CELLS", "")

cells <- tidyr::expand_grid(dgp = DGP_IDS, n = N_GRID)
if (nzchar(CELLS_GATE)) {
  wanted <- strsplit(CELLS_GATE, ",", fixed = TRUE)[[1]]
  cells <- cells |>
    dplyr::filter(paste(dgp, n, sep = ":") %in% trimws(wanted))
  if (nrow(cells) == 0L) {
    cli::cli_abort("PLC_CELLS={CELLS_GATE} matched no cell of {.val {paste(DGP_IDS, collapse='/')}} x {.val {N_GRID}}.")
  }
}

cli::cli_inform(c(
  "gates: PLC_REPS={N_REPS} PLC_CELLS={if (nzchar(CELLS_GATE)) CELLS_GATE else '<all>'}",
  "*" = "leaf_budget={LEAF_BUDGET} m_n={M_N} lambda_n=log(n)/n outcome_type={OUTCOME_TYPE}",
  "*" = "{nrow(cells)} cells x {N_REPS} reps x {length(ARMS)} arms = {nrow(cells) * N_REPS * length(ARMS)} estimate_att() calls"
))

# ---- 1 run ------------------------------------------------------------------

# on_error = "record" keeps the condition MESSAGE (never a plausible-looking NA
# estimate) so a failure RATE can be reported; section 2 then refuses to hand a
# contaminated run onward.
t_run <- system.time(
  results <- purrr::pmap_dfr(cells, function(dgp, n) {
    cli::cli_inform("cell {dgp} n={n}")
    run_cell(dgp, n, N_REPS, on_error = "record", heartbeat = 100L)
  })
)
cli::cli_inform("run took {round(t_run[['elapsed']] / 60, 2)} min")

# ---- 2 integrity checks ----------------------------------------------------

expected_rows <- nrow(cells) * N_REPS * length(ARMS)
if (nrow(results) != expected_rows) {
  cli::cli_abort("Expected {expected_rows} rows, got {nrow(results)}.")
}

# The paired design is the whole point: verify the two arms of each replicate
# really were fit to the same dataset.
pairing <- results |>
  dplyr::filter(is.na(error_message)) |>
  dplyr::summarise(n_hashes = dplyr::n_distinct(data_hash), .by = c(dgp, n, rep))
if (any(pairing$n_hashes != 1L)) {
  bad <- pairing[pairing$n_hashes != 1L, ]
  cli::cli_abort(c(
    "Paired design violated: {nrow(bad)} replicate(s) have arm-specific data.",
    i = "First: dgp={bad$dgp[1]} n={bad$n[1]} rep={bad$rep[1]}"
  ))
}

failures <- results |> dplyr::filter(!is.na(error_message))

# ---- 3 export ---------------------------------------------------------------

meta <- run_metadata(N_REPS, label = "full")
meta$cells <- cells
meta$elapsed_min <- t_run[["elapsed"]] / 60
out_path <- file.path(DIR_RESULTS, paste0("results_", meta$runid, ".rds"))
saveRDS(list(meta = meta, results = results), out_path)
cli::cli_inform(c("v" = "wrote {.file {out_path}}"))

# ---- 4 failure gate --------------------------------------------------------

# Spec decision rule, item 3: a non-trivial failure rate in either arm means
# this sim is not informative, so stop instead of summarising noise. The results
# file is written first so the failures are diagnosable.
if (nrow(failures) > 0L) {
  rate <- results |>
    dplyr::summarise(failure_rate = mean(!is.na(error_message)), .by = c(dgp, n, arm))
  print(as.data.frame(rate))
  cli::cli_abort(c(
    "{nrow(failures)} of {nrow(results)} arm-fits failed -- spec decision rule item 3.",
    x = "Distinct messages:",
    stats::setNames(unique(failures$error_message),
                    rep("*", length(unique(failures$error_message)))),
    i = "Results saved to {.file {out_path}} for diagnosis."
  ))
}
cli::cli_inform(c("v" = "0 failures in {nrow(results)} arm-fits"))
