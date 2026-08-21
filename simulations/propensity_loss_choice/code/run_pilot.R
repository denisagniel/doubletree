# ============================================================
# Propensity-tree loss ablation -- PILOT
#
# Stem:    2026-08-20_propensity-loss-choice
# Purpose: (a) verify empirically that the chosen leaf_budget is large enough
#              that structure selection is not starved on EITHER DGP
#              (spec reconciliation note: "state the chosen value explicitly");
#          (b) time 10 paired replicates per cell and extrapolate the full run.
# Inputs:  code/common.R (which sources the two existing DGPs)
# Outputs: results/pilot_<runid>.rds  (+ console report)
#
# Run from the doubletree package root:
#   Rscript simulations/propensity_loss_choice/code/run_pilot.R
# ============================================================

# ---- 0 setup ----------------------------------------------------------------

source(file.path("simulations", "propensity_loss_choice", "code", "common.R"))

set.seed(SEED_MASTER)

PILOT_REPS     <- as.integer(Sys.getenv("PLC_PILOT_REPS", "10"))
BUDGET_GRID    <- c(2L, 3L, 4L, 6L, 8L)   # budget-adequacy probe only
BUDGET_PROBE_N <- 2000L
BUDGET_PROBE_R <- 5L                       # prototype fits per (dgp, budget)
cli::cli_inform(c(
  "gates: PLC_PILOT_REPS={PILOT_REPS}",
  "*" = "leaf_budget under test = {LEAF_BUDGET}; m_n = {M_N}; lambda_n = log(n)/n"
))

# ---- 1 how many leaves does EXACT representation actually need? -------------

# theory.tex ass:sparsity holds with delta = 0 iff some tree of at most
# leaf_budget leaves reproduces the true nuisance exactly. For binary X that is
# a finite question, so answer it exactly instead of eyeballing it: enumerate
# the covariate cells of {0,1}^p, attach the TRUE nuisance value of each, and
# find the minimum-leaf tree that separates the level sets.
min_exact_leaves <- function(patterns, value) {
  features <- colnames(patterns)
  recurse <- function(rows, feats) {
    if (dplyr::n_distinct(value[rows]) == 1L) return(1L)
    if (length(feats) == 0L) return(NA_integer_)   # not tree-representable
    best <- Inf
    for (f in feats) {
      is0 <- rows[patterns[rows, f] == 0]
      is1 <- rows[patterns[rows, f] == 1]
      if (length(is0) == 0L || length(is1) == 0L) next
      rest <- setdiff(feats, f)
      lo <- recurse(is0, rest)
      hi <- recurse(is1, rest)
      if (is.na(lo) || is.na(hi)) next
      best <- min(best, lo + hi)
    }
    if (is.infinite(best)) NA_integer_ else as.integer(best)
  }
  recurse(seq_len(nrow(patterns)), features)
}

#' Structural fingerprint of a DGP's true nuisances
#'
#' Uses one large draw purely to enumerate the covariate cells and read off the
#' true nuisance value attached to each (the DGPs expose `e_true`/`true_e` and
#' `mu0_true`/`true_m0`, so no DGP logic is reimplemented here).
dgp_fingerprint <- function(dgp, n = 20000L) {
  d <- generate_data(dgp, n, rep = 0L)
  m0 <- if (!is.null(d$mu0_true)) d$mu0_true else d$true_m0
  X <- as.matrix(d$X)

  cell_key <- apply(X, 1, paste, collapse = "")
  first <- !duplicated(cell_key)
  patterns <- X[first, , drop = FALSE]

  tibble::tibble(
    dgp                = dgp,
    p                  = ncol(X),
    cells_observed     = nrow(patterns),
    cells_possible     = 2^ncol(X),
    distinct_e_values  = dplyr::n_distinct(round(d$e_true_vec, 12)),
    distinct_m0_values = dplyr::n_distinct(round(m0, 12)),
    min_leaves_e       = min_exact_leaves(patterns, round(d$e_true_vec[first], 12)),
    min_leaves_m0      = min_exact_leaves(patterns, round(m0[first], 12)),
    min_true_e         = min(d$e_true_vec),
    max_true_e         = max(d$e_true_vec),
    prop_e_extreme     = mean(d$e_true_vec < 0.1 | d$e_true_vec > 0.9)
  )
}

fingerprint <- purrr::map_dfr(DGP_IDS, dgp_fingerprint)

cli::cli_h1("1. Exact-representation requirement of the true nuisances")
print(as.data.frame(fingerprint), digits = 4)
budget_ok <- all(fingerprint$min_leaves_e <= LEAF_BUDGET,
                 fingerprint$min_leaves_m0 <= LEAF_BUDGET)
max_needed <- max(c(fingerprint$min_leaves_e, fingerprint$min_leaves_m0))
cli::cli_inform(stats::setNames(
  paste0("leaf_budget = ", LEAF_BUDGET, " ",
         if (budget_ok) "covers" else "DOES NOT cover",
         " the exact-representation requirement (max leaves needed = ",
         max_needed, ")"),
  if (budget_ok) "v" else "x"
))

# ---- 2 empirical budget adequacy -------------------------------------------

# Sweep leaf_budget on a handful of prototype fits per DGP. The signature of an
# ADEQUATE budget is: n_leaves_e / n_leaves_m0 stop growing with the budget
# (selection stops being budget-constrained) and land at the cardinality above.
# A budget that BINDS at its ceiling (n_leaves == leaf_budget for every fit,
# still rising when the ceiling rises) is a starved budget.
probe_budget <- function(dgp, leaf_budget, rep, arm) {
  d <- generate_data(dgp, BUDGET_PROBE_N, rep)
  fit <- doubletree::estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    leaf_budget     = leaf_budget,
    outcome_type    = OUTCOME_TYPE,
    lambda_n        = LAMBDA_N,
    m_n             = M_N,
    propensity_loss = ARMS[[arm]],
    verbose         = FALSE
  )
  tibble::tibble(
    dgp = dgp, leaf_budget = leaf_budget, rep = rep, arm = arm,
    n_leaves_e = fit$n_leaves_e, n_leaves_m0 = fit$n_leaves_m0,
    certified_e = fit$certified_e, certified_m0 = fit$certified_m0,
    used_search_e = fit$used_search_e, used_search_m0 = fit$used_search_m0,
    theta = fit$theta, true_att = d$true_att,
    max_e_hat = max(fit$nuisance_fits$propensity)
  )
}

budget_grid <- tidyr::expand_grid(
  dgp         = DGP_IDS,
  leaf_budget = BUDGET_GRID,
  rep         = seq_len(BUDGET_PROBE_R),
  arm         = names(ARMS)
)

t_budget <- system.time(
  budget_probe <- purrr::pmap_dfr(budget_grid, probe_budget)
)

budget_summary <- budget_probe |>
  dplyr::summarise(
    mean_leaves_e   = mean(n_leaves_e),
    max_leaves_e    = max(n_leaves_e),
    mean_leaves_m0  = mean(n_leaves_m0),
    max_leaves_m0   = max(n_leaves_m0),
    pct_certified_e = mean(certified_e),
    pct_certified_m0 = mean(certified_m0),
    pct_search_e    = mean(used_search_e),
    pct_search_m0   = mean(used_search_m0),
    mean_bias       = mean(theta - true_att),
    .by = c(dgp, arm, leaf_budget)
  ) |>
  dplyr::arrange(dgp, arm, leaf_budget)

cli::cli_h1("2. Budget-adequacy sweep (n = {BUDGET_PROBE_N}, {BUDGET_PROBE_R} fits per cell)")
print(as.data.frame(budget_summary), digits = 3)
cli::cli_inform("budget sweep took {round(t_budget[['elapsed']], 1)}s")

# ---- 3 paired-replicate prototype and timing -------------------------------

pilot_grid <- tidyr::expand_grid(dgp = DGP_IDS, n = N_GRID)

time_cell <- function(dgp, n) {
  el <- system.time(res <- run_cell(dgp, n, PILOT_REPS, on_error = "stop", heartbeat = 0L))
  list(res = res, elapsed = el[["elapsed"]])
}

cli::cli_h1("3. Paired prototype: {PILOT_REPS} reps x {nrow(pilot_grid)} cells x {length(ARMS)} arms")
timed <- purrr::pmap(pilot_grid, time_cell)
pilot <- dplyr::bind_rows(purrr::map(timed, "res"))

timing <- pilot_grid |>
  dplyr::mutate(
    elapsed_s     = purrr::map_dbl(timed, "elapsed"),
    s_per_rep     = elapsed_s / PILOT_REPS,
    s_per_arm_fit = elapsed_s / (PILOT_REPS * length(ARMS))
  )
print(as.data.frame(timing), digits = 3)

full_reps <- 300L
projected_s <- sum(timing$s_per_rep) * full_reps
cli::cli_inform(c(
  "v" = "prototype total: {round(sum(timing$elapsed_s), 1)}s for {PILOT_REPS * nrow(pilot_grid) * length(ARMS)} estimate_att() calls",
  "i" = "projected SEQUENTIAL full run at R={full_reps}: {round(projected_s / 60, 1)} min ({full_reps * nrow(pilot_grid) * length(ARMS)} calls)",
  "i" = "projected with 4 cells in parallel: {round(max(timing$s_per_rep) * full_reps / 60, 1)} min"
))

# ---- 4 prototype metrics (sanity, not a result) ----------------------------

pilot_summary <- pilot |>
  dplyr::summarise(
    reps        = dplyr::n(),
    failures    = sum(!is.na(error_message)),
    bias        = mean(theta - true_att),
    rmse        = sqrt(mean((theta - true_att)^2)),
    coverage    = mean(covered),
    ci_width    = mean(ci_width),
    mean_leaves_e = mean(n_leaves_e),
    max_e_hat   = max(max_e_hat),
    .by = c(dgp, n, arm)
  ) |>
  dplyr::arrange(dgp, n, arm)

same_partition <- pilot |>
  dplyr::select(dgp, n, rep, arm, hash_e) |>
  tidyr::pivot_wider(names_from = arm, values_from = hash_e) |>
  dplyr::summarise(same_partition_rate = mean(L == S), .by = c(dgp, n))

cli::cli_h1("4. Prototype metrics ({PILOT_REPS} reps -- indicative only)")
print(as.data.frame(pilot_summary), digits = 3)
print(as.data.frame(same_partition), digits = 3)

# ---- 5 export --------------------------------------------------------------

meta <- run_metadata(PILOT_REPS, label = "pilot")
out_path <- file.path(DIR_RESULTS, paste0("pilot_", meta$runid, ".rds"))
saveRDS(
  list(
    meta = meta, fingerprint = fingerprint, budget_probe = budget_probe,
    budget_summary = budget_summary, pilot = pilot,
    pilot_summary = pilot_summary, same_partition = same_partition,
    timing = timing
  ),
  out_path
)
cli::cli_inform(c("v" = "wrote {.file {out_path}}"))
