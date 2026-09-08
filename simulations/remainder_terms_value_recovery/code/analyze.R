# ---------------------------------------------------------------------------
# Analysis: turn results/replications.rds into the report tables.
#
# Usage: Rscript code/analyze.R          (writes results/tables.md)
# ---------------------------------------------------------------------------

library(readr)
library(fs)
library(dplyr)
library(tidyr)

RESULT_DIR <- Sys.getenv("RTVR_RESULT_DIR", "results")
source(fs::path("code", "dgp.R"))
source(fs::path("code", "bounds.R"))

reps <- readr::read_rds(fs::path(RESULT_DIR, "replications.rds"))
bounds <- readr::read_rds(fs::path(RESULT_DIR, "bounds.rds"))
calib <- readr::read_rds(fs::path(RESULT_DIR, "lambda_calibration.rds"))

stopifnot(nrow(reps) > 0, all(c("T_a", "T_eps", "D_w", "D_mu") %in% names(reps)))
stopifnot("refit identity violated" = max(reps$identity_gap) < IDENTITY_TOL)

N_FEAT <- P_COVARIATES * N_THRESH

# ---- per-replication scaled quantities ------------------------------------
scaled <- reps |>
  mutate(
    rn_T_eps = sqrt(n) * abs(T_eps),
    rn_T_a = sqrt(n) * abs(T_a),
    rn_T = sqrt(n) * abs(T_eps - T_a)
  )

q <- function(x, p) unname(stats::quantile(x, p, names = FALSE))

# ---- Table 1: fit quality and realised discrepancies ---------------------
tab_fit <- scaled |>
  group_by(n, L) |>
  summarise(
    reps = dplyr::n(),
    leaves_e = mean(leaves_e), leaves_mu = mean(leaves_mu),
    refits = mean(refits_e + refits_mu),
    D_w_mean = mean(D_w), D_w_sd = sd(D_w),
    D_mu_mean = mean(D_mu), D_mu_sd = sd(D_mu),
    a_norm_mean = mean(a_norm_mu),
    .groups = "drop"
  )

# ---- Table 2: the remainder terms ----------------------------------------
tab_terms <- scaled |>
  group_by(n, L) |>
  summarise(
    T_a_mean = mean(T_a), T_a_sd = sd(T_a),
    T_eps_mean = mean(T_eps), T_eps_sd = sd(T_eps),
    rn_Ta_mean = mean(rn_T_a), rn_Ta_sd = sd(rn_T_a),
    rn_Ta_p50 = q(rn_T_a, 0.50), rn_Ta_p95 = q(rn_T_a, 0.95),
    rn_Te_mean = mean(rn_T_eps), rn_Te_sd = sd(rn_T_eps),
    rn_Te_p50 = q(rn_T_eps, 0.50), rn_Te_p95 = q(rn_T_eps, 0.95),
    rn_T_mean = mean(rn_T), rn_T_p95 = q(rn_T, 0.95), rn_T_max = max(rn_T),
    .groups = "drop"
  )

# ---- Table 3: empirical vs the crude theoretical bound -------------------
tab_bound <- tab_terms |>
  left_join(bounds, by = c("n", "L")) |>
  left_join(tab_fit |> select(n, L, D_w_mean, D_mu_mean), by = c("n", "L")) |>
  rowwise() |>
  mutate(
    # lemma bound evaluated at the REALISED mean discrepancies
    bd_Te_realised = lemma_bound_scaled(n, rho_n, D_w_mean, D_mu_mean)$T_eps,
    bd_Ta_realised = lemma_bound_scaled(n, rho_n, D_w_mean, D_mu_mean)$T_a,
    # pure worst case: D ~ rho_n (the form the audit called vacuous)
    bd_Te_worst = lemma_bound_scaled(n, rho_n, rho_n, rho_n)$T_eps,
    bd_Ta_worst = lemma_bound_scaled(n, rho_n, rho_n, rho_n)$T_a,
    # same, with the exact log-cardinality radius r_n of the class searched
    bd_Ta_card = lemma_bound_scaled(n, r_n, D_w_mean, D_mu_mean)$T_a,
    bd_Te_card = lemma_bound_scaled(n, r_n, D_w_mean, D_mu_mean)$T_eps
  ) |>
  ungroup() |>
  mutate(
    ratio_Ta_realised = rn_Ta_mean / bd_Ta_realised,
    ratio_Ta_worst = rn_Ta_mean / bd_Ta_worst,
    ratio_Te_realised = rn_Te_mean / bd_Te_realised,
    ratio_Te_worst = rn_Te_mean / bd_Te_worst
  )

# ---- Table 4: ATT inference ----------------------------------------------
tab_att <- scaled |>
  group_by(n, L) |>
  summarise(
    bias = mean(theta_hat) - unique(TAU_ATT),
    sd_theta = sd(theta_hat),
    mean_se = mean(se),
    coverage = mean(covered),
    ci_width = mean(2 * stats::qnorm(0.975) * se),
    .groups = "drop"
  )

# ---- n-trend: does sqrt(n)|T| shrink as n grows at fixed Lbar? -----------
tab_trend <- tab_terms |>
  select(n, L, rn_T_mean) |>
  tidyr::pivot_wider(names_from = n, values_from = rn_T_mean,
                     names_prefix = "n_") |>
  arrange(L)

# ---- render ---------------------------------------------------------------
fmt <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) signif(x, digits))
  knitr_kable(df)
}
knitr_kable <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

cal <- vc_calibration(P_COVARIATES)
lines <- c(
  "# Remainder-term study: result tables",
  "",
  sprintf("Generated %s from `%s`.", format(Sys.time(), "%Y-%m-%d %H:%M"), RESULT_DIR),
  sprintf("Replications per cell: %d. Max refit-identity gap across all cells: %.1e.",
          min(tab_fit$reps), max(reps$identity_gap)),
  sprintf("Grid: p = %d, %d threshold(s)/coordinate, %d binary features, %d grid cells.",
          P_COVARIATES, N_THRESH, N_FEAT, N_GRID_CELL),
  sprintf("VC calibration: v(L) = %.3f * L * log2(%dL) + %.2f (anchors 185/427/686 at L = 10/20/30).",
          cal$c1, P_COVARIATES, cal$c0),
  "",
  "## Table 0. Lambda calibration and realised leaf counts",
  "",
  fmt(calib |> select(n, L, depth, lambda_e, leaves_e_cal, lambda_mu, leaves_mu_cal)),
  "",
  "## Table 1. Fit quality and exact population discrepancies",
  "",
  "`D_w = ||w_hat - w_0||_mu`, `D_mu = ||mu_hat - mu_0||_mu`, `a_norm = ||a_mu||_mu`;",
  "all computed exactly by summation over grid cells (no Monte-Carlo error).",
  "",
  fmt(tab_fit),
  "",
  "## Table 2. Remainder terms, raw and root-n scaled",
  "",
  "The theorem needs `sqrt(n)(T_eps - T_a) = o_p(1)`; `rn_*` columns are",
  "`sqrt(n) * |.|`.",
  "",
  fmt(tab_terms |> select(n, L, T_a_mean, T_a_sd, T_eps_mean, T_eps_sd)),
  "",
  fmt(tab_terms |> select(n, L, rn_Ta_mean, rn_Ta_sd, rn_Ta_p50, rn_Ta_p95,
                          rn_Te_mean, rn_Te_sd, rn_Te_p50, rn_Te_p95)),
  "",
  fmt(tab_terms |> select(n, L, rn_T_mean, rn_T_p95, rn_T_max)),
  "",
  "## Table 3. Empirical size versus the crude theoretical bound",
  "",
  "`rho_n = sqrt(kappa_n / n)` with `kappa_n = v * log(2 e n / v) + log(4n)`;",
  "`r_n = sqrt((1 + log|T|) / n)` uses the EXACT cardinality of the",
  "depth-restricted class actually searched. `bd_*_worst` sets D ~ rho_n",
  "(the form the audit found vacuous); `bd_*_realised` evaluates the same",
  "lemma expression at the realised mean D_w, D_mu.",
  "",
  fmt(tab_bound |> select(n, L, v_calibrated, kappa_n, rho_n, log_card, r_n)),
  "",
  fmt(tab_bound |> select(n, L, rn_Ta_mean, bd_Ta_realised, ratio_Ta_realised,
                          bd_Ta_worst, ratio_Ta_worst, bd_Ta_card)),
  "",
  fmt(tab_bound |> select(n, L, rn_Te_mean, bd_Te_realised, ratio_Te_realised,
                          bd_Te_worst, ratio_Te_worst, bd_Te_card)),
  "",
  "## Table 4. ATT point estimation and Wald coverage (nominal 95%)",
  "",
  fmt(tab_att),
  "",
  "## Table 5. sqrt(n)|T_eps - T_a| as n grows, at fixed leaf budget",
  "",
  fmt(tab_trend),
  ""
)

writeLines(lines, fs::path(RESULT_DIR, "tables.md"))
readr::write_rds(list(fit = tab_fit, terms = tab_terms, bound = tab_bound,
                      att = tab_att, trend = tab_trend),
                 fs::path(RESULT_DIR, "summaries.rds"))
cat("wrote", as.character(fs::path(RESULT_DIR, "tables.md")), "\n")
print(as.data.frame(tab_trend))
