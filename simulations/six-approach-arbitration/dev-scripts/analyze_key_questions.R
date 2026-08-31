#!/usr/bin/env Rscript
# =============================================================================
# analyze_key_questions.R -- manuscript-focused analysis of the arbitration run
# =============================================================================
# Targets the four questions the paper needs answered:
#
#   Q1. Does Alt A (single_tree) achieve nominal coverage? Under which DGPs / n?
#   Q2. Does Alt B (honest CI, all Rashomon methods) always restore coverage?
#   Q3. Is |delta/SE| a reliable diagnostic for when Alt A is adequate?
#   Q4. RMSE cost of interpretability: single tree vs crossfit baseline.
#
# Bonus: dt_averaged vs single_tree -- are they interchangeable?
#
# Usage:
#   Rscript dev-scripts/analyze_key_questions.R [path/to/run.rds]
#
# Schema (from run_one.R):
#   estimate, std_error, ci_lower, ci_upper, truth, error, covered
#   theta_crossfit, se_crossfit, delta, delta_over_se, covered_crossfit
#   intersection_nonempty, converged, rashomon_c_e, rashomon_c_m0
#   method, n, dgp, escalate
#
# Coverage target: 0.95. Truth: 0.15.
# Key methods:
#   crossfit        -- valid baseline (K separate trees)
#   single_tree     -- Alt A: 1 tree, honest CI from twin
#   dt_averaged     -- Alt B: 1 averaged tree, honest CI
#   msplit_averaged -- Alt B (modal): 1 averaged tree, honest CI
#   full            -- biased baseline
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
rds  <- if (length(args) >= 1) args[1] else {
  cand <- c(
    "20260722-114141_fffe596.rds",
    "results/20260722-114141_fffe596.rds"
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) stop("No run rds found; pass a path.", call. = FALSE)
  hit[1]
}
cat(sprintf("Reading: %s\n\n", rds))
r <- readRDS(rds)

# Source grid for INFEASIBLE_CELLS and DGP ordering
study_root <- {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa) == 1) dirname(dirname(normalizePath(sub("^--file=", "", fa)))) else "."
}
source(file.path(study_root, "config", "grid.R"))

# =============================================================================
# Helpers
# =============================================================================

cov_se <- function(x) {
  x <- x[is.finite(x)]
  p <- mean(x)
  sqrt(p * (1 - p) / length(x))
}

# DGP display order and labels
DGP_ORDER  <- c("simple", "moderate", "complex", "continuous")
DGP_LABELS <- c(simple = "Simple", moderate = "Moderate",
                complex = "Complex", continuous = "Continuous")

# Only non-escalation cells; drop feasibility-excluded rows
r_main <- r |>
  filter(!escalate) |>
  filter(is_feasible(data.frame(method = method, n = n, dgp = dgp,
                                escalate = escalate)))

# Minimum reps to report a cell (MC error on coverage ~ ±1.4pp at 1000 reps)
MIN_REPS <- 500L

# Cell summary: one row per (method, n, dgp)
cell_summary <- function(d, group_cols = c("method", "n", "dgp")) {
  d |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      reps          = sum(is.finite(estimate)),
      bias          = mean(error,         na.rm = TRUE),
      emp_sd        = sd(estimate,        na.rm = TRUE),
      rmse          = sqrt(mean(error^2,  na.rm = TRUE)),
      mean_se       = mean(std_error,     na.rm = TRUE),
      se_ratio      = mean_se / emp_sd,
      cover         = mean(covered,       na.rm = TRUE),   # reported CI (honest for Rashomon)
      cover_se      = cov_se(covered),
      cover_cf      = mean(covered_crossfit, na.rm = TRUE), # twin CI
      conv          = mean(converged,     na.rm = TRUE),
      int_ne        = mean(intersection_nonempty, na.rm = TRUE),
      abs_delta_med = median(abs(delta),  na.rm = TRUE),
      dse_med       = median(abs(delta_over_se), na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(reps >= MIN_REPS) |>
    mutate(
      dgp = factor(dgp, levels = DGP_ORDER),
      cover_fmt = sprintf("%.3f (±%.3f)", cover, 1.96 * cover_se)
    ) |>
    arrange(method, n, dgp)
}

fmt_num <- function(x, digits = 3) formatC(round(x, digits), format = "f", digits = digits)

# =============================================================================
# Q1 + Q2. Coverage: Alt A (single_tree) and Alt B (honest CI methods)
#           vs crossfit baseline, by DGP and n
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("Q1/Q2. COVERAGE: single_tree, dt_averaged, msplit_averaged vs crossfit\n")
cat("       Reported CI = honest bias-aware CI for Rashomon methods\n")
cat("=" |> strrep(72), "\n\n")

methods_cov <- c("crossfit", "single_tree", "dt_averaged",
                 "msplit_averaged", "full")
cs <- r_main |>
  filter(method %in% methods_cov) |>
  cell_summary()

# Wide table: one row per (n, dgp), columns per method
cov_wide <- cs |>
  select(method, n, dgp, reps, cover, cover_se) |>
  pivot_wider(
    names_from  = method,
    values_from = c(cover, cover_se, reps),
    names_glue  = "{method}_{.value}"
  )

cat("Coverage by (n, DGP) — target 0.950\n")
cat("Columns: crossfit | single_tree (Alt A) | dt_averaged (Alt B) | msplit_averaged (Alt B) | full (biased baseline)\n\n")

for (dgp_val in DGP_ORDER) {
  sub <- cov_wide |> filter(dgp == dgp_val)
  if (nrow(sub) == 0) {
    cat(sprintf("  %s: no adequate cells\n", DGP_LABELS[dgp_val]))
    next
  }
  cat(sprintf("DGP: %s\n", DGP_LABELS[dgp_val]))
  cat(sprintf("  %-6s  %-22s  %-22s  %-22s  %-22s  %-22s\n",
              "n",
              "crossfit",
              "single_tree (Alt A)",
              "dt_averaged (Alt B)",
              "msplit_avg  (Alt B)",
              "full (baseline)"))
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    fmt_cell <- function(meth) {
      cov <- row[[paste0(meth, "_cover")]]
      se  <- row[[paste0(meth, "_cover_se")]]
      rps <- row[[paste0(meth, "_reps")]]
      if (is.na(cov)) return(sprintf("%-22s", "—"))
      flag <- if (!is.na(cov) && abs(cov - 0.95) > 1.96 * se) " *" else "  "
      sprintf("%.3f ±%.3f (n=%d)%s", cov, 1.96 * se, as.integer(rps), flag)
    }
    cat(sprintf("  %-6d  %-22s  %-22s  %-22s  %-22s  %-22s\n",
                row$n,
                fmt_cell("crossfit"),
                fmt_cell("single_tree"),
                fmt_cell("dt_averaged"),
                fmt_cell("msplit_averaged"),
                fmt_cell("full")))
  }
  cat("\n")
}
cat("* = outside 95% CI for nominal 0.95 coverage (|cov - 0.95| > 1.96 * MC_SE)\n\n")

# =============================================================================
# Q1 continued: single_tree convergence rate (intersection non-empty %)
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("Q1b. SINGLE_TREE CONVERGENCE: intersection non-empty rate\n")
cat("     (when 0, method fell back to fold-specific; no single tree returned)\n")
cat("=" |> strrep(72), "\n\n")

conv_tab <- r_main |>
  filter(method == "single_tree") |>
  group_by(n, dgp) |>
  summarise(
    reps    = sum(is.finite(estimate)),
    int_ne  = mean(intersection_nonempty, na.rm = TRUE),
    conv    = mean(converged, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(reps >= MIN_REPS) |>
  mutate(dgp = factor(dgp, levels = DGP_ORDER)) |>
  arrange(dgp, n)

print(conv_tab |> rename(
  `n` = n, DGP = dgp,
  `reps` = reps,
  `P(intersection non-empty)` = int_ne,
  `P(converged)` = conv
), row.names = FALSE)
cat("\n")

# =============================================================================
# Q3. Fidelity diagnostic: is |delta/SE| a reliable signal?
#     For single_tree: bin reps by |delta/SE| quartile and check coverage.
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("Q3. FIDELITY DIAGNOSTIC: does |delta/SE| predict Alt A coverage?\n")
cat("    Hypothesis: when |delta/SE| is small, Alt A coverage ~ nominal;\n")
cat("    when large, Alt A undercovers (honest CI is doing the work).\n")
cat("=" |> strrep(72), "\n\n")

st <- r_main |>
  filter(method == "single_tree",
         is.finite(delta_over_se),
         is.finite(covered))

if (nrow(st) < 100) {
  cat("Insufficient single_tree rows with finite delta_over_se.\n\n")
} else {
  # Wald coverage (what Alt A without honest CI would give) vs honest CI coverage
  # Wald: estimate +/- 1.96 * se_crossfit (the twin SE used in the honest CI)
  z <- qnorm(0.975)
  st <- st |>
    mutate(
      wald_lower  = estimate - z * se_crossfit,
      wald_upper  = estimate + z * se_crossfit,
      covered_wald = as.integer(truth >= wald_lower & truth <= wald_upper),
      abs_dse     = abs(delta_over_se)
    )

  # Bin |delta/SE| into quartiles overall
  st <- st |>
    mutate(dse_quartile = cut(abs_dse,
                              breaks = quantile(abs_dse, probs = c(0, .25, .5, .75, 1),
                                               na.rm = TRUE),
                              include.lowest = TRUE,
                              labels = c("Q1 (small)", "Q2", "Q3", "Q4 (large)")))

  dse_tab <- st |>
    group_by(dse_quartile) |>
    summarise(
      n_reps         = n(),
      dse_range      = sprintf("[%.2f, %.2f]",
                               min(abs_dse, na.rm = TRUE),
                               max(abs_dse, na.rm = TRUE)),
      cover_honest   = round(mean(covered,      na.rm = TRUE), 3),
      cover_wald     = round(mean(covered_wald, na.rm = TRUE), 3),
      cover_cf_twin  = round(mean(covered_crossfit, na.rm = TRUE), 3),
      .groups = "drop"
    )

  cat("|delta/SE| quartile analysis (pooled across DGP and n)\n")
  cat("cover_honest  = Alt A with honest CI (reported)\n")
  cat("cover_wald    = Alt A with Wald CI using twin SE (counterfactual)\n")
  cat("cover_cf_twin = cross-fit twin Wald CI (what you'd get from Alt B inference)\n\n")
  print(dse_tab, row.names = FALSE)
  cat("\n")

  # Same analysis by DGP
  dse_dgp <- st |>
    group_by(dgp, dse_quartile) |>
    summarise(
      n_reps       = n(),
      cover_honest = round(mean(covered,      na.rm = TRUE), 3),
      cover_wald   = round(mean(covered_wald, na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    mutate(dgp = factor(dgp, levels = DGP_ORDER)) |>
    arrange(dgp, dse_quartile)

  cat("|delta/SE| vs coverage by DGP:\n")
  print(dse_dgp, row.names = FALSE)
  cat("\n")

  # Summary: what fraction of reps have |delta/SE| < 0.5, 1.0?
  cat("Fraction of single_tree reps with |delta/SE| below threshold:\n")
  for (thresh in c(0.25, 0.5, 1.0, 2.0)) {
    pct <- mean(st$abs_dse < thresh, na.rm = TRUE)
    cat(sprintf("  |delta/SE| < %.2f : %.1f%%\n", thresh, 100 * pct))
  }
  cat("\n")
}

# =============================================================================
# Q4. RMSE cost of interpretability
#     single_tree and dt_averaged vs crossfit benchmark, by DGP and n
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("Q4. RMSE COST OF INTERPRETABILITY: single_tree / dt_averaged vs crossfit\n")
cat("    RMSE ratio > 1 = cost; ratio = 1 = no cost\n")
cat("=" |> strrep(72), "\n\n")

rmse_methods <- c("crossfit", "single_tree", "dt_averaged", "msplit_averaged")
rmse_tab <- r_main |>
  filter(method %in% rmse_methods) |>
  cell_summary() |>
  select(method, n, dgp, reps, bias, rmse) |>
  pivot_wider(names_from = method, values_from = c(rmse, bias, reps),
              names_glue = "{method}_{.value}") |>
  mutate(
    rmse_ratio_single  = round(single_tree_rmse  / crossfit_rmse, 3),
    rmse_ratio_dtavg   = round(dt_averaged_rmse  / crossfit_rmse, 3),
    rmse_ratio_msplitavg = round(msplit_averaged_rmse / crossfit_rmse, 3)
  )

cat("RMSE ratios (vs crossfit baseline; 1.00 = no cost, >1 = worse)\n\n")
cat(sprintf("  %-6s  %-12s  %-12s  %-12s  %-12s\n",
            "n", "DGP", "single_tree", "dt_averaged", "msplit_avg"))
for (i in seq_len(nrow(rmse_tab))) {
  row <- rmse_tab[i, ]
  cat(sprintf("  %-6d  %-12s  %-12s  %-12s  %-12s\n",
              row$n, as.character(row$dgp),
              fmt_num(row$rmse_ratio_single),
              fmt_num(row$rmse_ratio_dtavg),
              fmt_num(row$rmse_ratio_msplitavg)))
}
cat("\n")

# =============================================================================
# BONUS: dt_averaged vs single_tree -- are they interchangeable?
#        Direct comparison on cells where both have >= MIN_REPS
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("BONUS. dt_averaged vs single_tree: interchangeable?\n")
cat("       Both use intersection structure + honest CI. Difference: leaf strategy.\n")
cat("       single_tree: leaves refit on all n (Alt A theory)\n")
cat("       dt_averaged: leaves averaged across K cross-fit refits\n")
cat("=" |> strrep(72), "\n\n")

comp_methods <- c("single_tree", "dt_averaged")
comp_tab <- r_main |>
  filter(method %in% comp_methods) |>
  cell_summary() |>
  select(method, n, dgp, reps, bias, rmse, cover, se_ratio) |>
  pivot_wider(names_from = method,
              values_from = c(bias, rmse, cover, se_ratio, reps),
              names_glue = "{method}_{.value}")

cat(sprintf("  %-6s  %-12s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n",
            "n", "DGP",
            "st_cov", "dtavg_cov",
            "st_rmse", "dtavg_rmse",
            "st_bias", "dtavg_bias"))
for (i in seq_len(nrow(comp_tab))) {
  row <- comp_tab[i, ]
  cat(sprintf("  %-6d  %-12s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n",
              row$n, as.character(row$dgp),
              fmt_num(row$single_tree_cover),
              fmt_num(row$dt_averaged_cover),
              fmt_num(row$single_tree_rmse),
              fmt_num(row$dt_averaged_rmse),
              fmt_num(row$single_tree_bias),
              fmt_num(row$dt_averaged_bias)))
}
cat("\n")

# =============================================================================
# ESCALATION SWEEP: coverage vs c for complex DGPs
# (only Rashomon methods with escalate = TRUE)
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("ESCALATION SWEEP: coverage vs Rashomon tolerance multiplier c\n")
cat("  epsilon_n = c * log(n)/n; c=1 is the theory value\n")
cat("  Does larger c (to force non-empty intersection) degrade coverage?\n")
cat("=" |> strrep(72), "\n\n")

r_esc <- r |>
  filter(escalate == TRUE,
         method %in% c("doubletree", "dt_averaged", "single_tree"),
         is.finite(estimate))

if (nrow(r_esc) < 100) {
  cat("Insufficient escalation rows.\n\n")
} else {
  # Bin by realized c (max of c_e, c_m0)
  r_esc <- r_esc |>
    mutate(c_max = pmax(rashomon_c_e, rashomon_c_m0, na.rm = TRUE),
           c_bin = cut(c_max,
                       breaks = c(0, 1, 2, 8, 32, 128, Inf),
                       labels = c("c=1", "c=2", "c=3-8", "c=9-32",
                                  "c=33-128", "c>128"),
                       include.lowest = TRUE))

  esc_tab <- r_esc |>
    group_by(method, n, dgp, c_bin) |>
    summarise(
      reps  = n(),
      cover = round(mean(covered, na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    filter(reps >= 50) |>   # lower threshold: escalation cells have fewer reps
    mutate(dgp = factor(dgp, levels = DGP_ORDER)) |>
    arrange(method, dgp, n, c_bin)

  print(esc_tab, row.names = FALSE)
  cat("\n")
}

# =============================================================================
# VERDICT SUMMARY
# =============================================================================
cat("=" |> strrep(72), "\n")
cat("VERDICT SUMMARY\n")
cat("=" |> strrep(72), "\n\n")

# Pull key numbers for the summary
cs_focus <- cs |>
  filter(method %in% c("crossfit", "single_tree", "dt_averaged")) |>
  select(method, n, dgp, reps, cover, cover_se, rmse, bias, int_ne)

# Q1: Does single_tree reach nominal?
st_cells <- cs_focus |> filter(method == "single_tree")
st_nominal <- st_cells |>
  mutate(nominal = abs(cover - 0.95) <= 1.96 * cover_se)

cat("Q1. single_tree covers nominally in",
    sum(st_nominal$nominal, na.rm = TRUE), "/",
    nrow(st_nominal), "adequately-sampled cells.\n")

worst <- st_cells |> slice_min(cover, n = 3)
cat("   Worst coverage cells:\n")
for (i in seq_len(nrow(worst))) {
  row <- worst[i, ]
  cat(sprintf("     n=%d dgp=%-12s cover=%.3f (reps=%d)\n",
              row$n, as.character(row$dgp), row$cover, row$reps))
}

# Q2: Does honest CI restore coverage vs Wald?
dtavg_cells <- cs_focus |> filter(method == "dt_averaged")
dtavg_nominal <- dtavg_cells |>
  mutate(nominal = abs(cover - 0.95) <= 1.96 * cover_se)
cat("\nQ2. dt_averaged (honest CI) covers nominally in",
    sum(dtavg_nominal$nominal, na.rm = TRUE), "/",
    nrow(dtavg_nominal), "adequately-sampled cells.\n")

# Q4: typical RMSE cost
if ("crossfit_rmse" %in% names(rmse_tab)) {
  st_cost <- rmse_tab |>
    filter(!is.na(rmse_ratio_single)) |>
    summarise(
      med_cost  = median(rmse_ratio_single, na.rm = TRUE),
      max_cost  = max(rmse_ratio_single,    na.rm = TRUE)
    )
  cat(sprintf("\nQ4. single_tree RMSE ratio vs crossfit: median=%.3f, max=%.3f\n",
              st_cost$med_cost, st_cost$max_cost))
}

cat("\nDone.\n")
