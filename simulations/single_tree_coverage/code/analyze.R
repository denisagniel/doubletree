# ============================================================
# analyze.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
# Purpose: Aggregate the per-rep results into per-(DGP, n) summaries -- structure
#   recovery, theta*_n coverage, theta_0 coverage, the fidelity diagnostic delta,
#   all with Monte-Carlo standard errors -- and render the two headline figures:
#   (1) recovery_rate vs n by DGP; (2) coverage (theta*_n vs theta_0) vs n by DGP.
# Inputs : results/results_<runid>.rds (latest, or path via arg) + theta_star_oracle.rds.
# Outputs: results/summary_<runid>.rds + .csv; figures/recovery_rate_vs_n.{pdf,png},
#          figures/coverage_star_vs_theta0.{pdf,png}.
# Run    : from repo root ->
#   Rscript doubletree/simulations/single_tree_coverage/code/analyze.R [results_file.rds]
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(fs)
})

STUDY_DIR <- fs::path("doubletree", "simulations", "single_tree_coverage")
source(fs::path(STUDY_DIR, "code", "oracle_theta_star.R"))   # THETA0, load_oracle

# RAND house style if available, else a clean fallback (r-code-conventions S4).
use_rand <- requireNamespace("randplot", quietly = TRUE)
if (use_rand) suppressPackageStartupMessages(library(randplot))
study_theme <- if (use_rand) randplot::theme_rand() else ggplot2::theme_minimal(base_size = 12)
dgp_palette <- if (use_rand) randplot::RandCatPal else
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

NOMINAL <- 0.95   # target coverage / recovery reference line

# --- locate the results file -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
results_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  cand <- fs::dir_ls(fs::path(STUDY_DIR, "results"), regexp = "results_\\d+.*\\.rds$")
  if (length(cand) == 0) {
    stop("No results_<runid>.rds found in ", fs::path(STUDY_DIR, "results"),
         ". Run run_coverage.R (or run_pilot.R) first.", call. = FALSE)
  }
  cand[which.max(file.info(cand)$mtime)]
}
cli::cli_inform("Analyzing {.path {results_path}}")

results    <- readr::read_rds(results_path)
oracle_tbl <- load_oracle(STUDY_DIR)

# Coverage/recovery are defined only over cells that actually RAN (status == "run").
run_rows <- dplyr::filter(results, status == "run")
if (nrow(run_rows) == 0) stop("No run rows (status == 'run') in results.", call. = FALSE)

# --- per-(DGP, n) summary with MC standard errors ----------------------------
# Recovery MC se uses R (all reps); coverage MC se uses n_converged (coverage is
# conditional on convergence). se = sqrt(p(1-p)/denominator).
mc_se <- function(p, denom) sqrt(p * (1 - p) / denom)

summary_tbl <- run_rows |>
  dplyr::summarise(
    reps          = dplyr::n(),
    n_converged   = sum(converged),
    recovery_rate = mean(converged),
    cov_star_honest   = mean(covers_star_honest[converged]),
    cov_star_single   = mean(covers_star_single[converged]),
    cov_star_cf       = mean(covers_star_cf[converged]),
    cov_theta0_honest = mean(covers_theta0_honest[converged]),
    cov_theta0_single = mean(covers_theta0_single[converged]),
    cov_theta0_cf     = mean(covers_theta0_cf[converged]),
    mean_delta         = mean(delta[converged]),
    sd_delta           = sd(delta[converged]),
    mean_delta_over_se = mean(delta_over_se[converged], na.rm = TRUE),
    mean_elapsed_s     = mean(elapsed_s),
    .by = c(dgp, n)
  ) |>
  dplyr::mutate(
    recovery_se       = mc_se(recovery_rate, reps),
    cov_star_honest_se   = mc_se(cov_star_honest, n_converged),
    cov_star_single_se   = mc_se(cov_star_single, n_converged),
    cov_star_cf_se       = mc_se(cov_star_cf, n_converged),
    cov_theta0_honest_se = mc_se(cov_theta0_honest, n_converged),
    cov_theta0_single_se = mc_se(cov_theta0_single, n_converged),
    cov_theta0_cf_se     = mc_se(cov_theta0_cf, n_converged)
  ) |>
  dplyr::left_join(
    dplyr::select(oracle_tbl, dgp, n, theta_star, theta0, theta0_realized, gap_vs_realized),
    by = c("dgp", "n")
  ) |>
  dplyr::arrange(dgp, n)

# --- save summaries ----------------------------------------------------------
runid <- sub(".*results_(\\d+_\\d+)\\.rds$", "\\1", fs::path_file(results_path))
if (identical(runid, fs::path_file(results_path))) runid <- format(Sys.time(), "%Y%m%d_%H%M%S")
readr::write_rds(summary_tbl, fs::path(STUDY_DIR, "results", paste0("summary_", runid, ".rds")))
readr::write_csv(summary_tbl, fs::path(STUDY_DIR, "results", paste0("summary_", runid, ".csv")))

cat("\n=== Per-(DGP, n) summary (feasible/run cells) ===\n")
summary_tbl |>
  dplyr::select(dgp, n, reps, recovery_rate, recovery_se,
                cov_star_honest, cov_star_honest_se,
                cov_theta0_honest, cov_theta0_honest_se, mean_delta) |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 4))) |>
  as.data.frame() |>
  print(row.names = FALSE)

# --- Figure 1: structure recovery vs n by DGP --------------------------------
fig_recovery <- summary_tbl |>
  ggplot2::ggplot(ggplot2::aes(n, recovery_rate, colour = dgp, group = dgp)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = recovery_rate - 1.96 * recovery_se,
                 ymax = pmin(1, recovery_rate + 1.96 * recovery_se)),
    width = 0.03
  ) +
  ggplot2::scale_x_log10(breaks = N_GRID) +
  ggplot2::scale_colour_manual(values = dgp_palette, name = "DGP") +
  ggplot2::labs(
    title = "Structure recovery (Rashomon intersection non-empty) vs n",
    subtitle = "Empirical proxy for the partition-fixing event E_n (Theorem 1 hypothesis)",
    x = "Sample size n (log scale)", y = "Recovery rate  (mean of `converged`)"
  ) +
  study_theme

# --- Figure 2: coverage of theta*_n vs theta_0 (honest CI) by DGP ------------
# Long frame built explicitly (one row per target x cell) for a two-line facet.
cov_long <- dplyr::bind_rows(
  summary_tbl |> dplyr::transmute(dgp, n, target = "theta*_n (Thm 1)",
                                  coverage = cov_star_honest, se = cov_star_honest_se),
  summary_tbl |> dplyr::transmute(dgp, n, target = "theta_0 (true ATT)",
                                  coverage = cov_theta0_honest, se = cov_theta0_honest_se)
)

fig_coverage <- cov_long |>
  ggplot2::ggplot(ggplot2::aes(n, coverage, colour = target, group = target)) +
  ggplot2::geom_hline(yintercept = NOMINAL, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = coverage - 1.96 * se, ymax = coverage + 1.96 * se),
    width = 0.03
  ) +
  ggplot2::facet_wrap(~ dgp) +
  ggplot2::scale_x_log10(breaks = N_GRID) +
  ggplot2::scale_colour_manual(values = dgp_palette, name = "Coverage target") +
  ggplot2::labs(
    title = "Honest-CI coverage: theta*_n (working model) vs theta_0 (true ATT)",
    subtitle = "Dashed line = 0.95 nominal; error bars +/- 1.96 MC se",
    x = "Sample size n (log scale)", y = "Coverage of honest 95% CI"
  ) +
  study_theme

# --- export figures (vector PDF for paper + PNG preview) ---------------------
fig_dir <- fs::path(STUDY_DIR, "figures")
fs::dir_create(fig_dir)
ggplot2::ggsave(fs::path(fig_dir, "recovery_rate_vs_n.pdf"), fig_recovery,
                width = 8, height = 5, device = "pdf")
ggplot2::ggsave(fs::path(fig_dir, "recovery_rate_vs_n.png"), fig_recovery,
                width = 8, height = 5, dpi = 150)
ggplot2::ggsave(fs::path(fig_dir, "coverage_star_vs_theta0.pdf"), fig_coverage,
                width = 9, height = 6, device = "pdf")
ggplot2::ggsave(fs::path(fig_dir, "coverage_star_vs_theta0.png"), fig_coverage,
                width = 9, height = 6, dpi = 150)

cli::cli_inform(c(
  "v" = "Wrote summary_{runid}.rds/.csv and 2 figures (pdf + png) to {.path {fig_dir}}.",
  "i" = "Theme: {if (use_rand) 'randplot::theme_rand()' else 'theme_minimal (randplot unavailable)'}."
))
