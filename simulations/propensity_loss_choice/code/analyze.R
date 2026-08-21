# ============================================================
# Propensity-tree loss ablation -- ANALYSIS
#
# Stem:    2026-08-20_propensity-loss-choice
# Purpose: Turn results_<runid>.rds into the spec's primary metrics (bias, RMSE,
#          95% coverage, mean CI width per arm), the diagnostic table (failure
#          rate, same-partition rate), the max(e_hat) distributions, and the
#          paired arm contrasts that the decision rule is stated in terms of.
# Inputs:  results/results_<runid>.rds (latest, or PLC_RESULTS=<path>)
# Outputs: tables/*.csv, tables/summary_<runid>.md,
#          figures/bias_rmse_coverage_<runid>.png
#
# Run from the doubletree package root:
#   Rscript simulations/propensity_loss_choice/code/analyze.R
# ============================================================

# ---- 0 setup ----------------------------------------------------------------

source(file.path("simulations", "propensity_loss_choice", "code", "common.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

results_path <- Sys.getenv("PLC_RESULTS", "")
if (!nzchar(results_path)) {
  candidates <- sort(list.files(DIR_RESULTS, "^results_.*\\.rds$", full.names = TRUE))
  if (length(candidates) == 0L) {
    cli::cli_abort("No {.file results_*.rds} in {.path {DIR_RESULTS}}; run {.file run_sim.R} first.")
  }
  results_path <- candidates[[length(candidates)]]
}
cli::cli_inform("gates: PLC_RESULTS={results_path}")

bundle  <- readRDS(results_path)
meta    <- bundle$meta
results <- bundle$results
runid   <- meta$runid

if (!identical(as.integer(unique(results$rep)) |> length(), as.integer(meta$n_reps))) {
  cli::cli_abort("Results hold {dplyr::n_distinct(results$rep)} distinct reps, metadata claims {meta$n_reps}.")
}

# A wide, one-row-per-replicate frame: the paired unit of analysis.
paired <- results |>
  dplyr::mutate(err = theta - true_att) |>
  dplyr::select(dgp, n, rep, arm, err, covered, ci_width, n_leaves_e, max_e_hat,
                hash_e, error_message) |>
  tidyr::pivot_wider(
    names_from  = arm,
    values_from = c(err, covered, ci_width, n_leaves_e, max_e_hat, hash_e, error_message)
  )

# ---- 1 primary metrics per arm ---------------------------------------------

primary <- results |>
  dplyr::summarise(
    reps        = dplyr::n(),
    bias        = mean(theta - true_att),
    bias_mcse   = stats::sd(theta - true_att) / sqrt(dplyr::n()),
    rmse        = sqrt(mean((theta - true_att)^2)),
    sd_theta    = stats::sd(theta),
    coverage    = mean(covered),
    cov_mcse    = sqrt(mean(covered) * (1 - mean(covered)) / dplyr::n()),
    # width_mcse BEFORE ci_width: dplyr::summarise() evaluates sequentially, so
    # naming the mean `ci_width` first would shadow the column and yield sd(scalar) = NA.
    width_mcse  = stats::sd(ci_width) / sqrt(dplyr::n()),
    ci_width    = mean(ci_width),
    mean_sigma  = mean(sigma),
    .by = c(dgp, n, arm)
  ) |>
  dplyr::relocate(ci_width, .before = width_mcse) |>
  dplyr::arrange(dgp, n, arm)

# ---- 2 paired arm contrasts (S - L) ----------------------------------------

# The decision rule asks whether S "matches L within MC noise". Because the
# design is paired, the right yardstick is the MC se of the WITHIN-REPLICATE
# difference, which is far tighter than the two arms' separate MC ses.
contrasts <- paired |>
  dplyr::summarise(
    reps            = dplyr::n(),
    d_bias          = mean(err_S - err_L),
    d_bias_mcse     = stats::sd(err_S - err_L) / sqrt(dplyr::n()),
    d_sqerr         = mean(err_S^2 - err_L^2),
    d_sqerr_mcse    = stats::sd(err_S^2 - err_L^2) / sqrt(dplyr::n()),
    rmse_ratio_S_L  = sqrt(mean(err_S^2)) / sqrt(mean(err_L^2)),
    d_coverage      = mean(covered_S) - mean(covered_L),
    d_cov_mcse      = stats::sd(as.numeric(covered_S) - as.numeric(covered_L)) / sqrt(dplyr::n()),
    d_width         = mean(ci_width_S - ci_width_L),
    d_width_mcse    = stats::sd(ci_width_S - ci_width_L) / sqrt(dplyr::n()),
    .by = c(dgp, n)
  ) |>
  dplyr::mutate(
    bias_z     = d_bias / d_bias_mcse,
    sqerr_z    = d_sqerr / d_sqerr_mcse,
    coverage_z = d_coverage / d_cov_mcse
  ) |>
  dplyr::arrange(dgp, n)

# ---- 3 diagnostics: failures, partition agreement, fitted propensity -------

diagnostics <- results |>
  dplyr::summarise(
    arm_fits     = dplyr::n(),
    failure_rate = mean(!is.na(error_message)),
    .by = c(dgp, n)
  ) |>
  dplyr::left_join(
    paired |>
      dplyr::summarise(
        same_partition_rate_e = mean(hash_e_L == hash_e_S),
        same_leafcount_rate_e = mean(n_leaves_e_L == n_leaves_e_S),
        mean_leaves_e_L       = mean(n_leaves_e_L),
        mean_leaves_e_S       = mean(n_leaves_e_S),
        .by = c(dgp, n)
      ),
    by = c("dgp", "n")
  ) |>
  dplyr::arrange(dgp, n)

e_hat_dist <- results |>
  dplyr::summarise(
    mean_max_e_hat = mean(max_e_hat),
    q50            = stats::quantile(max_e_hat, 0.50),
    q90            = stats::quantile(max_e_hat, 0.90),
    max_max_e_hat  = max(max_e_hat),
    pct_at_clip    = mean(max_e_hat >= 0.99 - 1e-12),
    max_true_e     = max(max_true_e),
    .by = c(dgp, n, arm)
  ) |>
  dplyr::arrange(dgp, n, arm)

leaf_dist <- results |>
  dplyr::summarise(count = dplyr::n(), .by = c(dgp, n, arm, n_leaves_e)) |>
  tidyr::pivot_wider(names_from = n_leaves_e, values_from = count,
                     names_prefix = "leaves_", values_fill = 0L) |>
  dplyr::arrange(dgp, n, arm)

# ---- 4 figure ---------------------------------------------------------------

fig_data <- primary |>
  dplyr::select(dgp, n, arm, bias, bias_mcse, rmse, coverage, cov_mcse,
                ci_width, width_mcse) |>
  dplyr::mutate(rmse_mcse = NA_real_) |>
  tidyr::pivot_longer(
    cols = c(bias, rmse, coverage, ci_width),
    names_to = "metric", values_to = "value"
  ) |>
  dplyr::mutate(
    mcse = dplyr::case_when(
      metric == "bias"     ~ bias_mcse,
      metric == "coverage" ~ cov_mcse,
      metric == "ci_width" ~ width_mcse,
      TRUE                 ~ NA_real_
    ),
    metric = factor(metric, levels = c("bias", "rmse", "coverage", "ci_width"),
                    labels = c("ATT bias", "ATT RMSE", "95% CI coverage",
                               "mean CI width")),
    dgp_label = factor(DGP_LABELS[dgp], levels = unname(DGP_LABELS)),
    arm_label = factor(ARM_LABELS[arm], levels = unname(ARM_LABELS)),
    n_label   = factor(n, levels = N_GRID, labels = paste0("n = ", N_GRID))
  )

reference_lines <- tibble::tibble(
  metric = factor(c("ATT bias", "95% CI coverage"),
                  levels = levels(fig_data$metric)),
  yintercept = c(0, 0.95)
)

fig <- ggplot(fig_data, aes(x = n_label, y = value, colour = arm_label,
                            shape = arm_label)) +
  geom_hline(data = reference_lines, aes(yintercept = yintercept),
             linetype = "dashed", colour = "grey45", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = value - 1.96 * mcse, ymax = value + 1.96 * mcse),
                width = 0.12, position = position_dodge(width = 0.45),
                na.rm = TRUE) +
  geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
  facet_grid(metric ~ dgp_label, scales = "free_y", switch = "y") +
  scale_colour_manual(values = c("#1B6CA8", "#C4471C")) +
  labs(
    title    = "Propensity-tree loss ablation: log-loss vs squared-error",
    subtitle = paste0(
      "Paired design, identical data per replicate; leaf_budget = ", meta$leaf_budget,
      ", m_n = ", meta$m_n, ", lambda_n = log(n)/n, R = ", meta$n_reps,
      " reps/cell.\nBars are +/- 1.96 MC se (RMSE has none)."
    ),
    x = NULL, y = NULL, colour = NULL, shape = NULL,
    caption = paste0("runid ", runid, " | doubletree ", meta$doubletree_sha,
                     " | optimaltrees ", meta$optimaltrees_sha)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    panel.grid.minor = element_blank()
  )

fig_path <- file.path(DIR_FIGURES, paste0("bias_rmse_coverage_", runid, ".png"))
ggsave(fig_path, fig, width = 9, height = 8, dpi = 200, bg = "white")

# ---- 5 export ---------------------------------------------------------------

readr::write_csv(primary,     file.path(DIR_TABLES, paste0("primary_",     runid, ".csv")))
readr::write_csv(contrasts,   file.path(DIR_TABLES, paste0("contrasts_",   runid, ".csv")))
readr::write_csv(diagnostics, file.path(DIR_TABLES, paste0("diagnostics_", runid, ".csv")))
readr::write_csv(e_hat_dist,  file.path(DIR_TABLES, paste0("e_hat_",       runid, ".csv")))
readr::write_csv(leaf_dist,   file.path(DIR_TABLES, paste0("leaves_e_",    runid, ".csv")))

md_path <- file.path(DIR_TABLES, paste0("summary_", runid, ".md"))
writeLines(c(
  paste0("# Propensity-loss ablation -- run ", runid),
  "",
  paste0("- spec: `", meta$spec, "`"),
  paste0("- leaf_budget = ", meta$leaf_budget, ", m_n = ", meta$m_n,
         ", lambda_n = ", meta$lambda_n, ", outcome_type = ", meta$outcome_type),
  paste0("- R = ", meta$n_reps, " reps per (dgp, n) cell; arms = ",
         paste(names(meta$arms), unname(meta$arms), sep = "=", collapse = ", ")),
  paste0("- seed_master = ", meta$seed_master,
         "; per-rep seed = f(seed_master, dgp, n, rep)"),
  paste0("- doubletree ", meta$doubletree_sha, ", optimaltrees ",
         meta$optimaltrees_sha, ", ", meta$r_version),
  "",
  "## Primary metrics per arm",
  "",
  knitr::kable(primary, digits = 4, format = "pipe"),
  "",
  "## Paired contrasts (Arm S minus Arm L, within replicate)",
  "",
  knitr::kable(contrasts, digits = 4, format = "pipe"),
  "",
  "## Diagnostics: failure rate and partition agreement",
  "",
  knitr::kable(diagnostics, digits = 4, format = "pipe"),
  "",
  "## max(e_hat) distribution per arm",
  "",
  knitr::kable(e_hat_dist, digits = 4, format = "pipe"),
  "",
  "## Selected propensity-tree leaf counts",
  "",
  knitr::kable(leaf_dist, format = "pipe")
), md_path)

# ---- 6 console report -------------------------------------------------------

cli::cli_h1("Primary metrics per arm")
print(as.data.frame(primary), digits = 4)
cli::cli_h1("Paired contrasts (S - L)")
print(as.data.frame(contrasts), digits = 3)
cli::cli_h1("Diagnostics (failure rate, same-partition rate)")
print(as.data.frame(diagnostics), digits = 4)
cli::cli_h1("max(e_hat)")
print(as.data.frame(e_hat_dist), digits = 4)
cli::cli_h1("Propensity-tree leaf counts")
print(as.data.frame(leaf_dist))

# Evidence for the spec's decision rule on the STRESS DGP. Reported, not acted
# on: changing estimate_att()'s default is a paper-facing call made by a human.
stress <- contrasts |> dplyr::filter(dgp == "weak_overlap")
cli::cli_h1("Decision-rule evidence (stress DGP only)")
for (i in seq_len(nrow(stress))) {
  s <- stress[i, ]
  cli::cli_inform(c(
    "*" = "n={s$n}: d_bias={round(s$d_bias, 4)} (z={round(s$bias_z, 2)}), d_sqerr={signif(s$d_sqerr, 3)} (z={round(s$sqerr_z, 2)}), RMSE ratio S/L={round(s$rmse_ratio_S_L, 3)}, d_coverage={round(s$d_coverage, 4)} (z={round(s$coverage_z, 2)})"
  ))
}
cli::cli_inform(c(
  i = "|z| < 2 on all three => Arm S matches Arm L within MC noise at that n.",
  i = "This script does NOT modify estimate_att()'s default propensity_loss."
))

cli::cli_inform(c("v" = "wrote {.file {fig_path}}", "v" = "wrote {.file {md_path}}"))
