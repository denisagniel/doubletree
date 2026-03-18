#' Analysis Script for β < d/2 Smoothness Regime Simulations
#'
#' Generates figures and tables characterizing performance across β regimes.
#' Outputs:
#' - Table: Coverage/RMSE/CI width by regime × method × n
#' - Figure 1: Coverage vs sample size (faceted by β regime)
#' - Figure 2: RMSE vs sample size with theoretical rates
#' - Figure 3: CI width by regime
#'
#' Three-way fidelity: Results feed manuscript Section 4 or Appendix C

library(dplyr)
library(ggplot2)
library(tidyr)

# Find most recent results directory
results_dirs <- list.dirs("results", full.names = TRUE, recursive = FALSE)
beta_dirs <- results_dirs[grepl("beta_study", results_dirs)]

if (length(beta_dirs) == 0) {
  stop("No beta_study results found. Run run_beta_study.R first.")
}

# Use most recent
results_dir <- beta_dirs[which.max(file.info(beta_dirs)$mtime)]
cat(sprintf("Analyzing results from: %s\n\n", results_dir))

# Load simulation results
results <- readRDS(file.path(results_dir, "simulation_results.rds"))

cat(sprintf("Loaded %d simulation results\n", nrow(results)))
cat(sprintf("Convergence rate: %.1f%%\n\n", 100 * mean(results$converged, na.rm = TRUE)))

# Filter to converged results for main analysis
results_converged <- results %>%
  filter(converged)

cat(sprintf("Analyzing %d converged results\n\n", nrow(results_converged)))

# Create figures directory
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

#------------------------------------------------------------------------------
# Figure 0: Tree Complexity Verification (s_n: Fitted vs Theoretical)
#------------------------------------------------------------------------------

cat("Generating Figure 0: Tree complexity verification (s_n)...\n")

# Extract tree complexity data (tree method only)
sn_data <- results_converged %>%
  filter(method == "tree") %>%
  group_by(dgp, beta, rate_regime, n) %>%
  summarize(
    mean_n_leaves_e = mean(n_leaves_e, na.rm = TRUE),
    sd_n_leaves_e = sd(n_leaves_e, na.rm = TRUE),
    mean_n_leaves_m0 = mean(n_leaves_m0, na.rm = TRUE),
    sd_n_leaves_m0 = sd(n_leaves_m0, na.rm = TRUE),
    theoretical_sn = first(theoretical_sn),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_e = sd_n_leaves_e / sqrt(n_reps),
    se_m0 = sd_n_leaves_m0 / sqrt(n_reps),
    beta_label = factor(
      sprintf("β = %d (%s)", beta, rate_regime),
      levels = c("β = 3 (valid)", "β = 2 (boundary)", "β = 1 (invalid)")
    )
  )

# Reshape for plotting (propensity and outcome separate)
sn_long <- sn_data %>%
  pivot_longer(
    cols = c(mean_n_leaves_e, mean_n_leaves_m0),
    names_to = "nuisance",
    values_to = "fitted_sn"
  ) %>%
  mutate(
    se = ifelse(nuisance == "mean_n_leaves_e", se_e, se_m0),
    nuisance_label = ifelse(nuisance == "mean_n_leaves_e",
                            "Propensity e(X)", "Outcome m0(X)")
  )

fig0 <- ggplot(sn_long, aes(x = n, y = fitted_sn, color = nuisance_label, shape = nuisance_label)) +
  geom_line(aes(y = theoretical_sn), color = "gray30", linetype = "dashed",
            size = 1.2, inherit.aes = FALSE) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = fitted_sn - 1.96*se, ymax = fitted_sn + 1.96*se),
                width = 50, alpha = 0.5) +
  facet_wrap(~ beta_label, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(400, 800, 1600)) +
  scale_y_continuous() +
  scale_color_manual(
    values = c("Propensity e(X)" = "#E41A1C", "Outcome m0(X)" = "#377EB8")
  ) +
  scale_shape_manual(
    values = c("Propensity e(X)" = 16, "Outcome m0(X)" = 17)
  ) +
  labs(
    x = "Sample Size",
    y = "Number of Leaves (Tree Complexity)",
    color = "Nuisance Function",
    shape = "Nuisance Function",
    title = "Tree Complexity: Fitted vs Theoretical s_n = n^(d/(2β+d))",
    subtitle = "Dashed line: theoretical prediction; Points: fitted trees (mean ± 95% CI)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig0_tree_complexity_verification.pdf"),
       fig0, width = 10, height = 4.5)
ggsave(file.path(fig_dir, "fig0_tree_complexity_verification.png"),
       fig0, width = 10, height = 4.5, dpi = 300)

cat(sprintf("✓ Figure 0 saved to: %s\n\n", fig_dir))

# Print s_n verification table
cat("\nTree Complexity Verification (n=800):\n")
cat(strrep("=", 70), "\n")
sn_verification <- sn_data %>%
  filter(n == 800) %>%
  select(beta, rate_regime, theoretical_sn, mean_n_leaves_e, mean_n_leaves_m0) %>%
  mutate(
    theoretical_sn = round(theoretical_sn, 1),
    mean_n_leaves_e = round(mean_n_leaves_e, 1),
    mean_n_leaves_m0 = round(mean_n_leaves_m0, 1),
    ratio_e = round(mean_n_leaves_e / theoretical_sn, 2),
    ratio_m0 = round(mean_n_leaves_m0 / theoretical_sn, 2)
  )
print(sn_verification)
cat("\nRatio interpretation: 1.0 = perfect match, <0.5 = over-regularized, >2.0 = under-regularized\n\n")

#------------------------------------------------------------------------------
# Table 1: Summary Statistics by β Regime
#------------------------------------------------------------------------------

cat("Generating Table 1: Summary statistics by β regime...\n")

summary_stats <- results_converged %>%
  group_by(dgp, beta, rate_regime, method, n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att, na.rm = TRUE),
    rmse = sqrt(mean((theta - true_att)^2, na.rm = TRUE)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    median_ci_width = median(ci_upper - ci_lower, na.rm = TRUE),
    mean_n_leaves_e = mean(n_leaves_e, na.rm = TRUE),     # NEW: tree complexity
    mean_n_leaves_m0 = mean(n_leaves_m0, na.rm = TRUE),   # NEW: tree complexity
    theoretical_sn = first(theoretical_sn),               # NEW: theoretical prediction
    .groups = "drop"
  )

# Format table for display
table1 <- summary_stats %>%
  mutate(
    beta_label = sprintf("β=%d (%s)", beta, rate_regime),
    coverage_pct = sprintf("%.1f%%", 100 * coverage),
    rmse_fmt = sprintf("%.4f", rmse),
    ci_width_fmt = sprintf("%.4f", mean_ci_width)
  ) %>%
  select(beta_label, method, n, coverage_pct, rmse_fmt, ci_width_fmt)

cat("\nTable 1: Performance by β Regime\n")
cat(strrep("=", 70), "\n")
print(table1, n = Inf)

# Save table
write.csv(table1, file.path(results_dir, "table1_performance_by_beta.csv"),
          row.names = FALSE)
cat(sprintf("\n✓ Table saved to: %s\n\n",
            file.path(results_dir, "table1_performance_by_beta.csv")))

#------------------------------------------------------------------------------
# Figure 1: Coverage vs Sample Size (by β regime)
#------------------------------------------------------------------------------

cat("Generating Figure 1: Coverage vs sample size by β regime...\n")

# Compute coverage with 95% CIs
coverage_data <- results_converged %>%
  group_by(dgp, beta, rate_regime, method, n) %>%
  summarize(
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att),
    n_reps = n(),
    se = sqrt(coverage * (1 - coverage) / n_reps),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = pmax(0, coverage - 1.96 * se),
    ci_upper = pmin(1, coverage + 1.96 * se),
    beta_label = factor(
      sprintf("β = %d (%s)", beta, rate_regime),
      levels = c("β = 3 (valid)", "β = 2 (boundary)", "β = 1 (invalid)")
    )
  )

fig1 <- ggplot(coverage_data, aes(x = n, y = coverage, color = method, shape = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40", size = 0.5) +
  geom_hline(yintercept = 0.90, linetype = "dotted", color = "gray60", size = 0.3) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 50, alpha = 0.5) +
  facet_wrap(~ beta_label, nrow = 1) +
  scale_x_continuous(breaks = c(400, 800, 1600)) +
  scale_y_continuous(limits = c(0.85, 1.00),
                     breaks = seq(0.85, 1.00, 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c("tree" = "#E41A1C", "forest" = "#377EB8", "linear" = "#4DAF4A"),
    labels = c("tree" = "Tree-DML", "forest" = "Forest-DML", "linear" = "Linear-DML")
  ) +
  scale_shape_manual(
    values = c("tree" = 16, "forest" = 17, "linear" = 15),
    labels = c("tree" = "Tree-DML", "forest" = "Forest-DML", "linear" = "Linear-DML")
  ) +
  labs(
    x = "Sample Size",
    y = "Coverage",
    color = "Method",
    shape = "Method",
    title = "Coverage by β Smoothness Regime",
    subtitle = "Target: 95% (dashed line); Acceptable: 90-98%"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig1_coverage_by_beta.pdf"),
       fig1, width = 10, height = 4)
ggsave(file.path(fig_dir, "fig1_coverage_by_beta.png"),
       fig1, width = 10, height = 4, dpi = 300)

cat(sprintf("✓ Figure 1 saved to: %s\n\n", fig_dir))

#------------------------------------------------------------------------------
# Figure 2: RMSE vs Sample Size with Theoretical Rates
#------------------------------------------------------------------------------

cat("Generating Figure 2: RMSE vs sample size with theoretical rates...\n")

rmse_data <- results_converged %>%
  group_by(dgp, beta, rate_regime, method, n) %>%
  summarize(
    rmse = sqrt(mean((theta - true_att)^2)),
    .groups = "drop"
  ) %>%
  mutate(
    beta_label = factor(
      sprintf("β = %d (%s)", beta, rate_regime),
      levels = c("β = 3 (valid)", "β = 2 (boundary)", "β = 1 (invalid)")
    )
  )

# Add theoretical rate lines (for tree method only, as reference)
# β=3: n^(-0.30), β=2: n^(-0.25), β=1: n^(-0.167)
n_seq <- seq(400, 1600, by = 10)
theoretical_rates <- expand.grid(
  n = n_seq,
  beta = c(1, 2, 3)
) %>%
  mutate(
    d = 4,
    rate_exponent = -beta / (2 * beta + d),
    # Scale to match approximate RMSE magnitude (calibrate constant)
    rmse_theoretical = 0.08 * (n / 800) ^ rate_exponent,
    beta_label = factor(
      sprintf("β = %d (%s)", beta,
              ifelse(beta == 3, "valid",
                     ifelse(beta == 2, "boundary", "invalid"))),
      levels = c("β = 3 (valid)", "β = 2 (boundary)", "β = 1 (invalid)")
    )
  )

fig2 <- ggplot(rmse_data, aes(x = n, y = rmse, color = method, shape = method)) +
  geom_line(data = theoretical_rates,
            aes(x = n, y = rmse_theoretical),
            color = "gray40", linetype = "dashed", size = 0.5,
            inherit.aes = FALSE) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ beta_label, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(400, 800, 1600)) +
  scale_y_log10(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(
    values = c("tree" = "#E41A1C", "forest" = "#377EB8", "linear" = "#4DAF4A"),
    labels = c("tree" = "Tree-DML", "forest" = "Forest-DML", "linear" = "Linear-DML")
  ) +
  scale_shape_manual(
    values = c("tree" = 16, "forest" = 17, "linear" = 15),
    labels = c("tree" = "Tree-DML", "forest" = "Forest-DML", "linear" = "Linear-DML")
  ) +
  labs(
    x = "Sample Size",
    y = "RMSE (log scale)",
    color = "Method",
    shape = "Method",
    title = "RMSE by β Smoothness Regime",
    subtitle = "Dashed lines: theoretical rates n^(-β/(2β+d))"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig2_rmse_by_beta.pdf"),
       fig2, width = 10, height = 4)
ggsave(file.path(fig_dir, "fig2_rmse_by_beta.png"),
       fig2, width = 10, height = 4, dpi = 300)

cat(sprintf("✓ Figure 2 saved to: %s\n\n", fig_dir))

#------------------------------------------------------------------------------
# Figure 3: CI Width by β Regime
#------------------------------------------------------------------------------

cat("Generating Figure 3: CI width distribution by β regime...\n")

ci_width_data <- results_converged %>%
  mutate(
    ci_width = ci_upper - ci_lower,
    beta_label = factor(
      sprintf("β = %d\n(%s)", beta, rate_regime),
      levels = c("β = 3\n(valid)", "β = 2\n(boundary)", "β = 1\n(invalid)")
    )
  ) %>%
  filter(method == "tree", n == 800)  # Focus on tree-DML at n=800

fig3 <- ggplot(ci_width_data, aes(x = beta_label, y = ci_width, fill = beta_label)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
  scale_fill_manual(
    values = c("β = 3\n(valid)" = "#4DAF4A",
               "β = 2\n(boundary)" = "#FFA500",
               "β = 1\n(invalid)" = "#E41A1C")
  ) +
  labs(
    x = "β Smoothness Regime",
    y = "Confidence Interval Width",
    title = "CI Width Distribution by β Regime (Tree-DML, n=800)",
    subtitle = "Wider CIs expected when β < d/2 (if method detects slower convergence)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig3_ci_width_by_beta.pdf"),
       fig3, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig3_ci_width_by_beta.png"),
       fig3, width = 7, height = 5, dpi = 300)

cat(sprintf("✓ Figure 3 saved to: %s\n\n", fig_dir))

#------------------------------------------------------------------------------
# Statistical Tests
#------------------------------------------------------------------------------

cat(strrep("=", 70), "\n")
cat("Statistical Tests\n")
cat(strrep("=", 70), "\n\n")

# Test 1: Is β=1 coverage significantly different from 95% (tree-DML, n=800)?
test1_data <- results_converged %>%
  filter(dgp == "beta_low", method == "tree", n == 800)

if (nrow(test1_data) > 0) {
  coverage_beta1 <- mean(test1_data$ci_lower <= test1_data$true_att &
                         test1_data$ci_upper >= test1_data$true_att)
  n_beta1 <- nrow(test1_data)

  z_stat <- (coverage_beta1 - 0.95) / sqrt(0.95 * 0.05 / n_beta1)
  p_value <- 2 * pnorm(-abs(z_stat))

  cat("Test 1: β=1 (low smoothness) vs 95% target\n")
  cat(sprintf("  Sample: Tree-DML at n=800 (%d replications)\n", n_beta1))
  cat(sprintf("  Observed coverage: %.3f (%.1f%%)\n", coverage_beta1, 100 * coverage_beta1))
  cat(sprintf("  Z-statistic: %.3f\n", z_stat))
  cat(sprintf("  P-value: %.5f\n", p_value))
  if (p_value < 0.05) {
    cat(sprintf("  ✓ Significantly different from 95%% (p < 0.05)\n"))
    cat(sprintf("    → β=1 < d/2 condition manifests empirically\n"))
  } else {
    cat(sprintf("  ✗ Not significantly different (p ≥ 0.05)\n"))
    cat(sprintf("    → β=1 < d/2 does not manifest at n=800\n"))
  }
  cat("\n")
}

# Test 2: β=3 vs β=1 coverage difference (tree-DML, n=800)
test2_data_high <- results_converged %>%
  filter(dgp == "beta_high", method == "tree", n == 800)
test2_data_low <- results_converged %>%
  filter(dgp == "beta_low", method == "tree", n == 800)

if (nrow(test2_data_high) > 0 && nrow(test2_data_low) > 0) {
  coverage_high <- mean(test2_data_high$ci_lower <= test2_data_high$true_att &
                        test2_data_high$ci_upper >= test2_data_high$true_att)
  coverage_low <- mean(test2_data_low$ci_lower <= test2_data_low$true_att &
                       test2_data_low$ci_upper >= test2_data_low$true_att)

  n_high <- nrow(test2_data_high)
  n_low <- nrow(test2_data_low)

  successes_high <- sum(test2_data_high$ci_lower <= test2_data_high$true_att &
                        test2_data_high$ci_upper >= test2_data_high$true_att)
  successes_low <- sum(test2_data_low$ci_lower <= test2_data_low$true_att &
                       test2_data_low$ci_upper >= test2_data_low$true_att)

  pooled_p <- (successes_high + successes_low) / (n_high + n_low)
  se_diff <- sqrt(pooled_p * (1 - pooled_p) * (1/n_high + 1/n_low))
  z_diff <- (coverage_high - coverage_low) / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))

  cat("Test 2: β=3 (high) vs β=1 (low) coverage difference\n")
  cat(sprintf("  Sample: Tree-DML at n=800\n"))
  cat(sprintf("  Coverage β=3: %.3f (%.1f%%)\n", coverage_high, 100 * coverage_high))
  cat(sprintf("  Coverage β=1: %.3f (%.1f%%)\n", coverage_low, 100 * coverage_low))
  cat(sprintf("  Difference: %.3f (%.1f percentage points)\n",
              coverage_high - coverage_low,
              100 * (coverage_high - coverage_low)))
  cat(sprintf("  Z-statistic: %.3f\n", z_diff))
  cat(sprintf("  P-value: %.5f\n", p_diff))
  if (p_diff < 0.05) {
    cat(sprintf("  ✓ Significantly different (p < 0.05)\n"))
    cat(sprintf("    → β regimes show measurable performance difference\n"))
  } else {
    cat(sprintf("  ✗ Not significantly different (p ≥ 0.05)\n"))
    cat(sprintf("    → Performance similar across β regimes at n=800\n"))
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Summary Report
#------------------------------------------------------------------------------

cat(strrep("=", 70), "\n")
cat("Analysis Summary\n")
cat(strrep("=", 70), "\n\n")

cat("Generated outputs:\n")
cat(sprintf("  • Figure 0 (s_n verification): %s\n",
            file.path(fig_dir, "fig0_tree_complexity_verification.pdf")))
cat(sprintf("  • Table 1: %s\n",
            file.path(results_dir, "table1_performance_by_beta.csv")))
cat(sprintf("  • Figure 1 (coverage): %s\n",
            file.path(fig_dir, "fig1_coverage_by_beta.pdf")))
cat(sprintf("  • Figure 2 (RMSE): %s\n",
            file.path(fig_dir, "fig2_rmse_by_beta.pdf")))
cat(sprintf("  • Figure 3 (CI width): %s\n",
            file.path(fig_dir, "fig3_ci_width_by_beta.pdf")))

cat("\n\nNext steps:\n")
cat("  1. Review figures and table for manuscript integration\n")
cat("  2. Update session notes with key findings\n")
cat("  3. Draft manuscript text for Section 4 or Appendix C\n")
cat("  4. Consider additional analyses if patterns suggest them\n\n")

cat("Analysis complete.\n")
