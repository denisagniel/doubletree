# analyze_results.R
# Analysis and visualization for extended Rashomon-DML simulations
#
# Generates:
#   - Figure 1: Intersection existence by n and epsilon
#   - Figure 2: Bias-interpretability tradeoff
#   - Figure 3: Coverage and CI width
#   - Figure 4: Rashomon overhead vs epsilon
#   - Table 1: Summary across all conditions

# Setup ------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 required for visualization")
}
library(ggplot2)

results_dir <- "simulations/results_extended"
figures_dir <- "simulations/figures"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Load results (from batch processing)
summary_table <- readRDS(file.path(results_dir, "simulation_summary.rds"))

message("Summary table: ", nrow(summary_table), " rows")

# Figure 1: Intersection existence by n and epsilon ---------------------------

# Extract pct_nonempty for Rashomon method only
rashomon_data <- summary_table[summary_table$method == "rashomon", ]

# Average across nuisances (e, m0, m1) for simplicity
if ("pct_nonempty_e" %in% names(rashomon_data)) {
  rashomon_data$pct_nonempty_avg <- rowMeans(
    rashomon_data[, c("pct_nonempty_e", "pct_nonempty_m0", "pct_nonempty_m1")],
    na.rm = TRUE
  )
} else {
  warning("pct_nonempty columns not found; skipping Figure 1")
  rashomon_data$pct_nonempty_avg <- NA
}

fig1 <- ggplot(rashomon_data, aes(x = factor(n), y = pct_nonempty_avg,
                                   fill = factor(epsilon), group = epsilon)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ dgp, ncol = 2, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_fill_brewer(palette = "Blues", name = expression(epsilon)) +
  labs(
    title = "Intersection Existence by Sample Size and Tolerance",
    x = "Sample Size (n)",
    y = "% Replications with Non-Empty Intersection",
    caption = "Averaged across e, m0, m1 nuisances"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "figure1_intersection_existence.pdf"),
       fig1, width = 10, height = 8)
message("Saved Figure 1: intersection existence")

# Figure 2: Bias-interpretability tradeoff ------------------------------------

# Compare bias across methods and epsilons
# Focus on n=400 for clarity (or facet by n if desired)
bias_data <- summary_table[summary_table$n == 400, ]

fig2 <- ggplot(bias_data, aes(x = method, y = abs(bias), fill = method)) +
  geom_boxplot(aes(group = interaction(method, epsilon))) +
  facet_grid(dgp ~ factor(epsilon), labeller = label_both, scales = "free_y") +
  scale_fill_manual(
    values = c("fold_specific" = "#1f77b4", "rashomon" = "#ff7f0e", "oracle" = "#2ca02c"),
    labels = c("Fold-Specific", "Rashomon", "Oracle")
  ) +
  labs(
    title = "Bias-Interpretability Tradeoff (n = 400)",
    x = "Method",
    y = "Absolute Bias",
    fill = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "figure2_bias_tradeoff.pdf"),
       fig2, width = 12, height = 10)
message("Saved Figure 2: bias-interpretability tradeoff")

# Figure 3: Coverage and CI width ----------------------------------------------

# Coverage (target 0.95)
coverage_data <- summary_table[summary_table$method != "oracle", ]  # Oracle has no inference

fig3a <- ggplot(coverage_data, aes(x = factor(n), y = coverage_95,
                                    color = method, group = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(dgp ~ epsilon, labeller = label_both) +
  scale_color_manual(
    values = c("fold_specific" = "#1f77b4", "rashomon" = "#ff7f0e"),
    labels = c("Fold-Specific", "Rashomon")
  ) +
  scale_y_continuous(limits = c(0.8, 1.0)) +
  labs(
    title = "95% CI Coverage by Method",
    x = "Sample Size (n)",
    y = "Coverage",
    color = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "figure3a_coverage.pdf"),
       fig3a, width = 12, height = 10)
message("Saved Figure 3a: coverage")

# CI width
fig3b <- ggplot(coverage_data, aes(x = factor(n), y = mean_ci_width,
                                    color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(dgp ~ epsilon, labeller = label_both, scales = "free_y") +
  scale_color_manual(
    values = c("fold_specific" = "#1f77b4", "rashomon" = "#ff7f0e"),
    labels = c("Fold-Specific", "Rashomon")
  ) +
  labs(
    title = "Mean 95% CI Width by Method",
    x = "Sample Size (n)",
    y = "Mean CI Width",
    color = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "figure3b_ci_width.pdf"),
       fig3b, width = 12, height = 10)
message("Saved Figure 3b: CI width")

# Figure 4: Rashomon overhead vs epsilon --------------------------------------

# Compute overhead as bias_rashomon - bias_fold_specific
# Group by dgp, n, epsilon
overhead_data <- merge(
  summary_table[summary_table$method == "rashomon", c("dgp", "n", "epsilon", "bias")],
  summary_table[summary_table$method == "fold_specific", c("dgp", "n", "epsilon", "bias")],
  by = c("dgp", "n", "epsilon"),
  suffixes = c("_rashomon", "_fold")
)
overhead_data$overhead <- overhead_data$bias_rashomon - overhead_data$bias_fold

fig4 <- ggplot(overhead_data, aes(x = epsilon, y = abs(overhead), color = factor(n))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ dgp, ncol = 2, labeller = label_both, scales = "free_y") +
  scale_color_brewer(palette = "Set1", name = "Sample Size") +
  labs(
    title = "Rashomon Overhead vs Tolerance",
    x = expression(paste("Rashomon Tolerance (", epsilon, ")")),
    y = "|Bias(Rashomon) - Bias(Fold-Specific)|",
    caption = "Excess bias from using intersection tree vs fold-specific optimum"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "figure4_rashomon_overhead.pdf"),
       fig4, width = 10, height = 8)
message("Saved Figure 4: Rashomon overhead")

# Table 1: Summary across all conditions --------------------------------------

# Aggregate by method and DGP (averaging over n and epsilon)
table1 <- aggregate(
  cbind(bias, rmse, coverage_95, mean_ci_width) ~ method + dgp,
  data = summary_table[summary_table$n == 400, ],  # Focus on n=400
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)

# Add pct_nonempty for Rashomon only
if ("pct_nonempty_avg" %in% names(rashomon_data)) {
  pct_nonempty_summary <- aggregate(
    pct_nonempty_avg ~ dgp,
    data = rashomon_data[rashomon_data$n == 400, ],
    FUN = function(x) round(mean(x, na.rm = TRUE), 3)
  )
  names(pct_nonempty_summary)[2] <- "pct_nonempty"

  table1 <- merge(table1, pct_nonempty_summary, by = "dgp", all.x = TRUE)
}

# Reorder columns
table1 <- table1[, c("dgp", "method", "bias", "rmse", "coverage_95",
                     "mean_ci_width", if("pct_nonempty" %in% names(table1)) "pct_nonempty")]

write.csv(table1, file.path(figures_dir, "table1_summary.csv"), row.names = FALSE)
message("Saved Table 1: summary")

print(table1)

# Key findings summary ---------------------------------------------------------

message("\n=== Key Findings ===")

# Intersection existence
if ("pct_nonempty_avg" %in% names(rashomon_data)) {
  smooth_dgps <- rashomon_data[rashomon_data$dgp %in% c("dgp1", "dgp2", "dgp3"), ]
  rough_dgp <- rashomon_data[rashomon_data$dgp == "dgp4", ]

  message("Intersection non-empty:")
  message("  Smooth DGPs (1-3): ",
          round(100 * mean(smooth_dgps$pct_nonempty_avg, na.rm = TRUE), 1), "%")
  message("  Rough DGP (4): ",
          round(100 * mean(rough_dgp$pct_nonempty_avg, na.rm = TRUE), 1), "%")
}

# Bias overhead
if (nrow(overhead_data) > 0) {
  message("\nRashomon overhead (mean |bias difference|):")
  overhead_summary <- aggregate(abs(overhead) ~ dgp, data = overhead_data,
                                FUN = function(x) round(mean(x, na.rm = TRUE), 4))
  for (i in 1:nrow(overhead_summary)) {
    message("  ", overhead_summary$dgp[i], ": ", overhead_summary$V1[i])
  }
}

# Coverage
coverage_rashomon <- summary_table[summary_table$method == "rashomon" & summary_table$n >= 400, ]
message("\nRashomon coverage (n >= 400):")
message("  Mean: ", round(mean(coverage_rashomon$coverage_95, na.rm = TRUE), 3))
message("  Min: ", round(min(coverage_rashomon$coverage_95, na.rm = TRUE), 3))
message("  Max: ", round(max(coverage_rashomon$coverage_95, na.rm = TRUE), 3))

message("\n=== Analysis Complete ===")
message("Figures saved to: ", figures_dir)
