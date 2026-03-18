# analyze_strategic_results.R
# Quick analysis of strategic subset for iteration
#
# Generates preliminary figures/tables to inform full grid decision

# Setup ------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 required")
}
library(ggplot2)
library(cli)

results_dir <- "simulations/results_extended"
figures_dir <- "simulations/figures_strategic"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Load and combine results -----------------------------------------------------

# Filter to strategic subset only: DGP 3-4, n in 200/400/800, eps in 0.05/0.1
all_files <- list.files(results_dir, pattern = "^result_.*\\.rds$", full.names = FALSE)
result_files <- grep("result_dgp[34]_n(200|400|800)_eps0\\.(050|100)\\.rds",
                     all_files, value = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in ", results_dir)
}

cli_alert_info("Found {length(result_files)} result files")

# Load results and annotate with filename metadata
all_results <- list()
for (fname in result_files) {
  # Parse filename: result_dgp3_n200_eps0.050.rds
  parts <- strsplit(fname, "_")[[1]]
  dgp <- sub("result", "", parts[1])
  if (dgp == "") dgp <- parts[2]  # handle result_dgp3 format
  dgp <- gsub("dgp", "", dgp)
  n_val <- as.integer(sub("n", "", parts[which(grepl("^n[0-9]", parts))]))
  eps_val <- as.numeric(sub("eps", "", sub("\\.rds", "", parts[which(grepl("^eps", parts))])))

  # Load results from this file
  file_results <- readRDS(file.path(results_dir, fname))

  # Annotate each result with metadata
  for (i in seq_along(file_results)) {
    if (!is.null(file_results[[i]])) {
      file_results[[i]]$dgp <- paste0("dgp", dgp)
      file_results[[i]]$n <- n_val
      if (file_results[[i]]$method == "rashomon") {
        file_results[[i]]$epsilon <- eps_val
      }
    }
  }

  all_results <- c(all_results, file_results)
}

cli_alert_info("Loaded {length(all_results)} individual results")

# Extract summary statistics ---------------------------------------------------

extract_summary <- function(result) {
  # Handle NULL or failed results
  if (is.null(result) || is.null(result$theta)) return(NULL)

  dgp <- if ("dgp" %in% names(result)) result$dgp else NA
  n <- if ("n" %in% names(result)) result$n else NA
  method <- if ("method" %in% names(result)) result$method else NA
  epsilon <- if ("epsilon" %in% names(result)) result$epsilon else NA

  # ATT estimation
  theta <- result$theta
  true_tau <- if ("true_tau" %in% names(result)) result$true_tau else NA
  bias <- theta - true_tau

  # CI
  ci_lower <- if (!is.null(result$ci_95)) result$ci_95[1] else NA
  ci_upper <- if (!is.null(result$ci_95)) result$ci_95[2] else NA
  ci_width <- ci_upper - ci_lower
  coverage <- as.integer(!is.na(ci_lower) && !is.na(ci_upper) &&
                         ci_lower <= true_tau && true_tau <= ci_upper)

  # Rashomon-specific
  pct_nonempty_e <- if ("pct_nonempty_e" %in% names(result)) result$pct_nonempty_e else NA
  pct_nonempty_m0 <- if ("pct_nonempty_m0" %in% names(result)) result$pct_nonempty_m0 else NA
  pct_nonempty_m1 <- if ("pct_nonempty_m1" %in% names(result)) result$pct_nonempty_m1 else NA
  n_intersecting_e <- if ("n_intersecting_e" %in% names(result)) result$n_intersecting_e else NA

  data.frame(
    dgp = dgp,
    n = n,
    method = method,
    epsilon = epsilon,
    theta = theta,
    true_tau = true_tau,
    bias = bias,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_width = ci_width,
    coverage = coverage,
    pct_nonempty_e = pct_nonempty_e,
    pct_nonempty_m0 = pct_nonempty_m0,
    pct_nonempty_m1 = pct_nonempty_m1,
    n_intersecting_e = n_intersecting_e,
    stringsAsFactors = FALSE
  )
}

summary_list <- lapply(all_results, extract_summary)
summary_table <- do.call(rbind, Filter(Negate(is.null), summary_list))

cli_alert_success("Summary table: {nrow(summary_table)} rows")

# Aggregate by config ----------------------------------------------------------

aggregate_by_config <- function(df) {
  # Group by dgp, n, method, epsilon
  configs <- unique(df[, c("dgp", "n", "method", "epsilon")])

  agg_list <- lapply(1:nrow(configs), function(i) {
    cfg <- configs[i, ]
    subset_df <- df[
      df$dgp == cfg$dgp &
      df$n == cfg$n &
      df$method == cfg$method &
      (is.na(cfg$epsilon) | is.na(df$epsilon) | df$epsilon == cfg$epsilon),
    ]

    if (nrow(subset_df) == 0) return(NULL)

    data.frame(
      dgp = cfg$dgp,
      n = cfg$n,
      method = cfg$method,
      epsilon = cfg$epsilon,
      n_reps = nrow(subset_df),
      bias = mean(subset_df$bias, na.rm = TRUE),
      rmse = sqrt(mean(subset_df$bias^2, na.rm = TRUE)),
      coverage = mean(subset_df$coverage, na.rm = TRUE),
      ci_width = mean(subset_df$ci_width, na.rm = TRUE),
      pct_nonempty_e = mean(subset_df$pct_nonempty_e, na.rm = TRUE),
      pct_nonempty_m0 = mean(subset_df$pct_nonempty_m0, na.rm = TRUE),
      pct_nonempty_m1 = mean(subset_df$pct_nonempty_m1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), agg_list))
}

agg_summary <- aggregate_by_config(summary_table)

cli_alert_success("Aggregated to {nrow(agg_summary)} configurations")

# Save summary
saveRDS(agg_summary, file.path(results_dir, "summary_strategic.rds"))
write.csv(agg_summary, file.path(results_dir, "summary_strategic.csv"), row.names = FALSE)

cli_alert_success("Saved summary tables")

# Print key findings -----------------------------------------------------------

cli_h1("Key Findings (Strategic Subset)")

cli_h2("Coverage by Method")
coverage_by_method <- aggregate(coverage ~ method, agg_summary, mean)
for (i in 1:nrow(coverage_by_method)) {
  cli_alert_info("{coverage_by_method$method[i]}: {round(coverage_by_method$coverage[i], 3)}")
}

cli_h2("Intersection Success (Rashomon only)")
rashomon_agg <- agg_summary[agg_summary$method == "rashomon", ]
if (nrow(rashomon_agg) > 0) {
  pct_nonempty_avg <- mean(c(
    mean(rashomon_agg$pct_nonempty_e, na.rm = TRUE),
    mean(rashomon_agg$pct_nonempty_m0, na.rm = TRUE),
    mean(rashomon_agg$pct_nonempty_m1, na.rm = TRUE)
  ), na.rm = TRUE)
  cli_alert_info("Average non-empty intersection: {round(pct_nonempty_avg * 100, 1)}%")
}

cli_h2("RMSE by Method")
rmse_by_method <- aggregate(rmse ~ method, agg_summary, mean)
for (i in 1:nrow(rmse_by_method)) {
  cli_alert_info("{rmse_by_method$method[i]}: {round(rmse_by_method$rmse[i], 4)}")
}

# Figure 1: Coverage by n and method -------------------------------------------

fig1 <- ggplot(agg_summary, aes(x = factor(n), y = coverage,
                                 fill = method, group = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~ dgp, ncol = 2, labeller = label_both) +
  scale_y_continuous(limits = c(0.8, 1), labels = scales::percent) +
  scale_fill_manual(
    values = c("fold_specific" = "#1f77b4", "rashomon" = "#ff7f0e", "oracle" = "#2ca02c"),
    labels = c("Fold-Specific", "Rashomon", "Oracle")
  ) +
  labs(
    title = "95% Confidence Interval Coverage",
    subtitle = "Strategic Subset: DGP 3 (smooth) and DGP 4 (rough)",
    x = "Sample Size (n)",
    y = "Empirical Coverage",
    fill = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "fig_coverage_strategic.pdf"), fig1, width = 10, height = 6)
cli_alert_success("Saved Figure 1: coverage")

# Figure 2: Intersection existence (Rashomon only) -----------------------------

rashomon_data <- agg_summary[agg_summary$method == "rashomon", ]

if (nrow(rashomon_data) > 0) {
  rashomon_data$pct_nonempty_avg <- rowMeans(
    rashomon_data[, c("pct_nonempty_e", "pct_nonempty_m0", "pct_nonempty_m1")],
    na.rm = TRUE
  )

  fig2 <- ggplot(rashomon_data, aes(x = factor(n), y = pct_nonempty_avg,
                                     fill = factor(epsilon))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ dgp, ncol = 2, labeller = label_both) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_fill_brewer(palette = "Blues", name = expression(epsilon)) +
    labs(
      title = "Rashomon Intersection Success Rate",
      subtitle = "% Replications with Non-Empty Intersection (avg across e, m0, m1)",
      x = "Sample Size (n)",
      y = "Success Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(file.path(figures_dir, "fig_intersection_strategic.pdf"), fig2, width = 10, height = 6)
  cli_alert_success("Saved Figure 2: intersection existence")
}

# Figure 3: Bias comparison ----------------------------------------------------

fig3 <- ggplot(agg_summary, aes(x = method, y = abs(bias), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dgp ~ factor(n), labeller = label_both, scales = "free_y") +
  scale_fill_manual(
    values = c("fold_specific" = "#1f77b4", "rashomon" = "#ff7f0e", "oracle" = "#2ca02c"),
    labels = c("Fold-Specific", "Rashomon", "Oracle")
  ) +
  labs(
    title = "Absolute Bias by Method",
    subtitle = "Strategic Subset",
    x = "Method",
    y = "|Bias|",
    fill = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_blank()
  )

ggsave(file.path(figures_dir, "fig_bias_strategic.pdf"), fig3, width = 10, height = 6)
cli_alert_success("Saved Figure 3: bias comparison")

# Summary output ---------------------------------------------------------------

cli_h1("Strategic Subset Complete!")
cli_alert_info("Figures saved to: {figures_dir}")
cli_alert_info("Summary tables saved to: {results_dir}")

cli_h2("Next Steps")
cli_alert_info("1. Review figures in {figures_dir}")
cli_alert_info("2. Check summary_strategic.csv for detailed numbers")
cli_alert_info("3. Decide: proceed with full grid or iterate on design?")
