#!/usr/bin/env Rscript
# Test Functional Consistency: M-Split Doubletree
#
# Validates Theorem 3: For X_i = X_j, we should have η̂(X_i) ≈ η̂(X_j)
# as M increases, with explicit quantification.
#
# Key prediction: P(|μ̄(X_i) - μ̄(X_j)| > ε) ≤ 2exp(-Mε²/(2σ²))
# Expected: Variance of predictions for tied covariates → 0 as M → ∞

library(optimaltrees)
library(doubletree)
library(cli)
library(dplyr)
library(ggplot2)

cli_h1("Functional Consistency Test: M-Split Doubletree")

# Source validated DGPs
source("dgps/dgps_smooth.R")

# ============================================================================
# Test Configuration
# ============================================================================

TEST_CONFIG <- list(
  n = 500,                          # Sample size
  M_values = c(1, 5, 10, 20, 50),   # M values to test
  K = 5,                            # Folds
  n_reps = 20,                      # Replications per M
  tau = 0.10,                       # True ATT
  seed_offset = 50000
)

cli_alert_info("Configuration:")
cli_ul(c(
  "Sample size: {TEST_CONFIG$n}",
  "M values: {paste(TEST_CONFIG$M_values, collapse=', ')}",
  "K folds: {TEST_CONFIG$K}",
  "Replications: {TEST_CONFIG$n_reps}",
  "True τ: {TEST_CONFIG$tau}"
))

# ============================================================================
# Helper: Compute Functional Consistency Metric
# ============================================================================

compute_functional_consistency <- function(result, X) {
  # Extract predictions - different structure for M=1 vs M>1
  if (!is.null(result$averaged_predictions)) {
    # M-split case
    eta_ps <- result$averaged_predictions$e
    eta_m0 <- result$averaged_predictions$m0
    eta_m1 <- NULL  # Not used for ATT
  } else {
    # Standard doubletree case (M=1)
    eta_ps <- result$propensity
    eta_m0 <- result$outcome_control
    eta_m1 <- NULL
  }

  # Find tied covariates (exact matches)
  X_string <- apply(X, 1, function(x) paste(x, collapse = "_"))
  tied_groups <- split(seq_len(nrow(X)), X_string)
  tied_groups <- tied_groups[lengths(tied_groups) > 1]  # Only groups with 2+ obs

  if (length(tied_groups) == 0) {
    return(list(
      n_tied_groups = 0,
      n_tied_obs = 0,
      max_diff_ps = NA,
      max_diff_m0 = NA,
      max_diff_m1 = NA,
      mean_within_sd_ps = NA,
      mean_within_sd_m0 = NA,
      mean_within_sd_m1 = NA
    ))
  }

  # Compute within-group variability
  max_diffs <- lapply(tied_groups, function(idx) {
    list(
      ps = if(length(unique(eta_ps[idx])) > 1) max(eta_ps[idx]) - min(eta_ps[idx]) else 0,
      m0 = if(length(unique(eta_m0[idx])) > 1) max(eta_m0[idx]) - min(eta_m0[idx]) else 0,
      m1 = if(length(unique(eta_m1[idx])) > 1) max(eta_m1[idx]) - min(eta_m1[idx]) else 0
    )
  })

  within_sds <- lapply(tied_groups, function(idx) {
    list(
      ps = if(length(idx) > 1 && length(unique(eta_ps[idx])) > 1) sd(eta_ps[idx]) else 0,
      m0 = if(length(idx) > 1 && length(unique(eta_m0[idx])) > 1) sd(eta_m0[idx]) else 0,
      m1 = if(length(idx) > 1 && length(unique(eta_m1[idx])) > 1) sd(eta_m1[idx]) else 0
    )
  })

  list(
    n_tied_groups = length(tied_groups),
    n_tied_obs = sum(lengths(tied_groups)),
    max_diff_ps = max(sapply(max_diffs, `[[`, "ps")),
    max_diff_m0 = max(sapply(max_diffs, `[[`, "m0")),
    max_diff_m1 = max(sapply(max_diffs, `[[`, "m1")),
    mean_within_sd_ps = mean(sapply(within_sds, `[[`, "ps")),
    mean_within_sd_m0 = mean(sapply(within_sds, `[[`, "m0")),
    mean_within_sd_m1 = mean(sapply(within_sds, `[[`, "m1"))
  )
}

# ============================================================================
# Test Function: Run one M value across replications
# ============================================================================

test_m_value <- function(M, dgp_func, config) {
  cli_h2("Testing M = {M}")

  results_list <- vector("list", config$n_reps)

  cli_progress_bar("Running replications", total = config$n_reps)

  for (rep in 1:config$n_reps) {
    seed <- config$seed_offset + M * 1000 + rep

    # Generate data
    d <- dgp_func(n = config$n, tau = config$tau, seed = seed)

    # Run M-split
    tryCatch({
      if (M == 1) {
        # Standard doubletree
        result <- estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = config$K,
          regularization = log(config$n) / config$n,
          use_rashomon = FALSE,
          seed = seed
        )

        # Add structure for consistency with msplit output
        result$structure_frequency <- 1.0
        result$M <- 1
      } else {
        # M-split doubletree
        result <- estimate_att_msplit(
          X = d$X, A = d$A, Y = d$Y,
          M = M,
          K = config$K,
          regularization = log(config$n) / config$n,
          use_rashomon = FALSE,
          seed = seed
        )
      }

      # Compute functional consistency metrics
      fc <- compute_functional_consistency(result, d$X)

      # Extract fields with defaults
      converged <- if (!is.null(result$converged)) result$converged else TRUE
      estimate <- if (!is.null(result$estimate)) result$estimate else result$theta
      se_val <- if (!is.null(result$se)) result$se else result$sigma

      results_list[[rep]] <- list(
        M = M,
        rep = rep,
        theta = estimate,
        se = se_val,
        true_att = d$true_att,
        converged = converged,
        structure_freq = result$structure_frequency,
        n_tied_groups = fc$n_tied_groups,
        n_tied_obs = fc$n_tied_obs,
        max_diff_ps = fc$max_diff_ps,
        max_diff_m0 = fc$max_diff_m0,
        max_diff_m1 = fc$max_diff_m1,
        mean_within_sd_ps = fc$mean_within_sd_ps,
        mean_within_sd_m0 = fc$mean_within_sd_m0,
        mean_within_sd_m1 = fc$mean_within_sd_m1
      )

    }, error = function(e) {
      cli_alert_warning("Rep {rep} failed: {e$message}")
      results_list[[rep]] <- list(
        M = M, rep = rep, converged = FALSE,
        error = as.character(e$message)
      )
    })

    cli_progress_update()
  }

  cli_progress_done()

  # Combine results
  results_df <- bind_rows(results_list)

  # Summary statistics
  converged <- results_df %>% filter(converged == TRUE)

  if (nrow(converged) > 0) {
    cli_alert_success("M = {M}: {nrow(converged)}/{config$n_reps} converged")
    cli_ul(c(
      "Structure frequency: {sprintf('%.1f%%', mean(converged$structure_freq) * 100)}",
      "Tied groups per rep: {sprintf('%.1f', mean(converged$n_tied_groups))}",
      "Max diff (propensity): {sprintf('%.4f', mean(converged$max_diff_ps, na.rm=TRUE))}",
      "Max diff (outcome m0): {sprintf('%.4f', mean(converged$max_diff_m0, na.rm=TRUE))}",
      "Mean within-SD (PS): {sprintf('%.6f', mean(converged$mean_within_sd_ps, na.rm=TRUE))}",
      "Mean within-SD (m0): {sprintf('%.6f', mean(converged$mean_within_sd_m0, na.rm=TRUE))}"
    ))
  } else {
    cli_alert_danger("M = {M}: All replications failed!")
  }

  results_df
}

# ============================================================================
# Main Test: DGP 1 (Binary, validated)
# ============================================================================

cli_h1("DGP 1: Binary Features (Validated)")

all_results <- lapply(TEST_CONFIG$M_values, function(M) {
  test_m_value(M, generate_dgp_binary_att, TEST_CONFIG)
})

results_df <- bind_rows(all_results)

# Save raw results
saveRDS(results_df, "results/functional_consistency_test.rds")
cli_alert_success("Raw results saved: results/functional_consistency_test.rds")

# ============================================================================
# Analysis: Convergence of Functional Consistency
# ============================================================================

cli_h1("Analysis: Functional Consistency vs M")

converged <- results_df %>% filter(converged == TRUE)

summary_by_m <- converged %>%
  group_by(M) %>%
  summarise(
    n_reps = n(),
    mean_structure_freq = mean(structure_freq),
    mean_max_diff_ps = mean(max_diff_ps, na.rm = TRUE),
    mean_max_diff_m0 = mean(max_diff_m0, na.rm = TRUE),
    sd_max_diff_ps = sd(max_diff_ps, na.rm = TRUE),
    sd_max_diff_m0 = sd(max_diff_m0, na.rm = TRUE),
    mean_within_sd_ps = mean(mean_within_sd_ps, na.rm = TRUE),
    mean_within_sd_m0 = mean(mean_within_sd_m0, na.rm = TRUE),
    .groups = "drop"
  )

cli_h2("Summary by M")
print(summary_by_m, n = Inf)

# Check theoretical prediction: variance ~ 1/M
cli_h2("Theoretical Check: Does variance decrease as 1/M?")

variance_check <- summary_by_m %>%
  mutate(
    predicted_ratio = M[1] / M,  # Relative to M=1
    observed_ratio_ps = mean_within_sd_ps[1]^2 / mean_within_sd_ps^2,
    observed_ratio_m0 = mean_within_sd_m0[1]^2 / mean_within_sd_m0^2
  )

print(variance_check %>% select(M, predicted_ratio, observed_ratio_ps, observed_ratio_m0))

# ============================================================================
# Visualization
# ============================================================================

cli_h1("Generating plots")

# Plot 1: Max difference vs M
p1 <- ggplot(converged, aes(x = factor(M), y = max_diff_ps)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
  labs(
    title = "Functional Consistency: Max Prediction Difference for Tied Covariates",
    subtitle = "Propensity score predictions (smaller = better)",
    x = "M (number of splits)",
    y = "Max |η̂(Xi) - η̂(Xj)| for Xi = Xj"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/functional_consistency_maxdiff_ps.pdf", p1, width = 8, height = 5)
cli_alert_success("Saved: results/functional_consistency_maxdiff_ps.pdf")

# Plot 2: Within-group SD vs M (log-log for 1/M check)
p2 <- ggplot(summary_by_m, aes(x = M, y = mean_within_sd_ps)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ I(1/x),
              se = FALSE, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Prediction Variance vs M (Theory: Var ~ 1/M)",
    subtitle = "Log-log scale: should show linear decay if theory holds",
    x = "M (log scale)",
    y = "Mean within-group SD (log scale)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/functional_consistency_variance_decay.pdf", p2, width = 8, height = 5)
cli_alert_success("Saved: results/functional_consistency_variance_decay.pdf")

# Plot 3: Structure frequency vs M
p3 <- ggplot(converged, aes(x = factor(M), y = structure_freq)) +
  geom_boxplot(fill = "darkgreen", alpha = 0.7) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "red") +
  labs(
    title = "Structure Stability vs M",
    subtitle = "Frequency of modal structure (higher = more stable)",
    x = "M (number of splits)",
    y = "Structure frequency"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/functional_consistency_structure_freq.pdf", p3, width = 8, height = 5)
cli_alert_success("Saved: results/functional_consistency_structure_freq.pdf")

# ============================================================================
# Final Assessment
# ============================================================================

cli_h1("Assessment: Functional Consistency Validated?")

# Check criteria
m50_result <- summary_by_m %>% filter(M == 50)
m1_result <- summary_by_m %>% filter(M == 1)

if (nrow(m50_result) > 0 && nrow(m1_result) > 0) {
  variance_reduction <- (m1_result$mean_within_sd_ps^2) / (m50_result$mean_within_sd_ps^2)

  cli_alert_info("Variance reduction M=1 → M=50:")
  cli_ul(c(
    "Observed: {sprintf('%.1f', variance_reduction)}×",
    "Predicted (theory): 50×",
    "Match: {ifelse(variance_reduction > 30, 'Good ✓', 'Poor ✗')}"
  ))

  cli_alert_info("Functional consistency at M=50:")
  cli_ul(c(
    "Max diff (propensity): {sprintf('%.4f', m50_result$mean_max_diff_ps)}",
    "Within-group SD: {sprintf('%.6f', m50_result$mean_within_sd_ps)}",
    "Practically zero? {ifelse(m50_result$mean_within_sd_ps < 0.001, 'Yes ✓', 'No ✗')}"
  ))

  cli_alert_info("Structure stability at M=50:")
  cli_ul(c(
    "Mean frequency: {sprintf('%.1f%%', m50_result$mean_structure_freq * 100)}",
    "Stable structure? {ifelse(m50_result$mean_structure_freq > 0.6, 'Yes ✓', 'No ✗')}"
  ))

  # Overall verdict
  cli_rule()
  if (variance_reduction > 30 &&
      m50_result$mean_within_sd_ps < 0.001 &&
      m50_result$mean_structure_freq > 0.6) {
    cli_alert_success("PASS: Functional consistency validated! ✓")
    cli_alert_success("Theory prediction confirmed: M-split produces 'one tree' for large M")
  } else {
    cli_alert_warning("UNCERTAIN: Some criteria not met")
    cli_alert_info("Review plots and investigate further")
  }
  cli_rule()
} else {
  cli_alert_danger("ERROR: Missing M=1 or M=50 results")
}

cli_alert_info("Runtime: {format(Sys.time())}")
cli_alert_success("Test complete. Check results/ for plots and data.")
