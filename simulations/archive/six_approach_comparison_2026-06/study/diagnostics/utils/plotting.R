# Plotting Utilities for Diagnostics
# Created: 2026-05-27

library(ggplot2)

#' Plot propensity score distribution
#'
#' @param e_hat Estimated propensity scores
#' @param e_true True propensity scores (optional)
#' @param A Treatment indicator (optional)
#' @return ggplot object
#' @export
plot_propensity_distribution <- function(e_hat, e_true = NULL, A = NULL) {
  df <- data.frame(e_hat = e_hat)

  if (!is.null(e_true)) {
    df$e_true <- e_true
  }

  if (!is.null(A)) {
    df$treatment <- factor(A, levels = c(0, 1), labels = c("Control", "Treated"))
  }

  p <- ggplot(df, aes(x = e_hat))

  if (!is.null(A)) {
    p <- p + geom_histogram(aes(fill = treatment), bins = 30, alpha = 0.6, position = "identity")
  } else {
    p <- p + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8)
  }

  if (!is.null(e_true)) {
    p <- p + geom_histogram(aes(x = e_true), bins = 30, fill = "red", alpha = 0.3)
  }

  p <- p +
    labs(
      title = "Propensity Score Distribution",
      x = "Propensity Score",
      y = "Count"
    ) +
    theme_minimal() +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", color = "darkred")

  return(p)
}

#' Plot calibration curve for propensity score
#'
#' @param e_hat Predicted propensity scores
#' @param A Observed treatment
#' @param n_bins Number of bins (default 10)
#' @return ggplot object
#' @export
plot_calibration <- function(e_hat, A, n_bins = 10) {
  bins <- cut(e_hat, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)

  calib_df <- data.frame(
    bin = bins,
    predicted = e_hat,
    observed = A
  )

  calib_summary <- aggregate(
    cbind(predicted, observed) ~ bin,
    data = calib_df,
    FUN = mean
  )

  # Add confidence intervals
  calib_summary$n <- aggregate(observed ~ bin, data = calib_df, FUN = length)$observed
  calib_summary$se <- sqrt(calib_summary$observed * (1 - calib_summary$observed) / calib_summary$n)
  calib_summary$ci_lower <- calib_summary$observed - 1.96 * calib_summary$se
  calib_summary$ci_upper <- calib_summary$observed + 1.96 * calib_summary$se

  p <- ggplot(calib_summary, aes(x = predicted, y = observed)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = "Propensity Score Calibration",
      x = "Predicted Propensity",
      y = "Observed Treatment Rate"
    ) +
    theme_minimal() +
    coord_fixed(ratio = 1) +
    lims(x = c(0, 1), y = c(0, 1))

  return(p)
}

#' Plot prediction error vs true values
#'
#' @param y_hat Predicted values
#' @param y_true True values
#' @param title Plot title
#' @return ggplot object
#' @export
plot_prediction_error <- function(y_hat, y_true, title = "Prediction Error") {
  df <- data.frame(
    y_true = y_true,
    y_hat = y_hat,
    error = y_hat - y_true
  )

  p <- ggplot(df, aes(x = y_true, y = error)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(
      title = title,
      x = "True Value",
      y = "Prediction Error (Predicted - True)"
    ) +
    theme_minimal()

  return(p)
}

#' Plot tree size distribution across folds
#'
#' @param n_leaves Vector of leaf counts
#' @param fold_ids Fold identifiers (optional)
#' @return ggplot object
#' @export
plot_tree_sizes <- function(n_leaves, fold_ids = NULL) {
  if (is.null(fold_ids)) {
    fold_ids <- seq_along(n_leaves)
  }

  df <- data.frame(
    fold = factor(fold_ids),
    n_leaves = n_leaves
  )

  p <- ggplot(df, aes(x = fold, y = n_leaves)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = mean(n_leaves), linetype = "dashed", color = "red") +
    labs(
      title = "Tree Complexity Across Folds",
      subtitle = sprintf("Mean: %.1f leaves (SD: %.1f)", mean(n_leaves), sd(n_leaves)),
      x = "Fold",
      y = "Number of Leaves"
    ) +
    theme_minimal()

  return(p)
}

#' Plot bias decomposition
#'
#' @param comp1_bias Bias from component 1
#' @param comp2_bias Bias from component 2
#' @param total_bias Total bias
#' @return ggplot object
#' @export
plot_bias_decomposition <- function(comp1_bias, comp2_bias, total_bias) {
  df <- data.frame(
    component = c("Outcome Model\n(Treated)", "Propensity Weighting\n(Control)", "Total"),
    bias = c(comp1_bias, comp2_bias, total_bias),
    type = c("Component", "Component", "Total")
  )

  p <- ggplot(df, aes(x = component, y = bias, fill = type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Component" = "steelblue", "Total" = "darkred")) +
    labs(
      title = "EIF Bias Decomposition",
      x = "",
      y = "Bias"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  return(p)
}

#' Plot lambda selection from CV
#'
#' @param lambda_grid Vector of lambda values tested
#' @param cv_losses Vector of CV losses
#' @param best_lambda Selected lambda value
#' @return ggplot object
#' @export
plot_lambda_selection <- function(lambda_grid, cv_losses, best_lambda) {
  df <- data.frame(
    lambda = lambda_grid,
    cv_loss = cv_losses
  )

  p <- ggplot(df, aes(x = lambda, y = cv_loss)) +
    geom_line(color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    geom_vline(xintercept = best_lambda, linetype = "dashed", color = "red") +
    scale_x_log10() +
    labs(
      title = "Cross-Validation: Lambda Selection",
      subtitle = sprintf("Selected λ = %.4f", best_lambda),
      x = "Lambda (log scale)",
      y = "CV Loss"
    ) +
    theme_minimal()

  return(p)
}

#' Create multi-panel diagnostic plot
#'
#' @param e_hat Estimated propensity scores
#' @param e_true True propensity scores
#' @param A Treatment indicator
#' @param mu0_hat Estimated outcome predictions
#' @param mu0_true True outcomes
#' @return Combined ggplot object (requires gridExtra)
#' @export
plot_diagnostics_grid <- function(e_hat, e_true, A, mu0_hat = NULL, mu0_true = NULL) {
  require(gridExtra)

  # Plot 1: Propensity distribution
  p1 <- plot_propensity_distribution(e_hat, e_true, A)

  # Plot 2: Calibration
  p2 <- plot_calibration(e_hat, A)

  # Plot 3: Propensity error
  p3 <- plot_prediction_error(e_hat, e_true, title = "Propensity Score Error")

  # Plot 4: Outcome error (if provided)
  if (!is.null(mu0_hat) && !is.null(mu0_true)) {
    p4 <- plot_prediction_error(mu0_hat, mu0_true, title = "Outcome Model Error")
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  } else {
    grid.arrange(p1, p2, p3, ncol = 2)
  }
}
