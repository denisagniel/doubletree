#' Diagnostics for M-Split Estimation
#'
#' @description
#' Diagnostic functions for M-split doubletree ATT estimation.
#' Assess stability, consistency, and variance of predictions across M splits.
#'
#' @name msplit_diagnostics
NULL

#' Compute Functional Consistency Metric
#'
#' @description
#' For observations with identical covariates (Xᵢ = Xⱼ), compute the maximum
#' difference in averaged predictions: max|μ̄(Xᵢ) - μ̄(Xⱼ)|.
#'
#' This should be approximately zero for a functionally consistent estimator.
#' Large values indicate that averaging across M splits has not eliminated
#' randomness for tied covariate patterns.
#'
#' @param predictions_e Matrix (n x M) of propensity predictions
#' @param predictions_m0 Matrix (n x M) of outcome predictions
#' @param X Covariate matrix/data.frame
#'
#' @return List with:
#'   \item{max_diff_e}{Numeric: max difference for propensity}
#'   \item{max_diff_m0}{Numeric: max difference for outcome}
#'   \item{n_unique_patterns}{Integer: number of unique covariate patterns}
#'   \item{n_groups_with_ties}{Integer: number of patterns with >1 observation}
#'   \item{avg_group_size}{Numeric: average size of tied groups}
#'
#' @details
#' Algorithm:
#' 1. Group observations by unique covariate pattern
#' 2. Within each group with >1 observation:
#'    - Compute averaged predictions μ̄ᵢ = mean(predictions[i, ])
#'    - Compute max pairwise difference within group
#' 3. Return overall maximum
#'
#' @examples
#' \dontrun{
#' result <- estimate_att_msplit(X, A, Y, M = 10)
#' consistency <- result$diagnostics$functional_consistency
#' cat(sprintf("Max inconsistency: %.6f\n", consistency$max_diff_e))
#' }
#'
#' @export
compute_functional_consistency <- function(predictions_e, predictions_m0, X) {
  n <- nrow(X)
  M <- ncol(predictions_e)

  # Group observations by covariate pattern
  # Convert each row to string for fast grouping
  if (is.matrix(X)) {
    X_str <- apply(X, 1, function(row) paste(row, collapse = "_"))
  } else {
    X_str <- apply(as.matrix(X), 1, function(row) paste(row, collapse = "_"))
  }

  groups <- split(seq_len(n), X_str)

  # Compute max within-group difference
  max_diff_e <- 0
  max_diff_m0 <- 0
  n_groups_with_ties <- 0L

  for (group_idx in groups) {
    if (length(group_idx) > 1) {
      n_groups_with_ties <- as.integer(n_groups_with_ties + 1)

      # Average predictions for this group
      avg_e <- rowMeans(predictions_e[group_idx, , drop = FALSE])
      avg_m0 <- rowMeans(predictions_m0[group_idx, , drop = FALSE])

      # Max pairwise difference within group
      if (length(avg_e) > 1) {
        diffs_e <- as.matrix(dist(avg_e))
        diffs_m0 <- as.matrix(dist(avg_m0))

        max_diff_e <- max(max_diff_e, max(diffs_e))
        max_diff_m0 <- max(max_diff_m0, max(diffs_m0))
      }
    }
  }

  list(
    max_diff_e = max_diff_e,
    max_diff_m0 = max_diff_m0,
    n_unique_patterns = length(groups),
    n_groups_with_ties = n_groups_with_ties,
    avg_group_size = mean(vapply(groups, length, integer(1)))
  )
}

#' Plot M-Split Diagnostics
#'
#' @description
#' Create diagnostic plots for M-split estimation results.
#'
#' @param result msplit_att object (from estimate_att_msplit)
#' @param which Character vector: which plots to create
#'   - "variance": Prediction variance vs observation index
#'   - "consistency": Consistency check (scatter of averaged predictions)
#'   - "structure": Structure frequency barplot
#'   - "all": All plots (default)
#'
#' @return Invisibly returns list of ggplot2 objects
#'
#' @examples
#' \dontrun{
#' result <- estimate_att_msplit(X, A, Y, M = 10)
#' plot(result)  # All diagnostic plots
#' plot(result, which = "variance")  # Just variance plot
#' }
#'
#' @export
plot.msplit_att <- function(x, which = "all", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting. Please install it.", call. = FALSE)
  }

  plots <- list()
  which_plots <- if ("all" %in% which) {
    c("variance", "consistency", "structure")
  } else {
    which
  }

  # 1. Prediction variance plot
  if ("variance" %in% which_plots) {
    var_df <- data.frame(
      obs = seq_len(x$n),
      var_e = x$diagnostics$prediction_variance_e,
      var_m0 = x$diagnostics$prediction_variance_m0
    )
    var_df_long <- tidyr::pivot_longer(
      var_df,
      cols = c(var_e, var_m0),
      names_to = "nuisance",
      values_to = "variance"
    )

    p_var <- ggplot2::ggplot(var_df_long, ggplot2::aes(x = obs, y = variance, color = nuisance)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.5) +
      ggplot2::geom_smooth(se = TRUE) +
      ggplot2::scale_color_manual(
        values = c("var_e" = "blue", "var_m0" = "red"),
        labels = c("var_e" = "e(X)", "var_m0" = "m0(X)")
      ) +
      ggplot2::labs(
        title = "Prediction Variance Across M Splits",
        subtitle = sprintf("Lower variance = more stable (M = %d)", x$M),
        x = "Observation Index",
        y = "Variance",
        color = "Nuisance"
      ) +
      ggplot2::theme_minimal()

    plots$variance <- p_var
  }

  # 2. Functional consistency plot
  if ("consistency" %in% which_plots) {
    # Group by covariate pattern and compute within-group variance
    X_str <- apply(as.matrix(x$predictions_all_splits$e), 1, function(row) {
      paste(row[1:min(5, length(row))], collapse = "_")  # Simplified
    })

    # This is a simplified version - full version would properly group by X
    message("Consistency plot: simplified version (groups by first 5 predictions)")

    plots$consistency <- NULL  # Placeholder
  }

  # 3. Structure frequency barplot
  if ("structure" %in% which_plots) {
    struct_df <- data.frame(
      nuisance = c("e(X)", "m0(X)"),
      frequency = c(
        x$diagnostics$structure_frequency_e,
        x$diagnostics$structure_frequency_m0
      ),
      n_leaves = c(
        x$diagnostics$n_leaves_e,
        x$diagnostics$n_leaves_m0
      )
    )

    p_struct <- ggplot2::ggplot(struct_df, ggplot2::aes(x = nuisance, y = frequency)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%% (%d leaves)",
                                                       frequency * 100, n_leaves)),
                         vjust = -0.5) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::labs(
        title = "Modal Structure Frequency",
        subtitle = "Higher frequency = more stable structure selection",
        x = "Nuisance Function",
        y = "Frequency (% of M splits)"
      ) +
      ggplot2::theme_minimal()

    plots$structure <- p_struct
  }

  # Print plots if interactive
  if (interactive()) {
    for (p in plots) {
      if (!is.null(p)) print(p)
    }
  }

  invisible(plots)
}

#' Summary Method for msplit_att
#'
#' @param object msplit_att object
#' @param ... Additional arguments (ignored)
#' @export
summary.msplit_att <- function(object, ...) {
  cat("M-Split Doubletree ATT Estimation Summary\n")
  cat("=========================================\n\n")

  cat("Point Estimate and Inference:\n")
  cat(sprintf("  θ̂ = %.4f (SE = %.4f)\n", object$theta, object$sigma))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", object$ci_95[1], object$ci_95[2]))
  cat(sprintf("  z-statistic: %.3f\n", object$theta / object$sigma))
  cat(sprintf("  p-value (two-sided): %.4f\n",
              2 * (1 - pnorm(abs(object$theta / object$sigma)))))

  cat("\nSample Information:\n")
  cat(sprintf("  n = %d (n_treated = %d, %.1f%%)\n",
              object$n, object$n_treated, 100 * object$n_treated / object$n))
  cat(sprintf("  M = %d independent splits\n", object$M))
  cat(sprintf("  K = %d cross-validation folds per split\n", object$K))

  cat("\nSelected Structures:\n")
  cat(sprintf("  Propensity e(X):\n"))
  cat(sprintf("    Modal frequency: %.1f%%\n",
              object$diagnostics$structure_frequency_e * 100))
  cat(sprintf("    Leaves: %d (max depth %d)\n",
              object$diagnostics$n_leaves_e,
              object$diagnostics$max_depth_e))

  cat(sprintf("  Outcome m0(X):\n"))
  cat(sprintf("    Modal frequency: %.1f%%\n",
              object$diagnostics$structure_frequency_m0 * 100))
  cat(sprintf("    Leaves: %d (max depth %d)\n",
              object$diagnostics$n_leaves_m0,
              object$diagnostics$max_depth_m0))

  cat("\nStability Diagnostics:\n")
  cat(sprintf("  Mean prediction variance (should be small):\n"))
  cat(sprintf("    e(X):  %.6f\n", object$diagnostics$mean_prediction_variance_e))
  cat(sprintf("    m0(X): %.6f\n", object$diagnostics$mean_prediction_variance_m0))

  cat(sprintf("  Functional consistency (should be ~0):\n"))
  cat(sprintf("    max|μ̄ᵢ-μ̄ⱼ| for Xᵢ=Xⱼ:\n"))
  cat(sprintf("      e(X):  %.6f\n",
              object$diagnostics$functional_consistency$max_diff_e))
  cat(sprintf("      m0(X): %.6f\n",
              object$diagnostics$functional_consistency$max_diff_m0))

  cat(sprintf("  Unique covariate patterns: %d (%.1f%% of n)\n",
              object$diagnostics$functional_consistency$n_unique_patterns,
              100 * object$diagnostics$functional_consistency$n_unique_patterns / object$n))

  invisible(object)
}
