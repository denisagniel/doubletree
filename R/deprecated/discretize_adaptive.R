# ============================================================================
# DEPRECATED: This file is deprecated as of 2025-03-11
#
# REASON: dmltree now uses treefarmr's built-in discretization with threshold
# encoding (k bins → k-1 features), which is more efficient than the one-hot
# encoding (k bins → k features) implemented here.
#
# REPLACEMENT: Pass discretize_method and discretize_bins directly to
# optimaltrees::fit_tree() or optimaltrees::cv_regularization()
#
# BACKGROUND: This manual discretization was needed before treefarmr supported
# continuous features. Now that treefarmr has built-in discretization with
# threshold encoding, this function is no longer needed and wastes features
# (33% more features than necessary).
#
# STATUS: Kept for backward compatibility testing only. Do not use in new code.
# ============================================================================

#' Adaptive Discretization for Continuous Features
#'
#' Discretizes continuous features into bins for tree-based methods.
#' Uses adaptive binning: b_n = max(2, ceiling(log(n)/3)) as suggested by
#' nonparametric theory for optimal bias-variance tradeoff.
#'
#' @param X Data.frame with continuous features
#' @param n_bins Integer or "adaptive". If "adaptive", uses b_n = max(2, ceiling(log(n)/3)).
#'   If integer, uses fixed number of bins.
#' @param method Character. "quantiles" (default) uses quantile-based thresholds.
#' @param breaks_list Optional list of pre-computed breaks (for test set discretization).
#'   If provided, uses these breaks instead of computing new ones.
#' @return List with:
#'   - X_discrete: Data.frame with discretized features (integer 0, 1, ..., b_n-1)
#'   - breaks_list: List of breaks used for each feature (for applying to test data)
#'   - n_bins: Actual number of bins used
#' @examples
#' X <- data.frame(x1 = runif(100), x2 = runif(100))
#' disc <- discretize_adaptive(X, n_bins = "adaptive")
#' # Apply same breaks to test data:
#' X_test <- data.frame(x1 = runif(10), x2 = runif(10))
#' disc_test <- discretize_adaptive(X_test, breaks_list = disc$breaks_list)
#' @noRd
discretize_adaptive <- function(X, n_bins = "adaptive", method = "quantiles", breaks_list = NULL) {

  .Deprecated(msg = paste(
    "discretize_adaptive() is deprecated.",
    "Use treefarmr's built-in discretization (pass discretize_method and discretize_bins",
    "to optimaltrees::fit_tree()) for more efficient threshold encoding."
  ))

  n <- nrow(X)
  p <- ncol(X)

  # Determine number of bins
  if (is.character(n_bins) && n_bins == "adaptive") {
    # Theory: b_n = max(2, ceiling(log(n) / 3))
    # At n=800: b_n = max(2, ceiling(6.68/3)) = max(2, 3) = 3
    n_bins_actual <- max(2, ceiling(log(n) / 3))
  } else if (is.numeric(n_bins) && n_bins >= 2) {
    n_bins_actual <- as.integer(n_bins)
  } else {
    stop("n_bins must be 'adaptive' or an integer >= 2")
  }

  # If breaks provided (test set), use them
  if (!is.null(breaks_list)) {
    if (length(breaks_list) != p) {
      stop("breaks_list must have same length as ncol(X)")
    }

    X_discrete <- X
    for (j in 1:p) {
      breaks <- breaks_list[[j]]
      if (is.null(breaks)) {
        # Feature was constant in training set
        X_discrete[[j]] <- 0L
      } else {
        # Cut using provided breaks
        X_discrete[[j]] <- as.integer(cut(X[[j]], breaks = breaks,
                                          labels = FALSE, include.lowest = TRUE)) - 1L
        # Handle values outside training range
        X_discrete[[j]][X[[j]] < min(breaks)] <- 0L
        X_discrete[[j]][X[[j]] > max(breaks)] <- length(breaks) - 2L
      }
    }

    # Convert to binary indicators (same as training data)
    X_binary <- data.frame(matrix(0L, nrow = nrow(X_discrete), ncol = 0))

    for (j in seq_along(X_discrete)) {
      col_name <- names(X_discrete)[j]
      unique_vals <- sort(unique(X_discrete[[j]]))

      # Create binary indicator for each unique value
      for (val in unique_vals) {
        new_col_name <- paste0(col_name, "_bin", val)
        X_binary[[new_col_name]] <- as.integer(X_discrete[[j]] == val)
      }
    }

    return(list(
      X_discrete = X_binary,  # Return binary indicators
      breaks_list = breaks_list,
      n_bins = length(breaks_list[[1]]) - 1
    ))
  }

  # Otherwise, compute breaks from training data
  breaks_list <- vector("list", p)
  names(breaks_list) <- names(X)

  X_discrete <- X

  for (j in 1:p) {
    x <- X[[j]]

    # Check if feature is constant
    if (length(unique(x)) == 1) {
      X_discrete[[j]] <- rep(0L, n)
      breaks_list[[j]] <- NULL
      next
    }

    # Check if feature is already binary (0/1)
    unique_vals <- unique(x)
    if (length(unique_vals) == 2 && all(sort(unique_vals) == c(0, 1))) {
      X_discrete[[j]] <- as.integer(x)
      breaks_list[[j]] <- c(-Inf, 0.5, Inf)
      next
    }

    # Discretize continuous feature
    if (method == "quantiles") {
      # Use quantiles to define breaks
      probs <- seq(0, 1, length.out = n_bins_actual + 1)
      breaks <- quantile(x, probs = probs, na.rm = TRUE)

      # Ensure breaks are unique (can happen with discrete or skewed data)
      breaks <- unique(breaks)

      # Ensure first and last breaks are -Inf and Inf for robust binning
      breaks[1] <- -Inf
      breaks[length(breaks)] <- Inf

    } else {
      stop("Only method='quantiles' is currently supported")
    }

    # Cut into bins, labeled 0, 1, ..., k-1
    X_discrete[[j]] <- as.integer(cut(x, breaks = breaks,
                                      labels = FALSE, include.lowest = TRUE)) - 1L
    breaks_list[[j]] <- breaks
  }

  # Ensure all columns are integer
  X_discrete <- as.data.frame(lapply(X_discrete, as.integer))

  # Convert to binary indicators (required by optimaltrees)
  # Each integer feature (0, 1, 2, ..., k-1) becomes k binary features
  X_binary <- data.frame(matrix(0L, nrow = nrow(X_discrete), ncol = 0))

  for (j in seq_along(X_discrete)) {
    col_name <- names(X_discrete)[j]
    unique_vals <- sort(unique(X_discrete[[j]]))

    # Create binary indicator for each unique value
    for (val in unique_vals) {
      new_col_name <- paste0(col_name, "_bin", val)
      X_binary[[new_col_name]] <- as.integer(X_discrete[[j]] == val)
    }
  }

  list(
    X_discrete = X_binary,  # Return binary indicators, not integers
    breaks_list = breaks_list,
    n_bins = n_bins_actual
  )
}
