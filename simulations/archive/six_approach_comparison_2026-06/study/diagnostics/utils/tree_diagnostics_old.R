# Tree Diagnostics Utilities
# Helper functions for analyzing tree structure and quality
# Created: 2026-05-27

#' Count number of leaves in a tree
#'
#' @param tree optimaltree object
#' @return Integer, number of terminal nodes
#' @export
count_leaves <- function(tree) {
  if (inherits(tree, "OptimalTreesModel")) {
    # Extract tree structure
    structure <- tree@tree_structure
    if (is.null(structure) || nrow(structure) == 0) {
      return(1L)  # Root only
    }
    # Count nodes without children (leaves)
    n_leaves <- sum(structure$left_is_leaf & structure$right_is_leaf) + 1
    return(n_leaves)
  } else {
    stop("Object is not an optimaltree")
  }
}

#' Compute maximum depth of tree
#'
#' @param tree optimaltree object
#' @return Integer, maximum depth (root = 0)
#' @export
max_depth <- function(tree) {
  if (inherits(tree, "OptimalTreesModel")) {
    structure <- tree@tree_structure
    if (is.null(structure) || nrow(structure) == 0) {
      return(0L)  # Root only
    }
    # Depth is length of longest path
    max_path_length <- max(sapply(structure$path, function(p) {
      if (p == "") return(0L)
      return(nchar(gsub("[^LR]", "", p)))
    }))
    return(max_path_length)
  } else {
    stop("Object is not an optimaltree")
  }
}

#' Extract tree structure as data frame
#'
#' @param tree optimaltree object
#' @return Data frame with split information
#' @export
extract_tree_structure <- function(tree) {
  if (inherits(tree, "OptimalTreesModel")) {
    structure <- tree@tree_structure
    if (is.null(structure) || nrow(structure) == 0) {
      return(data.frame(
        path = "",
        feature = NA_character_,
        threshold = NA_real_,
        depth = 0L,
        stringsAsFactors = FALSE
      ))
    }

    # Add depth column
    structure$depth <- sapply(structure$path, function(p) {
      if (p == "") return(0L)
      return(nchar(gsub("[^LR]", "", p)))
    })

    return(as.data.frame(structure))
  } else {
    stop("Object is not an optimaltree")
  }
}

#' Compute prediction accuracy metrics
#'
#' @param pred Predicted values (probabilities for classification)
#' @param true True values (binary for classification, continuous for regression)
#' @param type "classification" or "regression"
#' @return List of metrics
#' @export
compute_prediction_metrics <- function(pred, true, type = "classification") {
  if (type == "classification") {
    # Binary classification metrics
    # Assume pred is probability, true is binary

    # Log loss (negative log-likelihood)
    eps <- 1e-15  # Avoid log(0)
    pred_clipped <- pmax(eps, pmin(1 - eps, pred))
    log_loss <- -mean(true * log(pred_clipped) + (1 - true) * log(1 - pred_clipped))

    # Brier score (MSE of probabilities)
    brier_score <- mean((pred - true)^2)

    # Calibration: bin predictions and compute observed rate
    n_bins <- 10
    bins <- cut(pred, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
    calibration_df <- data.frame(
      bin = bins,
      predicted = pred,
      observed = true
    )

    calibration_summary <- aggregate(
      cbind(predicted, observed) ~ bin,
      data = calibration_df,
      FUN = function(x) c(mean = mean(x), n = length(x))
    )

    # Calibration slope (should be 1.0 if well-calibrated)
    if (nrow(calibration_summary) > 1) {
      calib_model <- lm(
        calibration_summary$observed[, "mean"] ~ calibration_summary$predicted[, "mean"]
      )
      calibration_slope <- coef(calib_model)[2]
      calibration_intercept <- coef(calib_model)[1]
    } else {
      calibration_slope <- NA_real_
      calibration_intercept <- NA_real_
    }

    return(list(
      log_loss = log_loss,
      brier_score = brier_score,
      calibration_slope = calibration_slope,
      calibration_intercept = calibration_intercept,
      calibration_data = calibration_summary
    ))

  } else if (type == "regression") {
    # Regression metrics
    bias <- mean(pred - true)
    rmse <- sqrt(mean((pred - true)^2))
    mae <- mean(abs(pred - true))
    max_error <- max(abs(pred - true))

    # R-squared
    ss_res <- sum((true - pred)^2)
    ss_tot <- sum((true - mean(true))^2)
    r_squared <- 1 - ss_res / ss_tot

    return(list(
      bias = bias,
      rmse = rmse,
      mae = mae,
      max_error = max_error,
      r_squared = r_squared
    ))
  } else {
    stop("type must be 'classification' or 'regression'")
  }
}

#' Analyze overlap quality from propensity scores
#'
#' @param e_hat Estimated propensity scores
#' @param A Treatment indicator
#' @return List of overlap diagnostics
#' @export
analyze_overlap <- function(e_hat, A) {
  # Overall distribution
  min_e <- min(e_hat)
  max_e <- max(e_hat)
  mean_e <- mean(e_hat)
  median_e <- median(e_hat)

  # Extreme values (poor overlap)
  extreme_low <- sum(e_hat < 0.05) / length(e_hat)
  extreme_high <- sum(e_hat > 0.95) / length(e_hat)
  extreme_total <- sum(e_hat < 0.05 | e_hat > 0.95) / length(e_hat)

  # Distribution by treatment group
  e_treated <- e_hat[A == 1]
  e_control <- e_hat[A == 0]

  # Effective sample size (Kish's ESS for weights)
  weights <- A / e_hat + (1 - A) / (1 - e_hat)
  ess <- sum(weights)^2 / sum(weights^2)
  ess_ratio <- ess / length(A)

  # Positivity violations
  positivity_violations <- sum(e_hat <= 0 | e_hat >= 1)

  return(list(
    min_e = min_e,
    max_e = max_e,
    mean_e = mean_e,
    median_e = median_e,
    extreme_low = extreme_low,
    extreme_high = extreme_high,
    extreme_total = extreme_total,
    ess = ess,
    ess_ratio = ess_ratio,
    positivity_violations = positivity_violations,
    n_treated = sum(A == 1),
    n_control = sum(A == 0),
    mean_e_treated = mean(e_treated),
    mean_e_control = mean(e_control)
  ))
}

#' Compare tree structures across folds
#'
#' @param tree_list List of optimaltree objects (one per fold)
#' @return List with structure comparison metrics
#' @export
compare_tree_structures <- function(tree_list) {
  n_trees <- length(tree_list)

  # Extract key features
  n_leaves <- sapply(tree_list, count_leaves)
  max_depths <- sapply(tree_list, max_depth)

  # Extract split features used
  features_used <- lapply(tree_list, function(tree) {
    structure <- extract_tree_structure(tree)
    unique(structure$feature[!is.na(structure$feature)])
  })

  # Common features across all folds
  common_features <- Reduce(intersect, features_used)

  # Feature usage frequency
  all_features <- unlist(features_used)
  feature_counts <- table(all_features)

  return(list(
    n_trees = n_trees,
    n_leaves = n_leaves,
    n_leaves_range = range(n_leaves),
    n_leaves_mean = mean(n_leaves),
    n_leaves_sd = sd(n_leaves),
    max_depths = max_depths,
    max_depth_range = range(max_depths),
    max_depth_mean = mean(max_depths),
    common_features = common_features,
    n_common_features = length(common_features),
    feature_counts = feature_counts,
    features_used = features_used
  ))
}

#' Compute oracle tree performance (fit tree to true function)
#'
#' @param X Covariate matrix
#' @param y_true True outcome (e.g., true propensity or true Y0)
#' @param loss Loss function to use
#' @param regularization Regularization parameter
#' @return List with oracle tree and metrics
#' @export
compute_oracle_performance <- function(X, y_true, loss = "squared_error",
                                       regularization = NULL) {
  # Fit tree to true function
  if (is.null(regularization)) {
    n <- nrow(X)
    regularization <- log(n) / n
  }

  oracle_tree <- optimaltrees::optimaltree(
    X = X,
    y = y_true,
    loss = loss,
    regularization = regularization
  )

  # Predict and compute error
  y_hat <- predict(oracle_tree, X)

  # Metrics depend on loss type
  if (loss %in% c("log_loss", "misclassification")) {
    metrics <- compute_prediction_metrics(y_hat, y_true, type = "classification")
  } else {
    metrics <- compute_prediction_metrics(y_hat, y_true, type = "regression")
  }

  # Tree complexity
  metrics$n_leaves <- count_leaves(oracle_tree)
  metrics$max_depth <- max_depth(oracle_tree)

  return(list(
    tree = oracle_tree,
    metrics = metrics
  ))
}
