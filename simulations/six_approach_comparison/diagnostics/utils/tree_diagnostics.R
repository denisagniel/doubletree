# Simplified Tree Diagnostics - Use optimaltrees package functions
# Created: 2026-05-27

#' Count number of leaves in a tree
#' Wrapper around optimaltrees::count_leaves_tree
count_leaves <- function(tree) {
  optimaltrees::count_leaves_tree(tree@trees[[1]])
}

#' Compute maximum depth of tree
#' Simple recursive traversal
max_depth <- function(tree_struct, depth = 0) {
  if (is.null(tree_struct) || !is.list(tree_struct)) {
    return(depth)
  }
  
  # Check if leaf
  if (!is.null(tree_struct$type) && tree_struct$type == "leaf") {
    return(depth)
  }
  
  # Recurse on children
  left_depth <- if (!is.null(tree_struct$true)) max_depth(tree_struct$true, depth + 1) else depth
  right_depth <- if (!is.null(tree_struct$false)) max_depth(tree_struct$false, depth + 1) else depth
  
  return(max(left_depth, right_depth))
}

#' Get tree structure for analysis
extract_tree_structure <- function(tree) {
  # Just return the raw tree structure
  tree@trees[[1]]
}

#' Compute prediction accuracy metrics
compute_prediction_metrics <- function(pred, true, type = "classification") {
  if (type == "classification") {
    # Binary classification metrics
    eps <- 1e-15
    pred_clipped <- pmax(eps, pmin(1 - eps, pred))
    
    # Log loss
    log_loss <- -mean(true * log(pred_clipped) + (1 - true) * log(1 - pred_clipped))
    
    # Brier score
    brier_score <- mean((pred - true)^2)
    
    # Calibration
    n_bins <- 10
    bins <- cut(pred, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
    calib_df <- data.frame(bin = bins, predicted = pred, observed = true)
    
    calib_summary <- aggregate(
      cbind(predicted, observed) ~ bin,
      data = calib_df,
      FUN = mean
    )
    
    # Calibration slope
    if (nrow(calib_summary) > 1) {
      calib_model <- lm(calib_summary$observed ~ calib_summary$predicted)
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
      calibration_data = calib_summary
    ))
    
  } else if (type == "regression") {
    # Regression metrics
    bias <- mean(pred - true)
    rmse <- sqrt(mean((pred - true)^2))
    mae <- mean(abs(pred - true))
    max_error <- max(abs(pred - true))
    
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
analyze_overlap <- function(e_hat, A) {
  min_e <- min(e_hat)
  max_e <- max(e_hat)
  mean_e <- mean(e_hat)
  median_e <- median(e_hat)
  
  extreme_low <- sum(e_hat < 0.05) / length(e_hat)
  extreme_high <- sum(e_hat > 0.95) / length(e_hat)
  extreme_total <- sum(e_hat < 0.05 | e_hat > 0.95) / length(e_hat)
  
  weights <- A / e_hat + (1 - A) / (1 - e_hat)
  ess <- sum(weights)^2 / sum(weights^2)
  ess_ratio <- ess / length(A)
  
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
    positivity_violations = positivity_violations
  ))
}

#' Compute oracle tree performance
compute_oracle_performance <- function(X, y_true, loss = "squared_error", regularization = NULL) {
  if (is.null(regularization)) {
    n <- nrow(X)
    regularization <- log(n) / n
  }
  
  oracle_tree <- optimaltrees::fit_tree(
    X = X,
    y = y_true,
    loss_function = loss,
    regularization = regularization
  )
  
  y_hat <- predict(oracle_tree, X)
  
  if (loss %in% c("log_loss", "misclassification")) {
    metrics <- compute_prediction_metrics(y_hat, y_true, type = "classification")
  } else {
    metrics <- compute_prediction_metrics(y_hat, y_true, type = "regression")
  }
  
  metrics$n_leaves <- count_leaves(oracle_tree)
  metrics$max_depth <- max_depth(oracle_tree@trees[[1]])
  
  return(list(
    tree = oracle_tree,
    metrics = metrics
  ))
}
