# Metrics Functions for Six-Approach Comparison
# Created: 2026-05-01
#
# Compute metrics for each replication:
# - Bias, RMSE (computed after combining results)
# - Coverage (computed per rep)
# - Structure similarity (requires pairwise comparison - do offline)

#' Compute coverage for a single replication
#'
#' @param theta_hat Point estimate
#' @param se Standard error
#' @param theta_true True parameter value
#' @param alpha Significance level (default 0.05)
#' @return Named vector with coverage indicators
compute_coverage <- function(theta_hat, se, theta_true, alpha = 0.05) {
  z <- qnorm(1 - alpha/2)
  ci_lower <- theta_hat - z * se
  ci_upper <- theta_hat + z * se

  covers <- (theta_true >= ci_lower) & (theta_true <= ci_upper)
  covers_lower <- (theta_true >= ci_lower)
  covers_upper <- (theta_true <= ci_upper)

  c(
    coverage = covers,
    coverage_lower = covers_lower,
    coverage_upper = covers_upper,
    ci_width = 2 * z * se,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

#' Compute bias-adjusted coverage
#'
#' For approaches with bias, test alternative CI constructions
#'
#' @param theta_hat Point estimate (possibly biased)
#' @param se Standard error
#' @param bias_estimate Estimated bias (e.g., theta_hat - theta_crossfit)
#' @param theta_true True parameter
#' @param method CI adjustment method: "standard", "additive", "conservative"
#' @return Named vector with coverage and CI width
compute_bias_adjusted_coverage <- function(theta_hat, se, bias_estimate,
                                          theta_true, method = "additive") {
  z <- qnorm(0.975)

  if (method == "standard") {
    # Naive CI (no correction)
    ci_lower <- theta_hat - z * se
    ci_upper <- theta_hat + z * se

  } else if (method == "additive") {
    # Add |bias| to margin
    ci_lower <- theta_hat - (z * se + abs(bias_estimate))
    ci_upper <- theta_hat + (z * se + abs(bias_estimate))

  } else if (method == "conservative") {
    # Add 2*|bias| to margin
    ci_lower <- theta_hat - (z * se + 2 * abs(bias_estimate))
    ci_upper <- theta_hat + (z * se + 2 * abs(bias_estimate))

  } else {
    stop("Unknown method: ", method)
  }

  covers <- (theta_true >= ci_lower) & (theta_true <= ci_upper)
  ci_width <- ci_upper - ci_lower

  c(
    coverage = covers,
    ci_width = ci_width,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

#' Extract tree structure summary
#'
#' Simple structure representation for comparison
#'
#' @param tree_structure Tree structure object
#' @return Character hash of structure
structure_hash <- function(tree_structure) {
  if (is.null(tree_structure)) return(NA_character_)

  # Simple hash: convert structure to string and hash
  # This is a placeholder - actual implementation depends on structure format
  structure_str <- paste(
    tree_structure$splits,
    tree_structure$features,
    collapse = "_"
  )
  digest::digest(structure_str, algo = "md5")
}

#' Compare two tree structures
#'
#' Check if two structures are identical
#'
#' @param struct1 First tree structure
#' @param struct2 Second tree structure
#' @return Logical: are structures identical?
structures_match <- function(struct1, struct2) {
  hash1 <- structure_hash(struct1)
  hash2 <- structure_hash(struct2)

  if (is.na(hash1) || is.na(hash2)) return(NA)

  hash1 == hash2
}

#' Compute leaf-level RMSE between predictions
#'
#' @param pred1 Predictions from first approach
#' @param pred2 Predictions from second approach
#' @return RMSE value
prediction_rmse <- function(pred1, pred2) {
  sqrt(mean((pred1 - pred2)^2, na.rm = TRUE))
}
