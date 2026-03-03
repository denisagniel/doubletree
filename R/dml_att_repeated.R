#' DML-ATT with Repeated Cross-Fitting
#'
#' Implements repeated sample splitting (Chernozhukov et al. 2018) to account for
#' fold-assignment randomness. Runs M independent cross-fits and combines:
#' - Point estimate: median or mean of M estimates
#' - Variance: accounts for both within-fold (influence function) and between-split variation
#'
#' @inheritParams dml_att
#' @param n_splits Integer. Number of independent cross-fit repetitions. Default 1 (no repetition).
#' @param aggregation Character. How to aggregate M estimates: "median" (default) or "mean".
#'
#' @return List with elements:
#'   - theta: aggregated point estimate
#'   - sigma: corrected standard error (accounts for fold randomness)
#'   - ci_95: 95% confidence interval with corrected SE
#'   - theta_splits: vector of M point estimates (one per split)
#'   - sigma_splits: vector of M within-fold SEs
#'   - between_var: variance across splits
#'   - within_var: mean within-fold variance
#'   - ... (other dml_att outputs)
#'
#' @references Chernozhukov et al. (2018), "Double/debiased machine learning for treatment
#'   and structural parameters", Econometrics Journal.
#'
#' @examples
#' \donttest{
#' library(treefarmr)
#' set.seed(42)
#' n <- 150
#' X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, plogis(0.5 * X$X1))
#' Y <- rbinom(n, 1, 0.4 + 0.2 * A + 0.1 * X$X1)
#'
#' # Standard DML (single split)
#' fit_single <- dml_att_repeated(X, A, Y, K = 5, n_splits = 1)
#' print(fit_single$theta)
#'
#' # Repeated cross-fitting (M=5 splits)
#' fit_repeated <- dml_att_repeated(X, A, Y, K = 5, n_splits = 5)
#' print(fit_repeated$theta)  # Median of 5 estimates
#' print(fit_repeated$sigma)  # Accounts for fold randomness
#' print(fit_repeated$theta_splits)  # See variation across splits
#' }
#'
#' @export
dml_att_repeated <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                              regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                              stratified = TRUE, seed = NULL,
                              verbose = FALSE, use_rashomon = FALSE,
                              rashomon_bound_multiplier = 0.05,
                              rashomon_bound_adder = 0,
                              max_leaves = NULL,
                              auto_tune_intersecting = FALSE,
                              n_splits = 1,
                              aggregation = c("median", "mean"),
                              ...) {

  aggregation <- match.arg(aggregation)

  if (n_splits == 1) {
    # No repetition - call standard dml_att
    return(dml_att(X, A, Y, K = K, outcome_type = outcome_type,
                   regularization = regularization,
                   cv_regularization = cv_regularization, cv_K = cv_K,
                   stratified = stratified,
                   seed = seed, verbose = verbose, use_rashomon = use_rashomon,
                   rashomon_bound_multiplier = rashomon_bound_multiplier,
                   rashomon_bound_adder = rashomon_bound_adder,
                   max_leaves = max_leaves,
                   auto_tune_intersecting = auto_tune_intersecting, ...))
  }

  # Run M independent cross-fits
  theta_splits <- numeric(n_splits)
  sigma_splits <- numeric(n_splits)
  results_list <- vector("list", n_splits)

  for (m in seq_len(n_splits)) {
    # Different seed for each split (if seed provided)
    split_seed <- if (!is.null(seed)) seed + m * 1000 else NULL

    result_m <- dml_att(
      X, A, Y, K = K, outcome_type = outcome_type,
      regularization = regularization,
      cv_regularization = cv_regularization, cv_K = cv_K,
      stratified = stratified,
      seed = split_seed, verbose = verbose, use_rashomon = use_rashomon,
      rashomon_bound_multiplier = rashomon_bound_multiplier,
      rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves,
      auto_tune_intersecting = auto_tune_intersecting, ...
    )

    theta_splits[m] <- result_m$theta
    sigma_splits[m] <- result_m$sigma
    results_list[[m]] <- result_m
  }

  # Aggregate point estimates
  theta <- if (aggregation == "median") {
    median(theta_splits)
  } else {
    mean(theta_splits)
  }

  # Variance estimation following Chernozhukov et al. (2018), Section 3.4
  #
  # Formula (3.13) for mean aggregation:
  #   σ²_mean = (1/S) * Σ[σ²_s + (θ_s - θ_mean)²]
  #
  # where σ²_s is Var[√n·θ̃_s] and (θ_s - θ_mean)² is on the θ scale.
  # To combine them, we need to scale the second term:
  #   σ²_mean = (1/S) * Σ[σ²_s + n·(θ_s - θ_mean)²]
  #
  # This accounts for both:
  # - Within-split variance (σ²_s from influence function)
  # - Between-split variance (from randomness in fold assignment)

  n <- nrow(X)

  # Compute split-specific total variances (on sqrt(n) scale)
  split_total_vars <- sigma_splits^2 + n * (theta_splits - theta)^2

  # Aggregated variance (mean or median)
  if (aggregation == "median") {
    # Median picks the split with median combined variance
    sigma_squared <- median(split_total_vars)
  } else {
    # Mean: average of split-specific total variances
    sigma_squared <- mean(split_total_vars)
  }

  sigma_scaled <- sqrt(sigma_squared)

  # Components for diagnostics
  within_var <- mean(sigma_splits^2)
  between_var <- mean(n * (theta_splits - theta)^2)

  # Debug output
  if (verbose) {
    message("Variance estimation (Chernozhukov et al. 2018):")
    message(sprintf("  Within variance (mean σ²): %.4f", within_var))
    message(sprintf("  Between variance (n·var(θ)): %.4f", between_var))
    message(sprintf("  Total variance: %.4f", sigma_squared))
    message(sprintf("  sigma_scaled (sqrt(n) scale): %.4f", sigma_scaled))
    message(sprintf("  SE on theta scale: %.4f", sigma_scaled / sqrt(n)))
  }

  # 95% CI: convert back to θ̂ scale
  ci_95 <- theta + c(-1, 1) * qnorm(0.975) * sigma_scaled / sqrt(n)

  # Return result with additional diagnostics
  list(
    theta = theta,
    sigma = sigma_scaled,  # On sqrt(n) scale per Chernozhukov et al. (2018)
    ci_95 = ci_95,
    theta_splits = theta_splits,
    sigma_splits = sigma_splits,
    between_var = between_var,  # On sqrt(n) scale
    within_var = within_var,    # On sqrt(n) scale
    within_var_frac = within_var / sigma_squared,
    between_var_frac = between_var / sigma_squared,
    n_splits = n_splits,
    aggregation = aggregation,
    results_list = results_list,  # Full results for diagnostics
    n = n,
    K = K
  )
}
