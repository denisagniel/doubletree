#' ATT Estimation with Repeated Cross-Fitting
#'
#' Implements repeated sample splitting (Chernozhukov et al. 2018) to account for
#' fold-assignment randomness. Runs M independent cross-fits and combines:
#' - Point estimate: median or mean of M estimates
#' - Variance: accounts for both within-fold (influence function) and between-split variation
#'
#' @inheritParams estimate_att
#' @param n_splits Integer. Number of independent cross-fit repetitions. Default 1 (no repetition).
#' @param aggregation Character. How to combine M point estimates:
#'   \itemize{
#'     \item \code{"median"} (default): Robust to outlier splits. Recommended when
#'       K is small (< 5), data is noisy, or model fitting is unstable.
#'     \item \code{"mean"}: Efficient under normality. Use when splits are well-behaved
#'       and you want to minimize variance. Theoretically aligned with Chernozhukov et al. (2018).
#'   }
#'   If unsure, use median (more robust). Mean is more efficient but sensitive to outliers.
#'
#' @return List with elements:
#'   - theta: aggregated point estimate
#'   - sigma: corrected standard error on the theta-hat scale (accounts for fold
#'     randomness); \eqn{\sqrt{\overline{\sigma_m^2 + (\hat\theta_m - \hat\theta)^2}}}
#'   - ci_95: 95\% confidence interval, \code{theta +/- 1.96 * sigma}
#'   - theta_splits: vector of M point estimates (one per split)
#'   - sigma_splits: vector of M within-fold SEs (theta-hat scale)
#'   - between_var: mean squared deviation of split estimates (theta-hat scale)
#'   - within_var: mean within-fold variance (theta-hat scale)
#'   - ... (other estimate_att outputs)
#'
#' @references Chernozhukov et al. (2018), "Double/debiased machine learning for treatment
#'   and structural parameters", Econometrics Journal.
#'
#' @examples
#' \donttest{
#' library(optimaltrees)
#' set.seed(42)
#' n <- 150
#' X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, plogis(0.5 * X$X1))
#' Y <- rbinom(n, 1, 0.4 + 0.2 * A + 0.1 * X$X1)
#'
#' # Standard estimation (single split)
#' fit_single <- att_repeated(X, A, Y, K = 5, n_splits = 1)
#' print(fit_single$theta)
#'
#' # Repeated cross-fitting (M=5 splits)
#' fit_repeated <- att_repeated(X, A, Y, K = 5, n_splits = 5)
#' print(fit_repeated$theta)  # Median of 5 estimates
#' print(fit_repeated$sigma)  # Accounts for fold randomness
#' print(fit_repeated$theta_splits)  # See variation across splits
#' }
#'
#' @export
att_repeated <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                              regularization = 0.1, cv_regularization = TRUE, cv_K = 5,
                              stratified = TRUE, seed = NULL,
                              verbose = FALSE, use_rashomon = FALSE,
                              rashomon_bound_multiplier = NULL,
                              rashomon_bound_adder = 0,
                              max_leaves = NULL,
                              auto_tune_intersecting = FALSE,
                              escalate_intersection = FALSE,
                              n_splits = 1,
                              aggregation = c("median", "mean"),
                              ...) {

  aggregation <- match.arg(aggregation)

  if (n_splits == 1) {
    # No repetition - call standard estimate_att
    return(estimate_att(X, A, Y, K = K, outcome_type = outcome_type,
                   regularization = regularization,
                   cv_regularization = cv_regularization, cv_K = cv_K,
                   stratified = stratified,
                   seed = seed, verbose = verbose, use_rashomon = use_rashomon,
                   rashomon_bound_multiplier = rashomon_bound_multiplier,
                   rashomon_bound_adder = rashomon_bound_adder,
                   max_leaves = max_leaves,
                   auto_tune_intersecting = auto_tune_intersecting,
                   escalate_intersection = escalate_intersection, ...))
  }

  # Run M independent cross-fits
  theta_splits <- numeric(n_splits)
  sigma_splits <- numeric(n_splits)
  results_list <- vector("list", n_splits)

  for (m in seq_len(n_splits)) {
    # Different seed for each split (if seed provided)
    split_seed <- if (!is.null(seed)) seed + m * 1000 else NULL

    result_m <- estimate_att(
      X, A, Y, K = K, outcome_type = outcome_type,
      regularization = regularization,
      cv_regularization = cv_regularization, cv_K = cv_K,
      stratified = stratified,
      seed = split_seed, verbose = verbose, use_rashomon = use_rashomon,
      rashomon_bound_multiplier = rashomon_bound_multiplier,
      rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves,
      auto_tune_intersecting = auto_tune_intersecting,
      escalate_intersection = escalate_intersection, ...
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
  # We account for both within-split variance (from the influence function) and
  # between-split variance (from fold-assignment randomness). BOTH components are
  # on the theta-hat scale:
  #   - sigma_splits[m] is the return of att_se(), which is ALREADY SE(theta_hat_m)
  #     (it includes the 1/n factor: sqrt(mean(psi^2)/n)). So sigma_splits^2 is a
  #     theta-hat-scale variance, O(1/n).
  #   - (theta_splits - theta)^2 is likewise a theta-hat-scale squared deviation.
  # The split-specific total variance is therefore their sum, with no n factor:
  #   var_m = sigma_splits[m]^2 + (theta_splits[m] - theta)^2
  # (The previous code multiplied the between term by n, assuming sigma was on the
  # sqrt(n) scale; that convention predates the att_se /n fix and made the two terms
  # differ by a factor of n. Corrected here.)

  n <- nrow(X)

  # Split-specific total variances (theta-hat scale)
  split_total_vars <- sigma_splits^2 + (theta_splits - theta)^2

  # Aggregated variance (mean or median), matching the point-estimate aggregation.
  if (aggregation == "median") {
    sigma_squared <- median(split_total_vars)
  } else {
    sigma_squared <- mean(split_total_vars)
  }

  # SE on the theta-hat scale (att_se convention). No further /sqrt(n).
  sigma_theta <- sqrt(sigma_squared)

  # Components for diagnostics (both on the theta-hat scale)
  within_var <- mean(sigma_splits^2)
  between_var <- mean((theta_splits - theta)^2)

  # Debug output
  if (verbose) {
    message("Variance estimation (Chernozhukov et al. 2018; theta-hat scale):")
    message(sprintf("  Within variance (mean sigma^2):        %.6f", within_var))
    message(sprintf("  Between variance (var(theta_splits)):  %.6f", between_var))
    message(sprintf("  Total variance:                        %.6f", sigma_squared))
    message(sprintf("  SE on theta scale:                     %.6f", sigma_theta))
  }

  # 95% CI on the theta-hat scale (sigma_theta is already SE(theta_hat)).
  ci_95 <- theta + c(-1, 1) * qnorm(0.975) * sigma_theta

  # Return result with additional diagnostics
  list(
    theta = theta,
    sigma = sigma_theta,  # SE on the theta-hat scale (att_se convention)
    ci_95 = ci_95,
    theta_splits = theta_splits,
    sigma_splits = sigma_splits,
    between_var = between_var,  # theta-hat scale
    within_var = within_var,    # theta-hat scale
    within_var_frac = within_var / sigma_squared,
    between_var_frac = between_var / sigma_squared,
    n_splits = n_splits,
    aggregation = aggregation,
    results_list = results_list,  # Full results for diagnostics
    n = n,
    K = K
  )
}
