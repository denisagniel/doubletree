#' Standard error of the ATT estimator
#'
#' Computes the standard error from the influence function (score) values
#' as sqrt((1/n) * sum_i psi_i^2). Also validates that the estimating equation
#' is satisfied (E[ψ(θ̂)] ≈ 0).
#'
#' @param score_values Numeric vector of length n (score at theta_hat).
#' @param n Optional. If provided, used as denominator; otherwise length(score_values).
#' @return Scalar (estimated standard error of the ATT estimator).
#' @export
att_se <- function(score_values, n = NULL) {
  if (is.null(n)) n <- length(score_values)

  # Validate that estimating equation is satisfied: E[ψ(θ̂)] ≈ 0
  # If mean is far from zero, θ̂ doesn't solve the estimating equation
  score_mean <- mean(score_values)
  score_sd <- sd(score_values)

  if (abs(score_mean) > 0.01 * score_sd) {
    warning("Score function mean (", round(score_mean, 4),
            ") is not close to zero (> 1% of SD). ",
            "This suggests the estimating equation may not be satisfied. ",
            "ATT estimate may be unreliable. ",
            "Check propensity and outcome models for misspecification.",
            call. = FALSE)
  }

  # Standard error: SE(θ̂) = sqrt(E[ψ²] / n) when E[ψ] ≈ 0
  sqrt(mean(score_values^2) / n)
}

#' Wald confidence interval for ATT
#'
#' @param theta Point estimate (scalar).
#' @param sigma Numeric. Estimated standard error (from att_se(...)).
#' @param n Sample size (for sqrt(n) scaling).
#' @param level Numeric. Confidence level (default 0.95).
#' @return Numeric vector of length 2 (lower, upper).
#' @export
att_ci <- function(theta, sigma, n, level = 0.95) {
  z <- qnorm(0.5 + level / 2)
  half <- z * sigma / sqrt(n)
  c(theta - half, theta + half)
}
