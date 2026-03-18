#' Standard error of the ATT estimator
#'
#' Computes the standard error from the influence function (score) values
#' as sqrt((1/n) * sum_i psi_i^2).
#'
#' @param score_values Numeric vector of length n (score at theta_hat).
#' @param n Optional. If provided, used as denominator; otherwise length(score_values).
#' @return Scalar (estimated standard error of the ATT estimator).
#' @export
att_se <- function(score_values, n = NULL) {
  if (is.null(n)) n <- length(score_values)
  sqrt(mean(score_values^2))
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
