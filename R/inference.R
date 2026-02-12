#' Variance of the ATT DML estimator
#'
#' Estimates sigma^2 as the empirical variance of the influence (score) values,
#' i.e. (1/n) * sum_i psi_i^2.
#'
#' @param score_values Numeric vector of length n (score at theta_hat).
#' @param n Optional. If provided, used as denominator; otherwise length(score_values).
#' @return Scalar (estimated variance of the influence function).
#' @export
dml_att_variance <- function(score_values, n = NULL) {
  if (is.null(n)) n <- length(score_values)
  mean(score_values^2)
}

#' Wald confidence interval for ATT
#'
#' @param theta Point estimate (scalar).
#' @param sigma Numeric. Estimated standard deviation (e.g. sqrt(dml_att_variance(...))).
#' @param n Sample size (for sqrt(n) scaling).
#' @param level Numeric. Confidence level (default 0.95).
#' @return Numeric vector of length 2 (lower, upper).
#' @export
dml_att_ci <- function(theta, sigma, n, level = 0.95) {
  z <- qnorm(0.5 + level / 2)
  half <- z * sigma / sqrt(n)
  c(theta - half, theta + half)
}
