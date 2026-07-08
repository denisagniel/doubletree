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
#' @param level Numeric. Confidence level (default 0.95).
#' @return Numeric vector of length 2 (lower, upper).
#' @export
att_ci <- function(theta, sigma, level = 0.95) {
  z <- qnorm(0.5 + level / 2)
  half <- z * sigma
  c(theta - half, theta + half)
}

#' Solve the ATT EIF for point estimate, SE, and 95% CI
#'
#' @description
#' Given plugged-in nuisance predictions \code{e_hat = e(X)} and
#' \code{m0_hat = m0(X)}, solve the efficient-influence-function estimating
#' equation for the ATT point estimate and return its standard error and Wald
#' interval. This is the single shared solver used by every \code{estimate_att*}
#' entry point (previously copy-pasted at each site, with minor drift between
#' \code{att_se()} and inline centered-variance forms).
#'
#' The estimate uses the closed form \eqn{\hat\theta = \sum_i \psi_i(0) / \sum_i (A_i/\hat\pi)},
#' where \eqn{\psi} is \code{\link{psi_att}}; because \eqn{\sum_i \psi_i(\hat\theta) = 0}
#' by construction, the centered and uncentered score variances coincide, so
#' \code{\link{att_se}} (uncentered, with an estimating-equation check) is used.
#'
#' @param Y Numeric outcome vector.
#' @param A Binary treatment vector (0/1).
#' @param e_hat Numeric vector of propensity predictions \eqn{e(X)}.
#' @param m0_hat Numeric vector of control-outcome predictions \eqn{m_0(X)}.
#' @param n Optional integer denominator for the SE (default \code{length(Y)}).
#' @return A list with \code{theta}, \code{sigma}, \code{ci_95},
#'   \code{score_values} (the score at \eqn{\hat\theta}), and \code{pi_hat}.
#' @seealso \code{\link{att_se}}, \code{\link{att_ci}}, \code{\link{psi_att}}
#' @keywords internal
eif_att_solve <- function(Y, A, e_hat, m0_hat, n = length(Y)) {
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  theta <- sum(score_at_zero) / sum(A / pi_hat)
  score_values <- psi_att(Y, A, theta, eta, pi_hat)
  sigma <- att_se(score_values, n)
  list(
    theta = theta,
    sigma = sigma,
    ci_95 = att_ci(theta, sigma, level = 0.95),
    score_values = score_values,
    pi_hat = pi_hat
  )
}
