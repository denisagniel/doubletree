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

#' Honest (bias-aware) critical value for a bias-to-SE ratio
#'
#' @description
#' Armstrong & Kolesar (2018) fixed-length critical value for a Wald interval when
#' the estimator carries a known bias bound. For a point estimate with standard
#' error \code{se} and a bias whose magnitude is bounded by \code{B}, the interval
#' \eqn{\hat\theta \pm cv(b)\cdot se} with \eqn{b = B/se} has coverage \code{level}
#' for the true parameter, where \eqn{cv(b)} is the \code{level} quantile of the
#' folded normal \eqn{|N(b, 1)|}. It solves
#' \eqn{P(-cv \le N(b,1) \le cv) = \Phi(cv - b) - \Phi(-cv - b) = level}.
#'
#' At \eqn{b = 0} this returns the ordinary \eqn{z_{(1+level)/2}} (e.g. 1.96 at 95\%);
#' it is increasing in \eqn{b} (more bias -> wider interval).
#'
#' @param b Non-negative bias-to-SE ratio \eqn{B / se}.
#' @param level Confidence level (default 0.95).
#' @return Scalar critical value \eqn{cv(b) \ge z_{(1+level)/2}}.
#' @references Armstrong, T. B. & Kolesar, M. (2018). Optimal inference in a class of
#'   regression models. \emph{Econometrica}, 86(2), 655-683.
#' @seealso \code{\link{honest_ci}}
#' @export
honest_cv <- function(b, level = 0.95) {
  if (!is.numeric(b) || length(b) != 1 || !is.finite(b) || b < 0) {
    stop("`b` must be a single finite non-negative number, got: ", b, call. = FALSE)
  }
  z0 <- qnorm(0.5 + level / 2)          # cv at b = 0 (e.g. 1.96)
  if (b == 0) return(z0)
  # coverage(cv) = Phi(cv - b) - Phi(-cv - b); increasing in cv, from <level at cv=z0
  # to 1 as cv -> Inf. Root is in (z0, b + z0]: at cv = b + z0, coverage >= Phi(z0) -
  # Phi(-2b - z0) >= level. uniroot on that bracket.
  f <- function(cv) pnorm(cv - b) - pnorm(-cv - b) - level
  uniroot(f, lower = z0, upper = b + z0, tol = .Machine$double.eps^0.5)$root
}

#' Honest (bias-aware) confidence interval centered at a biased estimate
#'
#' @description
#' Builds a confidence interval for the true parameter around a (possibly biased)
#' point estimate \code{theta_display}, using the standard error \code{se} of a valid
#' companion estimator and a conservative bias bound \eqn{B = |\delta| + z\cdot se_{\delta}}.
#' The half-width is \eqn{cv(B/se)\cdot se} (see \code{\link{honest_cv}}). Widening the
#' bias bound by \eqn{z\cdot se_\delta} guards against the sampling noise in the
#' plug-in bias estimate \code{delta}, giving a genuine (not merely nominal) coverage
#' guarantee for \eqn{\theta_0} when \eqn{B} dominates the true bias with high
#' probability.
#'
#' @param theta_display Point estimate to center the interval on (the biased display
#'   estimate, e.g. an averaged single tree).
#' @param se Standard error of the valid companion estimator (the cross-fit twin).
#' @param delta Plug-in bias estimate \eqn{\hat\theta_{display} - \hat\theta_{cf}}.
#' @param se_delta Standard error of \code{delta} (0 collapses to the raw \eqn{|\delta|} bound).
#' @param level Confidence level (default 0.95).
#' @return List with \code{ci} (length-2 numeric), the bias bound \code{B}, the
#'   critical value \code{cv}, and \code{half_width}.
#' @seealso \code{\link{honest_cv}}
#' @export
honest_ci <- function(theta_display, se, delta, se_delta = 0, level = 0.95) {
  if (!is.numeric(se) || length(se) != 1 || !is.finite(se) || se <= 0) {
    stop("`se` must be a single finite positive number, got: ", se, call. = FALSE)
  }
  z <- qnorm(0.5 + level / 2)
  se_delta <- if (is.finite(se_delta)) max(se_delta, 0) else 0
  B  <- abs(delta) + z * se_delta          # conservative bias bound
  cv <- honest_cv(B / se, level)
  half <- cv * se
  list(ci = theta_display + c(-1, 1) * half, B = B, cv = cv, half_width = half)
}

#' Solve the ATT EIF for point estimate, SE, and 95\% CI
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
