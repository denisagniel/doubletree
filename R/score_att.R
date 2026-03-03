#' ATT orthogonal score (psi)
#'
#' Neyman-orthogonal score for the Average Treatment effect on the Treated (ATT).
#' For observation i: psi_i = (A_i/pi)*(Y_i - m0_i - theta) - (1/pi)*(e_i*(1 - A_i)/(1 - e_i))*(Y_i - m0_i).
#' See manuscript equation (2).
#'
#' @param Y Numeric vector of length n (outcomes).
#' @param A Numeric vector of length n (treatment 0/1).
#' @param theta Scalar (candidate ATT value).
#' @param eta List with elements \code{e}, \code{m0}, each numeric vector of length n
#'   (propensity P(A=1|X), E[Y|A=0,X]). Propensities must already be
#'   clamped to avoid division by zero (done at prediction time in nuisance fitting).
#' @param pi_hat Scalar estimate of P(A=1) (e.g. mean(A)).
#' @param e_min,e_max Deprecated. Clamping is now done at prediction time. These
#'   parameters are kept for backward compatibility but have no effect.
#' @return Numeric vector of length n (score values).
#' @export
psi_att <- function(Y, A, theta, eta, pi_hat, e_min = 1e-6, e_max = 1 - 1e-6) {
  e <- eta$e  # Already clamped at prediction time
  m0 <- eta$m0
  n <- length(Y)
  if (length(A) != n || length(e) != n || length(m0) != n) {
    stop("Y, A, e, m0 must have the same length")
  }
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("pi_hat must be in (0, 1)")
  }
  # No clamping needed - already done at prediction time
  # Correct Chernozhukov (2018) ATT score
  term1 <- (A / pi_hat) * (Y - m0 - theta)
  term2 <- (1 / pi_hat) * (e * (1 - A) / (1 - e)) * (Y - m0)
  term2[!is.finite(term2)] <- 0  # Safety for any remaining edge cases
  term1 - term2
}
