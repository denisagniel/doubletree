#' ATT orthogonal score (psi)
#'
#' Neyman-orthogonal score for the Average Treatment effect on the Treated (ATT).
#' For observation i: psi_i = (A_i/pi)*(Y_i - m0_i - theta) - (1/pi)*((A_i - e_i)/(1 - e_i))*(m1_i - m0_i).
#' See manuscript equation (2).
#'
#' @param Y Numeric vector of length n (outcomes).
#' @param A Numeric vector of length n (treatment 0/1).
#' @param theta Scalar (candidate ATT value).
#' @param eta List with elements \code{e}, \code{m0}, \code{m1}, each numeric vector of length n
#'   (propensity P(A=1|X), E[Y|A=0,X], E[Y|A=1,X]).
#' @param pi_hat Scalar estimate of P(A=1) (e.g. mean(A)).
#' @param e_min Lower bound for clamping propensity (avoid division by zero). Default 1e-6.
#' @param e_max Upper bound for clamping propensity. Default 1 - 1e-6.
#' @return Numeric vector of length n (score values).
#' @export
psi_att <- function(Y, A, theta, eta, pi_hat, e_min = 1e-6, e_max = 1 - 1e-6) {
  e <- eta$e
  m0 <- eta$m0
  m1 <- eta$m1
  n <- length(Y)
  if (length(A) != n || length(e) != n || length(m0) != n || length(m1) != n) {
    stop("Y, A, e, m0, m1 must have the same length")
  }
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("pi_hat must be in (0, 1)")
  }
  e_clamp <- pmin(pmax(e, e_min), e_max)
  term1 <- (A / pi_hat) * (Y - m0 - theta)
  term2 <- (1 / pi_hat) * ((A - e_clamp) / (1 - e_clamp)) * (m1 - m0)
  term2[!is.finite(term2)] <- 0
  term1 - term2
}
