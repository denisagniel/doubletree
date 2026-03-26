# Propensity score bounds for numerical stability
# Following standard practice in causal inference (Crump et al. 2009, Li et al. 2018),
# clip to [0.01, 0.99] to prevent extreme inverse probability weights in the EIF.
# With these bounds, e/(1-e) is bounded by [0.01, 99], preventing numerical instability.
.PROPENSITY_LOWER_BOUND <- 0.01
.PROPENSITY_UPPER_BOUND <- 0.99

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
#' @param e_min,e_max DEPRECATED. These parameters are no longer used.
#'   Propensity clamping is now done at prediction time. Providing these
#'   arguments will trigger a warning.
#' @return Numeric vector of length n (score values).
#' @export
psi_att <- function(Y, A, theta, eta, pi_hat, e_min = NULL, e_max = NULL) {
  # Warn if deprecated parameters are provided
  if (!is.null(e_min) || !is.null(e_max)) {
    warning("e_min and e_max are deprecated and ignored. ",
            "Propensity clamping is now done at prediction time in nuisance fitting. ",
            "These parameters will be removed in a future version.",
            call. = FALSE)
  }

  e <- eta$e  # Already clamped at prediction time
  m0 <- eta$m0
  n <- length(Y)
  if (length(A) != n || length(e) != n || length(m0) != n) {
    stop("Y, A, e, m0 must have the same length")
  }
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("pi_hat must be in (0, 1)")
  }

  # Validate propensity scores are within acceptable bounds
  # With clipping to [0.01, 0.99], e/(1-e) is bounded by [0.01, 99]
  # Use strict inequalities since clamping may produce values exactly at bounds
  if (any(e > 0.99) || any(e < 0.01)) {
    stop("Propensity scores outside acceptable bounds [0.01, 0.99]. ",
         "This indicates numerical instability in the propensity model. ",
         "Check that propensity bounds are enforced at prediction time.",
         call. = FALSE)
  }

  # No clamping needed - already done at prediction time
  # Correct Chernozhukov (2018) ATT score
  term1 <- (A / pi_hat) * (Y - m0 - theta)
  term2 <- (1 / pi_hat) * (e * (1 - A) / (1 - e)) * (Y - m0)

  # If term2 has non-finite values after validation, that's a real error
  if (any(!is.finite(term2))) {
    stop("Non-finite values in score computation. This should not happen after ",
         "propensity validation. Please report this as a bug.",
         call. = FALSE)
  }

  term1 - term2
}
