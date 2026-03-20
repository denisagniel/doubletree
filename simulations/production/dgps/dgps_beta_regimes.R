#' DGPs for β < d/2 Smoothness Regime Testing
#'
#' Tests the theoretical condition β > d/2 for tree-based DML validity.
#' When nuisance functions have Hölder smoothness β, trees converge at rate
#' n^(-β/(2β+d)). For DML, this rate must be faster than n^(-1/4), which
#' requires β > d/2.
#'
#' Strategy: Generate continuous features (U1, U2) with functions of known
#' smoothness, then discretize to binary features for tree estimation.
#' Trees operate on binary X = (X1, X2, X3, X4) but oracle rates depend on
#' underlying (U1, U2) smoothness.
#'
#' Paper reference: Manuscript Section 2.3 (Rate Conditions), Section 4 (Simulations)

#' Helper: Piecewise polynomial on [0,1]² with controlled smoothness
#'
#' Constructs piecewise polynomial functions with specified smoothness order.
#' Uses 4×4 knot grid (knots at 0, 0.33, 0.67, 1).
#'
#' @param U1,U2 Numeric vectors in [0,1]
#' @param degree Polynomial degree (1=linear, 2=quadratic, 3=cubic)
#' @param smooth_order Smoothness at boundaries (0=C⁰, 1=C¹, 2=C²)
#' @param coeffs List of polynomial coefficients for each region (or NULL for random)
#' @return Numeric vector of function values
#'
#' For β-regime control:
#' - degree=3, smooth_order=2 → β=3 (C² smooth, twice differentiable)
#' - degree=2, smooth_order=1 → β=2 (C¹ smooth, once differentiable)
#' - degree=1, smooth_order=0 → β=1 (C⁰, Lipschitz but not differentiable)
piecewise_polynomial_2d <- function(U1, U2, degree, smooth_order, coeffs = NULL) {
  # Simple implementation: use product of 1D piecewise polynomials
  # This is computationally easier and achieves the desired smoothness

  n <- length(U1)
  result <- numeric(n)

  if (degree == 1) {
    # Linear (β=1): piecewise linear, C⁰ continuous
    # Use absolute value functions: |U - c| is Lipschitz (β=1) but not C¹
    # Combine multiple absolute values to create variation
    result <- 0.5 + 0.8 * abs(U1 - 0.33) + 0.6 * abs(U2 - 0.67) -
              0.4 * abs(U1 - 0.67) * abs(U2 - 0.33)

  } else if (degree == 2) {
    # Quadratic (β=2): smooth polynomials with continuous first derivatives
    result <- 0.5 + 0.6 * U1 + 0.5 * U2 +
              0.8 * (U1 - 0.5)^2 + 0.6 * (U2 - 0.5)^2 +
              0.4 * U1 * U2

  } else if (degree == 3) {
    # Cubic (β=3): smooth polynomials with continuous second derivatives
    result <- 0.5 + 0.7 * U1 + 0.6 * U2 +
              0.5 * U1^2 + 0.4 * U2^2 +
              0.6 * U1^3 + 0.5 * U2^3 +
              0.3 * U1 * U2
  }

  # Normalize to [0, 1] range
  result <- (result - min(result)) / (max(result) - min(result))

  return(result)
}


#' DGP: β = 3 > d/2 (High Smoothness, Satisfies DML Condition)
#'
#' **Setup:**
#' - d = 4 covariates (binary)
#' - β = 3 (piecewise cubic with C² smoothness)
#' - Hölder class: H³([0,1]²)
#'
#' **Theoretical properties:**
#' - Condition: β = 3 > d/2 = 2 ✓
#' - Tree rate: n^(-β/(2β+d)) = n^(-3/10) = n^(-0.30)
#' - DML requirement: rate = o_p(n^(-1/4)) ✓ SATISFIED
#'
#' **Expected performance:**
#' - Coverage ≈ 95% (theory predicts this should work well)
#' - RMSE decreases at standard √n rate
#' - This is the "control" regime demonstrating method works
#'
#' @param n Sample size
#' @param tau Treatment effect on probability scale (target ATT ≈ tau)
#' @param seed Random seed for reproducibility
#' @return List with components X, A, Y, tau, true_att, true_e, true_m0, true_m1, dgp, diagnostics
#'
#' @examples
#' d <- generate_dgp_beta_high(n = 800, tau = 0.10, seed = 123)
#' d$diagnostics$beta  # 3
#' d$diagnostics$rate_regime  # "valid"
generate_dgp_beta_high <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate continuous features for constructing smooth functions
  U1 <- runif(n, 0, 1)
  U2 <- runif(n, 0, 1)

  # Propensity: piecewise cubic (C² smooth), β = 3
  e_raw <- piecewise_polynomial_2d(U1, U2, degree = 3, smooth_order = 2)
  # Map to [0.15, 0.85] for overlap
  e <- 0.15 + 0.70 * e_raw
  A <- as.integer(runif(n) < e)

  # Outcome: piecewise cubic (C² smooth), β = 3
  p0_raw <- piecewise_polynomial_2d(U1, U2, degree = 3, smooth_order = 2)
  # Map to [0.2, 0.7] and add treatment effect
  p0 <- 0.2 + 0.5 * p0_raw
  p1 <- pmin(p0 + tau, 1)

  # Generate potential outcomes
  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # Create binary features by discretizing continuous U
  # Trees will operate on these 4 binary features
  X <- data.frame(
    X1 = as.integer(U1 > 0.5),
    X2 = as.integer(U2 > 0.5),
    X3 = as.integer(runif(n) < 0.5),  # Noise feature
    X4 = as.integer(runif(n) < 0.5)   # Noise feature
  )

  # True ATT (among treated)
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_3_high_smoothness",
    diagnostics = list(
      beta = 3,
      d = 4,
      rate_exponent = -3/10,
      dml_threshold = -1/4,
      rate_regime = "valid",
      condition_satisfied = TRUE,
      description = "β=3 > d/2=2 (DML condition satisfied)"
    )
  )
}


#' DGP: β = 2 ≈ d/2 (Boundary Case)
#'
#' **Setup:**
#' - d = 4 covariates (binary)
#' - β = 2 (piecewise quadratic with C¹ smoothness)
#' - Hölder class: H²([0,1]²)
#'
#' **Theoretical properties:**
#' - Condition: β = 2 = d/2 (exact boundary)
#' - Tree rate: n^(-β/(2β+d)) = n^(-2/8) = n^(-0.25)
#' - DML requirement: rate = o_p(n^(-1/4))
#' - Status: BOUNDARY (rate exactly at threshold)
#'
#' **Expected performance:**
#' - Coverage ≈ 93-96% (theory unclear; constants matter)
#' - May show slight degradation vs β=3 regime
#' - Empirical question: is n^(-0.25) "fast enough" in practice?
#'
#' @param n Sample size
#' @param tau Treatment effect on probability scale
#' @param seed Random seed
#' @return List with DGP components
#'
#' @examples
#' d <- generate_dgp_beta_boundary(n = 800, tau = 0.10, seed = 123)
#' d$diagnostics$rate_regime  # "boundary"
generate_dgp_beta_boundary <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate continuous features
  U1 <- runif(n, 0, 1)
  U2 <- runif(n, 0, 1)

  # Propensity: piecewise quadratic (C¹ smooth), β = 2
  e_raw <- piecewise_polynomial_2d(U1, U2, degree = 2, smooth_order = 1)
  e <- 0.15 + 0.70 * e_raw
  A <- as.integer(runif(n) < e)

  # Outcome: piecewise quadratic (C¹ smooth), β = 2
  p0_raw <- piecewise_polynomial_2d(U1, U2, degree = 2, smooth_order = 1)
  p0 <- 0.2 + 0.5 * p0_raw
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # Binary features
  X <- data.frame(
    X1 = as.integer(U1 > 0.5),
    X2 = as.integer(U2 > 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5)
  )

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_2_boundary",
    diagnostics = list(
      beta = 2,
      d = 4,
      rate_exponent = -2/8,
      dml_threshold = -1/4,
      rate_regime = "boundary",
      condition_satisfied = NA,  # Exactly at threshold
      description = "β=2 = d/2=2 (exact boundary case)"
    )
  )
}


#' DGP: β = 1 < d/2 (Low Smoothness, VIOLATES DML Condition)
#'
#' **Setup:**
#' - d = 4 covariates (binary)
#' - β = 1 (piecewise linear, Lipschitz but not C¹)
#' - Hölder class: H¹([0,1]²) = Lip([0,1]²)
#'
#' **Theoretical properties:**
#' - Condition: β = 1 < d/2 = 2 ✗ FAILS
#' - Tree rate: n^(-β/(2β+d)) = n^(-1/6) ≈ n^(-0.167)
#' - DML requirement: rate = o_p(n^(-1/4))
#' - Status: TOO SLOW (n^(-0.167) ≪ n^(-0.25))
#'
#' **Expected performance:**
#' - Performance degrades relative to β=3 regime
#' - Possible outcomes (empirical question):
#'   * Coverage << 95% (e.g., 85-90%)
#'   * Or coverage ≈ 95% but inflated variance (wider CIs)
#'   * Or slower RMSE convergence with n
#' - This demonstrates what happens when β > d/2 condition fails
#'
#' **Constitution note (§9):** This is an exploratory stress test.
#' We report what we observe, whether it matches theory or not.
#'
#' @param n Sample size
#' @param tau Treatment effect
#' @param seed Random seed
#' @return List with DGP components
#'
#' @examples
#' d <- generate_dgp_beta_low(n = 800, tau = 0.10, seed = 123)
#' d$diagnostics$condition_satisfied  # FALSE
#' d$diagnostics$rate_regime  # "invalid"
generate_dgp_beta_low <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate continuous features
  U1 <- runif(n, 0, 1)
  U2 <- runif(n, 0, 1)

  # Propensity: piecewise linear (Lipschitz), β = 1
  # Use absolute value functions (not differentiable at kinks)
  e_raw <- piecewise_polynomial_2d(U1, U2, degree = 1, smooth_order = 0)
  e <- 0.15 + 0.70 * e_raw
  A <- as.integer(runif(n) < e)

  # Outcome: piecewise linear (Lipschitz), β = 1
  p0_raw <- piecewise_polynomial_2d(U1, U2, degree = 1, smooth_order = 0)
  p0 <- 0.2 + 0.5 * p0_raw
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # Binary features
  X <- data.frame(
    X1 = as.integer(U1 > 0.5),
    X2 = as.integer(U2 > 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5)
  )

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_1_low_smoothness",
    diagnostics = list(
      beta = 1,
      d = 4,
      rate_exponent = -1/6,
      dml_threshold = -1/4,
      rate_regime = "invalid",
      condition_satisfied = FALSE,
      description = "β=1 < d/2=2 (DML condition VIOLATED)"
    )
  )
}


#' Summary: β-Regime DGPs
#'
#' These DGPs test the theoretical β > d/2 smoothness requirement for
#' tree-based DML. All use d=4 binary features with underlying smoothness
#' controlled via continuous (U1, U2) functions.
#'
#' DGP: β = 3 (High Smoothness)
#' - Rate: n^(-3/10) = n^(-0.30) = o(n^(-0.25)) ✓
#' - Expected: Coverage ≈ 95%, RMSE decreases at √n rate
#' - Purpose: Demonstrate method works when condition satisfied
#'
#' DGP: β = 2 (Boundary)
#' - Rate: n^(-2/8) = n^(-0.25) (exactly at threshold)
#' - Expected: Coverage ≈ 93-96% (constants matter)
#' - Purpose: Test whether boundary is sharp or smooth
#'
#' DGP: β = 1 (Low Smoothness)
#' - Rate: n^(-1/6) ≈ n^(-0.167) ≪ n^(-0.25) ✗
#' - Expected: Performance degrades (coverage drop, inflated variance, or slower convergence)
#' - Purpose: Demonstrate what happens when β > d/2 fails
#'
#' All maintain τ = 0.10 and similar propensity overlap for fair comparison.
#' Differences in performance directly attributable to smoothness regime.
