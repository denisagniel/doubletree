#' DGPs for β < d/2 Smoothness Regime Testing (Continuous Features)
#'
#' **KEY DESIGN PRINCIPLE:** Use continuous features throughout, let optimaltrees
#' handle discretization adaptively. This allows tree complexity s_n to grow
#' naturally as n increases, following theoretical s_n ~ n^(d/(2β+d)).
#'
#' **Why continuous features?**
#' Original approach (dgps_beta_regimes.R) manually discretized continuous U
#' to binary X, creating information bottleneck. Trees limited to 2^4=16 leaves
#' permanently, regardless of n. Smooth nuisances varied continuously within
#' discrete groups → catastrophic approximation error.
#'
#' **This approach:**
#' - Generate X ~ Uniform[0,1]^d (d=4 continuous features)
#' - Create smooth nuisances as polynomials/functions of continuous X
#' - Pass continuous X to optimaltrees with adaptive discretization:
#'   * discretize_method = "quantiles"
#'   * discretize_bins = "adaptive" → b_n = max(2, ⌈log(n)/3⌉)
#' - At n=800: b_n=3 bins × 4 features = up to 81 possible leaves
#' - Regularization prunes to s_n ~ n^(d/(2β+d)) leaves
#'
#' **Theoretical predictions (d=4, n=800):**
#' - β=3: s_n ~ n^(4/10) ≈ 16 leaves
#' - β=2: s_n ~ n^(4/8) ≈ 28 leaves
#' - β=1: s_n ~ n^(4/6) ≈ 70 leaves
#'
#' Tests the β > d/2 condition: When nuisance functions have Hölder smoothness
#' β, trees converge at rate n^(-β/(2β+d)). For valid DML inference, this rate
#' must be faster than n^(-1/4), requiring β > d/2.
#'
#' Paper reference: Manuscript Section 2.3 (Rate Conditions), Section 4 (Simulations)


#' DGP: β = 3 > d/2 (High Smoothness, Satisfies DML Condition)
#'
#' **Setup:**
#' - d = 4 continuous covariates X ~ Uniform[0,1]^4
#' - β = 3 (cubic polynomials - C² smooth)
#' - Hölder class: Functions with 3rd-order Hölder continuity
#'
#' **Theoretical properties:**
#' - Condition: β = 3 > d/2 = 2 ✓
#' - Theoretical s_n: n^(d/(2β+d)) = n^(4/10) = n^0.4 ≈ 16 at n=800
#' - Tree rate: n^(-β/(2β+d)) = n^(-3/10) = n^(-0.30)
#' - DML requirement: rate = o_p(n^(-1/4)) ✓ SATISFIED
#'
#' **Expected performance:**
#' - Coverage ≈ 95% (valid regime)
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
#' d$diagnostics$theoretical_sn  # ~16
#' sapply(d$X, function(x) length(unique(x)))  # Should be >> 100 (continuous)
generate_dgp_beta_high <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate 4 continuous features
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)
  X3 <- runif(n, 0, 1)
  X4 <- runif(n, 0, 1)
  X <- data.frame(X1, X2, X3, X4)

  # Propensity: cubic polynomial (β=3) on continuous X
  # Use X1, X2 as signal; X3, X4 as noise (to maintain consistency with primary sims)
  e_logit <- 0.5 + 1.2*X1 + 0.8*X2 - 0.6*X1^2 - 0.4*X2^2 +
             0.3*X1^3 + 0.2*X2^3 + 0.1*X1*X2 + 0.1*X1^2*X2 + 0.1*X1*X2^2

  # Map through logistic to ensure [0.1, 0.9] overlap
  e <- plogis(e_logit)
  e <- pmin(pmax(e, 0.1), 0.9)  # Clip for overlap

  # Generate treatment
  A <- as.integer(runif(n) < e)

  # Outcome: cubic polynomial (β=3) on continuous X
  m0_logit <- -0.2 + 0.9*X1 + 0.7*X2 - 0.5*X1^2 - 0.3*X2^2 +
              0.25*X1^3 + 0.15*X2^3 + 0.1*X1*X2

  # Map through logistic and add treatment effect
  p0 <- plogis(m0_logit)
  p1 <- pmin(p0 + tau, 1)  # Additive effect, capped at 1

  # Generate potential outcomes
  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT (among treated)
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  # Theoretical s_n prediction
  d <- 4
  beta <- 3
  theoretical_sn <- n^(d / (2*beta + d))  # n^(4/10) = n^0.4
  rate_exponent <- -beta / (2*beta + d)   # -3/10 = -0.30

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_high",
    diagnostics = list(
      beta = 3,
      d = 4,
      theoretical_sn = theoretical_sn,
      rate_exponent = rate_exponent,
      dml_threshold = -1/4,
      rate_regime = "valid",
      condition_satisfied = TRUE,
      description = "β=3 > d/2=2 (DML condition satisfied, cubic polynomials)"
    )
  )
}


#' DGP: β = 2 = d/2 (Boundary Case)
#'
#' **Setup:**
#' - d = 4 continuous covariates X ~ Uniform[0,1]^4
#' - β = 2 (quadratic polynomials - C¹ smooth)
#' - Hölder class: Functions with 2nd-order Hölder continuity
#'
#' **Theoretical properties:**
#' - Condition: β = 2 = d/2 (exact boundary)
#' - Theoretical s_n: n^(d/(2β+d)) = n^(4/8) = n^0.5 ≈ 28 at n=800
#' - Tree rate: n^(-β/(2β+d)) = n^(-2/8) = n^(-0.25)
#' - DML requirement: rate = o_p(n^(-1/4))
#' - Status: BOUNDARY (rate exactly at threshold)
#'
#' **Expected performance:**
#' - Coverage ≈ 92-96% (theory unclear; constants matter)
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
#' d$diagnostics$theoretical_sn  # ~28
#' d$diagnostics$rate_regime  # "boundary"
generate_dgp_beta_boundary <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate 4 continuous features
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)
  X3 <- runif(n, 0, 1)
  X4 <- runif(n, 0, 1)
  X <- data.frame(X1, X2, X3, X4)

  # Propensity: quadratic polynomial (β=2) on continuous X
  # Drop cubic terms from beta_high, keep up to degree 2
  e_logit <- 0.5 + 1.2*X1 + 0.8*X2 - 0.6*X1^2 - 0.4*X2^2 + 0.2*X1*X2

  e <- plogis(e_logit)
  e <- pmin(pmax(e, 0.1), 0.9)

  A <- as.integer(runif(n) < e)

  # Outcome: quadratic polynomial (β=2) on continuous X
  m0_logit <- -0.2 + 0.9*X1 + 0.7*X2 - 0.5*X1^2 - 0.3*X2^2 + 0.15*X1*X2

  p0 <- plogis(m0_logit)
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  # Theoretical s_n prediction
  d <- 4
  beta <- 2
  theoretical_sn <- n^(d / (2*beta + d))  # n^(4/8) = n^0.5
  rate_exponent <- -beta / (2*beta + d)   # -2/8 = -0.25

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_boundary",
    diagnostics = list(
      beta = 2,
      d = 4,
      theoretical_sn = theoretical_sn,
      rate_exponent = rate_exponent,
      dml_threshold = -1/4,
      rate_regime = "boundary",
      condition_satisfied = NA,  # Exactly at threshold
      description = "β=2 = d/2=2 (exact boundary, quadratic polynomials)"
    )
  )
}


#' DGP: β = 1 < d/2 (Low Smoothness, VIOLATES DML Condition)
#'
#' **Setup:**
#' - d = 4 continuous covariates X ~ Uniform[0,1]^4
#' - β = 1 (absolute value functions - Lipschitz but not differentiable)
#' - Hölder class: Lipschitz continuous functions
#'
#' **Theoretical properties:**
#' - Condition: β = 1 < d/2 = 2 ✗ FAILS
#' - Theoretical s_n: n^(d/(2β+d)) = n^(4/6) = n^(2/3) ≈ 70 at n=800
#' - Tree rate: n^(-β/(2β+d)) = n^(-1/6) ≈ n^(-0.167)
#' - DML requirement: rate = o_p(n^(-1/4))
#' - Status: TOO SLOW (n^(-0.167) < n^(-0.25))
#'
#' **Expected performance:**
#' - Performance degrades relative to β=3 regime
#' - Possible outcomes (empirical question):
#'   * Coverage < 95% (e.g., 85-92%)
#'   * Or coverage ≈ 95% but inflated variance (wider CIs)
#'   * Or slower RMSE convergence with n
#' - This demonstrates what happens when β > d/2 condition fails
#'
#' **Constitution note (§9):** Exploratory stress test. Report what we observe,
#' whether it matches theory or not.
#'
#' @param n Sample size
#' @param tau Treatment effect
#' @param seed Random seed
#' @return List with DGP components
#'
#' @examples
#' d <- generate_dgp_beta_low(n = 800, tau = 0.10, seed = 123)
#' d$diagnostics$theoretical_sn  # ~70
#' d$diagnostics$condition_satisfied  # FALSE
generate_dgp_beta_low <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate 4 continuous features
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)
  X3 <- runif(n, 0, 1)
  X4 <- runif(n, 0, 1)
  X <- data.frame(X1, X2, X3, X4)

  # Propensity: absolute value functions (β=1, Lipschitz but not differentiable)
  # |X - c| is Lipschitz continuous (β=1) but has non-differentiable kink at X=c
  e_logit <- 0.5 + 1.5*abs(X1 - 0.5) + 1.2*abs(X2 - 0.5) +
             0.3*abs(X1 - 0.3)*abs(X2 - 0.7) +
             0.2*abs(X1 - 0.7) + 0.2*abs(X2 - 0.3)

  e <- plogis(e_logit)
  e <- pmin(pmax(e, 0.1), 0.9)

  A <- as.integer(runif(n) < e)

  # Outcome: absolute value functions (β=1)
  m0_logit <- -0.2 + 1.0*abs(X1 - 0.4) + 0.8*abs(X2 - 0.6) +
              0.3*abs(X1 - 0.6)*abs(X2 - 0.4)

  p0 <- plogis(m0_logit)
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  # Theoretical s_n prediction
  d <- 4
  beta <- 1
  theoretical_sn <- n^(d / (2*beta + d))  # n^(4/6) = n^(2/3)
  rate_exponent <- -beta / (2*beta + d)   # -1/6 ≈ -0.167

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "beta_low",
    diagnostics = list(
      beta = 1,
      d = 4,
      theoretical_sn = theoretical_sn,
      rate_exponent = rate_exponent,
      dml_threshold = -1/4,
      rate_regime = "invalid",
      condition_satisfied = FALSE,
      description = "β=1 < d/2=2 (DML condition VIOLATED, absolute value functions)"
    )
  )
}


#' Summary: β-Regime DGPs with Continuous Features
#'
#' These DGPs test the theoretical β > d/2 smoothness requirement for
#' tree-based DML using continuous features throughout.
#'
#' **Key innovation vs original approach:**
#' - Original (dgps_beta_regimes.R): Generated continuous U, manually discretized
#'   to binary X → information bottleneck, trees limited to 16 leaves
#' - This approach: Generate continuous X, pass to optimaltrees with adaptive
#'   discretization → tree complexity s_n grows naturally with n
#'
#' **DGP: β = 3 (High Smoothness)**
#' - Nuisances: Cubic polynomials (C² smooth)
#' - Theoretical s_n: n^0.4 ≈ 16 leaves at n=800
#' - Rate: n^(-0.30) = o(n^(-0.25)) ✓
#' - Expected: Coverage ≈ 95%, RMSE decreases at √n rate
#'
#' **DGP: β = 2 (Boundary)**
#' - Nuisances: Quadratic polynomials (C¹ smooth)
#' - Theoretical s_n: n^0.5 ≈ 28 leaves at n=800
#' - Rate: n^(-0.25) (exactly at threshold)
#' - Expected: Coverage ≈ 92-96% (constants matter)
#'
#' **DGP: β = 1 (Low Smoothness)**
#' - Nuisances: Absolute value functions (Lipschitz, not differentiable)
#' - Theoretical s_n: n^(2/3) ≈ 70 leaves at n=800
#' - Rate: n^(-0.167) < n^(-0.25) ✗
#' - Expected: Performance degrades (coverage drop, wider CIs, or slower convergence)
#'
#' All maintain τ = 0.10 and similar propensity overlap for fair comparison.
#' Differences in performance directly attributable to smoothness regime.
