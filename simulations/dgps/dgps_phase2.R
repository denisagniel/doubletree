# dgps_phase2.R
# Phase 2 DGPs: Stress testing and cases where trees outperform linear
#
# DGP 7: Deep interactions (3-way) - Tree beats linear on complex interactions
# DGP 8: Double nonlinearity (sin/cos in BOTH e AND m0) - Tree beats linear
# DGP 9: Weak overlap - Stress test showing tree maintains coverage
#
# Design rationale:
# - DGPs 1-6 have linear nuisances → favor linear methods
# - DGP7-8: Show tree advantage when truth is nonlinear
# - DGP9: Stress test under extreme propensity scores
#
# Key insight (DGP8): EIF-based ATT is doubly robust, so misspecifying
# only m0(X) OR only e(X) doesn't hurt linear much. Must misspecify BOTH
# simultaneously to show tree advantage. Sin/cos functions achieve this.

# DGP 7: Deep 3-way interaction
# Tree should outperform linear here unless linear includes all interactions
generate_dgp_deep_interaction <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 binary features
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5)
  )

  # Propensity: moderate overlap, depends on 2-way interaction
  # e(X) has X1*X2 interaction
  logit_e <- -0.2 + 0.8 * X$X1 + 0.6 * X$X2 - 1.0 * X$X1 * X$X2
  e <- plogis(logit_e)
  A <- as.integer(runif(n) < e)

  # Outcome: 3-way interaction X1*X2*X3
  # This is where tree should beat linear (unless linear includes all interactions)
  # Use indicator for 3-way: all three = 1
  interaction_3way <- X$X1 * X$X2 * X$X3

  logit_p0 <- -0.5 + 0.5 * X$X1 + 0.4 * X$X2 + 0.3 * X$X3 +
              1.5 * interaction_3way  # Strong 3-way effect
  p0 <- plogis(logit_p0)
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

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "deep_interaction"
  )
}

# DGP 8: Double nonlinearity (e AND m0 both nonlinear)
# Key: Linear is misspecified on BOTH propensity AND outcome
# Double robustness can't save it - need both models right
# Tree should win by approximating both functions
generate_dgp_threshold <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 continuous features
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: STRONGLY NONLINEAR (sin/cos - fundamentally non-polynomial)
  # Linear regression cannot approximate trigonometric functions
  logit_e <- -0.3 + 1.5 * sin(2 * pi * X$X1) + 1.2 * cos(2 * pi * X$X2)
  e <- plogis(logit_e)
  A <- as.integer(runif(n) < e)

  # Outcome: ALSO STRONGLY NONLINEAR (different sin/cos + interaction)
  # Linear completely misspecified on both models
  mu0 <- 0.5 + 0.3 * sin(2 * pi * X$X1) + 0.25 * cos(2 * pi * X$X2) +
         0.2 * X$X1 * X$X2
  mu1 <- mu0 + tau

  sigma_y <- 0.15  # Moderate noise

  Y0 <- rnorm(n, mean = mu0, sd = sigma_y)
  Y1 <- rnorm(n, mean = mu1, sd = sigma_y)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(mu1[treated_idx] - mu0[treated_idx])
  } else {
    tau
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = mu0, true_m1 = mu1,
    dgp = "double_nonlinear"
  )
}

# DGP 9: Weak overlap (extreme propensity scores)
# Stress test: both methods may struggle, but tree's robustness should help
# Tests numerical stability of IPW weights
generate_dgp_weak_overlap <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 continuous features
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: STEEP function creates near 0 and near 1 probabilities
  # Most units are strongly treated or strongly control
  # Overlap region is narrow
  logit_e <- -4 + 8 * X$X1  # Very steep: e ≈ 0 when X1 < 0.5, e ≈ 1 when X1 > 0.5
  e <- plogis(logit_e)
  A <- as.integer(runif(n) < e)

  # Outcome: moderate complexity (smooth with interaction)
  logit_p0 <- -0.3 + 1.0 * X$X1 + 0.8 * X$X2 + 0.5 * X$X1 * X$X2
  p0 <- plogis(logit_p0)
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

  # Calculate overlap diagnostic
  overlap_region <- sum(e > 0.1 & e < 0.9)
  overlap_pct <- 100 * overlap_region / n

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "weak_overlap",
    overlap_pct = overlap_pct  # Diagnostic: % in [0.1, 0.9]
  )
}

# Alias functions for consistency with simulation infrastructure
generate_dgp7 <- generate_dgp_deep_interaction
generate_dgp8 <- generate_dgp_threshold
generate_dgp9 <- generate_dgp_weak_overlap
