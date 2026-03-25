# dgps_att_correct.R
# DGPs where tau directly controls the ATT (not just a logit shift)
#
# Strategy: Use additive effects on the probability scale
# p1 = min(p0 + tau, 1) ensures E[Y1 - Y0] ≈ tau

# DGP 1: Simple binary (4 features, controlled ATT)
generate_dgp_binary_att <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 binary features: X1, X2 are signal; X3, X4 are noise
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5)
  )

  # Propensity: depends on X1, X2
  e <- plogis(0.5 * X$X1 - 0.3 * X$X2)
  A <- as.integer(runif(n) < e)

  # Outcome: p0 depends on X1, X2; p1 = p0 + tau (additive on probability scale)
  p0 <- plogis(-0.3 + 0.5 * X$X1 + 0.4 * X$X2)
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
    tau  # Fallback (shouldn't happen with reasonable n)
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "binary_att_controlled"
  )
}

# DGP 2: Continuous features (4 features, controlled ATT)
generate_dgp_continuous_att <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 continuous features: X1, X2 are signal; X3, X4 are noise
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: smooth function of X1, X2
  e <- plogis(0.7 * X$X1 - 0.5 * X$X2)
  A <- as.integer(runif(n) < e)

  # Outcome: p0 depends on X1, X2; p1 = p0 + tau (additive)
  p0 <- plogis(-0.3 + 0.6 * X$X1 + 0.5 * X$X2)
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
    dgp = "continuous_att_controlled"
  )
}

# DGP 3: Moderate complexity (5 binary features, controlled ATT)
generate_dgp_moderate_att <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 5 binary features: X1, X2, X3 are signal; X4, X5 are noise
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5),
    X5 = as.integer(runif(n) < 0.5)
  )

  # Propensity: depends on X1, X2, X3
  e <- plogis(0.5 * X$X1 - 0.3 * X$X2 + 0.2 * X$X3)
  A <- as.integer(runif(n) < e)

  # Outcome: depends on X1, X2, X3
  p0 <- plogis(-0.2 + 0.4 * X$X1 + 0.3 * X$X2 + 0.2 * X$X3)
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
    dgp = "moderate_att_controlled"
  )
}
