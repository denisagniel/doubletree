# dgps_continuous.R
# DGPs with continuous and mixed features for validating tree-based DML-ATT
# beyond binary features
#
# DGP 4: Continuous features, binary outcome
# DGP 5: Continuous features, continuous outcome
# DGP 6: Mixed features (2 binary + 2 continuous), binary outcome

# DGP 4: Continuous features, binary outcome
generate_dgp_continuous_binary <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 continuous features: X1, X2 are signal; X3, X4 are noise
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: smooth function of X1, X2
  # Coefficients chosen to give reasonable overlap [0.25, 0.75]
  e <- plogis(-0.5 + 2 * X$X1 - 1.5 * X$X2)
  A <- as.integer(runif(n) < e)

  # Outcome: p0 depends on X1, X2; p1 = p0 + tau (additive on probability scale)
  # Nonlinear relationship to test tree flexibility
  p0 <- plogis(-0.5 + 1.5 * X$X1 + 1.2 * X$X2 + 0.8 * X$X1 * X$X2)
  p1 <- pmin(p0 + tau, 1)

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

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "continuous_binary"
  )
}

# DGP 5: Continuous features, continuous outcome
generate_dgp_continuous_continuous <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 continuous features: X1, X2 are signal; X3, X4 are noise
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: smooth function of X1, X2
  e <- plogis(-0.5 + 2 * X$X1 - 1.5 * X$X2)
  A <- as.integer(runif(n) < e)

  # Continuous outcome: Gaussian with mean depending on X1, X2
  # Scale chosen so tau = 0.10 is meaningful (outcome range ~[0, 1])
  mu0 <- 0.3 + 0.4 * X$X1 + 0.3 * X$X2
  mu1 <- mu0 + tau
  sigma_y <- 0.15  # Noise level

  Y0 <- rnorm(n, mean = mu0, sd = sigma_y)
  Y1 <- rnorm(n, mean = mu1, sd = sigma_y)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT (among treated)
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
    dgp = "continuous_continuous"
  )
}

# DGP 6: Mixed features (2 binary + 2 continuous), binary outcome
generate_dgp_mixed <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Mixed features: 2 binary + 2 continuous
  # Binary: X1, X2 (signal)
  # Continuous: X3, X4 (X3 signal, X4 noise)
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = runif(n, 0, 1),
    X4 = runif(n, 0, 1)
  )

  # Propensity: depends on both binary and continuous features
  e <- plogis(0.3 * X$X1 - 0.4 * X$X2 + 1.5 * X$X3)
  A <- as.integer(runif(n) < e)

  # Outcome: depends on mixed features with interaction
  p0 <- plogis(-0.3 + 0.5 * X$X1 + 0.4 * X$X2 + 1.0 * X$X3)
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
    dgp = "mixed_features"
  )
}
