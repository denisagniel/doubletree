# Data-Generating Processes for Six-Approach Comparison
# Created: 2026-05-01
# Updated: 2026-05-01 - Added continuous covariate DGP (dgp4)
#
# Four DGPs with increasing complexity:
# - Simple: Linear propensity, linear outcome (~2-3 optimal splits)
# - Moderate: 2-way interactions (~4-5 optimal splits)
# - Complex: 3-way interactions (~6-8 optimal splits)
# - Continuous: Mixed binary + continuous covariates (~4-6 optimal splits)
#
# All: True ATT = 0.15

expit <- function(x) 1 / (1 + exp(-x))

#' Generate Simple DGP
#'
#' Linear propensity and outcome functions, binary covariates only.
#' Optimal trees: ~2-3 splits per nuisance
#'
#' @param n Sample size
#' @return List with X (data.frame), A (vector), Y (vector), true_att (scalar)
#' @export
generate_dgp_simple <- function(n = 1000) {
  # Binary covariates
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = rbinom(n, 1, 0.5)
  )

  # Propensity score (function of x1, x2)
  e_true <- expit(-0.5 + 0.3 * X$x1 + 0.3 * X$x2)

  # Treatment
  A <- rbinom(n, 1, e_true)

  # Outcome function for controls (function of x1, x3)
  mu0_true <- 0.2 + 0.15 * X$x1 + 0.15 * X$x3

  # Potential outcomes
  Y0 <- rbinom(n, 1, mu0_true)
  Y1 <- rbinom(n, 1, pmin(mu0_true + 0.15, 1))  # ATT = 0.15

  # Observed outcome
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  true_att <- 0.15

  list(
    X = X,
    A = A,
    Y = Y,
    true_att = true_att,
    e_true = e_true,
    mu0_true = mu0_true
  )
}

#' Generate Moderate DGP
#'
#' 2-way interactions in propensity and outcome, binary covariates.
#' Optimal trees: ~4-5 splits per nuisance
#'
#' @param n Sample size
#' @return List with X, A, Y, true_att
#' @export
generate_dgp_moderate <- function(n = 1000) {
  # Binary covariates
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = rbinom(n, 1, 0.5),
    x4 = rbinom(n, 1, 0.5)
  )

  # Propensity with interaction
  e_true <- expit(-0.5 + 0.3 * X$x1 + 0.2 * X$x2 + 0.3 * X$x1 * X$x2)

  # Treatment
  A <- rbinom(n, 1, e_true)

  # Outcome function with interaction
  mu0_true <- 0.2 + 0.2 * X$x3 + 0.15 * X$x4 + 0.2 * X$x3 * X$x4

  # Potential outcomes
  Y0 <- rbinom(n, 1, mu0_true)
  Y1 <- rbinom(n, 1, pmin(mu0_true + 0.15, 1))

  # Observed outcome
  Y <- A * Y1 + (1 - A) * Y0

  true_att <- 0.15

  list(
    X = X,
    A = A,
    Y = Y,
    true_att = true_att,
    e_true = e_true,
    mu0_true = mu0_true
  )
}

#' Generate Complex DGP
#'
#' Multiple interactions, higher-order terms, binary covariates.
#' Optimal trees: ~6-8 splits per nuisance
#'
#' @param n Sample size
#' @return List with X, A, Y, true_att
#' @export
generate_dgp_complex <- function(n = 1000) {
  # Binary covariates
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = rbinom(n, 1, 0.5),
    x4 = rbinom(n, 1, 0.5),
    x5 = rbinom(n, 1, 0.5)
  )

  # Propensity with multiple interactions
  e_linear <- -0.5 + 0.2 * (X$x1 + X$x2 + X$x3)
  e_interact <- 0.3 * X$x1 * X$x2 + 0.2 * X$x2 * X$x3
  e_true <- expit(e_linear + e_interact)

  # Treatment
  A <- rbinom(n, 1, e_true)

  # Outcome function with multiple interactions
  mu0_linear <- 0.2 + 0.15 * (X$x3 + X$x4 + X$x5)
  mu0_interact <- 0.2 * X$x3 * X$x4 + 0.15 * X$x4 * X$x5
  mu0_true <- mu0_linear + mu0_interact

  # Potential outcomes
  Y0 <- rbinom(n, 1, mu0_true)
  Y1 <- rbinom(n, 1, pmin(mu0_true + 0.15, 1))

  # Observed outcome
  Y <- A * Y1 + (1 - A) * Y0

  true_att <- 0.15

  list(
    X = X,
    A = A,
    Y = Y,
    true_att = true_att,
    e_true = e_true,
    mu0_true = mu0_true
  )
}

#' Generate Continuous Covariate DGP
#'
#' Mixed binary and continuous covariates with interactions.
#' Tests whether tree methods handle discretization appropriately.
#' Optimal trees: ~4-6 splits per nuisance
#'
#' @param n Sample size
#' @return List with X, A, Y, true_att
#' @export
generate_dgp_continuous <- function(n = 1000) {
  # Mixed covariates: 2 binary, 2 continuous
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),                    # Binary
    x2 = rbinom(n, 1, 0.5),                    # Binary
    x3 = runif(n, min = -1, max = 1),          # Continuous [-1, 1]
    x4 = rnorm(n, mean = 0, sd = 1)            # Continuous N(0,1)
  )

  # Propensity score (function of x1, x3, x4)
  # Include both binary and continuous, plus interaction
  e_linear <- -0.5 + 0.3 * X$x1 + 0.4 * X$x3 + 0.2 * X$x4
  e_interact <- 0.2 * X$x1 * X$x3  # Binary × continuous interaction
  e_true <- expit(e_linear + e_interact)

  # Treatment
  A <- rbinom(n, 1, e_true)

  # Outcome function (function of x2, x3, x4)
  # Nonlinear in continuous covariates
  mu0_linear <- 0.2 + 0.15 * X$x2 + 0.2 * X$x3
  mu0_nonlinear <- 0.15 * (X$x4^2 / 2)  # Quadratic in x4
  mu0_interact <- 0.1 * X$x2 * X$x3     # Binary × continuous
  mu0_true <- mu0_linear + mu0_nonlinear + mu0_interact

  # Clip to [0, 1] for outcome probabilities
  mu0_true <- pmax(0.01, pmin(0.99, mu0_true))

  # Potential outcomes
  Y0 <- rbinom(n, 1, mu0_true)
  Y1 <- rbinom(n, 1, pmin(mu0_true + 0.15, 1))

  # Observed outcome
  Y <- A * Y1 + (1 - A) * Y0

  true_att <- 0.15

  list(
    X = X,
    A = A,
    Y = Y,
    true_att = true_att,
    e_true = e_true,
    mu0_true = mu0_true
  )
}
