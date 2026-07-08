# EIF Component Decomposition Utilities
# Helper functions for analyzing EIF components and their contribution to bias
# Created: 2026-05-27

#' Compute true propensity score for each DGP
#'
#' @param X Covariate data frame
#' @param dgp DGP number (1-4)
#' @return Vector of true propensity scores
#' @export
compute_true_propensity <- function(X, dgp) {
  expit <- function(x) 1 / (1 + exp(-x))

  if (dgp == 1) {
    # Simple DGP
    e_true <- expit(-0.5 + 0.3 * X$x1 + 0.3 * X$x2)
  } else if (dgp == 2) {
    # Moderate DGP
    e_true <- expit(-0.5 + 0.3 * X$x1 + 0.2 * X$x2 + 0.3 * X$x1 * X$x2)
  } else if (dgp == 3) {
    # Complex DGP
    e_linear <- -0.5 + 0.2 * (X$x1 + X$x2 + X$x3)
    e_interact <- 0.3 * X$x1 * X$x2 + 0.2 * X$x2 * X$x3
    e_true <- expit(e_linear + e_interact)
  } else if (dgp == 4) {
    # Continuous DGP
    e_linear <- -0.5 + 0.3 * X$x1 + 0.4 * X$x3 + 0.2 * X$x4
    e_interact <- 0.2 * X$x1 * X$x3
    e_true <- expit(e_linear + e_interact)
  } else {
    stop(sprintf("Unknown DGP: %d", dgp))
  }

  return(e_true)
}

#' Compute true outcome function for each DGP
#'
#' @param X Covariate data frame
#' @param dgp DGP number (1-4)
#' @return Vector of true E[Y(0)|X]
#' @export
compute_true_outcome <- function(X, dgp) {
  if (dgp == 1) {
    # Simple DGP
    mu0_true <- 0.2 + 0.15 * X$x1 + 0.15 * X$x3
  } else if (dgp == 2) {
    # Moderate DGP
    mu0_true <- 0.2 + 0.2 * X$x3 + 0.15 * X$x4 + 0.2 * X$x3 * X$x4
  } else if (dgp == 3) {
    # Complex DGP
    mu0_linear <- 0.2 + 0.15 * (X$x3 + X$x4 + X$x5)
    mu0_interact <- 0.2 * X$x3 * X$x4 + 0.15 * X$x4 * X$x5
    mu0_true <- mu0_linear + mu0_interact
  } else if (dgp == 4) {
    # Continuous DGP
    mu0_linear <- 0.2 + 0.15 * X$x2 + 0.2 * X$x3
    mu0_nonlinear <- 0.15 * (X$x4^2 / 2)
    mu0_interact <- 0.1 * X$x2 * X$x3
    mu0_true <- mu0_linear + mu0_nonlinear + mu0_interact
    # Clip to [0.01, 0.99]
    mu0_true <- pmax(0.01, pmin(0.99, mu0_true))
  } else {
    stop(sprintf("Unknown DGP: %d", dgp))
  }

  return(mu0_true)
}

#' Decompose EIF into components for bias analysis
#'
#' The EIF for ATT is:
#' ψ(O; θ, η) = [A/P(A=1)] × [Y - m0(X)]
#'              - [A/P(A=1) - 1] × [e(X)/(1-e(X))] × [Y - m0(X)]
#'              - θ
#'
#' This can be separated into:
#' - Component 1: Direct outcome model (treated units)
#' - Component 2: Propensity-weighted residuals (control units)
#'
#' @param X Covariate matrix/data.frame
#' @param A Treatment indicator
#' @param Y Observed outcome
#' @param e_hat Estimated propensity scores
#' @param mu0_hat Estimated outcome under control
#' @param e_true True propensity scores (optional)
#' @param mu0_true True outcome under control (optional)
#' @param theta_true True ATT (default 0.15)
#' @return List with component-wise diagnostics
#' @export
decompose_eif_components <- function(X, A, Y, e_hat, mu0_hat,
                                     e_true = NULL, mu0_true = NULL,
                                     theta_true = 0.15) {

  n <- length(A)
  treated <- A == 1
  control <- A == 0
  n_treated <- sum(treated)
  n_control <- sum(control)
  P_A1 <- mean(A)

  # ---------------------------
  # Component 1: Outcome model on treated units
  # ---------------------------
  # True component (if provided)
  if (!is.null(mu0_true)) {
    comp1_true <- mean(Y[treated] - mu0_true[treated])
  } else {
    comp1_true <- NA_real_
  }

  # Estimated component
  comp1_est <- mean(Y[treated] - mu0_hat[treated])

  # Bias in component 1
  comp1_bias <- if (!is.na(comp1_true)) comp1_est - comp1_true else NA_real_

  # Outcome model error on treated units
  if (!is.null(mu0_true)) {
    outcome_error_treated <- mu0_hat[treated] - mu0_true[treated]
    outcome_rmse_treated <- sqrt(mean(outcome_error_treated^2))
    outcome_mae_treated <- mean(abs(outcome_error_treated))
    outcome_bias_treated <- mean(outcome_error_treated)
  } else {
    outcome_rmse_treated <- NA_real_
    outcome_mae_treated <- NA_real_
    outcome_bias_treated <- NA_real_
  }

  # ---------------------------
  # Component 2: Propensity-weighted residuals on control units
  # ---------------------------
  # Propensity weights
  weights_est <- e_hat[control] / (1 - e_hat[control])

  if (!is.null(e_true)) {
    weights_true <- e_true[control] / (1 - e_true[control])
  } else {
    weights_true <- NULL
  }

  # Residuals (using estimated outcome model)
  residuals <- Y[control] - mu0_hat[control]

  # True component (if true propensity provided)
  if (!is.null(weights_true)) {
    comp2_true <- -mean(weights_true * residuals)
  } else {
    comp2_true <- NA_real_
  }

  # Estimated component
  comp2_est <- -mean(weights_est * residuals)

  # Bias in component 2
  comp2_bias <- if (!is.na(comp2_true)) comp2_est - comp2_true else NA_real_

  # Weight diagnostics
  mean_weight_est <- mean(weights_est)
  max_weight_est <- max(weights_est)
  min_weight_est <- min(weights_est)
  extreme_weights <- sum(weights_est > 10) / length(weights_est)

  if (!is.null(weights_true)) {
    mean_weight_true <- mean(weights_true)
    max_weight_true <- max(weights_true)
    weight_bias <- mean(weights_est - weights_true)
  } else {
    mean_weight_true <- NA_real_
    max_weight_true <- NA_real_
    weight_bias <- NA_real_
  }

  # Propensity score error on control units
  if (!is.null(e_true)) {
    ps_error_control <- e_hat[control] - e_true[control]
    ps_rmse_control <- sqrt(mean(ps_error_control^2))
    ps_mae_control <- mean(abs(ps_error_control))
    ps_bias_control <- mean(ps_error_control)
  } else {
    ps_rmse_control <- NA_real_
    ps_mae_control <- NA_real_
    ps_bias_control <- NA_real_
  }

  # Outcome model error on control units
  if (!is.null(mu0_true)) {
    outcome_error_control <- mu0_hat[control] - mu0_true[control]
    outcome_rmse_control <- sqrt(mean(outcome_error_control^2))
    outcome_mae_control <- mean(abs(outcome_error_control))
    outcome_bias_control <- mean(outcome_error_control)
  } else {
    outcome_rmse_control <- NA_real_
    outcome_mae_control <- NA_real_
    outcome_bias_control <- NA_real_
  }

  # ---------------------------
  # Overall ATT estimate
  # ---------------------------
  theta_hat <- comp1_est + comp2_est
  theta_bias <- theta_hat - theta_true

  # Total bias decomposition
  total_bias_explained <- if (!is.na(comp1_bias) && !is.na(comp2_bias)) {
    comp1_bias + comp2_bias
  } else {
    NA_real_
  }

  # ---------------------------
  # Return diagnostics
  # ---------------------------
  return(list(
    # Component 1: Outcome model on treated
    comp1_true = comp1_true,
    comp1_est = comp1_est,
    comp1_bias = comp1_bias,
    outcome_rmse_treated = outcome_rmse_treated,
    outcome_mae_treated = outcome_mae_treated,
    outcome_bias_treated = outcome_bias_treated,

    # Component 2: Propensity-weighted residuals on control
    comp2_true = comp2_true,
    comp2_est = comp2_est,
    comp2_bias = comp2_bias,
    ps_rmse_control = ps_rmse_control,
    ps_mae_control = ps_mae_control,
    ps_bias_control = ps_bias_control,
    outcome_rmse_control = outcome_rmse_control,
    outcome_mae_control = outcome_mae_control,
    outcome_bias_control = outcome_bias_control,

    # Weight diagnostics
    mean_weight_est = mean_weight_est,
    max_weight_est = max_weight_est,
    min_weight_est = min_weight_est,
    extreme_weights = extreme_weights,
    mean_weight_true = mean_weight_true,
    max_weight_true = max_weight_true,
    weight_bias = weight_bias,

    # Overall
    theta_hat = theta_hat,
    theta_true = theta_true,
    theta_bias = theta_bias,
    total_bias_explained = total_bias_explained,

    # Sample info
    n = n,
    n_treated = n_treated,
    n_control = n_control,
    P_A1 = P_A1
  ))
}

#' Compute influence function values for ATT
#'
#' @param X Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param e_hat Estimated propensity scores
#' @param mu0_hat Estimated outcome under control
#' @param theta_hat Estimated ATT
#' @return Vector of influence function values
#' @export
compute_eif_values <- function(X, A, Y, e_hat, mu0_hat, theta_hat) {
  n <- length(A)
  P_A1 <- mean(A)

  # EIF for ATT:
  # ψ = [A/P(A=1)] × [Y - m0(X)] - [A/P(A=1) - 1] × [e(X)/(1-e(X))] × [Y - m0(X)] - θ

  psi <- (A / P_A1) * (Y - mu0_hat) -
         (A / P_A1 - 1) * (e_hat / (1 - e_hat)) * (Y - mu0_hat) -
         theta_hat

  return(psi)
}

#' Compute standard error from EIF
#'
#' @param psi Vector of influence function values
#' @return Standard error estimate
#' @export
compute_eif_se <- function(psi) {
  n <- length(psi)
  # Correct SE formula: sqrt(var(psi) / n)
  se <- sqrt(mean(psi^2) / n)
  return(se)
}
