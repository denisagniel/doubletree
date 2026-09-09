#' ATT Estimation with ORACLE Nuisances (no estimation, no cross-fitting)
#'
#' Plugs the DGP's TRUE propensity e_0(X) and TRUE control outcome mu_0(X)
#' straight into the same ATT orthogonal score that att_linear() and
#' att_forest() use. Nothing is estimated, so there is nothing to cross-fit and
#' no fold seed to share: the only randomness is the draw itself.
#'
#' Paper reference: manuscript.tex `eq:score`,
#'   psi(O; theta, eta) = pi^{-1}[ A{Y - mu(X) - theta}
#'                                 - (1 - A) w(X){Y - mu(X)} ],  w = e/(1 - e).
#'
#' ROLE. This is the EFFICIENCY DENOMINATOR of claim C1 in
#' quality_reports/specs/2026-09-08_head-to-head-comparison.md (its §1 table):
#' the study reports RMSE RATIOS to this arm rather than raw "beats GLM"
#' comparisons. Without it the study collapses into a leaderboard, which the
#' paper makes no claim to win -- under standard regularity a
#' correctly-specified linear working model, a correctly-specified black-box
#' AIPW/DML estimator, and doubletree under exact structural sparsity all attain
#' the SAME semiparametric efficiency bound. "Matched the oracle" is a
#' statement the theory licenses; "beat the GLM" is not.
#'
#' Expected behavior: unbiased at every n up to Monte Carlo error, with
#' near-nominal Wald coverage. A coverage shortfall here is a harness bug or a
#' wrong theta_0, not a property of any estimator -- which is exactly why this
#' arm is worth running.
#'
#' @param X Data.frame of covariates. Not used in the score (the nuisances are
#'   supplied), but required and length-checked so the call site cannot silently
#'   pair nuisance vectors with the wrong dataset.
#' @param A Binary treatment (0/1).
#' @param Y Outcome (binary 0/1 or continuous).
#' @param e_true True propensity e_0(X_i), one value per row of \code{X}. Every
#'   DGP draw in this project returns this as \code{e_true}.
#' @param m0_true True control outcome mu_0(X_i), one per row. Returned by the
#'   DGP draws as \code{mu0_true}.
#' @param clip Length-2 propensity clip, \code{c(0.01, 0.99)} by default to
#'   match \code{att_linear()}/\code{att_forest()}'s convention. Clipping a TRUE
#'   propensity would bias the efficiency benchmark, so this arm asserts the
#'   clip is INERT rather than relying on it: see \code{require_no_clipping}.
#' @param require_no_clipping If \code{TRUE} (default), error when \code{clip}
#'   would actually bind. The clip is kept for convention-consistency with the
#'   other arms and is verifiably a no-op in every regime this arm is used in
#'   (R1: \eqn{e_0 \in [0.20, 0.75]}; R3: \eqn{e_0 \in [0.10, 0.85]}). Set
#'   \code{FALSE} only to deliberately study a clipped oracle, and say so.
#' @param verbose Print progress messages.
#'
#' @return List with components (same shape as \code{att_linear()}):
#'   - theta: Point estimate of ATT
#'   - sigma: Standard error
#'   - ci: 95% confidence interval (length 2 vector)
#'   - n_treated: Number of treated units
#'   - convergence: Always "converged" (nothing is fit)
#'   - method: "oracle_aipw"
#'   - hyperparams: list(clip, n_clipped, e_true_range)
#'
#' @details
#' There is no cross-fitting step, so the procedure is just:
#' 1. Take e_0, mu_0 as given (the DGP's own closed forms).
#' 2. Compute the `eq:score` orthogonal score.
#' 3. Solve sum_i psi_i(theta) = 0 for theta.
#' 4. Inference via the EIF-based variance estimator, identically to the other
#'    arms.
#'
#' The score is written in the SAME algebraic form as \code{att_linear()}'s and
#' \code{att_forest()}'s, term for term, deliberately: the point of this study
#' is a comparison across nuisance estimators, so the score must not be a
#' second source of difference between arms.
#'
#' @examples
#' source("simulations/head_to_head_comparison/code/common.R")
#' dgp <- make_dgp("A")
#' set.seed(1)
#' d <- dgp$draw(500)
#' fit <- att_oracle(d$X, d$A, d$Y, e_true = d$e_true, m0_true = d$mu0_true)
#' print(fit$theta)

att_oracle <- function(X, A, Y, e_true, m0_true,
                       clip = c(0.01, 0.99), require_no_clipping = TRUE,
                       verbose = FALSE) {

  # Input validation
  n <- nrow(X)
  if (length(A) != n || length(Y) != n) {
    stop("X, A, Y must have same number of rows")
  }
  if (length(e_true) != n || length(m0_true) != n) {
    stop("e_true and m0_true must have one value per row of X (got ",
         length(e_true), " and ", length(m0_true), " for n = ", n, ")")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary (0/1)")
  }
  if (!all(is.finite(e_true)) || any(e_true <= 0) || any(e_true >= 1)) {
    stop("e_true must be finite and strictly inside (0, 1); got range [",
         signif(min(e_true), 4), ", ", signif(max(e_true), 4), "]")
  }
  if (!all(is.finite(m0_true))) {
    stop("m0_true must be finite")
  }
  stopifnot(length(clip) == 2L, clip[[1]] < clip[[2]])

  n_treated <- sum(A == 1)
  if (n_treated == 0) {
    stop("No treated units")
  }
  if (n_treated == n) {
    stop("No control units")
  }

  # The clip is present for convention-parity with att_linear()/att_forest() and
  # is asserted to be inert rather than silently modifying the truth. A binding
  # clip here would mean this arm is no longer an ORACLE, and the RMSE ratios
  # every C1 statement is built on would be ratios to a biased denominator.
  n_clipped <- sum(e_true <= clip[[1]] | e_true >= clip[[2]])
  if (n_clipped > 0L && require_no_clipping) {
    stop("clip = (", clip[[1]], ", ", clip[[2]], ") binds on ", n_clipped,
         " of ", n, " TRUE propensities (range [", signif(min(e_true), 4), ", ",
         signif(max(e_true), 4), "]). A clipped oracle is not an oracle and ",
         "cannot serve as this study's efficiency denominator. Either use a DGP ",
         "whose e_0 is interior, or pass require_no_clipping = FALSE and say so ",
         "in the write-up.", call. = FALSE)
  }
  e_hat <- pmax(pmin(e_true, clip[[2]]), clip[[1]])
  m0_hat <- m0_true

  # --- Orthogonal Score and Point Estimate ---
  # ATT orthogonal score (manuscript.tex eq:score; Chernozhukov et al. 2018):
  # psi_i(theta) = (A_i/pi)*(Y_i - m0_i - theta) - (1/pi)*(e_i*(1 - A_i)/(1 - e_i))*(Y_i - m0_i)

  pi_hat <- mean(A)

  # Solve for theta: sum(psi_i(theta)) = 0
  score_at_zero_term1 <- (A / pi_hat) * (Y - m0_hat)
  score_at_zero_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)

  score_at_zero <- score_at_zero_term1 - score_at_zero_term2
  sum_a_over_pi <- sum(A / pi_hat)

  theta <- sum(score_at_zero) / sum_a_over_pi

  # --- Variance Estimation ---
  # Score at theta-hat (proper Neyman-orthogonal score)
  score_term1 <- (A / pi_hat) * (Y - m0_hat - theta)
  score_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)

  score_values <- score_term1 - score_term2

  # No `is.finite` scrubbing of the weight term here, unlike the estimated arms:
  # e_true is asserted interior above, so a non-finite w = e/(1 - e) is
  # impossible unless that assertion is wrong, and a silent zero would hide it.
  if (!all(is.finite(score_values))) {
    stop("Non-finite score values with an interior oracle propensity; ",
         "this is a bug, not a data condition.", call. = FALSE)
  }

  # Variance: Var(sqrt(n) * theta-hat) = E[psi(theta)^2]
  variance <- mean(score_values^2)
  sigma <- sqrt(variance / n)

  # 95% CI (normal approximation)
  ci <- theta + c(-1.96, 1.96) * sigma

  result <- list(
    theta = theta,
    sigma = sigma,
    ci = ci,
    n_treated = n_treated,
    convergence = "converged",
    method = "oracle_aipw",
    hyperparams = list(
      clip = clip,
      n_clipped = n_clipped,
      e_true_range = range(e_true)
    )
  )

  if (verbose) {
    cat(sprintf("\nOracle-AIPW ATT Results:\n"))
    cat(sprintf("  ATT estimate: %.4f\n", theta))
    cat(sprintf("  Std error:    %.4f\n", sigma))
    cat(sprintf("  95%% CI:       [%.4f, %.4f]\n", ci[1], ci[2]))
    cat(sprintf("  e_0 range:    [%.4f, %.4f], clipped: %d\n",
                min(e_true), max(e_true), n_clipped))
  }

  return(result)
}

# Alias for consistency with simulation scripts
method_oracle_aipw <- att_oracle
