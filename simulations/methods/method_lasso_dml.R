#' ATT Estimation with Cross-Fitted Lasso Nuisances (glmnet)
#'
#' Cross-fitted ATT estimator using the lasso over the EXPANDED
#' main-effects-plus-two-way-interactions basis for both propensity e(X) and
#' control outcome m0(X). Matches the structure of \code{att_linear()} and
#' \code{att_forest()} exactly -- same stratified K-fold construction, same
#' orthogonal score (manuscript.tex \code{eq:score}), same propensity clip
#' \eqn{[0.01, 0.99]} -- and differs from \code{att_linear(interactions = TRUE)}
#' in exactly one respect: the coefficients are L1-penalised with the penalty
#' chosen by cross-validation inside each outer fold.
#'
#' Paper reference: \code{sec:lasso-comparison}. That section argues about lasso
#' DML in prose; this is the empirical instance it obliges. Requested as
#' decision 3 of
#' quality_reports/specs/2026-09-08_head-to-head-comparison.md §8 and confirmed
#' by the user 2026-09-08.
#'
#' DEPENDENCY CONVENTION. \code{glmnet} is used ad hoc via
#' \code{requireNamespace()} and is deliberately NOT added to doubletree's
#' DESCRIPTION. That mirrors \code{method_forest.R}'s handling of \code{ranger},
#' which is likewise absent from DESCRIPTION (checked 2026-09-09): these files
#' live under \code{simulations/}, are never installed with the package, and are
#' not part of any package code path, so adding them to Imports/Suggests would
#' impose a dependency on package users for code they cannot reach.
#'
#' TUNING ASYMMETRY, STATED NOT HIDDEN. Per the study's tuning charter (spec §2),
#' this arm CV-tunes lambda, exactly as \code{doubletree::estimate_att_crossfit()}
#' CV-tunes its own regularisation by default. \code{att_linear()} and
#' \code{att_forest()} have no corresponding penalty to tune (a GLM has none; the
#' forest uses \code{ranger} defaults). So the arms split into a CV-tuned group
#' \{lasso-DML, estimate_att_crossfit()\} and an untuned group \{GLM x 2, RF\},
#' and the flagship \code{estimate_att()} sits in neither (its \code{leaf_budget}
#' is fixed in advance per DGP, matching the paper's own framing that the analyst
#' sets it). This is a real asymmetry in the comparison, not parity; it is
#' recorded here and in the study README rather than papered over.
#'
#' Expected behavior: on DGP-A (the shared-interaction gate-swapped DGP of spec
#' §3) the true nuisances ARE in this basis -- both are exactly representable as
#' main effects plus the single \eqn{X_1 \colon X_2} interaction -- so lasso-DML
#' should be approximately unbiased where \code{att_linear(interactions = FALSE)}
#' is not. That is the point of including it: it separates "the GLM arm lacked
#' the right basis" from "any linear-in-parameters method fails here".
#'
#' @param X Data.frame of covariates.
#' @param A Binary treatment (0/1).
#' @param Y Outcome (binary 0/1 or continuous).
#' @param K Number of folds for cross-fitting (default 5).
#' @param seed Random seed for fold creation. Pass the SAME value used for
#'   \code{att_linear()}/\code{att_forest()} in the same replication: the fold
#'   construction below is byte-identical to theirs, so an identical seed yields
#'   an identical partition and cross-fitting noise is not a nuisance confound
#'   between arms (spec §2, first bullet).
#' @param cv_nfolds Inner CV folds for \code{glmnet::cv.glmnet()} (default 5).
#' @param lambda_rule Which \code{cv.glmnet()} lambda to use: \code{"min"}
#'   (default, \code{lambda.min}) or \code{"1se"}. \code{"min"} is the
#'   prediction-optimal choice and therefore the right default for a NUISANCE
#'   estimate feeding an orthogonal score; \code{"1se"} deliberately
#'   under-fits, which for a nuisance is a bias source, not parsimony.
#' @param interactions Include two-way interactions in the basis (default
#'   \code{TRUE} -- this is the whole point of the arm). \code{FALSE} gives a
#'   main-effects lasso, retained only so the basis and the penalty can be
#'   varied one at a time if a diagnostic ever needs it.
#' @param verbose Print progress messages.
#'
#' @return List with components (same shape as \code{att_linear()}):
#'   - theta: Point estimate of ATT
#'   - sigma: Standard error
#'   - ci: 95% confidence interval (length 2 vector)
#'   - n_treated: Number of treated units
#'   - convergence: "converged"
#'   - method: "lasso_dml"
#'   - hyperparams: list(K, cv_nfolds, lambda_rule, interactions, n_basis,
#'       lambda_e, lambda_m0, nz_e, nz_m0)
#'
#' @details
#' Cross-fitting procedure (identical in shape to \code{att_linear()}):
#' 1. Split data into K folds, stratified by treatment.
#' 2. For each fold k:
#'    - CV-tune and fit lasso e(X) on folds != k over the expanded basis
#'      (binomial family); predict on fold k.
#'    - CV-tune and fit lasso m0(X) on CONTROL units in folds != k
#'      (binomial for binary Y, gaussian otherwise); predict on fold k.
#' 3. Clip e-hat to [0.01, 0.99] (att_linear()'s convention).
#' 4. Compute the eq:score orthogonal score; solve sum_i psi_i(theta) = 0.
#' 5. Inference via the EIF-based variance estimator.
#'
#' NO SILENT FALLBACKS. Unlike \code{att_linear()}, which falls back to a
#' main-effects \code{glm} when the interaction fit warns or errors, this arm has
#' no fallback: a lasso that cannot be fit is a condition the write-up needs to
#' know about, and a fallback would report main-effects numbers under the
#' "lasso-DML" label. Per this project's R conventions (an estimator failure must
#' surface, not be swallowed), any glmnet condition propagates.
#'
#' @examples
#' source("simulations/head_to_head_comparison/code/common.R")
#' dgp <- make_dgp("A")
#' set.seed(1)
#' d <- dgp$draw(1000)
#' fit <- att_lasso_dml(d$X, d$A, d$Y, K = 5, seed = 123)
#' print(fit$theta)

att_lasso_dml <- function(X, A, Y, K = 5, seed = NULL,
                          cv_nfolds = 5, lambda_rule = c("min", "1se"),
                          interactions = TRUE, verbose = FALSE) {

  lambda_rule <- match.arg(lambda_rule)

  # Check required package. Same posture as method_forest.R's ranger check:
  # ad hoc, not a DESCRIPTION dependency (see the DEPENDENCY CONVENTION note).
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' required for lasso-DML ATT estimation. Install with: install.packages('glmnet')")
  }

  # Input validation
  n <- nrow(X)
  if (length(A) != n || length(Y) != n) {
    stop("X, A, Y must have same number of rows")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary (0/1)")
  }

  # Detect outcome type
  outcome_type <- if (all(Y %in% c(0, 1))) "binary" else "continuous"
  if (verbose) {
    cat(sprintf("Detected %s outcome\n", outcome_type))
  }

  n_treated <- sum(A == 1)
  if (n_treated == 0) {
    stop("No treated units")
  }
  if (n_treated == n) {
    stop("No control units")
  }

  # --- Expanded basis, built ONCE on the full X -----------------------------
  # Built on all rows so the column set is identical across folds; any
  # per-fold-constant column is dropped by glmnet's own standardisation, not by
  # a differing model.matrix() contrast per fold (which would make the fold-wise
  # predictions incomparable).
  basis_formula <- if (interactions) ~ .^2 else ~ .
  B <- stats::model.matrix(basis_formula, data = X)
  B <- B[, colnames(B) != "(Intercept)", drop = FALSE]
  if (ncol(B) < 2L) {
    stop("Expanded basis has ", ncol(B), " column(s); glmnet needs at least 2. ",
         "X has ", ncol(X), " column(s).", call. = FALSE)
  }
  if (verbose) {
    cat(sprintf("Basis: %d columns (%s)\n", ncol(B),
                if (interactions) "main effects + two-way interactions" else "main effects"))
  }

  # Create folds (stratified by treatment).
  # BYTE-IDENTICAL to att_linear()/att_forest()'s fold code, so the same `seed`
  # yields the same partition across all three arms (spec §2).
  if (!is.null(seed)) set.seed(seed)

  treated_idx <- which(A == 1)
  control_idx <- which(A == 0)

  folds <- integer(n)
  folds[treated_idx] <- sample(rep(1:K, length.out = length(treated_idx)))
  folds[control_idx] <- sample(rep(1:K, length.out = length(control_idx)))

  # Initialize storage for cross-fitted predictions
  e_hat <- numeric(n)    # Propensity scores
  m0_hat <- numeric(n)   # Control outcome predictions

  lambda_e <- numeric(K)
  lambda_m0 <- numeric(K)
  nz_e <- integer(K)
  nz_m0 <- integer(K)

  # Cross-fitting loop
  for (k in 1:K) {
    if (verbose) {
      cat(sprintf("Fold %d/%d...\n", k, K))
    }

    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    # --- Propensity Score e(X) ---
    cv_e <- glmnet::cv.glmnet(
      x = B[train_idx, , drop = FALSE], y = A[train_idx],
      family = "binomial", alpha = 1, nfolds = cv_nfolds
    )
    lam_e <- if (lambda_rule == "min") cv_e$lambda.min else cv_e$lambda.1se
    e_hat[test_idx] <- as.numeric(stats::predict(
      cv_e, newx = B[test_idx, , drop = FALSE], s = lam_e, type = "response"
    ))
    lambda_e[[k]] <- lam_e
    nz_e[[k]] <- sum(as.numeric(stats::coef(cv_e, s = lam_e))[-1] != 0)

    # --- Control Outcome m0(X) ---
    # Train ONLY on control units in training folds
    control_train_idx <- train_idx[A[train_idx] == 0]

    if (length(control_train_idx) < 10) {
      warning(sprintf("Fold %d: Only %d control units in training set",
                      k, length(control_train_idx)))
    }

    fam_m0 <- if (outcome_type == "binary") "binomial" else "gaussian"
    cv_m0 <- glmnet::cv.glmnet(
      x = B[control_train_idx, , drop = FALSE], y = Y[control_train_idx],
      family = fam_m0, alpha = 1, nfolds = cv_nfolds
    )
    lam_m0 <- if (lambda_rule == "min") cv_m0$lambda.min else cv_m0$lambda.1se
    m0_hat[test_idx] <- as.numeric(stats::predict(
      cv_m0, newx = B[test_idx, , drop = FALSE], s = lam_m0, type = "response"
    ))
    lambda_m0[[k]] <- lam_m0
    nz_m0[[k]] <- sum(as.numeric(stats::coef(cv_m0, s = lam_m0))[-1] != 0)
  }

  # Clip propensity scores to avoid extreme weights (att_linear()'s convention)
  e_hat <- pmax(pmin(e_hat, 0.99), 0.01)

  # For binary outcomes, clip m0_hat probabilities
  if (outcome_type == "binary") {
    m0_hat <- pmax(pmin(m0_hat, 0.99), 0.01)
  }
  # For continuous outcomes, no clipping needed (m0_hat is on natural scale)

  # --- Orthogonal Score and Point Estimate ---
  # ATT orthogonal score (manuscript.tex eq:score; Chernozhukov et al. 2018):
  # psi_i(theta) = (A_i/pi)*(Y_i - m0_i - theta) - (1/pi)*(e_i*(1 - A_i)/(1 - e_i))*(Y_i - m0_i)

  pi_hat <- mean(A)

  # Solve for theta: sum(psi_i(theta)) = 0
  score_at_zero_term1 <- (A / pi_hat) * (Y - m0_hat)
  score_at_zero_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)
  score_at_zero_term2[!is.finite(score_at_zero_term2)] <- 0  # Safety

  score_at_zero <- score_at_zero_term1 - score_at_zero_term2
  sum_a_over_pi <- sum(A / pi_hat)

  theta <- sum(score_at_zero) / sum_a_over_pi

  # --- Variance Estimation ---
  # Score at theta-hat (proper Neyman-orthogonal score)
  score_term1 <- (A / pi_hat) * (Y - m0_hat - theta)
  score_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)
  score_term2[!is.finite(score_term2)] <- 0

  score_values <- score_term1 - score_term2

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
    method = "lasso_dml",
    hyperparams = list(
      K = K,
      cv_nfolds = cv_nfolds,
      lambda_rule = lambda_rule,
      interactions = interactions,
      n_basis = ncol(B),
      lambda_e = lambda_e,
      lambda_m0 = lambda_m0,
      nz_e = nz_e,
      nz_m0 = nz_m0
    )
  )

  if (verbose) {
    cat(sprintf("\nLasso-DML ATT Results:\n"))
    cat(sprintf("  ATT estimate: %.4f\n", theta))
    cat(sprintf("  Std error:    %.4f\n", sigma))
    cat(sprintf("  95%% CI:       [%.4f, %.4f]\n", ci[1], ci[2]))
    cat(sprintf("  mean nonzero: e = %.1f, m0 = %.1f of %d basis terms\n",
                mean(nz_e), mean(nz_m0), ncol(B)))
  }

  return(result)
}

# Alias for consistency with simulation scripts
method_lasso_dml <- att_lasso_dml
