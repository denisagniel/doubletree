# Instrumented version of estimate_att() to find hang location

dml_att_instrumented <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                   regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                   stratified = TRUE, seed = NULL, verbose = FALSE,
                   use_rashomon = FALSE, rashomon_bound_multiplier = 0.05,
                   rashomon_bound_adder = 0, max_leaves = NULL,
                   auto_tune_intersecting = FALSE,
                   discretize_method = "quantiles",
                   discretize_bins = "adaptive",
                   ...) {

  cat("[INSTRUMENT] Starting dml_att...\n"); flush.console()

  outcome_type <- match.arg(outcome_type)
  check_dml_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  cat("[INSTRUMENT] Validated data, n=", n, "\n"); flush.console()

  # Validate parameters
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("K must be a single integer >= 2, got: ", K, call. = FALSE)
  }

  n_treated <- sum(A == 1)
  n_control <- sum(A == 0)
  if (n_treated < K) {
    stop("Insufficient treated units for K-fold cross-fitting. ",
         "Need at least K=", K, " treated units, got: ", n_treated, call. = FALSE)
  }
  if (n_control < K) {
    stop("Insufficient control units for K-fold cross-fitting. ",
         "Need at least K=", K, " control units, got: ", n_control, call. = FALSE)
  }

  if (!cv_regularization && (!is.numeric(regularization) || length(regularization) != 1 || regularization <= 0)) {
    stop("regularization must be a single positive numeric value, got: ",
         regularization, call. = FALSE)
  }

  if (cv_regularization && (!is.numeric(cv_K) || length(cv_K) != 1 || cv_K < 2)) {
    stop("cv_K must be a single integer >= 2 when cv_regularization = TRUE, got: ",
         cv_K, call. = FALSE)
  }

  if (!is.numeric(rashomon_bound_multiplier) || length(rashomon_bound_multiplier) != 1 || rashomon_bound_multiplier < 0) {
    stop("rashomon_bound_multiplier must be a single non-negative numeric value, got: ",
         rashomon_bound_multiplier, call. = FALSE)
  }

  cat("[INSTRUMENT] Creating folds...\n"); flush.console()
  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)
  cat("[INSTRUMENT] Folds created\n"); flush.console()

  if (use_rashomon) {
    cat("[INSTRUMENT] Fitting Rashomon nuisances...\n"); flush.console()
    nuisance_fits <- fit_nuisances_rashomon(X, A, Y, fold_indices, outcome_type = outcome_type,
                                           regularization = regularization,
                                           cv_regularization = cv_regularization, cv_K = cv_K,
                                           verbose = verbose,
                                           rashomon_bound_multiplier = rashomon_bound_multiplier,
                                           rashomon_bound_adder = rashomon_bound_adder,
                                           max_leaves = max_leaves,
                                           auto_tune_intersecting = auto_tune_intersecting, ...)
    cat("[INSTRUMENT] Rashomon nuisances fitted\n"); flush.console()
    cat("[INSTRUMENT] Getting fold-specific eta (Rashomon)...\n"); flush.console()
    eta <- get_fold_specific_eta_rashomon(nuisance_fits, X, fold_indices)
    cat("[INSTRUMENT] Eta computed (Rashomon)\n"); flush.console()
  } else {
    cat("[INSTRUMENT] Fitting fold-specific nuisances (K=", K, ")...\n"); flush.console()
    nuisance_fits <- vector("list", K)
    for (k in seq_len(K)) {
      cat("[INSTRUMENT]   Fold ", k, "/", K, "...\n"); flush.console()
      nuisance_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices,
                                              outcome_type = outcome_type,
                                              regularization = regularization,
                                              cv_regularization = cv_regularization, cv_K = cv_K,
                                              verbose = verbose, ...)
      cat("[INSTRUMENT]   Fold ", k, " complete\n"); flush.console()
    }
    cat("[INSTRUMENT] All folds fitted\n"); flush.console()

    cat("[INSTRUMENT] Getting fold-specific eta...\n"); flush.console()
    eta <- get_fold_specific_eta(nuisance_fits, X, fold_indices)
    cat("[INSTRUMENT] Eta computed\n"); flush.console()
  }

  cat("[INSTRUMENT] Computing pi_hat...\n"); flush.console()
  pi_hat <- mean(A)
  cat("[INSTRUMENT] pi_hat=", pi_hat, "\n"); flush.console()

  cat("[INSTRUMENT] Computing sum_a_over_pi...\n"); flush.console()
  sum_a_over_pi <- sum(A / pi_hat)
  cat("[INSTRUMENT] sum_a_over_pi=", sum_a_over_pi, "\n"); flush.console()

  if (sum_a_over_pi < 1e-10) {
    stop("No treated units (sum(A) ~ 0) or pi_hat extremely small.", call. = FALSE)
  }

  cat("[INSTRUMENT] Computing score_at_zero...\n"); flush.console()
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  cat("[INSTRUMENT] score_at_zero computed, length=", length(score_at_zero), "\n"); flush.console()

  cat("[INSTRUMENT] Computing theta...\n"); flush.console()
  theta <- sum(score_at_zero) / sum_a_over_pi
  cat("[INSTRUMENT] theta=", theta, "\n"); flush.console()

  cat("[INSTRUMENT] Computing score_values...\n"); flush.console()
  score_values <- psi_att(Y, A, theta, eta, pi_hat)
  cat("[INSTRUMENT] score_values computed\n"); flush.console()

  cat("[INSTRUMENT] Computing variance...\n"); flush.console()
  sigma_sq <- dml_att_variance(score_values, n)
  cat("[INSTRUMENT] sigma_sq=", sigma_sq, "\n"); flush.console()

  sigma <- sqrt(sigma_sq)
  cat("[INSTRUMENT] sigma=", sigma, "\n"); flush.console()

  cat("[INSTRUMENT] Computing CI...\n"); flush.console()
  ci_95 <- dml_att_ci(theta, sigma, n, level = 0.95)
  cat("[INSTRUMENT] CI computed\n"); flush.console()

  cat("[INSTRUMENT] Adding predictions to nuisance_fits...\n"); flush.console()
  nuisance_fits$propensity <- eta$e
  nuisance_fits$outcome_control <- eta$m0
  cat("[INSTRUMENT] Predictions added\n"); flush.console()

  cat("[INSTRUMENT] Returning result...\n"); flush.console()
  list(
    theta = theta,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    nuisance_fits = nuisance_fits,
    fold_indices = fold_indices,
    n = n,
    K = K
  )
}
