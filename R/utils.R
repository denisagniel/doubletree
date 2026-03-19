#' Create cross-fitting fold indices
#'
#' Assigns each observation to one of K folds. Optionally stratifies by a
#' binary variable (e.g. treatment A) so that fold proportions are balanced.
#'
#' @param n Integer. Number of observations.
#' @param K Integer. Number of folds (>= 2).
#' @param strata Optional integer or numeric vector of length n (e.g. treatment).
#'   If provided, folds are stratified by strata so each fold has roughly
#'   the same proportion of each stratum.
#' @param seed Optional integer. Random seed for reproducibility.
#' @return Integer vector of length n with values in 1..K.
#' @examples
#' # Basic usage: 5-fold cross-fitting
#' folds <- create_folds(n = 100, K = 5, seed = 42)
#' table(folds)  # Each fold has ~20 observations
#'
#' # Stratified by treatment
#' A <- rbinom(100, 1, 0.3)
#' folds_stratified <- create_folds(n = 100, K = 5, strata = A, seed = 42)
#' # Check balance: each fold has similar treatment proportions
#' tapply(A, folds_stratified, mean)
#'
#' @export
create_folds <- function(n, K, strata = NULL, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("n must be a positive integer")
  }
  n <- as.integer(n)
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("K must be an integer >= 2")
  }
  K <- as.integer(min(K, n))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(strata) || length(strata) != n) {
    perm <- sample.int(n)
    return(as.integer(cut(perm, breaks = K, labels = seq_len(K))))
  }
  # Stratified: within each stratum assign to folds 1..K
  fold_vec <- integer(n)
  u <- unique(strata)
  for (s in u) {
    idx <- which(strata == s)
    n_s <- length(idx)
    if (n_s == 0) next
    perm <- sample.int(n_s)

    actual_K <- min(K, n_s)
    if (n_s < K) {
      warning("Stratum '", s, "' has only ", n_s, " observations (less than K=", K, "). ",
              "Using ", actual_K, " folds for this stratum instead. ",
              "This may affect cross-validation balance.",
              call. = FALSE, immediate. = TRUE)
    }

    fold_vec[idx] <- as.integer(cut(perm, breaks = actual_K, labels = seq_len(actual_K)))
  }
  fold_vec
}

#' Check (X, A, Y) inputs for ATT estimation
#'
#' Validates presence, dimensions, binary A, and Y according to outcome_type.
#' Called internally by estimate_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1 for optimaltrees).
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Integer or numeric vector of outcome; must be binary (0/1) if outcome_type is "binary", numeric if "continuous".
#' @param outcome_type Character. "binary" (default) or "continuous".
#' @param check_overlap Logical. If TRUE, warn when propensity range suggests weak overlap.
#' @return Invisible NULL; stops on error.
#' @noRd
check_att_data <- function(X, A, Y, outcome_type = "binary", check_overlap = TRUE) {
  if (is.null(X) || is.null(A) || is.null(Y)) {
    stop("X, A, and Y are required")
  }
  n <- NROW(X)
  if (length(A) != n || length(Y) != n) {
    stop("Length of A and Y must equal nrow(X)")
  }
  if (any(is.na(A)) || any(is.na(Y)) || any(is.na(X))) {
    stop("NA not allowed in X, A, or Y")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary (0/1)")
  }
  if (outcome_type == "binary") {
    if (!all(Y %in% c(0, 1))) {
      stop("Y must be binary (0/1) when outcome_type is \"binary\".")
    }
  } else if (outcome_type == "continuous") {
    if (!is.numeric(Y)) {
      stop("Y must be numeric when outcome_type is \"continuous\".")
    }
  } else {
    stop("outcome_type must be \"binary\" or \"continuous\".")
  }
  if (check_overlap && length(unique(A)) == 2L) {
    p1 <- mean(A)
    if (p1 < 0.02 || p1 > 0.98) {
      warning("Weak overlap detected: treatment proportion P(A=1) = ", round(p1, 3),
              "\nWhen < 2% or > 98%, ATT estimates may be unstable due to extrapolation.",
              "\nRecommendations:",
              "\n  - Collect more data to improve balance",
              "\n  - Consider trimming extreme propensity scores",
              "\n  - Report limited generalizability in results",
              "\n  - Check overlap plots before trusting estimates",
              call. = FALSE, immediate. = TRUE)
    }
  }
  invisible(NULL)
}
