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
    fold_vec[idx] <- as.integer(cut(perm, breaks = min(K, n_s), labels = seq_len(min(K, n_s))))
  }
  fold_vec
}

#' Check (X, A, Y) inputs for DML ATT
#'
#' Validates presence, dimensions, binary A, and Y according to outcome_type.
#' Called internally by dml_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1 for treefarmr).
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Integer or numeric vector of outcome; must be binary (0/1) if outcome_type is "binary", numeric if "continuous".
#' @param outcome_type Character. "binary" (default) or "continuous".
#' @param check_overlap Logical. If TRUE, warn when propensity range suggests weak overlap.
#' @return Invisible NULL; stops on error.
#' @noRd
check_dml_att_data <- function(X, A, Y, outcome_type = "binary", check_overlap = TRUE) {
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
      warning("Very few treated or controls; overlap may be weak.")
    }
  }
  invisible(NULL)
}
