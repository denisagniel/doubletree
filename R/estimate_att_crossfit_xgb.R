#' Cross-fitted ATT with gradient-boosted nuisances (fast anchor)
#'
#' @description
#' K-fold cross-fitted AIPW/EIF estimator of the ATT whose nuisances
#' \eqn{e(x)} and \eqn{m_0(x)} are fitted by \pkg{xgboost}. It is the fast,
#' flexible \emph{anchor} for the tree estimators: no optimal-tree search, no
#' structure selection, and a fixed, pre-specified boosting configuration with
#' early stopping on a validation split carved out of each fold's TRAINING
#' rows only (never the held-out fold), so each held-out prediction is
#' out-of-fold.
#'
#' Validity of the Wald interval rests on the usual cross-fitting product-rate
#' condition, \eqn{\|\hat e - e\|\,\|\hat m_0 - m_0\| = o_p(n^{-1/2})} for the
#' boosted fits; xgboost does not guarantee it, so it is an assumption of this
#' anchor, not a property the function verifies.
#'
#' The estimator shares \code{\link{create_folds}} and the EIF solve
#' (\code{eif_att_solve}) with \code{\link{estimate_att_crossfit}}. It clips
#' the propensity one-sidedly (\eqn{\hat e \le 1 - c}), as
#' \code{\link{estimate_att}} does. NOTE that
#' \code{estimate_att_crossfit()} clamps two-sidedly to [0.01, 0.99], so the
#' two cross-fit anchors differ in clipping as well as in the learner.
#'
#' Early stopping uses an explicit metric (\code{logloss} for the propensity and
#' a binary outcome, \code{rmse} otherwise), not the version-dependent default.
#' The caller's \code{.Random.seed} is restored on exit, so calling this inside a
#' simulation loop does not reset the loop's random stream. Validated against
#' xgboost 1.7 and 2.x; xgboost >= 3.0 changes the R interface and is refused
#' until it has been validated.
#'
#' @param X Covariates (data frame or matrix, numeric/binary; no NA).
#' @param A Binary treatment vector (0/1).
#' @param Y Outcome. Numeric for \code{outcome_type = "continuous"}, 0/1 for
#'   \code{"binary"}.
#' @param K Number of folds (default 5).
#' @param outcome_type \code{"continuous"} (default) or \code{"binary"}.
#' @param clip_c One-sided propensity clip constant (default 0.01).
#' @param stratified Stratify folds on \code{A} (default \code{TRUE}).
#' @param seed Integer seed for fold assignment, validation splits and
#'   xgboost. Required (no default): the anchor must be reproducible.
#' @param params Named list overriding the default xgboost parameters
#'   (\code{max_depth = 3, eta = 0.1, subsample = 0.8, colsample_bytree = 1,
#'   min_child_weight = 1, lambda = 1}). \code{objective}, \code{nthread},
#'   \code{seed} and \code{eval_metric} are set by this function and may not
#'   be overridden.
#' @param nrounds_max Maximum boosting rounds (default 500).
#' @param early_stopping_rounds Patience for early stopping (default 20).
#' @param valid_frac Fraction of each fold's training rows (controls only for
#'   \eqn{m_0}) held out for early stopping (default 0.2).
#' @param nthread xgboost threads (default 1: more threads can change results
#'   at the last digits, so reproducibility is the default).
#' @param verbose Logical; print xgboost's per-round evaluation log.
#' @return A list shaped like \code{\link{estimate_att_crossfit}}'s:
#'   \code{theta}, \code{sigma}, \code{ci_95}, \code{ci_95_wald},
#'   \code{score_values}, \code{nuisance_fits} (per-fold best iterations plus
#'   \code{propensity}, \code{outcome_control}), \code{fold_indices}, \code{n},
#'   \code{K}; plus \code{learner = "xgboost"}, \code{clip_c},
#'   \code{clip_rate_e}, \code{xgb_params} and \code{xgboost_version}.
#' @seealso \code{\link{estimate_att_crossfit}}, \code{\link{estimate_att}}
#' @export
estimate_att_crossfit_xgb <- function(X, A, Y, K = 5,
                                      outcome_type = c("continuous", "binary"),
                                      clip_c = 0.01, stratified = TRUE,
                                      seed, params = list(),
                                      nrounds_max = 500L,
                                      early_stopping_rounds = 20L,
                                      valid_frac = 0.2, nthread = 1L,
                                      verbose = FALSE) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("estimate_att_crossfit_xgb() needs the 'xgboost' package.", call. = FALSE)
  }
  outcome_type <- match.arg(outcome_type)
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max - 1100L) {
    stop("`seed` is required: a single integer (the anchor must be reproducible).",
         call. = FALSE)
  }
  xgb_version <- utils::packageVersion("xgboost")
  if (xgb_version >= "3.0.0") {
    stop("estimate_att_crossfit_xgb() is validated on xgboost < 3.0 (the 3.x R ",
         "interface changed best_iteration and iterationrange handling); installed: ",
         as.character(xgb_version), ".", call. = FALSE)
  }
  # Restore the caller's RNG stream on exit: this function reseeds internally.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
  }
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    stop("K must be a single integer >= 2, got: ", K, call. = FALSE)
  }
  if (!is.numeric(clip_c) || length(clip_c) != 1 || clip_c <= 0 || clip_c >= 0.5) {
    stop("clip_c must be a single number in (0, 0.5), got: ", clip_c, call. = FALSE)
  }
  if (!is.numeric(valid_frac) || length(valid_frac) != 1 ||
      valid_frac <= 0 || valid_frac >= 0.5) {
    stop("valid_frac must be a single number in (0, 0.5), got: ", valid_frac,
         call. = FALSE)
  }
  Xm <- as.matrix(X)
  storage.mode(Xm) <- "double"
  A <- as.integer(A)
  Y <- as.numeric(Y)
  n <- nrow(Xm)
  min_per_fold <- 2L
  if (sum(A == 1L) < K * min_per_fold || sum(A == 0L) < K * min_per_fold) {
    stop("Too few treated or control units for K = ", K, " fold cross-fitting.",
         call. = FALSE)
  }

  reserved <- c("objective", "nthread", "seed", "eval_metric")
  if (any(names(params) %in% reserved)) {
    stop("`params` may not set: ", paste(intersect(names(params), reserved), collapse = ", "),
         "; use the dedicated arguments.", call. = FALSE)
  }
  xgb_params <- utils::modifyList(
    list(max_depth = 3L, eta = 0.1, subsample = 0.8, colsample_bytree = 1,
         min_child_weight = 1, lambda = 1),
    params
  )

  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  e_hat <- numeric(n)
  m0_hat <- numeric(n)
  best_iter <- data.frame(fold = seq_len(K), e = NA_integer_, m0 = NA_integer_)

  for (k in seq_len(K)) {
    test_idx <- which(fold_indices == k)
    train_idx <- which(fold_indices != k)

    fit_e <- .xgb_fit_one(
      Xm[train_idx, , drop = FALSE], A[train_idx],
      objective = "binary:logistic", params = xgb_params,
      nrounds_max = nrounds_max, early_stopping_rounds = early_stopping_rounds,
      valid_frac = valid_frac, nthread = nthread,
      seed = seed + k, strata = A[train_idx], verbose = verbose
    )
    e_hat[test_idx] <- .xgb_predict(fit_e, Xm[test_idx, , drop = FALSE])
    best_iter$e[k] <- fit_e$best_iter

    ctrl <- train_idx[A[train_idx] == 0L]
    fit_m <- .xgb_fit_one(
      Xm[ctrl, , drop = FALSE], Y[ctrl],
      objective = if (outcome_type == "binary") "binary:logistic" else "reg:squarederror",
      params = xgb_params, nrounds_max = nrounds_max,
      early_stopping_rounds = early_stopping_rounds,
      valid_frac = valid_frac, nthread = nthread,
      seed = seed + 1000L + k, strata = NULL, verbose = verbose
    )
    m0_hat[test_idx] <- .xgb_predict(fit_m, Xm[test_idx, , drop = FALSE])
    best_iter$m0[k] <- fit_m$best_iter
  }

  hit_nrounds_max <- any(best_iter$e >= nrounds_max | best_iter$m0 >= nrounds_max)
  if (hit_nrounds_max) {
    warning("early stopping never fired for at least one fold's nuisance ",
            "(best iteration reached nrounds_max = ", nrounds_max, "); raise ",
            "`nrounds_max` or `eta`.", call. = FALSE)
  }

  # One-sided clip, as estimate_att(): no lower bound on e (theory.tex
  # ass:construct(c)); only e <= 1 - clip_c.
  clip_rate_e <- mean(e_hat > 1 - clip_c)
  e_hat <- pmin(1 - clip_c, e_hat)

  .att <- eif_att_solve(Y, A, e_hat, m0_hat, n)

  list(
    theta = .att$theta,
    sigma = .att$sigma,
    ci_95 = .att$ci_95,
    ci_95_wald = .att$ci_95,
    score_values = .att$score_values,
    nuisance_fits = list(best_iterations = best_iter,
                         propensity = e_hat, outcome_control = m0_hat),
    fold_indices = fold_indices,
    n = n,
    K = K,
    learner = "xgboost",
    clip_c = clip_c,
    clip_rate_e = clip_rate_e,
    xgb_params = xgb_params,
    hit_nrounds_max = hit_nrounds_max,
    xgboost_version = as.character(xgb_version)
  )
}

#' Fit one xgboost model with early stopping on a split of its own rows
#' @param strata Optional vector to stratify the validation split on.
#' @return list(model, best_iter)
#' @noRd
.xgb_fit_one <- function(Xtr, ytr, objective, params, nrounds_max,
                         early_stopping_rounds, valid_frac, nthread, seed,
                         strata = NULL, verbose = FALSE) {
  m <- nrow(Xtr)
  if (m < 20L) {
    stop("Too few training rows (", m, ") to fit and validate an xgboost ",
         "nuisance; reduce K.", call. = FALSE)
  }
  set.seed(seed)
  if (is.null(strata)) {
    is_val <- seq_len(m) %in% sample.int(m, max(2L, round(valid_frac * m)))
  } else {
    is_val <- logical(m)
    for (s in unique(strata)) {
      idx <- which(strata == s)
      nv <- max(1L, round(valid_frac * length(idx)))
      is_val[idx[sample.int(length(idx), nv)]] <- TRUE
    }
  }
  dtrain <- xgboost::xgb.DMatrix(Xtr[!is_val, , drop = FALSE], label = ytr[!is_val])
  dval <- xgboost::xgb.DMatrix(Xtr[is_val, , drop = FALSE], label = ytr[is_val])
  # xgboost's R interface ignores a `seed` parameter (it warns): randomness
  # (subsampling) is driven by R's RNG, so seed it explicitly right before
  # training.
  p <- c(params, list(objective = objective, nthread = nthread,
                      eval_metric = if (objective == "binary:logistic") "logloss" else "rmse"))
  set.seed(seed)
  # `watchlist` is named `evals` in newer xgb.train() signatures.
  eval_arg <- if ("evals" %in% names(formals(xgboost::xgb.train))) "evals" else "watchlist"
  args <- list(params = p, data = dtrain, nrounds = nrounds_max,
               early_stopping_rounds = early_stopping_rounds,
               verbose = if (verbose) 1L else 0L)
  args[[eval_arg]] <- list(val = dval)
  model <- do.call(xgboost::xgb.train, args)
  best <- model$best_iteration
  if (is.null(best)) {
    stop("xgboost returned no best_iteration; early stopping did not run.", call. = FALSE)
  }
  list(model = model, best_iter = as.integer(best))
}

#' Predict from .xgb_fit_one() at the early-stopped iteration
#' @noRd
.xgb_predict <- function(fit, Xnew) {
  as.numeric(stats::predict(fit$model, xgboost::xgb.DMatrix(Xnew),
                            iterationrange = c(1L, fit$best_iter + 1L)))
}
