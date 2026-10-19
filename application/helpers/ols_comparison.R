## ============================================================================
## application/helpers/ols_comparison.R
##
## Parametric comparators for the tree estimators: logistic e(X) and OLS m0(X)
## (main effects, or main effects + all pairwise interactions), cross-fitted on
## the SAME folds as the anchor, plus the shared EIF solve and out-of-fold
## nuisance metrics. Pure functions; 07b_compare_ols.R wires them to the data.
##
## NOT simulations/methods/method_linear.R::att_linear(): that one clips
## two-sided, swallows glm/lm warnings and zeroes non-finite scores. Here the
## clip is one-sided (as estimate_att()), separation and non-convergence are
## errors, and aliased OLS coefficients are COUNTED and reported. Setting an
## aliased coefficient to 0 makes predictions independent of which column lm()
## dropped ONLY for rows inside the training column space; a pairwise product
## that is identically 0 (or collinear) among the training controls but varies
## on a held-out or treated row gets coefficient 0 there. That is a property of
## the data, so the aliased count is reported alongside every OLS fit.
## ============================================================================

#' Model matrix for the parametric comparators.
#' @param X Data frame of 0/1 covariates.
#' @param interactions Logical: add all pairwise products.
#' @return Numeric matrix, no intercept column.
ols_design <- function(X, interactions = FALSE) {
  M <- as.matrix(X)
  storage.mode(M) <- "double"
  if (!interactions) return(M)
  p <- ncol(M)
  idx <- utils::combn(p, 2L)
  P <- M[, idx[1, ], drop = FALSE] * M[, idx[2, ], drop = FALSE]
  colnames(P) <- paste(colnames(M)[idx[1, ]], colnames(M)[idx[2, ]], sep = ":")
  cbind(M, P)
}

#' Fit OLS on `train` rows and predict `test` rows.
#' @return list(pred, n_aliased, n_coef)
ols_fit_predict <- function(M, y, train, test) {
  d <- data.frame(y = y[train], M[train, , drop = FALSE], check.names = FALSE)
  fit <- stats::lm(y ~ ., data = d)
  aliased <- is.na(stats::coef(fit))
  b <- stats::coef(fit)
  b[aliased] <- 0
  pred <- as.numeric(cbind(1, M[test, , drop = FALSE]) %*% b)
  list(pred = pred, n_aliased = sum(aliased), n_coef = length(b))
}

#' Fit logistic regression on `train` rows and predict `test` rows.
#' Errors (no silent fallback) on non-convergence or (quasi-)separation.
logit_fit_predict <- function(M, a, train, test) {
  d <- data.frame(a = a[train], M[train, , drop = FALSE], check.names = FALSE)
  fit <- withCallingHandlers(
    stats::glm(a ~ ., data = d, family = stats::binomial()),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1|did not converge", conditionMessage(w))) {
        stop("logistic propensity: ", conditionMessage(w), call. = FALSE)
      }
    }
  )
  if (!fit$converged) stop("logistic propensity did not converge.", call. = FALSE)
  ## glm() only warns when fitted probabilities hit machine zero/one, which a
  ## separated fit can sit just short of; check explicitly.
  if (any(stats::fitted(fit) < 1e-6 | stats::fitted(fit) > 1 - 1e-6)) {
    stop("logistic propensity: (quasi-)separation -- training fitted ",
         "probabilities within 1e-6 of 0 or 1.", call. = FALSE)
  }
  b <- stats::coef(fit)
  b[is.na(b)] <- 0
  as.numeric(stats::plogis(cbind(1, M[test, , drop = FALSE]) %*% b))
}

#' Cross-fitted parametric nuisances on the given folds.
#' @param fold_indices Integer fold id per row (1..K).
#' @return list(e, m0_main, m0_inter, n_aliased) with out-of-fold predictions.
crossfit_parametric <- function(X, A, Y, fold_indices) {
  n <- nrow(X)
  M1 <- ols_design(X, FALSE)
  M2 <- ols_design(X, TRUE)
  out <- list(e = numeric(n), m0_main = numeric(n), m0_inter = numeric(n))
  aliased <- data.frame(fold = sort(unique(fold_indices)), main = NA_integer_, inter = NA_integer_)
  for (k in aliased$fold) {
    te <- which(fold_indices == k)
    tr <- which(fold_indices != k)
    trc <- tr[A[tr] == 0L]
    out$e[te] <- logit_fit_predict(M1, A, tr, te)
    f1 <- ols_fit_predict(M1, Y, trc, te)
    f2 <- ols_fit_predict(M2, Y, trc, te)
    out$m0_main[te] <- f1$pred
    out$m0_inter[te] <- f2$pred
    aliased$main[aliased$fold == k] <- f1$n_aliased
    aliased$inter[aliased$fold == k] <- f2$n_aliased
  }
  out$n_aliased <- aliased
  out
}

#' Full-sample parametric nuisances (matched to the flagship's no-splitting scheme).
full_sample_parametric <- function(X, A, Y) {
  n <- nrow(X)
  all_rows <- seq_len(n)
  ctrl <- which(A == 0L)
  M1 <- ols_design(X, FALSE)
  M2 <- ols_design(X, TRUE)
  f1 <- ols_fit_predict(M1, Y, ctrl, all_rows)
  f2 <- ols_fit_predict(M2, Y, ctrl, all_rows)
  list(e = logit_fit_predict(M1, A, all_rows, all_rows),
       m0_main = f1$pred, m0_inter = f2$pred,
       n_aliased = c(main = f1$n_aliased, inter = f2$n_aliased))
}

#' ATT by the shared EIF solve, with the flagship's one-sided propensity clip.
att_from_nuisances <- function(Y, A, e, m0, clip_c) {
  fit <- doubletree:::eif_att_solve(Y, A, pmin(1 - clip_c, e), m0, length(Y))
  data.frame(theta = fit$theta, sigma = fit$sigma,
             ci_lo = fit$ci_95[[1]], ci_hi = fit$ci_95[[2]])
}

#' Out-of-fold nuisance metrics with fold-wise SEs.
#'
#' m0: squared and absolute error on CONTROL rows. e: Brier score and log-loss
#' on all rows. Each metric is computed per fold, then averaged; the SE is the
#' standard error of the K fold values.
#' @param fold_indices Fold id per row.
#' @param e_list,m0_list Named lists of out-of-fold prediction vectors.
nuisance_metrics <- function(A, Y, fold_indices, e_list, m0_list, clip_c) {
  folds <- sort(unique(fold_indices))
  summarise <- function(vals) c(mean = mean(vals), se = stats::sd(vals) / sqrt(length(vals)))
  rows <- list()
  for (nm in names(m0_list)) {
    mse <- mae <- numeric(length(folds))
    for (i in seq_along(folds)) {
      r <- which(fold_indices == folds[i] & A == 0L)
      err <- Y[r] - m0_list[[nm]][r]
      mse[i] <- mean(err^2); mae[i] <- mean(abs(err))
    }
    rows[[length(rows) + 1L]] <- data.frame(nuisance = "m0", learner = nm, metric = "mse",
                                           t(summarise(mse)))
    rows[[length(rows) + 1L]] <- data.frame(nuisance = "m0", learner = nm, metric = "mae",
                                           t(summarise(mae)))
  }
  for (nm in names(e_list)) {
    p <- pmin(1 - clip_c, pmax(clip_c, e_list[[nm]]))
    brier <- logloss <- numeric(length(folds))
    for (i in seq_along(folds)) {
      r <- which(fold_indices == folds[i])
      brier[i] <- mean((A[r] - p[r])^2)
      logloss[i] <- -mean(A[r] * log(p[r]) + (1 - A[r]) * log(1 - p[r]))
    }
    rows[[length(rows) + 1L]] <- data.frame(nuisance = "e", learner = nm, metric = "brier",
                                           t(summarise(brier)))
    rows[[length(rows) + 1L]] <- data.frame(nuisance = "e", learner = nm, metric = "logloss",
                                           t(summarise(logloss)))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Per-fold held-out MSE of m0 on CONTROL rows, for a set of learners.
#' @param m0_list Named list of out-of-fold prediction vectors (length n).
#' @return data.frame(rep, fold, learner, mse, n_control).
fold_mse_m0 <- function(A, Y, fold_indices, m0_list, rep_id = 1L) {
  folds <- sort(unique(fold_indices))
  for (k in folds) {
    if (!any(fold_indices == k & A == 0L)) {
      stop("fold_mse_m0: fold ", k, " has no control rows.", call. = FALSE)
    }
  }
  do.call(rbind, lapply(names(m0_list), function(nm) {
    do.call(rbind, lapply(folds, function(k) {
      r <- which(fold_indices == k & A == 0L)
      data.frame(rep = rep_id, fold = k, learner = nm,
                 mse = mean((Y[r] - m0_list[[nm]][r])^2), n_control = length(r))
    }))
  }))
}

#' Paired comparison of a reference learner (the tree) with the others.
#'
#' For each other learner: mean held-out MSE of both, the paired difference
#' (reference - other; negative favours the reference; `diff_rep1` is the
#' repetition-1 mean that `diff_se_rep1` goes with), and a skill score
#' 1 - MSE / MSE(constant baseline). The paired SE uses the K fold differences
#' of repetition 1 only (the folds of one repetition are the independent-ish
#' units; folds from different repetitions reuse the same patients, so pooling
#' them would understate the SE). Means are over all repetitions.
#' @param fm Output of rbind()-ing fold_mse_m0() over repetitions; must include
#'   the `baseline` learner.
paired_mse_summary <- function(fm, reference = "tree_cv", baseline = "constant") {
  if (!all(c(reference, baseline) %in% fm$learner)) {
    stop("paired_mse_summary: `reference` and `baseline` must both be in fm$learner.",
         call. = FALSE)
  }
  mean_mse <- tapply(fm$mse, fm$learner, mean)
  others <- setdiff(unique(fm$learner), c(reference, baseline))
  if (length(others) == 0L) stop("paired_mse_summary: no learners to compare.", call. = FALSE)
  rep1 <- fm[fm$rep == min(fm$rep), ]
  do.call(rbind, lapply(others, function(o) {
    a <- rep1[rep1$learner == reference, ]
    b <- rep1[rep1$learner == o, ]
    a <- a[order(a$fold), ]; b <- b[order(b$fold), ]
    d <- a$mse - b$mse
    data.frame(
      reference = reference, other = o,
      mse_reference = unname(mean_mse[reference]), mse_other = unname(mean_mse[o]),
      mse_baseline = unname(mean_mse[baseline]),
      diff_ref_minus_other = unname(mean_mse[reference] - mean_mse[o]),
      diff_rep1 = mean(d),
      diff_se_rep1 = stats::sd(d) / sqrt(length(d)),
      skill_reference = 1 - unname(mean_mse[reference] / mean_mse[baseline]),
      skill_other = 1 - unname(mean_mse[o] / mean_mse[baseline]),
      n_reps = length(unique(fm$rep)), n_folds = length(unique(fm$fold))
    )
  }))
}
