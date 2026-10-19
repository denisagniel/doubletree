## ============================================================================
## application/tests/test-ols-comparison.R
## Parametric comparators: hand-checkable cases on tiny data.
## ============================================================================

source(file.path(app_dir_for_tests, "helpers", "ols_comparison.R"))

test_that("ols_design adds all pairwise products with stable names", {
  X <- data.frame(a = c(1, 0, 1), b = c(1, 1, 0), c = c(0, 1, 1))
  M <- ols_design(X, interactions = TRUE)
  expect_equal(ncol(M), 3L + choose(3L, 2L))
  expect_equal(M[, "a:b"], c(1, 0, 0))
  expect_equal(M[, "b:c"], c(0, 1, 0))
  expect_equal(ncol(ols_design(X, FALSE)), 3L)
})

test_that("ols_fit_predict reproduces an exact linear signal and counts aliased columns", {
  set.seed(1)
  n <- 200
  X <- data.frame(x1 = rbinom(n, 1, .5), x2 = rbinom(n, 1, .5))
  X$x3 <- X$x1                     # perfectly collinear with x1
  y <- 3 + 2 * X$x1 - 1 * X$x2
  M <- ols_design(X)
  r <- ols_fit_predict(M, y, train = 1:150, test = 151:200)
  expect_equal(r$pred, y[151:200], tolerance = 1e-8)
  expect_equal(r$n_aliased, 1L)
})

test_that("logit_fit_predict errors on perfect separation instead of returning garbage", {
  X <- data.frame(x = c(rep(0, 30), rep(1, 30)))
  a <- c(rep(0L, 30), rep(1L, 30))
  M <- ols_design(X)
  expect_error(logit_fit_predict(M, a, seq_len(60), seq_len(60)), "logistic propensity")
})

test_that("nuisance_metrics matches hand-computed fold-wise means", {
  ## Two folds, controls only for m0. Fold 1 controls: errors 1, 3 -> mse 5, mae 2.
  ## Fold 2 controls: errors 2, 2 -> mse 4, mae 2.
  A <- c(0, 0, 1, 0, 0, 1)
  fold <- c(1, 1, 1, 2, 2, 2)
  Y <- c(11, 13, 5, 12, 12, 7)
  m0 <- list(const10 = rep(10, 6))
  e <- list(half = rep(0.5, 6))
  mt <- nuisance_metrics(A, Y, fold, e, m0, clip_c = 0.01)
  mse <- mt[mt$learner == "const10" & mt$metric == "mse", ]
  expect_equal(mse$mean, mean(c(5, 4)))
  expect_equal(mse$se, stats::sd(c(5, 4)) / sqrt(2))
  expect_equal(mt[mt$learner == "const10" & mt$metric == "mae", "mean"], 2)
  brier <- mt[mt$learner == "half" & mt$metric == "brier", ]
  expect_equal(brier$mean, 0.25)
  expect_equal(mt[mt$learner == "half" & mt$metric == "logloss", "mean"], log(2))
})

test_that("crossfit_parametric returns out-of-fold predictions for every row", {
  set.seed(2)
  n <- 400
  X <- as.data.frame(matrix(rbinom(n * 4, 1, .5), n, 4)); names(X) <- paste0("x", 1:4)
  A <- rbinom(n, 1, plogis(-.2 + .8 * X$x1))
  Y <- 1 + X$x1 + rnorm(n)
  fold <- rep(1:4, length.out = n)
  r <- crossfit_parametric(X, A, Y, fold)
  expect_length(r$e, n); expect_length(r$m0_main, n); expect_length(r$m0_inter, n)
  expect_true(all(r$e > 0 & r$e < 1))
  expect_equal(nrow(r$n_aliased), 4L)
})

test_that("fold_mse_m0 and paired_mse_summary match hand computation", {
  ## Two folds. Controls: rows 1,2 (fold 1) and 4,5 (fold 2).
  A <- c(0, 0, 1, 0, 0, 1)
  fold <- c(1, 1, 1, 2, 2, 2)
  Y <- c(10, 20, 5, 10, 30, 7)
  tree <- c(11, 19, 0, 12, 28, 0)       # errors 1,1 | 2,2 -> mse 1, 4
  ols  <- c(12, 22, 0, 10, 30, 0)       # errors 2,2 | 0,0 -> mse 4, 0
  const <- rep(15, 6)                   # errors 5,5 | 5,15 -> mse 25, 125
  fm <- fold_mse_m0(A, Y, fold, list(tree_cv = tree, ols_main = ols, constant = const))
  expect_equal(fm$mse[fm$learner == "tree_cv"], c(1, 4))
  expect_equal(fm$mse[fm$learner == "ols_main"], c(4, 0))
  s <- paired_mse_summary(fm)
  expect_equal(s$mse_reference, 2.5)
  expect_equal(s$mse_other, 2)
  expect_equal(s$diff_ref_minus_other, 0.5)
  ## fold diffs (tree - ols) = -3, 4 -> sd = sqrt(24.5), se = sd / sqrt(2)
  expect_equal(s$diff_se_rep1, stats::sd(c(-3, 4)) / sqrt(2))
  expect_equal(s$skill_reference, 1 - 2.5 / 75)
  expect_equal(s$skill_other, 1 - 2 / 75)
})

test_that("held-out OLS predictions do not depend on held-out outcomes", {
  set.seed(5)
  n <- 300
  X <- as.data.frame(matrix(rbinom(n * 3, 1, .5), n, 3)); names(X) <- paste0("x", 1:3)
  A <- rbinom(n, 1, .5)
  Y <- 2 + X$x1 + rnorm(n)
  fold <- rep(1:3, length.out = n)
  r1 <- crossfit_parametric(X, A, Y, fold)
  Y2 <- Y; Y2[fold == 1] <- Y2[fold == 1] + 1e3
  r2 <- crossfit_parametric(X, A, Y2, fold)
  expect_identical(r1$m0_main[fold == 1], r2$m0_main[fold == 1])
  expect_identical(r1$e[fold == 1], r2$e[fold == 1])
})

test_that("att_from_nuisances applies the one-sided clip and matches a hand computation", {
  skip_if_not_installed("doubletree")
  set.seed(6)
  n <- 500
  A <- rbinom(n, 1, 0.4); Y <- rnorm(n, 10 + 2 * A)
  m0 <- rep(10, n)
  ## Propensity above 1 - clip_c is clipped; below clip_c is NOT (one-sided).
  A[2] <- 0L; Y[2] <- 1000     # a control with a huge residual and a tiny e
  e <- rep(0.4, n); e[1] <- 0.9999; e[2] <- 0.0005
  a <- att_from_nuisances(Y, A, e, m0, clip_c = 0.01)
  e_clip <- pmin(0.99, e)
  pi_hat <- mean(A)
  psi0 <- (A - (1 - A) * e_clip / (1 - e_clip)) * (Y - m0) / pi_hat
  expect_equal(a$theta, sum(psi0) / sum(A / pi_hat), tolerance = 1e-8)
  ## A two-sided clip would give a different answer when e[2] is tiny.
  e_two <- pmin(0.99, pmax(0.01, e))
  psi_two <- (A - (1 - A) * e_two / (1 - e_two)) * (Y - m0) / pi_hat
  expect_false(isTRUE(all.equal(sum(psi0), sum(psi_two))))
})

test_that("paired_mse_summary and fold_mse_m0 validate their inputs", {
  A <- c(1, 1, 0, 0); fold <- c(1, 1, 2, 2)
  expect_error(fold_mse_m0(A, c(1, 2, 3, 4), fold, list(a = rep(0, 4))), "no control rows")
  fm <- data.frame(rep = 1, fold = 1:2, learner = "tree_cv", mse = 1:2, n_control = 5)
  expect_error(paired_mse_summary(fm), "both be in")
})
