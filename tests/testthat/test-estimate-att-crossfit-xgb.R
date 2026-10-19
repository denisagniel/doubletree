# estimate_att_crossfit_xgb(): the fast xgboost anchor.

make_xgb_dgp <- function(n = 2000, seed = 1) {
  set.seed(seed)
  X <- as.data.frame(matrix(rbinom(n * 5, 1, 0.5), n, 5))
  names(X) <- paste0("x", 1:5)
  e <- plogis(-0.4 + 0.9 * X$x1 - 0.6 * X$x2)
  A <- rbinom(n, 1, e)
  Y0 <- 5 + 2 * X$x1 + 1.5 * X$x2 * X$x3 + rnorm(n)
  # Constant effect of 2 on everyone, so ATT = 2.
  list(X = X, A = A, Y = Y0 + 2 * A)
}

test_that("estimate_att_crossfit_xgb recovers a known ATT and returns the crossfit shape", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp()
  fit <- estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 5, seed = 20260919)
  expect_true(all(c("theta", "sigma", "ci_95", "ci_95_wald", "score_values",
                    "nuisance_fits", "fold_indices", "n", "K") %in% names(fit)))
  expect_equal(fit$learner, "xgboost")
  expect_length(fit$score_values, nrow(d$X))
  expect_lt(abs(fit$theta - 2), 4 * fit$sigma)
  expect_true(fit$ci_95[1] < fit$theta && fit$theta < fit$ci_95[2])
  expect_true(all(fit$nuisance_fits$propensity <= 1 - fit$clip_c))
})

test_that("estimate_att_crossfit_xgb is reproducible for a fixed seed", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 800)
  f1 <- estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 3, seed = 7)
  f2 <- estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 3, seed = 7)
  expect_identical(f1$theta, f2$theta)
  expect_identical(f1$fold_indices, f2$fold_indices)
})

test_that("estimate_att_crossfit_xgb validates its inputs", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 400)
  expect_error(estimate_att_crossfit_xgb(d$X, d$A, d$Y), "seed")
  expect_error(estimate_att_crossfit_xgb(d$X, d$A, d$Y, seed = 1, K = 1), "K must")
  expect_error(estimate_att_crossfit_xgb(d$X, d$A, d$Y, seed = 1, clip_c = 0.7), "clip_c")
  expect_error(estimate_att_crossfit_xgb(d$X, d$A, d$Y, seed = 1, outcome_type = "binary"),
               "binary")
})

test_that("estimate_att_crossfit_xgb handles a binary outcome", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 1500, seed = 3)
  Yb <- rbinom(nrow(d$X), 1, plogis(-0.5 + d$X$x1 + 0.4 * d$A))
  fit <- estimate_att_crossfit_xgb(d$X, d$A, Yb, K = 4, outcome_type = "binary", seed = 11)
  expect_true(is.finite(fit$theta) && is.finite(fit$sigma))
})

test_that("held-out predictions do not depend on held-out outcomes (no leakage)", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 800)
  f1 <- estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 3, seed = 7)
  te <- f1$fold_indices == 1
  Y2 <- d$Y
  Y2[te] <- Y2[te] + 1e4
  f2 <- estimate_att_crossfit_xgb(d$X, d$A, Y2, K = 3, seed = 7)
  expect_identical(f1$nuisance_fits$outcome_control[te], f2$nuisance_fits$outcome_control[te])
  expect_identical(f1$nuisance_fits$propensity[te], f2$nuisance_fits$propensity[te])
})

test_that("early stopping fires and the propensity tracks the truth", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 2000)
  f <- estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 3, seed = 7, nrounds_max = 2000L)
  expect_false(f$hit_nrounds_max)
  expect_true(all(unlist(f$nuisance_fits$best_iterations[, c("e", "m0")]) < 2000L))
  e_true <- plogis(-0.4 + 0.9 * d$X$x1 - 0.6 * d$X$x2)
  expect_gt(stats::cor(f$nuisance_fits$propensity, e_true), 0.9)
})

test_that("the caller's RNG stream is preserved", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 600)
  set.seed(99)
  before <- .Random.seed
  estimate_att_crossfit_xgb(d$X, d$A, d$Y, K = 3, seed = 7)
  expect_identical(.Random.seed, before)
})

test_that("reserved xgboost params are rejected", {
  skip_if_not_installed("xgboost")
  d <- make_xgb_dgp(n = 400)
  expect_error(estimate_att_crossfit_xgb(d$X, d$A, d$Y, seed = 1, params = list(nthread = 8)),
               "may not set")
})
