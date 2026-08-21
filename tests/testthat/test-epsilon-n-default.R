# Tests for the theory-aligned Rashomon tolerance default:
# rashomon_bound_multiplier = NULL resolves to optimaltrees::select_epsilon_n(n)
# = log(n)/n (o(n^{-1/2})), and data-adaptive auto_tune_intersecting warns.

make_binary_dgp <- function(n, seed = 42) {
  set.seed(seed)
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
  Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

test_that("NULL rashomon_bound_multiplier resolves to theory epsilon_n = log(n)/n", {
  d <- make_binary_dgp(150)
  fit <- estimate_att_rashomon(d$X, d$A, d$Y, K = 3,
                      rashomon_bound_multiplier = NULL, verbose = FALSE)
  expect_equal(fit$epsilon_n, optimaltrees::select_epsilon_n(150))
  expect_equal(fit$epsilon_n, log(150) / 150)
})

test_that("explicit rashomon_bound_multiplier is honored", {
  d <- make_binary_dgp(150)
  fit <- estimate_att_rashomon(d$X, d$A, d$Y, K = 3,
                      rashomon_bound_multiplier = 0.05, verbose = FALSE)
  expect_equal(fit$epsilon_n, 0.05)
})

test_that("resolved epsilon_n is much smaller than the old 0.05 default", {
  # Sanity: theory value at n=150 is well below the legacy constant.
  expect_lt(optimaltrees::select_epsilon_n(150), 0.05)
})

test_that("auto_tune_intersecting = TRUE warns it voids valid inference", {
  d <- make_binary_dgp(150)
  expect_warning(
    estimate_att_rashomon(d$X, d$A, d$Y, K = 3,
                 auto_tune_intersecting = TRUE, verbose = FALSE),
    "post-selection|valid-inference|o\\(n"
  )
})

test_that("estimate_att_crossfit has no Rashomon tolerance knobs and never warns", {
  d <- make_binary_dgp(120)
  # The post-selection warning is Rashomon-specific: the cross-fit path exposes no
  # auto_tune_intersecting / rashomon_bound_multiplier at all, so there is nothing to
  # warn about. (Pre-split, this test passed auto_tune_intersecting = TRUE alongside
  # use_rashomon = FALSE on the combined estimator, which silently ignored the flag.)
  expect_false("auto_tune_intersecting" %in% names(formals(estimate_att_crossfit)))
  expect_false("rashomon_bound_multiplier" %in% names(formals(estimate_att_crossfit)))
  expect_no_warning(
    estimate_att_crossfit(d$X, d$A, d$Y, K = 3, verbose = FALSE)
  )
})
