# Smoke tests for att_repeated (exported repeated-cross-fitting wrapper with
# Chernozhukov et al. (2018) variance combination). Nothing else in the package or
# the arbitration study exercises it, so these guard against silent rot.

test_that("att_repeated with n_splits = 1 delegates to estimate_att", {
  skip_if_not_installed("optimaltrees")

  set.seed(101)
  n <- 150
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.15 * A + 0.1 * X$x1)

  rep1 <- att_repeated(X, A, Y, K = 3, n_splits = 1, seed = 7, verbose = FALSE)
  direct <- estimate_att(X, A, Y, K = 3, seed = 7, verbose = FALSE)

  # n_splits = 1 is a pass-through to estimate_att with the same args/seed.
  expect_equal(rep1$theta, direct$theta)
  expect_equal(rep1$sigma, direct$sigma)
})

test_that("att_repeated with n_splits > 1 runs and returns variance-combination fields", {
  skip_if_not_installed("optimaltrees")

  set.seed(202)
  n <- 200
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A + 0.1 * X$x1)

  res <- att_repeated(X, A, Y, K = 3, n_splits = 3, seed = 11, verbose = FALSE)

  # Core estimate
  expect_type(res$theta, "double")
  expect_true(is.finite(res$theta))
  expect_true(is.finite(res$sigma) && res$sigma >= 0)
  expect_length(res$ci_95, 2)
  expect_true(res$ci_95[1] <= res$theta && res$theta <= res$ci_95[2])

  # Per-split bookkeeping
  expect_length(res$theta_splits, 3)
  expect_length(res$sigma_splits, 3)
  expect_equal(res$n_splits, 3)

  # Chernozhukov variance decomposition components are finite and non-negative.
  # (within_var/between_var are means of the per-split within/between terms; their
  # fractions sum to 1 only under mean aggregation, since the default median
  # aggregation divides by median(total_var) != mean(total_var) -- checked below.)
  expect_true(is.finite(res$within_var) && res$within_var >= 0)
  expect_true(is.finite(res$between_var) && res$between_var >= 0)
})

test_that("att_repeated variance fractions sum to 1 under mean aggregation", {
  skip_if_not_installed("optimaltrees")

  set.seed(404)
  n <- 200
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A + 0.1 * X$x1)

  res <- att_repeated(X, A, Y, K = 3, n_splits = 3, seed = 13,
                      aggregation = "mean", verbose = FALSE)

  # Under mean aggregation sigma^2 = mean(within_s + between_s), so the within and
  # between fractions partition the total variance exactly.
  expect_equal(res$within_var_frac + res$between_var_frac, 1, tolerance = 1e-8)
})

test_that("att_repeated sigma is on the theta-hat scale (not sqrt(n)-inflated)", {
  skip_if_not_installed("optimaltrees")

  # Regression guard for F1: att_se() returns SE on the theta-hat scale
  # (sqrt(mean(psi^2)/n)), so the combined repeated-splitting sigma must ALSO be
  # theta-hat scale, i.e. the same order of magnitude as a single split's SE. The
  # old code multiplied the between-split term by n, making sigma ~sqrt(n) too big.
  set.seed(505)
  n <- 300
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A + 0.1 * X$x1)

  res <- att_repeated(X, A, Y, K = 3, n_splits = 5, seed = 21, aggregation = "mean")

  # The combined SE cannot be smaller than the smallest within-split SE, and must
  # stay the same ORDER as the (mean) within-split SE -- never sqrt(n) larger.
  mean_within_se <- sqrt(mean(res$sigma_splits^2))
  expect_gte(res$sigma, min(res$sigma_splits))
  # Total variance = within + between, so sigma >= mean_within_se and, unless the
  # between-split term dominates pathologically, within a small factor of it.
  expect_gte(res$sigma, mean_within_se - 1e-8)
  expect_lt(res$sigma, 3 * mean_within_se)
  # Hard scale check: sigma is O(n^{-1/2}). The buggy sqrt(n)-inflated value would
  # be ~sqrt(300) ~ 17x larger and blow past this bound.
  expect_lt(res$sigma, 5 * mean_within_se)

  # within_var and between_var are on the same (theta-hat) scale: the between term
  # is a variance of estimates each ~att_se in spread, so it cannot exceed the
  # within term by orders of magnitude for a homogeneous DGP.
  expect_lt(res$between_var, 100 * res$within_var)

  # Total variance decomposition is exact under mean aggregation.
  expect_equal(res$sigma^2, res$within_var + res$between_var, tolerance = 1e-8)
})

test_that("att_repeated respects the aggregation argument", {
  skip_if_not_installed("optimaltrees")

  set.seed(303)
  n <- 160
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, 0.4)
  Y <- rbinom(n, 1, 0.3 + 0.2 * A)

  res_med <- att_repeated(X, A, Y, K = 3, n_splits = 3, seed = 5, aggregation = "median")
  res_mean <- att_repeated(X, A, Y, K = 3, n_splits = 3, seed = 5, aggregation = "mean")

  expect_equal(res_med$aggregation, "median")
  expect_equal(res_mean$aggregation, "mean")
  # Same seed -> same per-split thetas; only the aggregation differs.
  expect_equal(res_med$theta_splits, res_mean$theta_splits)
  expect_equal(res_mean$theta, mean(res_mean$theta_splits))
  expect_equal(res_med$theta, median(res_med$theta_splits))
})
