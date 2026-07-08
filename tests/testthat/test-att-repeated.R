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
