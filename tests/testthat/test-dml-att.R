# DML ATT and score tests
# Requires treefarmr to be installed (dmltree Imports treefarmr).

test_that("psi_att returns vector of length n and is linear in theta", {
  n <- 5
  Y <- c(1, 0, 1, 0, 1)
  A <- c(1, 0, 1, 0, 1)
  eta <- list(e = rep(0.5, n), m0 = rep(0.4, n), m1 = rep(0.6, n))
  pi_hat <- 0.6
  s0 <- psi_att(Y, A, theta = 0, eta, pi_hat)
  s1 <- psi_att(Y, A, theta = 1, eta, pi_hat)
  expect_length(s0, n)
  expect_length(s1, n)
  # psi(theta) = psi(0) - theta * (A/pi)
  expect_equal(s1, s0 - 1 * (A / pi_hat), tolerance = 1e-10)
})

test_that("dml_att returns list with theta, sigma, ci_95 and runs with binary data", {
  skip_if_not_installed("treefarmr")
  set.seed(42)
  n <- 120
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
  Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
  fit <- dml_att(X, A, Y, K = 3)
  expect_type(fit$theta, "double")
  expect_length(fit$theta, 1)
  expect_type(fit$sigma, "double")
  expect_length(fit$ci_95, 2)
  expect_true(fit$ci_95[1] < fit$theta)
  expect_true(fit$ci_95[2] > fit$theta)
  expect_equal(fit$n, n)
  expect_equal(fit$K, 3)
})

test_that("create_folds returns integer vector 1..K", {
  f <- create_folds(100, K = 5)
  expect_length(f, 100)
  expect_true(all(f >= 1 & f <= 5))
  expect_type(f, "integer")
})

test_that("dml_att_variance and dml_att_ci work", {
  scores <- rnorm(100, 0, 1)
  v <- dml_att_variance(scores)
  expect_true(v > 0)
  ci <- dml_att_ci(0.5, sqrt(v), 100, level = 0.95)
  expect_length(ci, 2)
  expect_true(ci[1] < 0.5)
  expect_true(ci[2] > 0.5)
})
