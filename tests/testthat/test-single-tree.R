# Tests for estimate_att_single_tree (Alternative A): one honest tree per
# nuisance (goal i) + cross-fit twin and delta fidelity diagnostic (goal ii).

make_binary_dgp <- function(n, seed = 20260707) {
  set.seed(seed)
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                  X3 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(-0.5 + 0.3 * X$X1 + 0.3 * X$X2))
  Y <- rbinom(n, 1, 0.2 + 0.15 * X$X1 + 0.15 * X$X3 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

test_that("returns one tree per nuisance plus both estimators and the diagnostic", {
  d <- make_binary_dgp(400)
  fit <- estimate_att_single_tree(d$X, d$A, d$Y, K = 3, verbose = FALSE)

  # Single tree per nuisance (goal i)
  expect_true(is.list(fit$tree_e))
  expect_true(is.list(fit$tree_m0))
  expect_true(fit$converged)

  # Both estimators present (goal ii)
  expect_true(is.numeric(fit$theta_single))
  expect_true(is.numeric(fit$theta_crossfit))
  expect_length(fit$ci_95_single, 2)
  expect_length(fit$ci_95_crossfit, 2)

  # Diagnostic: delta = single - crossfit
  expect_equal(fit$delta, fit$theta_single - fit$theta_crossfit)
  expect_equal(fit$delta_over_se, fit$delta / fit$sigma_crossfit)

  # theory epsilon_n
  expect_equal(fit$epsilon_n, optimaltrees::select_epsilon_n(400))
})

test_that("inference = 'single' vs 'crossfit' selects the reported target", {
  d <- make_binary_dgp(400)
  # Fix seed so both runs use identical folds -> identical underlying estimates.
  fit_s <- estimate_att_single_tree(d$X, d$A, d$Y, K = 3, inference = "single", seed = 1)
  fit_c <- estimate_att_single_tree(d$X, d$A, d$Y, K = 3, inference = "crossfit", seed = 1)

  expect_equal(fit_s$theta, fit_s$theta_single)
  expect_equal(fit_c$theta, fit_c$theta_crossfit)
  # Same seed -> same folds -> same underlying estimates regardless of report choice
  expect_equal(fit_s$theta_single, fit_c$theta_single)
  expect_equal(fit_s$theta_crossfit, fit_c$theta_crossfit)
})

test_that("single tree faithfully tracks the cross-fit twin on a clean DGP", {
  # On a well-separated binary DGP the margin holds, so the single tree should
  # closely match its cross-fit twin: |delta| well within one SE.
  d <- make_binary_dgp(500)
  fit <- estimate_att_single_tree(d$X, d$A, d$Y, K = 3)
  expect_lt(abs(fit$delta_over_se), 0.5)
})
