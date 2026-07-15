# Phase B (2026-07-15): every shared-Rashomon-structure display estimator reports an
# HONEST bias-aware CI built from the FULLY fold-specific twin (per-fold structure AND
# leaves), not the shared-intersection twin (which shares the display tree's selection
# variance and undercovers -- Phase-A diagnostic). Pins:
#   - estimate_att(use_rashomon=TRUE) now returns a twin + honest ci_95 (was Wald-only);
#   - estimate_att(use_rashomon=FALSE) is unchanged (Wald; it IS the fully-fold-specific
#     estimator, so no honesty correction and no twin);
#   - estimate_att_msplit / estimate_att_single_tree carry the twin + honest CI;
#   - the honest CI is at least as wide as the twin Wald and centered at the display point;
#   - se_delta = 0 everywhere (the chosen conservative-but-tightest bound).

make_binary_dgp <- function(n, seed = 20260715) {
  set.seed(seed)
  X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                  X3 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(-0.5 + 0.3 * X$X1 + 0.3 * X$X2))
  Y <- rbinom(n, 1, 0.2 + 0.15 * X$X1 + 0.15 * X$X3 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

# Helper: honest CI centered at display, at least as wide as the twin Wald.
expect_honest_ci <- function(fit) {
  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$theta_crossfit))          # twin present
  expect_true(is.finite(fit$sigma_crossfit) && fit$sigma_crossfit > 0)
  expect_equal(fit$delta, fit$theta - fit$theta_crossfit, tolerance = 1e-10)
  expect_equal(fit$se_delta, 0)                        # Phase B: tightest bound
  center <- mean(fit$ci_95)
  expect_equal(center, fit$theta, tolerance = 1e-8)    # centered at display estimate
  naive_half  <- qnorm(0.975) * fit$sigma_crossfit
  honest_half <- diff(fit$ci_95) / 2
  expect_gte(honest_half, naive_half - 1e-8)           # honest >= twin Wald
}

test_that("estimate_att(use_rashomon=TRUE) reports a fully-fold-specific twin + honest CI", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(500)
  fit <- estimate_att(d$X, d$A, d$Y, K = 3, use_rashomon = TRUE, verbose = FALSE)
  expect_honest_ci(fit)
  # ci_95_wald is the plain (undercovering) Wald interval, distinct from the honest one
  # whenever delta != 0.
  expect_false(is.null(fit$ci_95_wald))
})

test_that("estimate_att(use_rashomon=FALSE) is the fully-fold-specific estimator (Wald, no twin)", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(500)
  fit <- estimate_att(d$X, d$A, d$Y, K = 3, use_rashomon = FALSE, verbose = FALSE)
  # It IS the twin, so no honesty correction: twin/delta fields are NA and ci_95 == Wald.
  expect_true(is.na(fit$theta_crossfit))
  expect_true(is.na(fit$delta))
  expect_equal(fit$ci_95, fit$ci_95_wald, tolerance = 1e-12)
})

test_that("estimate_att_msplit reports a fully-fold-specific twin + honest CI", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(500)
  fit <- estimate_att_msplit(d$X, d$A, d$Y, M = 4, K = 3, verbose = FALSE)
  expect_honest_ci(fit)
  expect_true(is.finite(fit$sigma_wald) && fit$sigma_wald > 0)  # modal Wald kept for reference
})

test_that("estimate_att_single_tree reports an honest CI from the fully-fold-specific twin", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(500)
  fit <- estimate_att_single_tree(d$X, d$A, d$Y, K = 3, inference = "single", verbose = FALSE)
  # Reported point = single tree; reported CI = honest interval; both twin + honest present.
  expect_equal(fit$theta, fit$theta_single, tolerance = 1e-12)
  expect_equal(fit$ci_95, fit$ci_95_honest, tolerance = 1e-12)
  expect_honest_ci(fit)
})
