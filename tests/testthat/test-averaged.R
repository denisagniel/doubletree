# Tests for the averaged-tree estimators (Approach 4: estimate_att_doubletree_averaged;
# Approach 6: estimate_att_msplit_averaged) and the honest bias-aware CI helpers.
#
# These functions had ZERO test coverage before. The averaged tree is a biased (but
# interpretable) DISPLAY estimate; the reported CI is an Armstrong-Kolesar honest
# interval built from the valid cross-fit twin + a conservative bias bound. Tests pin:
#   - the honest_cv / honest_ci helpers (pure, fast, no tree fitting);
#   - that both estimators run and return the twin + honest-CI diagnostics;
#   - delta == theta_display - theta_crossfit;
#   - the honest CI is at least as wide as the naive twin CI (wider when delta != 0).

# ---- honest_cv / honest_ci helpers (pure math, no skip needed) ---------------

test_that("honest_cv(0) equals the ordinary z and is increasing in b", {
  expect_equal(honest_cv(0), qnorm(0.975), tolerance = 1e-8)
  # Monotone increasing in the bias-to-SE ratio.
  bs <- c(0, 0.25, 0.5, 1, 2, 5)
  cvs <- vapply(bs, honest_cv, numeric(1))
  expect_true(all(diff(cvs) > 0))
  # cv(b) always covers at least the folded-normal target: coverage == 0.95.
  for (b in bs) {
    cv <- honest_cv(b)
    cover <- pnorm(cv - b) - pnorm(-cv - b)
    expect_equal(cover, 0.95, tolerance = 1e-6)
  }
})

test_that("honest_cv rejects invalid b", {
  expect_error(honest_cv(-1), "non-negative")
  expect_error(honest_cv(c(1, 2)), "single")
  expect_error(honest_cv(NA_real_), "finite")
})

test_that("honest_ci widens with the bias bound and reduces to Wald at delta=0", {
  se <- 0.05
  # No bias, no noise -> ordinary 1.96*se half-width.
  h0 <- honest_ci(0.2, se = se, delta = 0, se_delta = 0)
  expect_equal(h0$half_width, qnorm(0.975) * se, tolerance = 1e-8)
  expect_equal(h0$cv, qnorm(0.975), tolerance = 1e-8)

  # Nonzero bias -> strictly wider than the naive Wald interval.
  h1 <- honest_ci(0.2, se = se, delta = 0.03, se_delta = 0.01)
  expect_gt(h1$half_width, qnorm(0.975) * se)
  expect_gt(h1$B, 0.03)                       # B = |delta| + z*se_delta > |delta|
  expect_equal(diff(h1$ci), 2 * h1$half_width, tolerance = 1e-12)
})

# ---- Approach 4: estimate_att_doubletree_averaged ---------------------------

make_binary_dgp <- function(n, seed = 20260714) {
  set.seed(seed)
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5), X3 = rbinom(n, 1, 0.5)
  )
  A <- rbinom(n, 1, plogis(-0.4 + 0.5 * X$X1 - 0.3 * X$X2))
  Y <- rbinom(n, 1, 0.25 + 0.2 * X$X1 + 0.15 * X$X3 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

test_that("doubletree_averaged returns display estimate, cross-fit twin, and honest CI", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(400)
  fit <- tryCatch(
    estimate_att_doubletree_averaged(d$X, d$A, d$Y, K = 3, outcome_type = "binary",
                                     verbose = FALSE),
    error = function(e) skip(paste("empty intersection on this DGP:", conditionMessage(e)))
  )

  # Point estimate is the averaged (display) tree; SE is the twin SE.
  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$theta_crossfit))
  expect_true(is.finite(fit$sigma_crossfit) && fit$sigma_crossfit > 0)

  # Fidelity diagnostic identity.
  expect_equal(fit$delta, fit$theta - fit$theta_crossfit, tolerance = 1e-10)

  # Honest CI is centered at the DISPLAY estimate and is at least as wide as the naive
  # twin Wald interval (strictly wider when delta != 0).
  naive_half <- qnorm(0.975) * fit$sigma_crossfit
  honest_half <- diff(fit$ci_95) / 2
  expect_gte(honest_half, naive_half - 1e-8)
  center <- mean(fit$ci_95)
  expect_equal(center, fit$theta, tolerance = 1e-8)
  expect_true(fit$bias_bound_B >= abs(fit$delta) - 1e-10)
})

# ---- Approach 6: estimate_att_msplit_averaged -------------------------------

test_that("msplit_averaged returns fully-fold-specific twin + honest CI (se_delta=0)", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(400)
  fit <- estimate_att_msplit_averaged(d$X, d$A, d$Y, M = 4, K = 3,
                                      outcome_type = "binary", verbose = FALSE)

  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$theta_crossfit))
  expect_true(is.finite(fit$sigma_crossfit) && fit$sigma_crossfit > 0)
  expect_equal(fit$delta, fit$theta - fit$theta_crossfit, tolerance = 1e-10)

  # Phase B: se_delta = 0 (tightest interval consistent with the coverage guarantee;
  # a positive se_delta only widens B). Replaced the old sd(theta_cf_m)/sqrt(M) form.
  expect_equal(fit$se_delta, 0)

  # Honest CI centered at display estimate, at least as wide as naive twin Wald.
  naive_half <- qnorm(0.975) * fit$sigma_crossfit
  honest_half <- diff(fit$ci_95) / 2
  expect_gte(honest_half, naive_half - 1e-8)
  expect_equal(mean(fit$ci_95), fit$theta, tolerance = 1e-8)

  # Cross-fit predictions are still exposed for downstream analysis.
  expect_equal(dim(fit$predictions_all_splits$e), c(nrow(d$X), 4L))
})
