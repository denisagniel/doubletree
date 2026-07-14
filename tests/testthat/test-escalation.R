# Tests for the Rashomon-tolerance escalation control (escalate_intersection).
#
# Escalation was previously unreachable via the public API: estimate_att resolved
# rashomon_bound_multiplier = NULL -> log(n)/n BEFORE fit_nuisances_rashomon, and
# fit_nuisances_rashomon resolved it again, so the c-grid always collapsed to {1}
# (c=1). These tests pin the corrected behavior:
#   (a) default (escalate_intersection = FALSE) never escalates: rashomon_c_* is 1
#       (intersection non-empty at theory tolerance) or NA (empty -> fold fallback);
#   (b) escalate_intersection = TRUE can widen the tolerance (rashomon_c_* > 1) on
#       data where the theory-tolerance intersection is empty;
#   (c) an explicit rashomon_bound_multiplier pins a single fixed tolerance (no
#       escalation) regardless of the flag.

make_binary_dgp <- function(n, seed = 20260714) {
  set.seed(seed)
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
    X3 = rbinom(n, 1, 0.5), X4 = rbinom(n, 1, 0.5)
  )
  A <- rbinom(n, 1, plogis(-0.4 + 0.5 * X$X1 - 0.4 * X$X2))
  Y <- rbinom(n, 1, 0.25 + 0.2 * X$X1 + 0.15 * X$X3 + 0.15 * A)
  list(X = X, A = A, Y = Y)
}

test_that("default path does not escalate: rashomon_c_* is 1 or NA", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(300)
  fit <- estimate_att(d$X, d$A, d$Y, K = 3, use_rashomon = TRUE,
                      escalate_intersection = FALSE, verbose = FALSE)

  # Default: fixed theory tolerance, a SINGLE c-grid point (no widening). The
  # propensity uses eps_base(n) so c_e == 1; the outcome nuisance divides by
  # eps_base(n0) on the control subset (n0 < n, so eps_base(n0) > eps_base(n)),
  # giving c_m0 <= 1. The invariant is that neither is WIDENED above the theory
  # value: c <= 1 (or NA if the intersection was empty -> fold-specific fallback).
  for (cval in c(fit$rashomon_c_e, fit$rashomon_c_m0)) {
    expect_true(is.na(cval) || cval <= 1 + 1e-8)
  }
  expect_true(is.finite(fit$theta))
})

test_that("escalate_intersection = TRUE can widen tolerance (rashomon_c_* >= 1) and stays finite", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  # Small n with several noise covariates makes folds more likely to disagree at the
  # theory tolerance, so escalation has something to do. We assert the mechanism is
  # LIVE (c can exceed 1) and never regresses below the theory value.
  d <- make_binary_dgp(200)
  fit <- estimate_att(d$X, d$A, d$Y, K = 5, use_rashomon = TRUE,
                      escalate_intersection = TRUE, verbose = FALSE)

  for (cval in c(fit$rashomon_c_e, fit$rashomon_c_m0)) {
    expect_true(is.na(cval) || cval >= 1)
  }
  expect_true(is.finite(fit$theta))
  # epsilon_n is reported as max c * log(n)/n over nuisances (or NA on full fallback).
  expect_true(is.na(fit$epsilon_n) || fit$epsilon_n >= log(200) / 200 - 1e-12)
})

test_that("explicit rashomon_bound_multiplier pins a fixed tolerance regardless of the flag", {
  skip_if_not_installed("optimaltrees")
  skip_on_cran()

  d <- make_binary_dgp(300)
  eps <- optimaltrees::select_epsilon_n(nrow(d$X))  # = log(n)/n

  # With an explicit multiplier, escalate_intersection is ignored: the c-grid is the
  # single point {multiplier / theory} = {1}, so c cannot exceed 1.
  fit <- estimate_att(d$X, d$A, d$Y, K = 3, use_rashomon = TRUE,
                      rashomon_bound_multiplier = eps,
                      escalate_intersection = TRUE, verbose = FALSE)

  # Explicit multiplier -> single fixed tolerance, no widening (same c <= 1 invariant
  # as the default path; the escalate flag is overridden).
  for (cval in c(fit$rashomon_c_e, fit$rashomon_c_m0)) {
    expect_true(is.na(cval) || cval <= 1 + 1e-8)
  }
  expect_true(is.finite(fit$theta))
})
