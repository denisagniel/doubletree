# =============================================================================
# dgp.R -- data-generating processes for the "six-approach-arbitration" study
# =============================================================================
# Ported from simulations/six_approach_comparison/code/dgps.R (2026-05-01).
# Four DGPs of increasing complexity; all have true ATT = 0.15, binary outcome.
#   simple     : linear propensity + outcome, binary X       (~2-3 splits)
#   moderate   : 2-way interactions, binary X                (~4-5 splits)
#   complex    : multiple/3-way interactions, binary X       (~6-8 splits)
#   continuous : mixed binary + continuous X, x4^2 nonlinearity  (STRESS regime:
#                the quadratic requires many thresholds; tree discretization
#                struggles -- Constitution Section 9 stress regime)
#
# The RNG seed is set by run_one() BEFORE generate_data() is called, so do NOT
# call set.seed() here.
# =============================================================================

expit <- function(x) 1 / (1 + exp(-x))

# Dispatch on config$dgp. One branch per regime in GRID$dgp.
generate_data <- function(config) {
  switch(
    config$dgp,
    simple     = .dgp_simple(config$n),
    moderate   = .dgp_moderate(config$n),
    complex    = .dgp_complex(config$n),
    continuous = .dgp_continuous(config$n),
    stop(sprintf("Unknown dgp regime: '%s'", config$dgp))  # no silent fallback
  )
}

# --- simple: linear propensity (x1,x2) + linear outcome (x1,x3), binary X ----
.dgp_simple <- function(n) {
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                  x3 = rbinom(n, 1, 0.5))
  e_true <- expit(-0.5 + 0.3 * X$x1 + 0.3 * X$x2)
  A <- rbinom(n, 1, e_true)
  mu0 <- 0.2 + 0.15 * X$x1 + 0.15 * X$x3
  Y0 <- rbinom(n, 1, mu0)
  Y1 <- rbinom(n, 1, pmin(mu0 + 0.15, 1))          # ATT = 0.15
  list(X = X, A = A, Y = A * Y1 + (1 - A) * Y0)
}

# --- moderate: 2-way interactions, binary X ----------------------------------
.dgp_moderate <- function(n) {
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                  x3 = rbinom(n, 1, 0.5), x4 = rbinom(n, 1, 0.5))
  e_true <- expit(-0.5 + 0.3 * X$x1 + 0.2 * X$x2 + 0.3 * X$x1 * X$x2)
  A <- rbinom(n, 1, e_true)
  mu0 <- 0.2 + 0.2 * X$x3 + 0.15 * X$x4 + 0.2 * X$x3 * X$x4
  Y0 <- rbinom(n, 1, mu0)
  Y1 <- rbinom(n, 1, pmin(mu0 + 0.15, 1))
  list(X = X, A = A, Y = A * Y1 + (1 - A) * Y0)
}

# --- complex: multiple interactions, binary X --------------------------------
.dgp_complex <- function(n) {
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                  x3 = rbinom(n, 1, 0.5), x4 = rbinom(n, 1, 0.5),
                  x5 = rbinom(n, 1, 0.5))
  e_true <- expit(-0.5 + 0.2 * (X$x1 + X$x2 + X$x3) +
                    0.3 * X$x1 * X$x2 + 0.2 * X$x2 * X$x3)
  A <- rbinom(n, 1, e_true)
  # Intercept 0.05 keeps max(mu0) = 0.85 so pmin(mu0 + 0.15, 1) never clips.
  mu0 <- 0.05 + 0.15 * (X$x3 + X$x4 + X$x5) +
    0.2 * X$x3 * X$x4 + 0.15 * X$x4 * X$x5
  Y0 <- rbinom(n, 1, mu0)
  Y1 <- rbinom(n, 1, pmin(mu0 + 0.15, 1))
  list(X = X, A = A, Y = A * Y1 + (1 - A) * Y0)
}

# --- continuous (STRESS): mixed binary + continuous, x4^2 nonlinearity -------
.dgp_continuous <- function(n) {
  X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                  x3 = runif(n, -1, 1), x4 = rnorm(n, 0, 1))
  e_true <- expit(-0.5 + 0.3 * X$x1 + 0.4 * X$x3 + 0.2 * X$x4 +
                    0.2 * X$x1 * X$x3)
  A <- rbinom(n, 1, e_true)
  mu0 <- 0.2 + 0.15 * X$x2 + 0.2 * X$x3 + 0.15 * (X$x4^2 / 2) +
    0.1 * X$x2 * X$x3
  mu0 <- pmax(0.01, pmin(0.99, mu0))
  Y0 <- rbinom(n, 1, mu0)
  Y1 <- rbinom(n, 1, pmin(mu0 + 0.15, 1))
  list(X = X, A = A, Y = A * Y1 + (1 - A) * Y0)
}

# The estimand: true ATT is 0.15 in every DGP (by construction above).
true_value <- function(config) {
  0.15
}
