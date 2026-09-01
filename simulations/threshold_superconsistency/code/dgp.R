# ============================================================
# dgp.R
# Study: threshold_superconsistency  (doubletree)
#
# One DGP: a regression function that is EXACTLY a depth-2 tree with a genuine
# JUMP of size kappa at every split boundary. This is the object the paper's
# split-location theory is about -- superconsistency (|t_hat - t*| = O_p(1/n))
# is claimed to come from the jump, not from grid refinement, so the DGP must
# have a real discontinuity and not a smooth ramp.
#
# Why a single regression tree and not the full ATT machinery: the claim under
# test is purely about recovering ONE tree's split locations. Simulating a
# regression tree directly isolates that and removes every confound coming from
# the propensity/outcome/EIF plumbing.
#
# Sourced by run_cell.R and analyze.R. Not run directly.
# ============================================================

# True tree (depth 2, 3 leaves). X3 is a pure noise coordinate.
#
#                     X1 <= T1_STAR ?
#                   /                 \
#                yes                   no
#                 |                     |
#          X2 <= T2_STAR ?           mu = 0.6   (leaf C)
#          /            \
#        yes             no
#         |               |
#    mu = 0.3        mu = 0.0
#    (leaf A)        (leaf B)
#
# Leaf values are spaced exactly KAPPA apart, so the MINIMUM jump across each
# split boundary is exactly KAPPA:
#   across split 1 (X1 = T1_STAR): |0.3 - 0.6| = KAPPA on the {X2 <= T2_STAR}
#     face, |0.0 - 0.6| = 2*KAPPA on the other face  -> min = KAPPA
#   across split 2 (X2 = T2_STAR): |0.3 - 0.0| = KAPPA
KAPPA   <- 0.3    # jump size in the regression function at each split boundary
T1_STAR <- 0.45   # true root threshold on X1
T2_STAR <- 0.62   # true second-split threshold on X2 (inside {X1 <= T1_STAR})
SIGMA   <- 0.25   # residual sd
P_COORD <- 3L     # number of covariates (X3 is noise)

# T1_STAR / T2_STAR are deliberately NOT round quantiles of Uniform(0,1).
# If a true threshold coincided with a grid point of the quantile
# discretization, Stage A alone would already localise it at the n^{-1/2} rate
# of a sample quantile and the diagnostic could not distinguish "Stage B beat
# the grid" from "the grid happened to be right". Both values sit strictly
# between grid points for every bin count used here.

#' True regression function of the jump-tree DGP
#'
#' @param X Data frame with numeric columns X1, X2, X3.
#' @return Numeric vector of conditional means, taking values in {0, KAPPA, 2*KAPPA}.
mu_true <- function(X) {
  ifelse(
    X$X1 > T1_STAR,
    2 * KAPPA,                                  # leaf C
    ifelse(X$X2 <= T2_STAR, KAPPA, 0)           # leaf A / leaf B
  )
}

#' Generate one replicate from the jump-tree DGP
#'
#' Covariates are iid Uniform(0,1) so the density is bounded away from zero at
#' both split boundaries -- the design condition the O_p(1/n) split-location
#' rate needs.
#'
#' @param n Sample size.
#' @return List with X (data frame, n x 3), Y (numeric), mu (numeric, the truth),
#'   and the true thresholds.
generate_jump_tree <- function(n) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 100)
  X <- data.frame(
    X1 = stats::runif(n),
    X2 = stats::runif(n),
    X3 = stats::runif(n)   # noise coordinate: never enters mu
  )
  mu <- mu_true(X)
  list(
    X = X,
    Y = mu + stats::rnorm(n, sd = SIGMA),
    mu = mu,
    t1_star = T1_STAR,
    t2_star = T2_STAR,
    kappa = KAPPA
  )
}

#' Population risk reduction from each of the two true splits
#'
#' Used to sanity-check that the regularization used at each n is small enough
#' for the OPTIMAL tree to actually contain both splits. The second split is the
#' binding one: its population gain is constant in n while lambda = log(n)/n
#' shrinks, so the true topology becomes easier to recover as n grows, and the
#' check below tells us the smallest n at which it is recoverable at all.
#'
#' @return List with the two population gains (on the mean-squared-error scale)
#'   and the total variance of mu.
population_split_gains <- function() {
  p_left  <- T1_STAR                 # P(X1 <= t1*)
  p_a     <- T2_STAR                 # P(X2 <= t2* | X1 <= t1*)
  # var of mu overall
  probs <- c(1 - p_left, p_left * p_a, p_left * (1 - p_a))
  vals  <- c(2 * KAPPA, KAPPA, 0)
  m1    <- sum(probs * vals)
  var_mu <- sum(probs * vals^2) - m1^2
  # after the root split, residual within-node variance (right node is constant)
  var_within_left <- p_a * (1 - p_a) * KAPPA^2
  resid_after_root <- p_left * var_within_left
  list(
    var_mu = var_mu,
    gain_split1 = var_mu - resid_after_root,
    gain_split2 = resid_after_root,   # the second split removes all of it
    resid_after_root = resid_after_root
  )
}
