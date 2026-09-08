# ---------------------------------------------------------------------------
# Fitting + exact evaluation of the remainder terms T_a and T_epsilon.
#
# Everything the analyst can fit is piecewise constant on the 243 grid cells, so
# every population integral below is an EXACT finite sum over cells -- there is
# no Monte-Carlo error in D_w, D_mu, ||a||_mu, or the population projection
# Pi^nu. Only the realised T_a / T_epsilon are sample quantities, as they must be.
#
# Definitions (doubletree/inst/paper/refining-grid-standalone.tex, proof of the
# full-sample CLT, Step 2; identical display in
# _refining-grid-inferential-draft.tex lines 660-731):
#
#   w_0      = e_0 / (1 - e_0),           w_hat = e_hat / (1 - e_hat)
#   Delta_w  = w_hat_{tau_e} - w_0
#   b_i      = Delta_w(Z_i) - mean{ Delta_w(Z_k) : A_k = 0, leaf_mu(k) = leaf_mu(i) }
#   T        = n^{-1} sum_{i:A_i=0} b_i { Y_i - mu_hat_{tau_mu}(Z_i) }
#   eps_i    = Y_i - mu_0(Z_i)
#   a_mu     = Pi^nu_{tau_mu} mu_0 - mu_0,   dnu(z) = (1 - e_0(z)) dP_Z(z)
#   T_eps    = n^{-1} sum_{i:A_i=0} b_i eps_i
#   T_a      = n^{-1} sum_{i:A_i=0} b_i a_mu(Z_i)
#   T        = T_eps - T_a                                      (exact identity)
#
#   ||f||_mu^2 = E[(1 - e_0(Z)) f(Z)^2]
#   D_w = ||w_hat - w_0||_mu ,  D_mu = ||mu_hat - mu_0||_mu
# ---------------------------------------------------------------------------

E_CLIP <- c(0.01, 0.99)      # doubletree/R/score_att.R convention
IDENTITY_TOL <- 1e-8         # T == T_eps - T_a must hold to machine precision


#' Leaf id per row of a binarised design matrix, for a fitted optimaltrees model
leaf_ids <- function(fit, Xbin_df) {
  optimaltrees:::assign_leaf_ids(fit, Xbin_df)
}


#' Fit one leaf-constant nuisance at a given regularization
#'
#' Leaf values are recomputed as within-leaf empirical means of the fitting
#' sample. For log_loss (propensity) and squared_error (outcome) the leaf-optimal
#' constant IS that mean, so this reproduces the fitted tree's own leaf values
#' while guaranteeing the exact refit identity the theory's Step 2 relies on
#' (sum_{i in leaf, A_i = 0} (Y_i - mu_hat_i) = 0).
#'
#' @return list(fit, leaf_train, leaf_value, n_leaves)
fit_leaf_constant <- function(Xbin_df, y, loss_function, lambda, max_depth,
                              time_limit = 300) {
  fit <- suppressWarnings(optimaltrees::fit_tree(
    Xbin_df, y,
    loss_function = loss_function,
    regularization = lambda,
    max_depth = max_depth,
    time_limit = time_limit,
    verbose = FALSE
  ))
  leaf_train <- leaf_ids(fit, Xbin_df)
  leaf_value <- tapply(y, leaf_train, mean)
  list(fit = fit, leaf_train = leaf_train,
       leaf_value = leaf_value, n_leaves = length(leaf_value))
}


#' Search lambda so the realised leaf count lands as close to `target` as
#' possible without exceeding it
#'
#' The package's bisect_lambda_to_budget() only certifies leaves <= budget and
#' stops at its starting lambda when that already holds, which here returns trees
#' far under budget (2 leaves at budget 30). We want the budget to BIND, since
#' the whole question is behaviour at leaf count Lbar.
calibrate_lambda <- function(Xbin_df, y, loss_function, target, max_depth,
                             lambda_hi = 0.05, lambda_lo = 1e-7,
                             n_iter = 24L, time_limit = 300) {
  leaves_at <- function(lam) {
    fit_leaf_constant(Xbin_df, y, loss_function, lam, max_depth, time_limit)$n_leaves
  }
  hi <- lambda_hi                        # large lambda  -> few leaves
  lo <- lambda_lo                        # small lambda  -> many leaves
  n_lo <- leaves_at(lo)
  if (n_lo <= target) {
    # cannot reach the budget even with (almost) no penalty: the data simply do
    # not support Lbar leaves. Report it rather than silently accepting fewer.
    return(list(lambda = lo, n_leaves = n_lo, binding = FALSE))
  }
  best <- list(lambda = hi, n_leaves = leaves_at(hi), binding = TRUE)
  for (i in seq_len(n_iter)) {
    mid <- sqrt(hi * lo)                 # geometric bisection
    n_mid <- leaves_at(mid)
    if (n_mid <= target) {
      hi <- mid
      if (n_mid > best$n_leaves) best <- list(lambda = mid, n_leaves = n_mid,
                                              binding = TRUE)
    } else {
      lo <- mid
    }
    if (best$n_leaves == target) break
  }
  best
}


#' Exact population projection of mu_0 onto a fitted partition, under nu
#'
#' @param cell_leaf integer leaf id of each grid cell under tau_mu
#' @return list(cell_a, a_norm_mu) with cell_a = Pi^nu mu_0 - mu_0 per cell
project_mu0 <- function(cell_leaf, truth) {
  wt <- truth$cell_mass * (1 - truth$cell_e0)          # dnu
  num <- tapply(wt * truth$cell_mu0, cell_leaf, sum)
  den <- tapply(wt, cell_leaf, sum)
  stopifnot("project_mu0: leaf with zero control mass" = all(den > 0))
  proj <- (num / den)[as.character(cell_leaf)]
  cell_a <- as.numeric(proj) - truth$cell_mu0
  list(cell_a = cell_a, a_norm_mu = sqrt(sum(wt * cell_a^2)))
}


#' Exact ||f_hat - f_0||_mu for a cell-constant fitted function
norm_mu <- function(cell_fit, cell_true, truth) {
  wt <- truth$cell_mass * (1 - truth$cell_e0)
  sqrt(sum(wt * (cell_fit - cell_true)^2))
}


#' Fit at the cell's calibrated lambda, enforcing the leaf budget on THIS draw
#'
#' lambda is calibrated once per (n, Lbar) cell, so an individual draw can
#' overshoot the budget. The theory reads the estimator inside the class
#' T_{Lbar,n} of trees with at most Lbar leaves, so an overshoot would be outside
#' the object being studied. Raise lambda geometrically until the budget holds.
fit_within_budget <- function(Xbin_df, y, loss_function, lambda, budget,
                              max_depth, time_limit = 300, max_tries = 12L) {
  lam <- lambda
  for (k in seq_len(max_tries)) {
    out <- fit_leaf_constant(Xbin_df, y, loss_function, lam, max_depth, time_limit)
    if (out$n_leaves <= budget) {
      out$lambda_used <- lam
      out$n_refits <- k - 1L
      return(out)
    }
    lam <- lam * 1.6
  }
  stop("fit_within_budget: leaf budget ", budget, " not met after ", max_tries,
       " lambda increases (last lambda ", signif(lam, 3), ").", call. = FALSE)
}


#' One replication: fit both nuisances, compute every reported quantity
#'
#' @param lambda_e,lambda_mu regularizations calibrated for this (n, L) cell
#' @return one-row data.frame
run_replication <- function(n, truth, lambda_e, lambda_mu, time_limit = 300) {
  dat <- simulate_data(n, truth)
  X_all <- as.data.frame(dat$Xbin)
  ctrl <- dat$A == 0L
  X_ctrl <- as.data.frame(dat$Xbin[ctrl, , drop = FALSE])
  Xcells <- as.data.frame(binarise(truth$cell_centres))

  # ---- propensity tree (all units, log_loss) ------------------------------
  e_fit <- fit_within_budget(X_all, dat$A, "log_loss", lambda_e, truth$L,
                             truth$depth_e, time_limit)
  cell_leaf_e <- leaf_ids(e_fit$fit, Xcells)
  cell_e_hat <- pmin(E_CLIP[2], pmax(E_CLIP[1],
                                     as.numeric(e_fit$leaf_value[as.character(cell_leaf_e)])))
  cell_w_hat <- cell_e_hat / (1 - cell_e_hat)

  # ---- control-outcome tree (controls only, squared_error) ----------------
  m_fit <- fit_within_budget(X_ctrl, dat$Y[ctrl], "squared_error", lambda_mu,
                             truth$L, truth$depth_mu, time_limit)
  cell_leaf_mu <- leaf_ids(m_fit$fit, Xcells)
  cell_mu_hat <- as.numeric(m_fit$leaf_value[as.character(cell_leaf_mu)])

  # ---- exact population discrepancies ------------------------------------
  D_w <- norm_mu(cell_w_hat, truth$cell_w0, truth)
  D_mu <- norm_mu(cell_mu_hat, truth$cell_mu0, truth)
  proj <- project_mu0(cell_leaf_mu, truth)

  # ---- per-unit fitted values, via each unit's grid cell ------------------
  w_hat_i <- cell_w_hat[dat$cell]
  mu_hat_i <- cell_mu_hat[dat$cell]
  a_i <- proj$cell_a[dat$cell]
  leaf_mu_i <- cell_leaf_mu[dat$cell]

  # ---- realised remainder terms ------------------------------------------
  delta_w <- w_hat_i - dat$w0
  dw_c <- delta_w[ctrl]
  leaf_c <- leaf_mu_i[ctrl]
  leaf_mean_dw <- tapply(dw_c, leaf_c, mean)
  b_c <- dw_c - as.numeric(leaf_mean_dw[as.character(leaf_c)])

  T_eps <- sum(b_c * dat$eps[ctrl]) / n
  T_a <- sum(b_c * a_i[ctrl]) / n
  T_direct <- sum(delta_w[ctrl] * (dat$Y[ctrl] - mu_hat_i[ctrl])) / n
  identity_gap <- abs(T_direct - (T_eps - T_a))

  # ---- ATT point estimate, plug-in variance, Wald interval ---------------
  resid <- dat$Y - mu_hat_i
  pi_hat <- mean(dat$A)
  theta_hat <- (mean(dat$A * resid) - mean((1 - dat$A) * w_hat_i * resid)) / pi_hat
  g <- dat$A * (resid - theta_hat) - (1 - dat$A) * w_hat_i * resid
  V_hat <- mean(g^2) / pi_hat^2
  se <- sqrt(V_hat / n)
  covered <- abs(theta_hat - truth$theta0) <= stats::qnorm(0.975) * se

  data.frame(
    n = n, L = truth$L,
    leaves_e = e_fit$n_leaves, leaves_mu = m_fit$n_leaves,
    refits_e = e_fit$n_refits, refits_mu = m_fit$n_refits,
    D_w = D_w, D_mu = D_mu, a_norm_mu = proj$a_norm_mu,
    T_eps = T_eps, T_a = T_a, T_total = T_eps - T_a,
    identity_gap = identity_gap,
    theta_hat = theta_hat, se = se, covered = covered,
    pi_hat = pi_hat, n_control = sum(ctrl)
  )
}
