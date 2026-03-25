#' Stress-Test DGPs for DML-ATT
#'
#' These DGPs deliberately violate or stretch key assumptions to test
#' method robustness. Research constitution §9 requires "regimes where
#' method should struggle" to avoid quiet favoritism.
#'
#' Each DGP documents:
#' - What assumption is violated/stretched
#' - Expected failure mode
#' - When method recovers (if ever)
#'
#' Paper reference: Manuscript Section 4, Table 2 (stress-test results)

#' DGP 4: Weak Overlap
#'
#' **Assumption stretched:** Positivity (propensity scores near boundaries)
#'
#' **Design:**
#' - Propensity scores pushed toward 0.05 and 0.95 (borderline positivity)
#' - Some covariate patterns have very few treated or control units
#' - Treatment effect remains well-defined (τ = 0.10)
#'
#' **Expected failure mode:**
#' - Valid estimates but LARGE variance (CI width 2-3× larger than DGP 1)
#' - Coverage maintains 95% but efficiency loss is severe
#' - Extreme weights in DML score (some obs have weight ≈ 1/0.05 = 20)
#'
#' **Recovery:** Does not "recover" with larger n (weak overlap persists)
#' but CIs narrow at √n rate, just from a higher baseline.
#'
#' @examples
#' d <- generate_dgp_weak_overlap(800, tau = 0.10, seed = 123)
#' mean(d$true_e < 0.1 | d$true_e > 0.9)  # % with extreme propensity

generate_dgp_weak_overlap <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 4 binary features (same as DGP 1 but more extreme propensity)
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5)
  )

  # Propensity: MORE EXTREME coefficients push toward boundaries
  # Range: approximately [0.08, 0.92] (borderline positivity)
  linear_pred <- -2.5 + 4.0 * X$X1 + 2.5 * X$X2
  e <- plogis(linear_pred)
  A <- as.integer(runif(n) < e)

  # Outcome: same additive structure as DGP 1
  p0 <- plogis(-0.3 + 0.5 * X$X1 + 0.4 * X$X2)
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT (among treated)
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  # Diagnostics
  prop_extreme <- mean(e < 0.1 | e > 0.9)

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "weak_overlap",
    diagnostics = list(
      propensity_range = range(e),
      prop_extreme_propensity = prop_extreme
    )
  )
}


#' DGP 5: Non-Smooth (Piecewise Constant)
#'
#' **Assumption stretched:** Smoothness (nuisances are step functions)
#'
#' **Design:**
#' - Propensity and outcome are piecewise constant (4 regions)
#' - Sharp discontinuities at region boundaries
#' - Trees naturally handle this (splits align with discontinuities)
#' - Forests okay (ensembles approximate steps)
#' - Linear models MISSPECIFIED (assumes smooth functions)
#'
#' **Expected failure mode:**
#' - Tree-DML: Good performance (natural fit)
#' - Forest-DML: Good performance (can approximate steps)
#' - Linear-DML: LARGE BIAS (2-3× worse than trees)
#'
#' **Recovery:** Linear never recovers (fundamentally misspecified).
#' Trees and forests perform well at all sample sizes.
#'
#' @examples
#' d <- generate_dgp_piecewise(800, tau = 0.10, seed = 123)
#' # Visualize regions (2D binary covariates)
#' table(X1 = d$X$X1, X2 = d$X$X2)

generate_dgp_piecewise <- function(n, tau = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 2 continuous features (easier to visualize regions)
  # Extended with 2 noise features for consistency
  X <- data.frame(
    X1 = runif(n, 0, 1),
    X2 = runif(n, 0, 1),
    X3 = as.integer(runif(n) < 0.5),  # Noise
    X4 = as.integer(runif(n) < 0.5)   # Noise
  )

  # Propensity: PIECEWISE CONSTANT (4 regions)
  # Region 1: X1 < 0.5, X2 < 0.5 → e = 0.25
  # Region 2: X1 >= 0.5, X2 < 0.5 → e = 0.45
  # Region 3: X1 < 0.5, X2 >= 0.5 → e = 0.60
  # Region 4: X1 >= 0.5, X2 >= 0.5 → e = 0.75
  e <- ifelse(X$X1 < 0.5 & X$X2 < 0.5, 0.25,
       ifelse(X$X1 >= 0.5 & X$X2 < 0.5, 0.45,
       ifelse(X$X1 < 0.5 & X$X2 >= 0.5, 0.60, 0.75)))
  A <- as.integer(runif(n) < e)

  # Outcome: PIECEWISE CONSTANT (4 regions, different boundaries)
  # Region 1: X1 < 0.3 → p0 = 0.20
  # Region 2: 0.3 <= X1 < 0.7 → p0 = 0.40
  # Region 3: X1 >= 0.7 → p0 = 0.60
  p0 <- ifelse(X$X1 < 0.3, 0.20,
        ifelse(X$X1 < 0.7, 0.40, 0.60))
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "piecewise_nonsmooth",
    diagnostics = list(
      n_propensity_regions = length(unique(e)),
      n_outcome_regions = length(unique(p0))
    )
  )
}


#' DGP 6: High-Dimensional Sparse Signal
#'
#' **Assumption stretched:** Dimensionality (p = 8 features, n/p = 100)
#'
#' **Design:**
#' - 8 binary features → 256 possible covariate patterns
#' - ONLY 3 features matter (X1, X2, X3) → sparse signal
#' - Other 5 features are pure noise
#' - At n = 800: only 800/256 ≈ 3 observations per pattern (on average)
#' - Pattern-based methods (trees) may struggle with data sparsity
#'
#' **Expected failure mode:**
#' - At n = 400: Coverage << 95% (underpowered, high variance)
#' - At n = 800: Coverage ≈ 90-92% (borderline)
#' - At n = 1600: Coverage recovers to 95% (sufficient data per pattern)
#' - RMSE higher than DGPs 1-3 (slower convergence rate)
#'
#' **Recovery:** Recovers at n ≥ 1600. Demonstrates curse of dimensionality.
#'
#' **Constitution note:** This stress test is NOT about breaking the method,
#' but showing where sample size requirements increase. Method is correctly
#' specified; just needs more data.
#'
#' @examples
#' d <- generate_dgp_high_dim(1600, tau = 0.10, seed = 123)
#' ncol(d$X)  # 8 features
#' # Only X1, X2, X3 used in d$true_e and d$true_m0

generate_dgp_high_dim <- function(n, tau = 0.10, p = 8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # p binary features (default 8)
  X <- as.data.frame(matrix(
    as.integer(runif(n * p) < 0.5),
    nrow = n, ncol = p
  ))
  names(X) <- paste0("X", 1:p)

  # Signal ONLY in first 3 features (sparse signal)
  # Other p-3 features are pure noise
  linear_pred_e <- 0.5 * X$X1 - 0.3 * X$X2 + 0.2 * X$X3
  e <- plogis(linear_pred_e)
  A <- as.integer(runif(n) < e)

  # Outcome: signal only in X1, X2, X3
  linear_pred_m0 <- -0.3 + 0.4 * X$X1 + 0.3 * X$X2 + 0.2 * X$X3
  p0 <- plogis(linear_pred_m0)
  p1 <- pmin(p0 + tau, 1)

  Y0 <- as.integer(runif(n) < p0)
  Y1 <- as.integer(runif(n) < p1)
  Y <- A * Y1 + (1 - A) * Y0

  # True ATT
  treated_idx <- which(A == 1)
  true_att <- if (length(treated_idx) > 0) {
    mean(p1[treated_idx] - p0[treated_idx])
  } else {
    tau
  }

  # Diagnostics
  n_patterns <- nrow(unique(X))  # Actual unique patterns observed

  list(
    X = X, A = A, Y = Y,
    tau = tau,
    true_att = true_att,
    true_e = e, true_m0 = p0, true_m1 = p1,
    dgp = "high_dim_sparse",
    diagnostics = list(
      p = p,
      n_possible_patterns = 2^p,
      n_observed_patterns = n_patterns,
      avg_obs_per_pattern = n / 2^p
    )
  )
}


#' Summary of Stress-Test DGPs
#'
#' DGP 4 (Weak Overlap):
#' - Assumption: Positivity stretched (e near 0.05, 0.95)
#' - Failure: Large variance (CI width 2-3× baseline)
#' - Recovery: No (inherent to design)
#' - Expected: Valid inference but inefficient
#'
#' DGP 5 (Piecewise):
#' - Assumption: Smoothness violated (step functions)
#' - Failure: Linear-DML has large bias (2-3× trees)
#' - Recovery: Linear never (misspecified); trees always good
#' - Expected: Trees excel, demonstrate flexible modeling
#'
#' DGP 6 (High-Dim):
#' - Assumption: Low dimension stretched (p=8 vs p=4)
#' - Failure: Coverage <95% at n=400, n=800
#' - Recovery: Yes, at n ≥ 1600
#' - Expected: Curse of dimensionality, sample size matters
#'
#' All DGPs maintain τ = 0.10 on probability scale for comparability.
