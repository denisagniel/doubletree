# test_feature_spaces.R
# Quick test of different feature spaces to make CV regularization informative

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

message("=== Testing Feature Spaces for CV Regularization ===\n")

n <- 400
tau <- 0.15
K <- 5
n_reps <- 5  # Very quick test

# DGP 1: Original (2 binary features) - BASELINE
generate_data_2binary <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.5 * X1 - 0.2 * X2)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X1 + 0.1 * X2))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X1 + 0.1 * X2 + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau)
}

# DGP 2: 6 binary features (64 covariate patterns)
generate_data_6binary <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5),
    X5 = as.integer(runif(n) < 0.5),
    X6 = as.integer(runif(n) < 0.5)
  )
  # True function depends only on X1, X2, X3 (others are noise)
  e <- plogis(0.5 * X$X1 - 0.2 * X$X2 + 0.1 * X$X3)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X$X1 + 0.1 * X$X2))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X$X1 + 0.1 * X$X2 + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau)
}

# DGP 3: 2 continuous features (discretized to 4 bins each = 16 patterns)
generate_data_2continuous <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1_cont <- runif(n)
  X2_cont <- runif(n)
  # Discretize to 4 bins for tree fitting
  X <- data.frame(
    X1 = as.integer(cut(X1_cont, breaks = 4, labels = FALSE)) - 1L,
    X2 = as.integer(cut(X2_cont, breaks = 4, labels = FALSE)) - 1L
  )
  # Convert to binary representation (2 bits each = 4 binary features)
  X <- data.frame(
    X1_bit1 = as.integer(X$X1 %% 2),
    X1_bit2 = as.integer(X$X1 %/% 2),
    X2_bit1 = as.integer(X$X2 %% 2),
    X2_bit2 = as.integer(X$X2 %/% 2)
  )
  # Use original continuous values for propensity
  e <- plogis(0.5 * X1_cont - 0.2 * X2_cont)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X1_cont + 0.1 * X2_cont))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X1_cont + 0.1 * X2_cont + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau)
}

# DGP 4: 8 binary features (256 patterns) - stress test
generate_data_8binary <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- data.frame(
    X1 = as.integer(runif(n) < 0.5),
    X2 = as.integer(runif(n) < 0.5),
    X3 = as.integer(runif(n) < 0.5),
    X4 = as.integer(runif(n) < 0.5),
    X5 = as.integer(runif(n) < 0.5),
    X6 = as.integer(runif(n) < 0.5),
    X7 = as.integer(runif(n) < 0.5),
    X8 = as.integer(runif(n) < 0.5)
  )
  # True function uses X1, X2, X3 (others are noise)
  e <- plogis(0.5 * X$X1 - 0.2 * X$X2 + 0.1 * X$X3)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X$X1 + 0.1 * X$X2))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X$X1 + 0.1 * X$X2 + tau))
  Y <- A * Y1 + (1 - A) * Y0
  list(X = X, A = A, Y = Y, tau = tau)
}

test_dgp <- function(dgp_fn, dgp_name, n, tau, n_reps) {
  message("\n=== ", dgp_name, " ===")

  # Test 1: Check CV loss variability
  d <- dgp_fn(n, tau, seed = 123)
  message("  Testing CV loss variability:")
  cv_e <- optimaltrees::cv_regularization(d$X, d$A, loss_function = "log_loss",
                                       K = 5, refit = FALSE, verbose = FALSE)
  cv_range <- max(cv_e$cv_loss) - min(cv_e$cv_loss)
  cv_rel_range <- cv_range / mean(cv_e$cv_loss)

  message(sprintf("    Lambda grid: %s", paste(round(cv_e$lambda_grid, 5), collapse=", ")))
  message(sprintf("    CV loss range: %.4f (%.1f%% of mean)", cv_range, cv_rel_range * 100))
  message(sprintf("    Selected: %.5f, Fixed: %.5f", cv_e$best_lambda, log(n)/n))

  is_informative <- cv_rel_range > 0.01  # >1% variation
  message(sprintf("    CV informative? %s", ifelse(is_informative, "YES", "NO (flat)")))

  if (!is_informative) {
    message("    Skipping coverage test - CV not informative")
    return(list(
      dgp = dgp_name,
      cv_informative = FALSE,
      cv_rel_range = cv_rel_range,
      coverage_fixed = NA,
      coverage_cv = NA,
      n_features = ncol(d$X)
    ))
  }

  # Test 2: Quick coverage test (5 reps)
  message(sprintf("  Running %d reps for coverage test...", n_reps))

  results_fixed <- list()
  results_cv <- list()

  for (i in seq_len(n_reps)) {
    d <- dgp_fn(n, tau, seed = i * 100)

    results_fixed[[i]] <- estimate_att(
      d$X, d$A, d$Y, K = K,
      use_rashomon = FALSE,
      regularization = log(n)/n,
      verbose = FALSE,
      seed = i * 100
    )

    results_cv[[i]] <- estimate_att(
      d$X, d$A, d$Y, K = K,
      use_rashomon = FALSE,
      cv_regularization = TRUE,
      cv_K = 5,
      verbose = FALSE,
      seed = i * 100
    )
  }

  compute_coverage <- function(results, true_att) {
    ci_low <- sapply(results, function(r) r$ci_95[1])
    ci_high <- sapply(results, function(r) r$ci_95[2])
    mean(ci_low <= true_att & true_att <= ci_high)
  }

  cov_fixed <- compute_coverage(results_fixed, tau)
  cov_cv <- compute_coverage(results_cv, tau)

  message(sprintf("    Coverage (fixed): %.2f", cov_fixed))
  message(sprintf("    Coverage (CV):    %.2f", cov_cv))

  list(
    dgp = dgp_name,
    cv_informative = TRUE,
    cv_rel_range = cv_rel_range,
    coverage_fixed = cov_fixed,
    coverage_cv = cov_cv,
    n_features = ncol(d$X)
  )
}

# Run tests
results <- list(
  test_dgp(generate_data_2binary, "2 binary features (4 patterns)", n, tau, n_reps),
  test_dgp(generate_data_6binary, "6 binary features (64 patterns)", n, tau, n_reps),
  test_dgp(generate_data_2continuous, "2 continuous → 4 binary (16 patterns)", n, tau, n_reps),
  test_dgp(generate_data_8binary, "8 binary features (256 patterns)", n, tau, n_reps)
)

message("\n\n=== SUMMARY ===\n")
for (r in results) {
  message(sprintf("%s:", r$dgp))
  message(sprintf("  Features: %d", r$n_features))
  message(sprintf("  CV informative: %s (%.1f%% variation)",
                  ifelse(r$cv_informative, "YES", "NO"), r$cv_rel_range * 100))
  if (r$cv_informative) {
    message(sprintf("  Coverage - Fixed: %.2f, CV: %.2f", r$coverage_fixed, r$coverage_cv))
  }
  message("")
}

message("=== Recommendation ===")
best <- which.max(sapply(results, function(r) {
  if (!r$cv_informative) return(-Inf)
  # Score by: CV works AND coverage is good
  r$cv_rel_range * 10 + ifelse(!is.na(r$coverage_cv), r$coverage_cv, 0)
}))

if (best > 0 && results[[best]]$cv_informative) {
  message(sprintf("Best option: %s", results[[best]]$dgp))
  message("CV is informative and coverage is reasonable")
} else {
  message("None of these feature spaces make CV meaningfully informative")
  message("Recommend: Use fixed regularization with theory-driven value")
}
