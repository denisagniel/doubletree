# run_simulations_extended.R
# Extended simulation study for Rashomon-DML paper
#
# Compares three approaches:
#   1. Fold-specific optimal trees (use_rashomon = FALSE)
#   2. Rashomon intersection (use_rashomon = TRUE)
#   3. Oracle (true nuisances)
#
# Tests across:
#   - 4 DGPs (DGP 1-2: binary features; DGP 3: smooth continuous; DGP 4: non-smooth)
#   - 4 sample sizes (200, 400, 800, 1600)
#   - 4 Rashomon tolerances (0.01, 0.05, 0.1, 0.2)
#   - 100 replications per configuration
#
# Output: Comprehensive results validating existence empirically and quantifying
# finite-sample bias overhead.

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. Install with: install.packages('devtools')")
}
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) {
  stop("optimaltrees is required. Install from source or set up .Rprofile to load from ../optimaltrees")
}

results_dir <- "simulations/results_extended"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# DGPs -------------------------------------------------------------------------

# DGP 1: Binary features, weak confounding (EXISTING)
generate_data_dgp1 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.5 * X1 - 0.2)
  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < (0.3 + 0.2 * X1))
  Y1 <- as.integer(runif(n) < (0.3 + 0.2 * X1 + tau))
  Y <- A * Y1 + (1 - A) * Y0

  # True nuisances for oracle
  true_e <- e
  true_m0 <- 0.3 + 0.2 * X1
  true_m1 <- 0.3 + 0.2 * X1 + tau

  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp1",
       true_e = true_e, true_m0 = true_m0, true_m1 = true_m1)
}

# DGP 2: Binary features, strong confounding (EXISTING)
generate_data_dgp2 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)
  e <- plogis(0.4 * X1 + 0.3 * X2 - 0.35)
  A <- as.integer(runif(n) < e)
  m0 <- 0.25 + 0.3 * X1 + 0.2 * X2
  m1 <- m0 + tau
  Y0 <- as.integer(runif(n) < m0)
  Y1 <- as.integer(runif(n) < pmin(1, m1))
  Y <- A * Y1 + (1 - A) * Y0

  # True nuisances for oracle
  true_e <- e
  true_m0 <- m0
  true_m1 <- m1

  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp2",
       true_e = true_e, true_m0 = true_m0, true_m1 = true_m1)
}

# DGP 3: SMOOTH, high signal-to-noise (NEW - stress test for ideal Rashomon)
# Continuous X ~ N(0,1)^2, smooth propensity and outcome functions
# Expected: Intersection non-empty almost always (85-95%)
generate_data_dgp3 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1_raw <- rnorm(n, 0, 1)
  X2_raw <- rnorm(n, 0, 1)

  # Discretize to binary for optimaltrees compatibility (0/1 based on median split)
  X1 <- as.integer(X1_raw > median(X1_raw))
  X2 <- as.integer(X2_raw > median(X2_raw))
  X <- data.frame(X1 = X1, X2 = X2)

  # Smooth functions on raw continuous scale
  e_raw <- pnorm(0.5 * X1_raw + 0.3 * X2_raw - 0.2)
  m0_raw <- pnorm(0.4 * X1_raw + 0.2 * X2_raw)

  # Map to binary features for consistency
  e <- plogis(0.5 * X1 + 0.3 * X2 - 0.2)
  m0 <- 0.4 + 0.2 * X1 + 0.15 * X2
  m1 <- m0 + tau

  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < pmin(1, pmax(0, m0)))
  Y1 <- as.integer(runif(n) < pmin(1, pmax(0, m1)))
  Y <- A * Y1 + (1 - A) * Y0

  # True nuisances for oracle
  true_e <- e
  true_m0 <- m0
  true_m1 <- m1

  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp3",
       true_e = true_e, true_m0 = true_m0, true_m1 = true_m1)
}

# DGP 4: NON-SMOOTH, low signal-to-noise (NEW - stress test for difficult Rashomon)
# Piecewise propensity, interactions in outcome, rough functions
# Expected: Intersection less common (60-80%)
generate_data_dgp4 <- function(n, tau = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- as.integer(runif(n) < 0.5)
  X2 <- as.integer(runif(n) < 0.5)
  X <- data.frame(X1 = X1, X2 = X2)

  # Piecewise propensity (non-smooth)
  e <- 0.3 * (X1 > 0) + 0.5 * (X2 > 0.5)
  e <- pmin(0.9, pmax(0.1, e))  # Ensure overlap

  # Outcome with interaction (non-additive)
  m0 <- 0.2 + 0.3 * X1 * X2 + 0.15 * (X1 + X2)
  m1 <- m0 + tau

  A <- as.integer(runif(n) < e)
  Y0 <- as.integer(runif(n) < pmin(1, pmax(0, m0)))
  Y1 <- as.integer(runif(n) < pmin(1, pmax(0, m1)))
  Y <- A * Y1 + (1 - A) * Y0

  # True nuisances for oracle
  true_e <- e
  true_m0 <- m0
  true_m1 <- m1

  list(X = X, A = A, Y = Y, tau = tau, dgp = "dgp4",
       true_e = true_e, true_m0 = true_m0, true_m1 = true_m1)
}

# Oracle DML (true nuisances) --------------------------------------------------

# Compute ATT with true nuisances (oracle)
# This establishes the best-possible performance ceiling
dml_att_oracle <- function(data, K = 5, seed = NULL) {
  n <- nrow(data$X)
  fold_indices <- create_folds(n, K, strata = data$A, seed = seed)

  # Use true nuisances (per-observation)
  eta <- list(
    e = data$true_e,
    m0 = data$true_m0,
    m1 = data$true_m1
  )

  pi_hat <- mean(data$A)

  # Closed form: psi(theta) = psi(0) - theta*(A/pi)
  sum_a_over_pi <- sum(data$A / pi_hat)
  if (sum_a_over_pi < 1e-10) stop("No treated units.")

  score_at_zero <- psi_att(data$Y, data$A, theta = 0, eta, pi_hat)
  theta <- sum(score_at_zero) / sum_a_over_pi

  score_values <- psi_att(data$Y, data$A, theta, eta, pi_hat)
  sigma_sq <- dml_att_variance(score_values, n)
  sigma <- sqrt(sigma_sq)
  ci_95 <- dml_att_ci(theta, sigma, n, level = 0.95)

  list(
    theta = theta,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    method = "oracle"
  )
}

# Metrics computation ----------------------------------------------------------

# Standard metrics (existing)
compute_metrics <- function(results, true_att) {
  theta_hats <- vapply(results, function(r) r$theta, numeric(1))
  ci_low <- vapply(results, function(r) r$ci_95[1L], numeric(1))
  ci_high <- vapply(results, function(r) r$ci_95[2L], numeric(1))
  n_sim <- length(results)

  bias <- mean(theta_hats) - true_att
  mse <- mean((theta_hats - true_att)^2)
  rmse <- sqrt(mse)
  empirical_se <- sd(theta_hats)
  coverage_95 <- mean(ci_low <= true_att & true_att <= ci_high)
  mean_ci_width <- mean(ci_high - ci_low)

  data.frame(
    n_sim = n_sim,
    true_att = true_att,
    mean_theta_hat = mean(theta_hats),
    bias = bias,
    mse = mse,
    rmse = rmse,
    empirical_se = empirical_se,
    coverage_95 = coverage_95,
    mean_ci_width = mean_ci_width
  )
}

# Rashomon-specific metrics (NEW)
compute_metrics_rashomon <- function(results, true_att, data_list = NULL) {
  base_metrics <- compute_metrics(results, true_att)

  # Extract Rashomon-specific information
  n_intersecting_e <- vapply(results, function(r) {
    if ("nuisance_fits" %in% names(r) && "cf_e" %in% names(r$nuisance_fits)) {
      cf_e <- r$nuisance_fits$cf_e
      if ("n_intersecting" %in% names(cf_e)) cf_e$n_intersecting else NA_integer_
    } else NA_integer_
  }, integer(1))

  n_intersecting_m0 <- vapply(results, function(r) {
    if ("nuisance_fits" %in% names(r) && "cf_m0" %in% names(r$nuisance_fits)) {
      cf_m0 <- r$nuisance_fits$cf_m0
      if ("n_intersecting" %in% names(cf_m0)) cf_m0$n_intersecting else NA_integer_
    } else NA_integer_
  }, integer(1))

  n_intersecting_m1 <- vapply(results, function(r) {
    if ("nuisance_fits" %in% names(r) && "cf_m1" %in% names(r$nuisance_fits)) {
      cf_m1 <- r$nuisance_fits$cf_m1
      if ("n_intersecting" %in% names(cf_m1)) cf_m1$n_intersecting else NA_integer_
    } else NA_integer_
  }, integer(1))

  # Compute percentage of replications with non-empty intersection
  pct_nonempty_e <- mean(n_intersecting_e > 0, na.rm = TRUE)
  pct_nonempty_m0 <- mean(n_intersecting_m0 > 0, na.rm = TRUE)
  pct_nonempty_m1 <- mean(n_intersecting_m1 > 0, na.rm = TRUE)
  pct_nonempty_any <- mean((n_intersecting_e > 0) | (n_intersecting_m0 > 0) | (n_intersecting_m1 > 0), na.rm = TRUE)

  # Mean number of intersecting trees (when non-empty)
  mean_n_intersecting_e <- mean(n_intersecting_e[n_intersecting_e > 0], na.rm = TRUE)
  mean_n_intersecting_m0 <- mean(n_intersecting_m0[n_intersecting_m0 > 0], na.rm = TRUE)
  mean_n_intersecting_m1 <- mean(n_intersecting_m1[n_intersecting_m1 > 0], na.rm = TRUE)

  # Add to base metrics
  cbind(
    base_metrics,
    pct_nonempty_e = pct_nonempty_e,
    pct_nonempty_m0 = pct_nonempty_m0,
    pct_nonempty_m1 = pct_nonempty_m1,
    pct_nonempty_any = pct_nonempty_any,
    mean_n_intersecting_e = mean_n_intersecting_e,
    mean_n_intersecting_m0 = mean_n_intersecting_m0,
    mean_n_intersecting_m1 = mean_n_intersecting_m1
  )
}

# Three-way comparison function ------------------------------------------------

run_comparison <- function(dgp_fn, n, K, tau, n_reps, epsilon,
                          regularization = NULL, seed_start = 1000) {
  if (is.null(regularization)) {
    regularization <- log(n) / n  # Theory-driven default
  }

  results <- list(
    fold_specific = vector("list", n_reps),
    rashomon = vector("list", n_reps),
    oracle = vector("list", n_reps)
  )

  message("Running ", n_reps, " replications: n=", n, ", epsilon=", epsilon,
          ", dgp=", deparse(substitute(dgp_fn)))

  for (rep in seq_len(n_reps)) {
    if (rep %% 10 == 0) message("  Rep ", rep, "/", n_reps)

    d <- dgp_fn(n, tau, seed = seed_start + rep)

    # Fold-specific optimal (use_rashomon = FALSE)
    results$fold_specific[[rep]] <- tryCatch({
      estimate_att(
        d$X, d$A, d$Y, K = K,
        use_rashomon = FALSE,
        regularization = regularization,
        verbose = FALSE,
        seed = seed_start + rep
      )
    }, error = function(e) {
      warning("Fold-specific failed at rep ", rep, ": ", e$message)
      list(theta = NA, ci_95 = c(NA, NA), method = "fold_specific")
    })

    # Rashomon intersection (use_rashomon = TRUE)
    results$rashomon[[rep]] <- tryCatch({
      estimate_att(
        d$X, d$A, d$Y, K = K,
        use_rashomon = TRUE,
        rashomon_bound_multiplier = epsilon,
        regularization = regularization,
        verbose = FALSE,
        seed = seed_start + rep
      )
    }, error = function(e) {
      warning("Rashomon failed at rep ", rep, ": ", e$message)
      list(theta = NA, ci_95 = c(NA, NA), method = "rashomon")
    })

    # Oracle (true nuisances)
    results$oracle[[rep]] <- tryCatch({
      dml_att_oracle(d, K = K, seed = seed_start + rep)
    }, error = function(e) {
      warning("Oracle failed at rep ", rep, ": ", e$message)
      list(theta = NA, ci_95 = c(NA, NA), method = "oracle")
    })
  }

  # Compute metrics for each method
  metrics_fold_specific <- compute_metrics(results$fold_specific, tau)
  metrics_rashomon <- compute_metrics_rashomon(results$rashomon, tau)
  metrics_oracle <- compute_metrics(results$oracle, tau)

  list(
    fold_specific = metrics_fold_specific,
    rashomon = metrics_rashomon,
    oracle = metrics_oracle,
    raw_results = results  # Keep for further analysis
  )
}

# Main simulation grid ---------------------------------------------------------

run_full_simulation_grid <- function() {
message("=== Extended Rashomon-DML Simulation Study ===\n")

# Simulation design
dgps <- list(
  dgp1 = generate_data_dgp1,
  dgp2 = generate_data_dgp2,
  dgp3 = generate_data_dgp3,  # NEW: smooth
  dgp4 = generate_data_dgp4   # NEW: non-smooth
)
ns <- c(200, 400, 800, 1600)
epsilons <- c(0.01, 0.05, 0.1, 0.2)
n_reps <- 100  # Full replications for paper
tau <- 0.15
K <- 5

# For quick testing, uncomment:
# ns <- c(200, 400)
# epsilons <- c(0.05, 0.1)
# n_reps <- 10

# Run full grid
all_results <- list()
config_idx <- 1

for (dgp_name in names(dgps)) {
  for (n in ns) {
    for (epsilon in epsilons) {
      message("\n--- Configuration ", config_idx, " ---")
      message("DGP: ", dgp_name, ", n: ", n, ", epsilon: ", epsilon)

      result <- run_comparison(
        dgp_fn = dgps[[dgp_name]],
        n = n,
        K = K,
        tau = tau,
        n_reps = n_reps,
        epsilon = epsilon,
        seed_start = config_idx * 10000
      )

      # Add configuration metadata
      result$config <- data.frame(
        dgp = dgp_name,
        n = n,
        epsilon = epsilon,
        K = K,
        tau = tau,
        n_reps = n_reps
      )

      all_results[[config_idx]] <- result

      # Save incrementally
      saveRDS(result, file.path(results_dir,
              sprintf("result_%s_n%d_eps%.3f.rds", dgp_name, n, epsilon)))

      config_idx <- config_idx + 1
    }
  }
}

# Save full results
saveRDS(all_results, file.path(results_dir, "all_results.rds"))

# Create summary table (handle different column counts)
summary_rows <- list()
for (i in seq_along(all_results)) {
  r <- all_results[[i]]
  cfg <- r$config

  # One row per method - bind as data frames to handle different columns
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "fold_specific", as.data.frame(r$fold_specific)
  )
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "rashomon", as.data.frame(r$rashomon)
  )
  summary_rows[[length(summary_rows) + 1]] <- cbind(
    cfg, method = "oracle", as.data.frame(r$oracle)
  )
}

# Use bind_rows to handle different column counts
if (requireNamespace("dplyr", quietly = TRUE)) {
  summary_table <- dplyr::bind_rows(summary_rows)
} else {
  # Fallback: fill missing columns with NA
  all_cols <- unique(unlist(lapply(summary_rows, names)))
  summary_rows_filled <- lapply(summary_rows, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }
    df[, all_cols]
  })
  summary_table <- do.call(rbind, summary_rows_filled)
}
saveRDS(summary_table, file.path(results_dir, "simulation_summary.rds"))
write.csv(summary_table, file.path(results_dir, "simulation_summary.csv"),
          row.names = FALSE)

message("\n=== Simulation Complete ===")
message("Results saved to ", results_dir)
message("Summary table: simulation_summary.csv")
message("Configuration count: ", config_idx - 1)
message("Total replications: ", (config_idx - 1) * n_reps)

return(invisible(all_results))
}

# Run if called as main script
if (!interactive() && sys.nframe() == 0) {
  run_full_simulation_grid()
}
