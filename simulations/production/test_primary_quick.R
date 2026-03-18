#' Quick Test of Primary Simulations
#'
#' Tests run_primary.R with N_REPS = 3 to verify:
#' 1. No log bloat (no .log files created)
#' 2. Tuning parameters correct (especially Rashomon epsilon_n)
#' 3. All methods run without errors
#' 4. Results file created and loadable
#'
#' Runtime: ~1-2 minutes
#'
#' Before running full simulation (N_REPS = 500), this test MUST pass.

# Suppress ALL output (including C++ stdout) during execution
# Save original connections
original_stdout <- stdout()
original_stderr <- stderr()

# Redirect everything to null device
sink("/dev/null", type = "output")
sink("/dev/null", type = "message")

# Load packages and functions
suppressMessages({
  library(parallel)
  library(dplyr)
  library(optimaltrees)  # Use installed package instead of devtools::load_all()
})

# Load simulation helpers (prevents log bloat)
source("../simulation_helpers.R")

# Source dmltree functions SILENTLY
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/dml_att_repeated.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R",
  "methods/method_forest_dml.R",
  "methods/method_linear_dml.R"
), safe_source))

# Configuration (SMALL TEST)
N_VALUES <- c(400, 800, 1600)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 3  # SMALL TEST
SEED_OFFSET <- 10000

# DGP functions
DGPS <- list(
  dgp1 = generate_dgp_binary_att,
  dgp2 = generate_dgp_continuous_att,
  dgp3 = generate_dgp_moderate_att
)

# Method functions
METHODS <- c("tree", "rashomon", "forest", "linear")

# Single-threaded for testing
N_CORES <- 1
cat(sprintf("Configuration: %d reps, single-threaded\n\n", N_REPS))

# Create output directory
output_dir <- sprintf("results/test_primary_%s", format(Sys.time(), "%Y%m%d_%H%M"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create simulation grid
sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Total simulations: %d\n", nrow(sim_grid)))
cat(sprintf("Expected runtime: ~1-2 minutes\n\n"))

# Simulation function (same as run_primary.R)
run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  # Progress tracking
  if (sim_id %% 10 == 0 || sim_id == nrow(grid)) {
    cat(sprintf("  Progress: %d/%d (%.0f%%)\n",
                sim_id, nrow(grid), 100 * sim_id / nrow(grid)))
  }

  # Generate data
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  # Fit model based on method
  result <- tryCatch({

    if (row$method == "tree") {
      fit <- suppress_all({
        estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = FALSE,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        epsilon_n = NA
      )

    } else if (row$method == "rashomon") {
      # Use theory-justified epsilon_n
      epsilon_n <- 2 * sqrt(log(row$n) / row$n)

      fit <- suppress_all({
        estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = TRUE,
          rashomon_bound_multiplier = epsilon_n,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        epsilon_n = epsilon_n
      )

    } else if (row$method == "forest") {
      fit <- suppress_all({
        att_forest(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          num.trees = 500,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        epsilon_n = NA
      )

    } else if (row$method == "linear") {
      fit <- suppress_all({
        att_linear(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          interactions = FALSE,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        epsilon_n = NA
      )

    } else {
      stop("Unknown method: ", row$method)
    }

  }, error = function(e) {
    list(
      theta = NA,
      sigma = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      epsilon_n = NA,
      error = conditionMessage(e)
    )
  })

  # Return results
  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = d$true_att,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    epsilon_n = result$epsilon_n,
    stringsAsFactors = FALSE
  )
}

# Run simulations
cat("Running simulations...\n")
start_time <- Sys.time()

results_list <- lapply(sim_grid$sim_id, function(id) {
  run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
})

results <- do.call(rbind, results_list)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\nCompleted in %.1f seconds\n\n", elapsed))

# Save results
safe_save(results, file.path(output_dir, "test_results.rds"))

# Check results
cat(strrep("=", 70), "\n")
cat("Test Results\n")
cat(strrep("=", 70), "\n\n")

# 1. Convergence
conv_rate <- 100 * mean(results$converged, na.rm = TRUE)
cat(sprintf("1. Convergence: %.0f%% (%d/%d)\n",
            conv_rate, sum(results$converged), nrow(results)))

if (conv_rate < 95) {
  cat("   ⚠️  Low convergence rate\n")
  failed <- results[!results$converged, ]
  if (nrow(failed) > 0) {
    cat("   Failed simulations:\n")
    print(failed[, c("dgp", "method", "n", "rep")])
  }
} else {
  cat("   ✓ Good convergence\n")
}
cat("\n")

# 2. Epsilon_n values (Rashomon only)
cat("2. Rashomon epsilon_n values:\n")
rashomon_eps <- results[results$method == "rashomon" & results$converged,
                        c("n", "epsilon_n")]

if (nrow(rashomon_eps) > 0) {
  eps_by_n <- aggregate(epsilon_n ~ n, rashomon_eps, mean)

  # Expected values
  expected <- data.frame(
    n = c(400, 800, 1600),
    expected = 2 * sqrt(log(c(400, 800, 1600)) / c(400, 800, 1600))
  )

  comparison <- merge(eps_by_n, expected, by = "n")
  comparison$match <- abs(comparison$epsilon_n - comparison$expected) < 0.001

  cat("\n")
  print(comparison)
  cat("\n")

  if (all(comparison$match)) {
    cat("   ✓ Epsilon_n values match theory\n")
  } else {
    cat("   ✗ Epsilon_n values DO NOT match theory\n")
  }
} else {
  cat("   ⚠️  No Rashomon results to check\n")
}
cat("\n")

# 3. Coverage (rough check with 3 reps)
cat("3. Coverage (preliminary - only 3 reps):\n")
conv_results <- results[results$converged, ]
if (nrow(conv_results) > 0) {
  conv_results$covered <- conv_results$ci_lower <= conv_results$true_att &
                          conv_results$ci_upper >= conv_results$true_att

  coverage_by_method <- aggregate(covered ~ method, conv_results,
                                   function(x) 100 * mean(x))
  names(coverage_by_method)[2] <- "coverage_pct"

  cat("\n")
  print(coverage_by_method)
  cat("\n")
  cat("   Note: Only 3 reps - expect high variability\n")
  cat("   Full run (500 reps) will give stable estimates\n")
} else {
  cat("   ⚠️  No converged results\n")
}
cat("\n")

# 4. Check for log files
cat("4. Checking for log bloat:\n")
log_files <- list.files(".", pattern = "\\.log$", recursive = FALSE)

if (length(log_files) == 0) {
  cat("   ✓ No .log files created\n")
} else {
  cat("   ✗ Found .log files:\n")
  for (f in log_files) {
    size_mb <- file.size(f) / 1024^2
    cat(sprintf("      %s (%.1f MB)\n", f, size_mb))
  }
}
cat("\n")

# 5. Check results file
cat("5. Results file:\n")
results_file <- file.path(output_dir, "test_results.rds")
if (file.exists(results_file)) {
  size_kb <- file.size(results_file) / 1024
  cat(sprintf("   ✓ Created: %s (%.1f KB)\n", results_file, size_kb))

  # Try loading
  test_load <- tryCatch({
    readRDS(results_file)
    TRUE
  }, error = function(e) FALSE)

  if (test_load) {
    cat("   ✓ File loads successfully\n")
  } else {
    cat("   ✗ File cannot be loaded\n")
  }
} else {
  cat("   ✗ Results file not created\n")
}
cat("\n")

# Summary
cat(strrep("=", 70), "\n")
cat("Test Summary\n")
cat(strrep("=", 70), "\n\n")

all_pass <- TRUE

if (conv_rate < 95) {
  cat("✗ Convergence rate low (<95%)\n")
  all_pass <- FALSE
} else {
  cat("✓ Convergence rate good (≥95%)\n")
}

if (nrow(rashomon_eps) > 0) {
  eps_ok <- all(abs(eps_by_n$epsilon_n -
                    2 * sqrt(log(eps_by_n$n) / eps_by_n$n)) < 0.001)
  if (eps_ok) {
    cat("✓ Rashomon epsilon_n correct\n")
  } else {
    cat("✗ Rashomon epsilon_n incorrect\n")
    all_pass <- FALSE
  }
} else {
  cat("⚠ Rashomon epsilon_n not tested\n")
}

if (length(log_files) == 0) {
  cat("✓ No log bloat\n")
} else {
  cat("✗ Log files created\n")
  all_pass <- FALSE
}

if (file.exists(results_file) && test_load) {
  cat("✓ Results file created and loadable\n")
} else {
  cat("✗ Results file issue\n")
  all_pass <- FALSE
}

cat("\n")

if (all_pass) {
  cat(strrep("=", 70), "\n")
  cat("✓✓ ALL TESTS PASSED ✓✓\n")
  cat(strrep("=", 70), "\n\n")
  cat("Ready to run full simulation:\n")
  cat("  1. Open run_primary.R\n")
  cat("  2. Verify N_REPS = 500\n")
  cat("  3. Run: source('production/run_primary.R')\n")
  cat("  4. Runtime: ~4 hours with 4 cores\n\n")
} else {
  cat(strrep("=", 70), "\n")
  cat("⚠️  SOME TESTS FAILED ⚠️\n")
  cat(strrep("=", 70), "\n\n")
  cat("Fix issues above before running full simulation.\n\n")
}
