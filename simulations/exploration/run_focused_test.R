# run_focused_test.R
# Focused test: Large n, smooth DGP, theory-guided epsilon
#
# Goal: Get theory working with favorable conditions
# - DGP 3 only (smooth, high SNR - best case)
# - Large n: {800, 1600, 3200} (easier for theory)
# - Epsilon: c * sqrt(log(n)/n) with c ∈ {1, 2, 3} (theory-guided)
#
# Total: 1 DGP × 3 n × 3 ε × 100 reps = 900 runs (~45 min)

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. Install with: install.packages('devtools')")
}
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) {
  stop("optimaltrees is required. Install from source or set up .Rprofile")
}

library(cli)

results_dir <- "simulations/results_focused"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Load DGPs from extended script -----------------------------------------------

source("simulations/run_simulations_extended.R", local = TRUE)

# Theory-guided epsilon function -----------------------------------------------

epsilon_theory <- function(n, c = 2) {
  # ε_n = c * sqrt(log(n)/n)
  # Default c = 2 (moderate choice from {1, 2, 3})
  c * sqrt(log(n) / n)
}

# Simulation grid --------------------------------------------------------------

dgp_name <- "dgp3"  # Smooth, high SNR only
sample_sizes <- c(800, 1600, 3200)
c_values <- c(1, 2, 3)  # Multipliers for epsilon
n_reps <- 100

# Compute epsilon values
grid <- expand.grid(n = sample_sizes, c = c_values)
grid$epsilon <- mapply(epsilon_theory, grid$n, grid$c)

# Total runs: 3 n × 3 c × 100 reps = 900
total_configs <- nrow(grid)
total_runs <- total_configs * n_reps

cli_h1("Focused Test: Large n, Smooth DGP, Theory-Guided Epsilon")
cli_alert_info("Total runs: {total_runs}")
cli_alert_info("DGP: {dgp_name} (smooth, high SNR)")
cli_alert_info("Sample sizes: {paste(sample_sizes, collapse = ', ')}")
cli_alert_info("Epsilon formula: c * sqrt(log(n)/n), c ∈ {{1, 2, 3}}")
cli_text("")
cli_text("Epsilon grid:")
for (i in 1:nrow(grid)) {
  cli_text("  n={grid$n[i]}, c={grid$c[i]} → ε={round(grid$epsilon[i], 4)}")
}

start_time <- Sys.time()

# Run simulations --------------------------------------------------------------

for (i in 1:nrow(grid)) {
  n <- grid$n[i]
  c_val <- grid$c[i]
  eps <- grid$epsilon[i]

  pct_complete <- round(100 * i / total_configs)

  cli_h2("Config {i}/{total_configs}: n={n}, c={c_val}, ε={round(eps, 4)} ({pct_complete}% complete)")

  # Storage for this config
  results_fold_specific <- vector("list", n_reps)
  results_rashomon <- vector("list", n_reps)
  results_oracle <- vector("list", n_reps)

  cli_progress_bar("Replications", total = n_reps)

  for (rep in 1:n_reps) {
    seed <- 1000 * i + rep

    # Generate data
    data <- generate_data_dgp3(n = n, seed = seed)

    # Method 1: Fold-specific optimal trees
    tryCatch({
      results_fold_specific[[rep]] <- estimate_att(
        X = data$X,
        A = data$A,
        Y = data$Y,
        use_rashomon = FALSE,
        K = 5,
        seed = seed
      )
      results_fold_specific[[rep]]$method <- "fold_specific"
      results_fold_specific[[rep]]$dgp <- dgp_name
      results_fold_specific[[rep]]$n <- n
      results_fold_specific[[rep]]$true_tau <- data$tau
    }, error = function(e) {
      cli_alert_warning("Fold-specific failed (rep {rep}): {e$message}")
      results_fold_specific[[rep]] <- NULL
    })

    # Method 2: Rashomon intersection with theory-guided epsilon
    tryCatch({
      results_rashomon[[rep]] <- estimate_att(
        X = data$X,
        A = data$A,
        Y = data$Y,
        use_rashomon = TRUE,
        rashomon_epsilon = eps,
        K = 5,
        seed = seed
      )
      results_rashomon[[rep]]$method <- "rashomon"
      results_rashomon[[rep]]$dgp <- dgp_name
      results_rashomon[[rep]]$n <- n
      results_rashomon[[rep]]$epsilon <- eps
      results_rashomon[[rep]]$c_value <- c_val
      results_rashomon[[rep]]$true_tau <- data$tau
    }, error = function(e) {
      cli_alert_warning("Rashomon failed (rep {rep}): {e$message}")
      results_rashomon[[rep]] <- NULL
    })

    # Method 3: Oracle (true nuisances)
    tryCatch({
      results_oracle[[rep]] <- dml_att_oracle(data, K = 5, seed = seed)
      results_oracle[[rep]]$method <- "oracle"
      results_oracle[[rep]]$dgp <- dgp_name
      results_oracle[[rep]]$n <- n
      results_oracle[[rep]]$true_tau <- data$tau
    }, error = function(e) {
      cli_alert_warning("Oracle failed (rep {rep}): {e$message}")
      results_oracle[[rep]] <- NULL
    })

    cli_progress_update()
  }

  cli_progress_done()

  # Combine results for this config
  combined <- c(
    Filter(Negate(is.null), results_fold_specific),
    Filter(Negate(is.null), results_rashomon),
    Filter(Negate(is.null), results_oracle)
  )

  # Save
  filename <- sprintf("result_%s_n%d_c%d_eps%.4f.rds",
                      dgp_name, n, c_val, eps)
  saveRDS(combined, file.path(results_dir, filename))
  cli_alert_success("Saved {filename}")
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cli_h1("Focused Test Complete!")
cli_alert_success("Total time: {round(elapsed, 1)} minutes")
cli_alert_info("Results saved to: {results_dir}")
cli_alert_info("Next: source('simulations/analyze_focused_results.R')")

# Quick preview of results -----------------------------------------------------

cli_h2("Quick Preview")

all_files <- list.files(results_dir, pattern = "^result_.*\\.rds$", full.names = TRUE)
if (length(all_files) > 0) {
  # Load first file as sample
  sample_result <- readRDS(all_files[1])

  # Check if Rashomon metadata is present
  rashomon_results <- Filter(function(x) !is.null(x) && x$method == "rashomon", sample_result)
  if (length(rashomon_results) > 0) {
    r <- rashomon_results[[1]]
    cli_text("Sample Rashomon result structure:")
    cli_text("  Has pct_nonempty_e: {!is.null(r$pct_nonempty_e)}")
    cli_text("  Has n_intersecting_e: {!is.null(r$n_intersecting_e)}")
    if (!is.null(r$pct_nonempty_e)) {
      cli_text("  pct_nonempty_e value: {r$pct_nonempty_e}")
    }
  }
}
