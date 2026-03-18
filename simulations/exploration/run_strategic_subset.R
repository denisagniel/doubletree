# run_strategic_subset.R
# Strategic subset for fast iteration (Option B)
#
# Scope: 2 DGPs × 3 n × 2 ε × 100 reps = 1,200 runs (~1 hour)
# - DGP 3 (smooth, high SNR) - ideal case
# - DGP 4 (non-smooth, low SNR) - stress test
# - n ∈ {200, 400, 800}
# - ε ∈ {0.05, 0.1}
#
# Purpose: Quick validation of core claims before committing to full grid

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required. Install with: install.packages('devtools')")
}
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) {
  stop("optimaltrees is required. Install from source or set up .Rprofile")
}

library(cli)

results_dir <- "simulations/results_strategic"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Load DGPs from extended script -----------------------------------------------

source("simulations/run_simulations_extended.R", local = TRUE)

# Simulation grid (strategic subset) -------------------------------------------

dgp_list <- c("dgp3", "dgp4")
sample_sizes <- c(200, 400, 800)
epsilon_vals <- c(0.05, 0.1)
n_reps <- 100

# Total runs: 2 × 3 × 2 × 100 = 1,200
total_runs <- length(dgp_list) * length(sample_sizes) * length(epsilon_vals) * n_reps
cli_alert_info("Strategic subset: {total_runs} runs")
cli_alert_info("DGPs: {paste(dgp_list, collapse = ', ')}")
cli_alert_info("Sample sizes: {paste(sample_sizes, collapse = ', ')}")
cli_alert_info("Tolerances: {paste(epsilon_vals, collapse = ', ')}")
cli_alert_info("Replications: {n_reps}")

start_time <- Sys.time()

# Run simulations --------------------------------------------------------------

run_counter <- 0

for (dgp_name in dgp_list) {
  for (n in sample_sizes) {
    for (eps in epsilon_vals) {

      run_counter <- run_counter + 1
      pct_complete <- round(100 * run_counter / (length(dgp_list) * length(sample_sizes) * length(epsilon_vals)))

      cli_h2("Config {run_counter}/{length(dgp_list) * length(sample_sizes) * length(epsilon_vals)}: {dgp_name}, n={n}, eps={eps} ({pct_complete}% complete)")

      # Select DGP function
      dgp_fn <- switch(dgp_name,
                       dgp3 = generate_data_dgp3,
                       dgp4 = generate_data_dgp4,
                       stop("Unknown DGP: ", dgp_name))

      # Storage for this config
      results_fold_specific <- vector("list", n_reps)
      results_rashomon <- vector("list", n_reps)
      results_oracle <- vector("list", n_reps)

      cli_progress_bar("Replications", total = n_reps)

      for (rep in 1:n_reps) {
        seed <- 1000 * run_counter + rep

        # Generate data
        data <- dgp_fn(n = n, seed = seed)

        # Method 1: Fold-specific optimal trees (baseline)
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
          results_fold_specific[[rep]]$true_tau <- data$tau
        }, error = function(e) {
          cli_alert_warning("Fold-specific failed (rep {rep}): {e$message}")
          results_fold_specific[[rep]] <- NULL
        })

        # Method 2: Rashomon intersection
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
          results_rashomon[[rep]]$epsilon <- eps
          results_rashomon[[rep]]$true_tau <- data$tau
        }, error = function(e) {
          cli_alert_warning("Rashomon failed (rep {rep}): {e$message}")
          results_rashomon[[rep]] <- NULL
        })

        # Method 3: Oracle (true nuisances)
        tryCatch({
          results_oracle[[rep]] <- dml_att_oracle(data, K = 5, seed = seed)
          results_oracle[[rep]]$method <- "oracle"
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
      filename <- sprintf("result_%s_n%d_eps%.3f.rds", dgp_name, n, eps)
      saveRDS(combined, file.path(results_dir, filename))
      cli_alert_success("Saved {filename}")
    }
  }
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cli_alert_success("Strategic subset complete!")
cli_alert_info("Total time: {round(elapsed, 1)} minutes")
cli_alert_info("Results saved to: {results_dir}")
cli_alert_info("Next: source('simulations/analyze_strategic_results.R')")
