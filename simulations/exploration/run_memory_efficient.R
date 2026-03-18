# run_memory_efficient.R
# Memory-efficient simulation: don't store raw results, only metrics
#
# Key changes:
# 1. Don't store raw_results (saves ~90% memory)
# 2. Force garbage collection after each config
# 3. Process and discard results immediately
# 4. Only keep summary metrics in memory

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/run_simulations_extended.R", local = TRUE)

results_dir <- "simulations/results_extended"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

message("=== Memory-Efficient Rashomon-DML Simulation ===\n")

# Simulation design
dgps <- list(
  dgp1 = generate_data_dgp1,
  dgp2 = generate_data_dgp2,
  dgp3 = generate_data_dgp3,
  dgp4 = generate_data_dgp4
)
ns <- c(200, 400, 800, 1600)
epsilons <- c(0.01, 0.05, 0.1, 0.2)
n_reps <- 100
tau <- 0.15
K <- 5

# For testing, use smaller grid:
# dgps <- list(dgp1 = generate_data_dgp1, dgp3 = generate_data_dgp3)
# ns <- c(200, 400)
# epsilons <- c(0.05, 0.1)
# n_reps <- 10

total_configs <- length(dgps) * length(ns) * length(epsilons)
message("Total configurations: ", total_configs)
message("Replications per config: ", n_reps)
message("Total runs: ", total_configs * n_reps, "\n")

# Track memory
initial_mem <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024
message("Initial memory: ", round(initial_mem, 1), " MB\n")

# Run simulations with memory efficiency
summary_rows <- list()
config_idx <- 1
start_time <- Sys.time()

for (dgp_name in names(dgps)) {
  for (n in ns) {
    for (epsilon in epsilons) {
      message("\n--- Configuration ", config_idx, "/", total_configs, " ---")
      message("DGP: ", dgp_name, ", n: ", n, ", epsilon: ", epsilon)

      config_start <- Sys.time()

      # Run comparison (returns metrics + raw_results)
      result <- run_comparison(
        dgp_fn = dgps[[dgp_name]],
        n = n,
        K = K,
        tau = tau,
        n_reps = n_reps,
        epsilon = epsilon,
        seed_start = config_idx * 10000
      )

      config_elapsed <- as.numeric(difftime(Sys.time(), config_start, units = "secs"))

      # Extract config metadata
      cfg <- data.frame(
        dgp = dgp_name,
        n = n,
        epsilon = epsilon,
        K = K,
        tau = tau,
        n_reps = n_reps
      )

      # Save ONLY metrics (not raw results) to conserve memory
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        cfg, method = "fold_specific", as.data.frame(result$fold_specific)
      )
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        cfg, method = "rashomon", as.data.frame(result$rashomon)
      )
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        cfg, method = "oracle", as.data.frame(result$oracle)
      )

      # Save individual config results to disk (in case of crash)
      # But save only metrics, not raw_results
      result_to_save <- list(
        config = cfg,
        fold_specific = result$fold_specific,
        rashomon = result$rashomon,
        oracle = result$oracle
        # NOTE: raw_results NOT saved to reduce file size
      )
      saveRDS(result_to_save, file.path(results_dir,
              sprintf("result_%s_n%d_eps%.3f.rds", dgp_name, n, epsilon)))

      # CRITICAL: Clear large objects and force garbage collection
      rm(result, result_to_save)
      gc(verbose = FALSE)

      # Report progress
      elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      configs_remaining <- total_configs - config_idx
      estimated_remaining <- (elapsed_total / config_idx) * configs_remaining

      current_mem <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024

      message(sprintf("  Config time: %.1f sec", config_elapsed))
      message(sprintf("  Current memory: %.1f MB (delta: %.1f MB)",
                     current_mem, current_mem - initial_mem))
      message(sprintf("  Progress: %d/%d (%.1f%%), elapsed: %.1f min, ETA: %.1f min",
                     config_idx, total_configs,
                     100 * config_idx / total_configs,
                     elapsed_total, estimated_remaining))

      config_idx <- config_idx + 1
    }
  }
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
final_mem <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024

message("\n=== Simulation Complete ===")
message("Total time: ", round(elapsed, 2), " minutes (", round(elapsed/60, 2), " hours)")
message("Initial memory: ", round(initial_mem, 1), " MB")
message("Final memory: ", round(final_mem, 1), " MB")
message("Memory increase: ", round(final_mem - initial_mem, 1), " MB")

# Create summary table
message("\nCreating summary table...")
if (requireNamespace("dplyr", quietly = TRUE)) {
  summary_table <- dplyr::bind_rows(summary_rows)
} else {
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

message("Results saved to: ", results_dir)
message("\nRun analysis with:")
message("  Rscript simulations/analyze_results.R")
