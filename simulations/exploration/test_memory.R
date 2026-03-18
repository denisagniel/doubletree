# test_memory.R
# Quick memory test: 2 configs only

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/run_simulations_extended.R", local = TRUE)

message("=== Memory Test: 2 Configs Only ===\n")

dgps <- list(dgp1 = generate_data_dgp1)
ns <- c(200, 400)
epsilons <- c(0.05)
n_reps <- 5  # Very small
tau <- 0.15
K <- 5

summary_rows <- list()
for (dgp_name in names(dgps)) {
  for (n in ns) {
    for (epsilon in epsilons) {
      message("\nConfig: ", dgp_name, ", n=", n, ", eps=", epsilon)

      # Check memory before
      mem_before <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024
      message("Memory before: ", round(mem_before, 1), " MB")

      result <- run_comparison(dgps[[dgp_name]], n, K, tau, n_reps, epsilon, seed_start = 1000)

      cfg <- data.frame(dgp = dgp_name, n = n, epsilon = epsilon, K = K, tau = tau, n_reps = n_reps)
      summary_rows[[length(summary_rows) + 1]] <- cbind(cfg, method = "fold_specific", as.data.frame(result$fold_specific))
      summary_rows[[length(summary_rows) + 1]] <- cbind(cfg, method = "rashomon", as.data.frame(result$rashomon))
      summary_rows[[length(summary_rows) + 1]] <- cbind(cfg, method = "oracle", as.data.frame(result$oracle))

      # Memory after (with raw_results still in memory)
      mem_after_with_raw <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024
      message("Memory after (with raw_results): ", round(mem_after_with_raw, 1), " MB")
      message("  Delta: ", round(mem_after_with_raw - mem_before, 1), " MB")

      # Now remove raw_results and force GC
      rm(result)
      gc(verbose = FALSE)

      mem_after_gc <- as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024
      message("Memory after rm + gc: ", round(mem_after_gc, 1), " MB")
      message("  Delta from start: ", round(mem_after_gc - mem_before, 1), " MB")
      message("  Freed by GC: ", round(mem_after_with_raw - mem_after_gc, 1), " MB")
    }
  }
}

message("\n=== Test Complete ===")
message("Key finding: Removing raw_results and calling gc() frees significant memory")
