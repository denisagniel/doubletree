# run_all_batches.R
# Execute all batch scripts sequentially
#
# Each batch runs in a fresh R process, preventing memory accumulation

batch_dir <- "simulations/batches"

# Get all batch files in order
batch_files <- list.files(batch_dir, pattern = "^batch_\\d+.*\\.R$", full.names = TRUE)
batch_files <- sort(batch_files)

if (length(batch_files) == 0) {
  stop("No batch files found. Run create_batches.R first.")
}

message("=== Running All Simulation Batches ===\n")
message("Total batches: ", length(batch_files))
message("Each batch runs in fresh R process (prevents memory accumulation)\n")

total_start <- Sys.time()

for (i in seq_along(batch_files)) {
  batch_file <- batch_files[i]
  batch_name <- basename(batch_file)

  message("\n===============================================")
  message("STARTING BATCH ", i, "/", length(batch_files), ": ", batch_name)
  message("===============================================\n")

  batch_start <- Sys.time()

  # Run batch in fresh R process
  cmd <- paste("Rscript", shQuote(batch_file))
  exit_code <- system(cmd)

  batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))

  if (exit_code != 0) {
    warning("Batch ", i, " (", batch_name, ") failed with exit code ", exit_code)
    message("\nContinuing to next batch...\n")
  } else {
    message("\nBatch ", i, " completed successfully in ", round(batch_elapsed, 2), " minutes")
  }

  # Brief pause to let system settle
  Sys.sleep(2)
}

total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))

message("\n\n=== ALL BATCHES COMPLETE ===")
message("Total time: ", round(total_elapsed, 2), " minutes (", round(total_elapsed/60, 2), " hours)")
message("Successful batches: ", sum(file.exists(list.files("simulations/results_extended",
                                                           pattern = "^result_.*\\.rds$",
                                                           full.names = TRUE))))

message("\nCombine results with:")
message("  Rscript simulations/combine_batch_results.R")
