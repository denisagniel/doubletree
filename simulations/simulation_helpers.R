#' Simulation Helper Functions
#'
#' Reusable utilities for running large-scale simulations without log bloat.
#' Source this file at the beginning of simulation scripts.
#'
#' See: LOGGING_PROTOCOL.md for full documentation

#' Get null device path (cross-platform)
#'
#' @return Path to null device ("/dev/null" on Unix, "NUL" on Windows)
nullfile <- function() {
  if (.Platform$OS.type == "windows") {
    "NUL"
  } else {
    "/dev/null"
  }
}

#' Complete output suppression for simulation workers
#'
#' Suppresses messages, warnings, stdout, and stderr. Use this to wrap
#' the entire body of parallel worker functions to prevent log file bloat.
#'
#' @param expr Expression to evaluate silently
#' @return Result of expr (no output side effects)
#'
#' @examples
#' run_single_sim <- function(sim_id, ...) {
#'   result <- suppress_all({
#'     # All simulation code here
#'     fit <- estimate_att(X, A, Y, verbose = FALSE)
#'     list(theta = fit$theta, sigma = fit$sigma)
#'   })
#'   data.frame(sim_id = sim_id, theta = result$theta, ...)
#' }
suppress_all <- function(expr) {
  null_device <- nullfile()

  suppressMessages(suppressWarnings({
    result <- invisible(capture.output({
      force(expr)
    }, file = null_device))
  }))

  # Return the expression result, not capture.output result
  expr
}

#' Progress message (console only, infrequent)
#'
#' Prints progress messages ONLY to console (never to log files) and
#' only when running interactively. Use this instead of cat() in loops.
#'
#' @param current Current iteration number
#' @param total Total number of iterations
#' @param every Print every N iterations (default: 50)
#' @param label Optional label for the progress message
#'
#' @examples
#' for (i in 1:500) {
#'   progress_msg(i, 500, every = 50, label = "Simulations")
#'   # ... run simulation ...
#' }
progress_msg <- function(current, total, every = 50, label = "Progress") {
  if (current %% every == 0 && interactive()) {
    pct <- round(100 * current / total, 1)
    cat(sprintf("[%s] %s: %d/%d (%.1f%%)\n",
                format(Sys.time(), "%H:%M:%S"), label, current, total, pct))
  }
}

#' Check for large files in directory
#'
#' Scans directory for files exceeding size threshold. Use at end of
#' simulations to detect unexpected log bloat.
#'
#' @param path Directory to scan (default: current directory)
#' @param min_mb Minimum file size in MB to report (default: 10)
#' @param recursive Scan subdirectories (default: TRUE)
#'
#' @examples
#' # At end of simulation script:
#' check_large_files(".", min_mb = 10)
check_large_files <- function(path = ".", min_mb = 10, recursive = TRUE) {
  files <- list.files(path, recursive = recursive, full.names = TRUE)
  sizes <- file.size(files)
  large <- files[!is.na(sizes) & sizes > min_mb * 1024^2]

  if (length(large) > 0) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("⚠️  WARNING: Large files detected\n")
    cat(strrep("=", 70), "\n\n")

    cat(sprintf("Found %d files > %d MB:\n\n", length(large), min_mb))

    for (f in large) {
      size_mb <- file.size(f) / 1024^2
      cat(sprintf("  %s: %.1f MB\n", basename(f), size_mb))
    }

    cat("\nCheck LOGGING_PROTOCOL.md to prevent log bloat.\n")
    cat(strrep("=", 70), "\n\n")
  }
}

#' Safe source (suppress output)
#'
#' Source R files without printing messages.
#'
#' @param file Path to R file
#' @param ... Additional arguments passed to source()
safe_source <- function(file, ...) {
  invisible(suppressMessages(source(file, verbose = FALSE, ...)))
}

#' Memory monitoring with automatic GC
#'
#' Monitors memory usage and forces garbage collection if approaching limit.
#' Use inside simulation loops to prevent memory exhaustion.
#'
#' @param threshold_mb Memory threshold in MB for warning (default: 12000)
#' @param force_gc Force garbage collection (default: TRUE if over threshold)
#' @param verbose Print memory status (default: FALSE, console only if TRUE)
#'
#' @examples
#' for (rep in 1:500) {
#'   # ... run simulation ...
#'   if (rep %% 10 == 0) {
#'     monitor_memory(threshold_mb = 12000, verbose = interactive())
#'   }
#' }
monitor_memory <- function(threshold_mb = 12000, force_gc = NULL, verbose = FALSE) {
  mem_info <- gc(verbose = FALSE, full = FALSE)
  used_mb <- sum(mem_info[, "used"]) * 0.001  # Convert Kb to MB

  if (verbose && interactive()) {
    cat(sprintf("  [Memory: %.0f MB]\n", used_mb))
  }

  # Auto GC if over threshold
  if (used_mb > threshold_mb) {
    if (verbose && interactive()) {
      cat(sprintf("  ⚠️  Memory usage high (%.0f MB) - forcing GC\n", used_mb))
    }

    gc(verbose = FALSE, full = TRUE)

    mem_after <- sum(gc(verbose = FALSE, full = FALSE)[, "used"]) * 0.001
    if (verbose && interactive()) {
      cat(sprintf("  [Memory after GC: %.0f MB]\n", mem_after))
    }
  } else if (!is.null(force_gc) && force_gc) {
    gc(verbose = FALSE, full = FALSE)
  }

  invisible(used_mb)
}

#' Batch-safe result saving
#'
#' Saves results with atomic write (write to temp, then rename) to prevent
#' corruption if process is killed mid-write.
#'
#' @param object R object to save
#' @param file Output file path
#'
#' @examples
#' safe_save(results, "results/batch_01.rds")
safe_save <- function(object, file) {
  temp_file <- paste0(file, ".tmp")

  tryCatch({
    saveRDS(object, temp_file)
    file.rename(temp_file, file)

    if (interactive()) {
      cat(sprintf("✓ Saved: %s (%.1f MB)\n",
                  basename(file), file.size(file) / 1024^2))
    }
  }, error = function(e) {
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
    stop("Failed to save results: ", conditionMessage(e))
  })

  invisible(file)
}
