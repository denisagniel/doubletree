# Simulation Logging Protocol

**Purpose:** Prevent massive log files (>10 GB) from simulation runs.

**Date:** 2026-03-12

---

## The Problem

Large-scale simulations with progress tracking and verbose output can create log files that:
- Exceed 10 GB in size
- Crash file sync services (OneDrive, Dropbox)
- Make git repositories unusable
- Fill up disk space

**Root causes:**
1. Progress messages in loops (e.g., "Rep 1/500...") × thousands of runs
2. Tree JSON output from `optimaltrees` (can be MB per tree)
3. Verbose output from `ranger`, `glmnet`, etc.
4. Parallel processing multiplying output across cores
5. Writing to log files instead of console only

---

## Rules

### Rule 1: NEVER redirect console output to log files

**Bad:**
```r
# This creates massive files
sink("simulation.log")
# ... run simulations ...
sink()
```

**Good:**
```r
# Let output go to console only (user can redirect if needed)
cat("Starting simulations...\n")
```

### Rule 2: Suppress ALL output in parallel simulation loops

**Required pattern:**
```r
run_single_sim <- function(...) {
  # Suppress EVERYTHING inside the worker function
  suppressMessages(suppressWarnings({
    invisible(capture.output({

      # All simulation code here
      fit <- dml_att(...)

    }, file = nullfile()))
  }))

  # Return results (no printing)
  return(results_df)
}

# Parallel execution
results_list <- mclapply(
  sim_grid$sim_id,
  run_single_sim,
  mc.cores = N_CORES
)
```

### Rule 3: Progress tracking to console ONLY, infrequently

**Bad:**
```r
for (i in 1:500) {
  cat(sprintf("Rep %d/%d...\n", i, 500))  # 500 lines of output
}
```

**Good:**
```r
for (i in 1:500) {
  # Progress ONLY to console, infrequently
  if (i %% 50 == 0 && interactive()) {
    cat(sprintf("[%s] Progress: %d/%d (%.0f%%)\n",
                Sys.time(), i, 500, 100*i/500))
  }
}
```

### Rule 4: Use `verbose = FALSE` for all package calls

**Always:**
```r
dml_att(..., verbose = FALSE)
ranger(..., verbose = FALSE)
cv.glmnet(..., verbose = FALSE)
```

### Rule 5: Use null device for unwanted output

**Pattern:**
```r
# Cross-platform null device
nullfile <- function() {
  if (.Platform$OS.type == "windows") {
    "NUL"
  } else {
    "/dev/null"
  }
}

# Use it
capture.output({
  # noisy code
}, file = nullfile())
```

---

## Standard Suppression Helper

Add to beginning of every simulation script:

```r
#' Complete output suppression for simulation workers
#'
#' Use this to wrap the entire body of parallel worker functions.
#' Suppresses messages, warnings, stdout, and stderr.
#'
#' @param expr Expression to evaluate silently
#' @return Result of expr (no side effects printed)
suppress_all <- function(expr) {
  nullfile <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"

  suppressMessages(suppressWarnings({
    invisible(capture.output({
      result <- expr
    }, file = nullfile))
  }))

  result
}

# Usage in worker function:
run_single_sim <- function(...) {
  suppress_all({
    # ALL simulation code here
    fit <- dml_att(...)
    # ...
  })

  # Return results
  data.frame(...)
}
```

---

## Checklist for Every Simulation Script

Before running production simulations, verify:

- [ ] No `sink()` calls that write to files
- [ ] No log files created explicitly
- [ ] All package calls have `verbose = FALSE`
- [ ] Worker functions wrap all code in `suppress_all()`
- [ ] Progress messages use `interactive()` check
- [ ] Progress messages are infrequent (e.g., every 50-100 reps)
- [ ] `capture.output()` uses `nullfile()`, not `tempfile()`
- [ ] Test with 3-5 reps first to verify no output bloat

---

## File Size Monitoring

Add to long-running simulations:

```r
# Before starting
cat("Starting simulations at:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n\n")

# Check for unexpected log files after completion
list_large_files <- function(path = ".", min_mb = 10) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  sizes <- file.size(files)
  large <- files[!is.na(sizes) & sizes > min_mb * 1024^2]

  if (length(large) > 0) {
    cat("\n⚠️  Large files detected (>", min_mb, "MB):\n")
    for (f in large) {
      cat(sprintf("  %s: %.1f MB\n", f, file.size(f) / 1024^2))
    }
  }
}

# At end of script
list_large_files(".", min_mb = 10)
```

---

## Emergency: Clean Up Existing Large Files

```bash
# Find large log files
find . -name "*.log" -size +100M -ls

# Remove them (after verifying they're not needed)
find . -name "*.log" -size +100M -delete

# Check .Rout files too
find . -name "*.Rout" -size +50M -delete
```

---

## .gitignore Additions

Ensure these patterns are in `.gitignore`:

```
*.log
*.Rout
*_log.txt
simulation_output*.txt
temp_*.txt
.Rproj.user/
.Rhistory
.RData
```

Already present in:
- `optimaltrees/.gitignore`
- `doubletree/.gitignore`
- `.gitignore` (root)

---

## Testing Protocol

**Phase 1: Micro test (1 minute)**
```r
N_REPS <- 3
# Run and verify:
# - No log files created
# - Minimal console output
# - Results correct
```

**Phase 2: Check file sizes**
```bash
# Before starting full simulation
du -sh .

# After completion
du -sh .
# Should be < 50 MB increase (results only)
```

**Phase 3: Production run**
```r
N_REPS <- 500
# Monitor first 10 reps:
# - No files growing unexpectedly
# - Console output minimal
```

---

## Updated Simulation Template

```r
# Load packages SILENTLY
suppressMessages({
  library(parallel)
  library(dplyr)
  devtools::load_all("../../../optimaltrees")
})

# Source functions SILENTLY
invisible(sapply(c(
  "../../R/dml_att.R",
  "../../R/inference.R",
  # ... other files
), source, verbose = FALSE))

# Null device helper
nullfile <- function() {
  if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
}

# Suppression helper
suppress_all <- function(expr) {
  suppressMessages(suppressWarnings({
    invisible(capture.output(expr, file = nullfile()))
  }))
}

# Worker function with COMPLETE suppression
run_single_sim <- function(sim_id, ...) {

  # Progress: console only, infrequent
  if (sim_id %% 50 == 0 && interactive()) {
    cat(sprintf("[%s] %d simulations complete\n", Sys.time(), sim_id))
  }

  # ALL simulation code wrapped
  result <- suppress_all({
    # Generate data
    d <- generate_dgp(...)

    # Fit model (NO VERBOSE OUTPUT)
    fit <- dml_att(
      X = d$X, A = d$A, Y = d$Y,
      verbose = FALSE,  # CRITICAL
      ...
    )

    # Extract results
    list(theta = fit$theta, sigma = fit$sigma, ...)
  })

  # Return clean data frame (no printing)
  data.frame(
    sim_id = sim_id,
    theta = result$theta,
    ...
  )
}

# Run simulations
results <- mclapply(
  sim_grid$sim_id,
  run_single_sim,
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Combine and save
results_df <- do.call(rbind, results)
saveRDS(results_df, "results/simulation_results.rds")

# Summary to console ONLY
cat(sprintf("\n✓ Complete: %d simulations\n", nrow(results_df)))
cat(sprintf("Convergence: %.1f%%\n", 100 * mean(results_df$converged)))
```

---

## Recovery After Bloat

If you've already created large files:

1. **Stop the simulation** immediately
2. **Check file sizes**: `find . -name "*.log" -size +10M`
3. **Delete large logs**: `rm simulation.log` (after verifying not needed)
4. **Fix the script** using patterns above
5. **Test with 3 reps** before restarting full run
6. **Resume** from last saved batch if using batch processing

---

## Next Steps

1. **Audit existing scripts** - Check all `doubletree/simulations/**/*.R` files
2. **Add `suppress_all()` helper** to all simulation scripts
3. **Remove all `sink()` calls** that write to files
4. **Test each script** with N_REPS = 3 before production runs
5. **Update `MEMORY_SAFE_SIMULATIONS.md`** to reference this protocol

---

**Last updated:** 2026-03-12 after log file bloat incident
