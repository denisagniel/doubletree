#!/usr/bin/env Rscript
# =============================================================================
# probe_tier2_boundary.R -- map the feasibility boundary for the Rashomon methods
# =============================================================================
# The three Rashomon-intersection methods (doubletree, dt_averaged, single_tree)
# blow up at the stress cell (n=2000, dgp=continuous): escalate_rashomon_intersection
# widens epsilon_n (c = 1..1000) until the K-fold intersection is non-empty, and the
# Rashomon set explodes (~123k trees) en route, OOMing the process.
#
# This probe finds WHERE each method crosses feasible -> infeasible, and WHY, by
# running ONE unit per (method x n) at the continuous DGP inside a memory-capped
# killable subprocess (same guard as profile_per_method.R). For each cell it reports:
#   status  : ok | killed(mem) | error/empty-intersection | slow
#   secs    : wall time
#   peak_gb : max resident memory observed (incl. C++/Rcpp GOSDT, via ps)
#   c_e/c_m0: the escalated Rashomon multiplier per nuisance (1 = theory value;
#             large = intersection empty at theory tolerance, had to widen)
#   theta / converged : sanity of the estimate when it does run
#
# A cell that is `ok` at low peak_gb is submittable as-is. A cell that is `killed`
# or errors (empty intersection) is the boundary: exclude+document it, or run it in a
# big-mem/long-wall array. Result is DATA about the method under stress, not a bug.
#
# RUN ON O2 (matches run hardware), on a COMPUTE node with headroom, under tmux:
#   srun --pty -p interactive -c 4 --mem=96G -t 0-03:00 bash
#   hostname                        # confirm compute-*  (NOT login*)
#   Rscript slurm/probe_tier2_boundary.R --mem-cap-gb 64 [--n 500,1000,2000] \
#           [--methods doubletree,dt_averaged,single_tree] [--study-dir .]
# =============================================================================

suppressPackageStartupMessages(library(optparse))

for (pkg in c("callr", "ps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' required for the memory guard; install it first.", pkg),
         call. = FALSE)
  }
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--methods", type = "character", dest = "methods",
              default = "doubletree,dt_averaged,single_tree",
              help = "Comma-separated methods to probe [default %default]"),
  make_option("--n", type = "character", dest = "ns", default = "500,1000,2000",
              help = "Comma-separated sample sizes [default %default]"),
  make_option("--dgp", type = "character", dest = "dgp", default = "continuous",
              help = "DGP regime (the stress cell) [default %default]"),
  make_option("--mem-cap-gb", type = "double", dest = "mem_cap_gb", default = 64,
              help = "Per-unit RSS ceiling in GB; set BELOW your salloc --mem [default %default]"),
  make_option("--poll-secs", type = "double", dest = "poll_secs", default = 0.5,
              help = "RSS poll interval, seconds [default %default]"),
  make_option("--study-dir", type = "character", dest = "study_dir", default = NA_character_,
              help = "Study dir [default: auto-detected from this script's path]")
)))

# Resolve study dir from the script path if not given (works from any cwd).
if (is.na(opt$study_dir)) {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
  opt$study_dir <- if (length(file_arg) == 1 && nzchar(file_arg))
    dirname(dirname(normalizePath(file_arg))) else "."
}
if (!file.exists(file.path(opt$study_dir, "config", "grid.R")))
  stop(sprintf("Cannot find config/grid.R under '%s'. Pass --study-dir.", opt$study_dir),
       call. = FALSE)

source(file.path(opt$study_dir, "config", "grid.R"))
source(file.path(opt$study_dir, "R", "dgp.R"))
source(file.path(opt$study_dir, "R", "estimators.R"))
source(file.path(opt$study_dir, "R", "run_one.R"))
suppressPackageStartupMessages({ library(doubletree); library(optimaltrees) })

methods <- strsplit(opt$methods, ",")[[1]]
ns      <- as.integer(strsplit(opt$ns, ",")[[1]])

# Run one unit in a killable callr subprocess; poll RSS (captures GOSDT's C++
# allocations, invisible to gc()); kill if it exceeds the cap. Returns timing,
# peak memory, and the run_one result row (or the error) -- never throws for a
# unit-level failure; a blow-up is DATA.
run_unit_capped <- function(unit_row, study_dir, cap_gb, poll_secs, libpath) {
  cap_bytes <- cap_gb * 1024^3
  study_dir_abs <- normalizePath(study_dir)
  proc <- callr::r_bg(
    function(unit_row, study_dir) {
      source(file.path(study_dir, "config", "grid.R"))
      source(file.path(study_dir, "R", "dgp.R"))
      source(file.path(study_dir, "R", "estimators.R"))
      source(file.path(study_dir, "R", "run_one.R"))
      suppressPackageStartupMessages({ library(doubletree); library(optimaltrees) })
      run_one(unit_row)
    },
    args = list(unit_row = unit_row, study_dir = study_dir_abs),
    libpath = libpath, supervise = TRUE
  )
  t0 <- proc.time()[["elapsed"]]; peak <- 0; killed <- FALSE
  h <- tryCatch(ps::ps_handle(proc$get_pid()), error = function(e) NULL)
  while (proc$is_alive()) {
    if (!is.null(h)) {
      rss <- tryCatch(ps::ps_memory_info(h)[["rss"]], error = function(e) NA_real_)
      if (is.finite(rss)) { peak <- max(peak, rss)
        if (rss > cap_bytes) { proc$kill_tree(); killed <- TRUE; break } }
    }
    Sys.sleep(poll_secs)
  }
  proc$wait(timeout = 5000)
  elapsed <- proc.time()[["elapsed"]] - t0
  res <- if (killed) NULL else tryCatch(proc$get_result(), error = function(e) e)
  list(elapsed = elapsed, peak_gb = peak / 1024^3, killed = killed, res = res)
}

ut <- unit_table()
cat(sprintf("Tier-2 boundary probe: dgp=%s, cap=%.0f GB, methods={%s}, n={%s}\n\n",
            opt$dgp, opt$mem_cap_gb, paste(methods, collapse=","),
            paste(ns, collapse=",")))
cat(sprintf("%-14s %5s %8s %8s %6s %6s %8s %s\n",
            "method","n","status","secs","peakGB","c_e","c_m0","theta"))

rows <- list()
for (meth in methods) for (nn in ns) {
  cell <- ut[ut$method==meth & ut$dgp==opt$dgp & ut$n==nn, , drop=FALSE]
  if (nrow(cell) == 0) next
  r <- run_unit_capped(cell[1, , drop=FALSE], opt$study_dir,
                       opt$mem_cap_gb, opt$poll_secs, .libPaths())
  status <- "ok"; c_e <- NA; c_m0 <- NA; theta <- NA
  if (r$killed) {
    status <- "KILLED(mem)"
  } else if (inherits(r$res, "error") || inherits(r$res, "condition")) {
    status <- "ERROR"
  } else if (is.data.frame(r$res)) {
    row <- r$res
    theta <- row$estimate
    if ("rashomon_c_e" %in% names(row)) c_e <- row$rashomon_c_e
    if ("rashomon_c_m0" %in% names(row)) c_m0 <- row$rashomon_c_m0
    # An all-NA estimate with a converged==FALSE => empty intersection, not OOM.
    if (is.na(theta)) status <- "empty-intersect"
  }
  cat(sprintf("%-14s %5d %8s %8.1f %6.2f %6s %8s %s\n",
              meth, nn, status, r$elapsed, r$peak_gb,
              ifelse(is.na(c_e),"-",format(c_e)), ifelse(is.na(c_m0),"-",format(c_m0)),
              ifelse(is.na(theta),"-",sprintf("%.4f",theta))))
  rows[[length(rows)+1]] <- data.frame(method=meth, n=nn, status=status,
    secs=r$elapsed, peak_gb=r$peak_gb, c_e=c_e, c_m0=c_m0, theta=theta)
}

cat("\nInterpretation:\n")
cat("  ok, low peakGB           -> submittable as-is (Tier 2, sized to peak)\n")
cat("  KILLED(mem)              -> Rashomon set OOMs; big-mem array or exclude+document\n")
cat("  empty-intersect / ERROR  -> intersection empty even at c_max; method breaks here\n")
cat("  large c_e/c_m0           -> theory tolerance was empty; escalation widened it\n")
