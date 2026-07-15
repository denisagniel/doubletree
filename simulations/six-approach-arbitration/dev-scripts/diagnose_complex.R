# =============================================================================
# diagnose_complex.R -- Phase A diagnostic for complex-DGP undercoverage
# =============================================================================
# THROWAWAY diagnostic (not part of the study runner). Investigates why the
# shared-Rashomon-tree ATT estimators undercover on the complex DGP, with
# coverage WORSENING as n grows (see plan abundant-churning-scone.md).
#
# Per rep on the complex DGP it computes THREE estimators of the same ATT:
#   oracle  : plug the TRUE e(X), m0(X) into the EIF solve      -> validity gate
#   shared  : estimate_att(use_rashomon=TRUE)  (shared structure, the suspect)
#   foldspec: estimate_att(use_rashomon=FALSE) (fully fold-specific, the honest twin)
#
# Then it evaluates, BY n:
#   (1) coverage of each estimator's Wald CI          (oracle & foldspec ~0.95 expected)
#   (2) shared-tree structural bias E[theta_shared-0.15] vs se_shared (bias/SE growth)
#   (3) whether an HONEST CI for the shared estimate, using the FOLD-SPECIFIC estimator
#       as the twin (delta = theta_shared - theta_foldspec), restores >=0.90 coverage.
#       This directly tests the proposed inference fix (user idea 2026-07-15).
#
# -----------------------------------------------------------------------------
# MEMORY SAFETY (2026-07-15 rewrite): the previous version fanned out
# `parallel::mclapply` over (detectCores()-1) workers, each running a full
# Rashomon `estimate_att` up to n=2000. Nine concurrent GOSDT fits exhausted
# 16 GB and HARD-CRASHED the host (3rd such OOM: 07-01, 07-07, 07-15). This
# version runs STRICTLY SERIALLY, each unit in a KILLABLE callr subprocess whose
# resident memory is polled via `ps` (captures GOSDT's C++ allocations, which
# gc() cannot see) and killed if it exceeds --mem-cap-gb. A per-unit wall clock
# (--max-secs) kills pathological tails. Completed rows are SAVED INCREMENTALLY
# after every unit, so an interruption never loses finished work.
#
# Usage:
#   Rscript dev-scripts/diagnose_complex.R [n_reps] \
#       [--mem-cap-gb 6] [--max-secs 900] [--poll-secs 0.5]
# Run time: serial, so long. Prototype with a small n_reps (e.g. 5) first;
# a per-unit ETA is printed so you can extrapolate before committing to 100.
# =============================================================================

suppressPackageStartupMessages({
  library(doubletree)   # parent needs honest_ci; children load their own copy
})

# ---- locate study root & source DGP (reuse smoke_run's resolver pattern) -----
.find_study_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) == 1) {
    here <- dirname(normalizePath(sub("^--file=", "", fa)))
    cand <- normalizePath(file.path(here, ".."), mustWork = FALSE)
    if (dir.exists(file.path(cand, "R"))) return(cand)
  }
  for (cand in c("simulations/six-approach-arbitration",
                 "doubletree/simulations/six-approach-arbitration", ".")) {
    if (dir.exists(file.path(cand, "R"))) return(normalizePath(cand))
  }
  stop("Could not locate the six-approach-arbitration study root.", call. = FALSE)
}
study_root <- .find_study_root()
source(file.path(study_root, "R", "dgp.R"))   # .dgp_complex, expit, true_value

# ---- config -----------------------------------------------------------------
# Minimal flag parser (no optparse dep): positional n_reps + --key value flags.
raw_args <- commandArgs(trailingOnly = TRUE)
.get_flag <- function(name, default) {
  hit <- which(raw_args == paste0("--", name))
  if (length(hit) == 1 && hit < length(raw_args)) as.numeric(raw_args[hit + 1L]) else default
}
positional <- raw_args[!grepl("^--", raw_args) &
                         !(seq_along(raw_args) %in%
                             (which(grepl("^--", raw_args)) + 1L))]
N_REPS   <- if (length(positional) >= 1) as.integer(positional[1]) else 100L
MEM_CAP  <- .get_flag("mem-cap-gb", 6)     # per-unit RSS ceiling (GB); < host RAM
MAX_SECS <- .get_flag("max-secs", 900)     # per-unit wall-clock kill (seconds)
POLL     <- .get_flag("poll-secs", 0.5)    # RSS poll interval (seconds)

N_GRID <- c(500L, 1000L, 2000L)
K      <- 5L
TRUTH  <- 0.15
Z      <- qnorm(0.975)
OUT_RDS <- file.path(study_root, "dev-scripts", "diagnose_complex_results.rds")

covered <- function(theta, se, truth = TRUTH) {
  as.integer(truth >= theta - Z * se & truth <= theta + Z * se)
}

# One (n, rep) work unit, computed INSIDE a fresh callr subprocess so its memory
# can be capped and reclaimed. Self-contained: sources dgp.R and defines its own
# true-nuisance helper (child does not inherit the parent's environment). Returns
# a one-row data.frame of RAW estimates; coverage/honest-CI are computed in the
# parent. Seed is deterministic per (n, rep) for reproducibility across runs.
.unit_body <- function(n, rep_id, study_root, K, TRUTH) {
  suppressPackageStartupMessages(library(doubletree))
  source(file.path(study_root, "R", "dgp.R"))   # .dgp_complex, expit

  # True nuisances for the complex DGP (mirrors .dgp_complex in dgp.R exactly).
  true_eta <- function(X) {
    e <- expit(-0.9 + 0.6 * (X$x1 + X$x2 + X$x3) +
                 0.8 * X$x1 * X$x2 + 0.6 * X$x2 * X$x3)
    m0 <- 0.05 + 0.15 * (X$x3 + X$x4 + X$x5) +
      0.2 * X$x3 * X$x4 + 0.15 * X$x4 * X$x5
    list(e = e, m0 = m0)
  }

  set.seed(20260715L + as.integer(n) * 100003L + rep_id)
  d <- .dgp_complex(n)
  X <- d$X; A <- d$A; Y <- d$Y

  # (0) oracle: true nuisances into the shared EIF solver.
  # eif_att_solve is internal (@keywords internal, not exported) -> use ::: .
  tru <- true_eta(X)
  orc <- tryCatch(doubletree:::eif_att_solve(Y, A, tru$e, tru$m0, n),
                  error = function(e) NULL)
  # (shared) Rashomon shared-structure estimator (the suspect)
  sh <- tryCatch(
    estimate_att(X, A, Y, K = K, use_rashomon = TRUE,
                 rashomon_bound_multiplier = NULL, auto_tune_intersecting = FALSE,
                 max_depth = 4L, verbose = FALSE),
    error = function(e) NULL)
  # (foldspec) fully fold-specific estimator (the honest twin)
  fs <- tryCatch(
    estimate_att(X, A, Y, K = K, use_rashomon = FALSE, max_depth = 4L, verbose = FALSE),
    error = function(e) NULL)

  if (is.null(sh) || is.null(fs) || is.null(orc)) {
    return(data.frame(n = n, rep_id = rep_id, ok = FALSE,
      theta_oracle = NA_real_, se_oracle = NA_real_,
      theta_shared = NA_real_, se_shared = NA_real_, intersect = NA_integer_,
      theta_fs = NA_real_, se_fs = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(n = n, rep_id = rep_id, ok = TRUE,
    theta_oracle = orc$theta, se_oracle = orc$sigma,
    theta_shared = sh$theta, se_shared = sh$sigma,
    intersect = as.integer(isTRUE(sh$converged)),
    theta_fs = fs$theta, se_fs = fs$sigma, stringsAsFactors = FALSE)
}

# NA row for a killed/errored/timed-out unit (keeps the incremental table shape).
.na_row <- function(n, rep_id) {
  data.frame(n = n, rep_id = rep_id, ok = FALSE,
    theta_oracle = NA_real_, se_oracle = NA_real_,
    theta_shared = NA_real_, se_shared = NA_real_, intersect = NA_integer_,
    theta_fs = NA_real_, se_fs = NA_real_, stringsAsFactors = FALSE)
}

# Run one unit in a killable subprocess; poll RSS (captures GOSDT's C++ memory,
# invisible to gc()); kill on RSS cap OR wall-clock timeout. Never throws: a
# blow-up is DATA, recorded as an ok=FALSE row.
run_unit_capped <- function(n, rep_id, cap_gb, max_secs, poll_secs, libpath) {
  cap_bytes <- cap_gb * 1024^3
  proc <- callr::r_bg(
    .unit_body,
    args = list(n = n, rep_id = rep_id, study_root = study_root, K = K, TRUTH = TRUTH),
    libpath = libpath, supervise = TRUE)

  t0 <- proc.time()[["elapsed"]]; peak <- 0; killed <- FALSE; reason <- "ok"
  h <- tryCatch(ps::ps_handle(proc$get_pid()), error = function(e) NULL)
  while (proc$is_alive()) {
    if (!is.null(h)) {
      rss <- tryCatch(ps::ps_memory_info(h)[["rss"]], error = function(e) NA_real_)
      if (is.finite(rss)) {
        peak <- max(peak, rss)
        if (rss > cap_bytes) { proc$kill_tree(); killed <- TRUE; reason <- "mem cap"; break }
      }
    }
    if (proc.time()[["elapsed"]] - t0 > max_secs) {
      proc$kill_tree(); killed <- TRUE; reason <- "wall timeout"; break
    }
    Sys.sleep(poll_secs)
  }
  proc$wait(timeout = 5000)
  elapsed <- proc.time()[["elapsed"]] - t0

  if (killed) {
    row <- .na_row(n, rep_id)
  } else {
    res <- tryCatch(proc$get_result(), error = function(e) e)
    if (inherits(res, "condition")) { row <- .na_row(n, rep_id); reason <- conditionMessage(res) }
    else row <- res
  }
  list(row = row, elapsed = elapsed, peak_gb = peak / 1024^3,
       killed = killed, reason = reason)
}

# ---- serial run with incremental save ----------------------------------------
grid <- expand.grid(rep_id = seq_len(N_REPS), n = N_GRID)
grid <- grid[order(grid$n, grid$rep_id), , drop = FALSE]   # n-blocks, ascending
cat(sprintf(paste0("Serial run: %d units (%d reps x %d n). Per-unit cap %.1f GB, ",
                   "wall %.0fs, poll %.1fs.\n"),
            nrow(grid), N_REPS, length(N_GRID), MEM_CAP, MAX_SECS, POLL))

rows <- vector("list", nrow(grid))
n_killed <- 0L; t_start <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(grid))) {
  r <- run_unit_capped(grid$n[i], grid$rep_id[i], MEM_CAP, MAX_SECS, POLL, .libPaths())
  rows[[i]] <- r$row
  if (r$killed) n_killed <- n_killed + 1L
  # Incremental save: flush accumulated rows after EVERY unit (crash-proof).
  saveRDS(do.call(rbind, rows[seq_len(i)]), OUT_RDS)
  if (i %% 10L == 0L || i == nrow(grid) || r$killed) {
    avg <- (proc.time()[["elapsed"]] - t_start) / i
    cat(sprintf("  [%3d/%3d] n=%4d rep=%3d  %5.1fs  peak %.2fGB  %-13s | ETA %.0fmin\n",
                i, nrow(grid), grid$n[i], grid$rep_id[i], r$elapsed, r$peak_gb,
                r$reason, avg * (nrow(grid) - i) / 60))
  }
}
res <- do.call(rbind, rows)
res <- res[isTRUE(res$ok) | res$ok %in% TRUE, , drop = FALSE]
res <- res[!is.na(res$ok) & res$ok, , drop = FALSE]

if (nrow(res) == 0L) {
  stop(sprintf("All %d units failed (killed=%d). Nothing to summarize. See %s.",
               nrow(grid), n_killed, OUT_RDS), call. = FALSE)
}

# ---- honest CI using the FOLD-SPECIFIC estimator as the twin -----------------
# delta = theta_shared - theta_fs (fold-specific is the valid, orthogonal twin).
# se_delta unknown per-rep here; use se_delta=0 (raw |delta| bound) as the
# conservative-but-not-noise-inflated baseline. honest_ci widens by cv(|delta|/se)*se.
res$delta_fs <- res$theta_shared - res$theta_fs
res$cov_shared_wald <- covered(res$theta_shared, res$se_shared)
res$cov_fs_wald     <- covered(res$theta_fs,     res$se_fs)
res$cov_oracle      <- covered(res$theta_oracle, res$se_oracle)
res$cov_shared_honest <- vapply(seq_len(nrow(res)), function(i) {
  hc <- honest_ci(res$theta_shared[i], res$se_fs[i], res$delta_fs[i],
                  se_delta = 0, level = 0.95)
  as.integer(TRUTH >= hc$ci[1] & TRUTH <= hc$ci[2])
}, integer(1))

# ---- summary by n ------------------------------------------------------------
summ <- function(df) {
  data.frame(
    n = df$n[1], reps = nrow(df),
    # validity gates
    cov_oracle = round(mean(df$cov_oracle), 3),
    cov_foldspec = round(mean(df$cov_fs_wald), 3),
    # the bug
    cov_shared = round(mean(df$cov_shared_wald), 3),
    # bias/SE growth (shared)
    bias_shared = round(mean(df$theta_shared) - TRUTH, 4),
    se_shared = round(mean(df$se_shared), 4),
    bias_over_se = round((mean(df$theta_shared) - TRUTH) / mean(df$se_shared), 2),
    se_ratio_shared = round(mean(df$se_shared) / sd(df$theta_shared), 2),
    # proposed inference fix
    mean_abs_delta = round(mean(abs(df$delta_fs)), 4),
    cov_shared_honest = round(mean(df$cov_shared_honest), 3),
    intersect_rate = round(mean(df$intersect), 2),
    stringsAsFactors = FALSE)
}
by_n <- do.call(rbind, lapply(split(res, res$n), summ))
rownames(by_n) <- NULL

cat("\n===================== PHASE A DIAGNOSTIC (complex DGP) =====================\n")
cat(sprintf("reps/n requested = %d, usable = %d, killed = %d, K = %d\n\n",
            N_REPS, nrow(res), n_killed, K))
print(by_n, row.names = FALSE)

cat("\nREAD:\n")
cat(" * cov_oracle ~0.95 => EIF/SE machinery is correct (gate).\n")
cat(" * cov_foldspec ~0.95 & improving => fully fold-specific estimator is the valid twin.\n")
cat(" * cov_shared falling as n grows => reproduces the shared-tree bug.\n")
cat(" * bias_over_se rising with n => structural bias vanishes slower than SE.\n")
cat(" * cov_shared_honest >=0.90 => the fold-specific-twin honest CI FIXES it (inference lever).\n")

saveRDS(res, OUT_RDS)   # final save: enriched with coverage/honest columns
cat(sprintf("\nSaved per-rep results (incrementally + final): %s\n", OUT_RDS))
