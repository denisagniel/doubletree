# =============================================================================
# verify_phaseB_coverage.R -- post-fix regression: does the PACKAGED honest CI cover?
# =============================================================================
# THROWAWAY verification (Phase B, 2026-07-15). Unlike diagnose_complex.R (which
# hand-built the honest CI from a fold-specific twin), this checks the CI the PACKAGE
# now reports directly -- estimate_att(use_rashomon=TRUE)$ci_95, msplit$ci_95,
# single_tree$ci_95 -- so it verifies the shipped code path, not a reconstruction.
#
# For each (n, rep) on the complex DGP it records coverage + width + power (excl 0) for:
#   doubletree  = estimate_att(use_rashomon=TRUE)   (now honest)
#   crossfit    = estimate_att(use_rashomon=FALSE)  (valid Wald baseline)
#   single_tree = estimate_att_single_tree(inference="single")  (honest)
#   msplit      = estimate_att_msplit               (now honest)
#
# MEMORY-SAFE: strictly serial, each unit in a killable callr subprocess (RSS cap +
# wall timeout), incremental save. Mirrors diagnose_complex.R's guard.
#   Rscript dev-scripts/verify_phaseB_coverage.R [n_reps] [--mem-cap-gb 6] [--max-secs 900]
# =============================================================================

suppressPackageStartupMessages(library(doubletree))

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
  stop("Could not locate study root.", call. = FALSE)
}
study_root <- .find_study_root()
source(file.path(study_root, "R", "dgp.R"))

raw_args <- commandArgs(trailingOnly = TRUE)
.get_flag <- function(name, default) {
  hit <- which(raw_args == paste0("--", name))
  if (length(hit) == 1 && hit < length(raw_args)) as.numeric(raw_args[hit + 1L]) else default
}
positional <- raw_args[!grepl("^--", raw_args) &
                         !(seq_along(raw_args) %in% (which(grepl("^--", raw_args)) + 1L))]
N_REPS   <- if (length(positional) >= 1) as.integer(positional[1]) else 60L
MEM_CAP  <- .get_flag("mem-cap-gb", 6)
MAX_SECS <- .get_flag("max-secs", 900)
POLL     <- .get_flag("poll-secs", 0.5)
N_GRID <- c(500L, 1000L, 2000L)
TRUTH  <- 0.15; Z <- qnorm(0.975)
OUT_RDS <- file.path(study_root, "dev-scripts", "verify_phaseB_results.rds")

.unit_body <- function(n, rep_id, study_root, TRUTH) {
  suppressPackageStartupMessages(library(doubletree))
  source(file.path(study_root, "R", "dgp.R"))
  set.seed(20260716L + as.integer(n) * 100003L + rep_id)
  d <- .dgp_complex(n); X <- d$X; A <- d$A; Y <- d$Y
  ci_row <- function(ci) if (is.null(ci) || any(!is.finite(ci))) c(NA, NA) else ci
  grab <- function(expr) tryCatch(expr, error = function(e) NULL)

  dt <- grab(estimate_att(X, A, Y, K = 5, use_rashomon = TRUE, max_depth = 4L, verbose = FALSE))
  cf <- grab(estimate_att(X, A, Y, K = 5, use_rashomon = FALSE, max_depth = 4L, verbose = FALSE))
  st <- grab(estimate_att_single_tree(X, A, Y, K = 5, inference = "single", verbose = FALSE))
  ms <- grab(estimate_att_msplit(X, A, Y, M = 5, K = 5, verbose = FALSE))

  one <- function(fit) {
    if (is.null(fit)) return(c(theta = NA, lo = NA, hi = NA))
    ci <- ci_row(fit$ci_95); c(theta = fit$theta, lo = ci[1], hi = ci[2])
  }
  d_dt <- one(dt); d_cf <- one(cf); d_st <- one(st); d_ms <- one(ms)
  data.frame(n = n, rep_id = rep_id,
    dt_theta = d_dt["theta"], dt_lo = d_dt["lo"], dt_hi = d_dt["hi"],
    cf_theta = d_cf["theta"], cf_lo = d_cf["lo"], cf_hi = d_cf["hi"],
    st_theta = d_st["theta"], st_lo = d_st["lo"], st_hi = d_st["hi"],
    ms_theta = d_ms["theta"], ms_lo = d_ms["lo"], ms_hi = d_ms["hi"],
    row.names = NULL, stringsAsFactors = FALSE)
}
.na_row <- function(n, rep_id) {
  z <- .unit_body; d <- data.frame(n = n, rep_id = rep_id)
  cols <- c("dt","cf","st","ms")
  for (p in cols) for (s in c("theta","lo","hi")) d[[paste0(p,"_",s)]] <- NA_real_
  d
}

run_unit_capped <- function(n, rep_id, cap_gb, max_secs, poll_secs, libpath) {
  cap_bytes <- cap_gb * 1024^3
  proc <- callr::r_bg(.unit_body,
    args = list(n = n, rep_id = rep_id, study_root = study_root, TRUTH = TRUTH),
    libpath = libpath, supervise = TRUE)
  t0 <- proc.time()[["elapsed"]]; peak <- 0; killed <- FALSE; reason <- "ok"
  h <- tryCatch(ps::ps_handle(proc$get_pid()), error = function(e) NULL)
  while (proc$is_alive()) {
    if (!is.null(h)) {
      rss <- tryCatch(ps::ps_memory_info(h)[["rss"]], error = function(e) NA_real_)
      if (is.finite(rss)) { peak <- max(peak, rss)
        if (rss > cap_bytes) { proc$kill_tree(); killed <- TRUE; reason <- "mem cap"; break } }
    }
    if (proc.time()[["elapsed"]] - t0 > max_secs) {
      proc$kill_tree(); killed <- TRUE; reason <- "wall timeout"; break }
    Sys.sleep(poll_secs)
  }
  proc$wait(timeout = 5000)
  row <- if (killed) .na_row(n, rep_id) else {
    res <- tryCatch(proc$get_result(), error = function(e) e)
    if (inherits(res, "condition")) { reason <- conditionMessage(res); .na_row(n, rep_id) } else res
  }
  list(row = row, elapsed = proc.time()[["elapsed"]] - t0, peak_gb = peak / 1024^3,
       killed = killed, reason = reason)
}

grid <- expand.grid(rep_id = seq_len(N_REPS), n = N_GRID)
grid <- grid[order(grid$n, grid$rep_id), , drop = FALSE]
cat(sprintf("Phase B coverage check: %d units, serial, cap %.1fGB, wall %.0fs.\n",
            nrow(grid), MEM_CAP, MAX_SECS))
rows <- vector("list", nrow(grid)); n_killed <- 0L; t_start <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(grid))) {
  r <- run_unit_capped(grid$n[i], grid$rep_id[i], MEM_CAP, MAX_SECS, POLL, .libPaths())
  rows[[i]] <- r$row
  if (r$killed) n_killed <- n_killed + 1L
  saveRDS(do.call(rbind, rows[seq_len(i)]), OUT_RDS)
  if (i %% 10L == 0L || i == nrow(grid) || r$killed) {
    avg <- (proc.time()[["elapsed"]] - t_start) / i
    cat(sprintf("  [%3d/%3d] n=%4d rep=%3d  %5.1fs  peak %.2fGB  %-12s | ETA %.0fmin\n",
                i, nrow(grid), grid$n[i], grid$rep_id[i], r$elapsed, r$peak_gb, r$reason,
                avg * (nrow(grid) - i) / 60))
  }
}
res <- do.call(rbind, rows)

cov <- function(lo, hi) mean(TRUTH >= lo & TRUTH <= hi, na.rm = TRUE)
pow <- function(lo, hi) mean(lo > 0 | hi < 0, na.rm = TRUE)
wid <- function(lo, hi) mean(hi - lo, na.rm = TRUE)
summ <- function(df) data.frame(n = df$n[1], reps = sum(is.finite(df$dt_theta)),
  cov_dt = round(cov(df$dt_lo, df$dt_hi), 3), cov_cf = round(cov(df$cf_lo, df$cf_hi), 3),
  cov_st = round(cov(df$st_lo, df$st_hi), 3), cov_ms = round(cov(df$ms_lo, df$ms_hi), 3),
  w_dt = round(wid(df$dt_lo, df$dt_hi), 3), pow_dt = round(pow(df$dt_lo, df$dt_hi), 2),
  pow_cf = round(pow(df$cf_lo, df$cf_hi), 2))
by_n <- do.call(rbind, lapply(split(res, res$n), summ)); rownames(by_n) <- NULL
cat(sprintf("\n===== PHASE B PACKAGED-CI COVERAGE (complex DGP) =====\nkilled=%d\n\n", n_killed))
print(by_n, row.names = FALSE)
cat("\nEXPECT: cov_dt/cov_st/cov_ms >= 0.90 at all n (honest); cov_cf ~0.95 (Wald baseline).\n")
saveRDS(res, OUT_RDS)
cat(sprintf("Saved: %s\n", OUT_RDS))
