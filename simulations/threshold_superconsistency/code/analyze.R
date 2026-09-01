# ============================================================
# analyze.R
# Study: threshold_superconsistency  (doubletree)
#
# Reads the per-cell checkpoints, computes mean |t_hat - t*| per sample size, and
# regresses log(mean error) on log(n) to read off the convergence rate.
#
# WHAT THE SLOPE MEANS
#   ~ -1    : the jump-driven superconsistency mechanism is biting. The threshold
#             is being located at the resolution of the data spacing (1/n), which
#             is what the paper's O_p(1/n) split-location claim needs and what
#             makes a pre-specified grid unnecessary.
#   ~ -1/3  : it is NOT biting. -1/3 is the mesh of the Stage A quantile grid
#             (ceiling(n^(1/3)) bins), so a slope there says Stage B failed to
#             improve on the grid and the error is grid-limited, i.e. exactly the
#             regime the rewrite is supposed to escape.
#   ~ -1/2  : neither -- would indicate the estimate is behaving like a sample
#             quantile / smooth-argmin rather than a jump change-point.
#
# Stage A is reported alongside Stage B as an internal control: Stage A's slope
# SHOULD sit near -1/3 by construction, so if it does not, the diagnostic itself
# is suspect and the Stage B number should not be trusted either.
#
# Run directly:  Rscript code/analyze.R
# ============================================================

suppressPackageStartupMessages({ library(cli) })

STUDY <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/threshold_superconsistency"

#' Load every cell checkpoint in results/
load_cells <- function(dir = file.path(STUDY, "results")) {
  files <- sort(list.files(dir, pattern = "^cell_n[0-9]+_r[0-9]+\\.rds$", full.names = TRUE))
  if (length(files) == 0L) cli::cli_abort("No cell checkpoints found in {dir}.")
  payloads <- lapply(files, readRDS)
  list(
    files   = files,
    results = do.call(rbind, lapply(payloads, `[[`, "results")),
    memory  = do.call(rbind, lapply(payloads, `[[`, "memory")),
    meta    = lapply(payloads, `[[`, "meta")
  )
}

#' Per-cell mean absolute threshold error
#'
#' Reported three ways on purpose:
#'   * mean over all topology-recovered reps -- the headline the design asked for;
#'   * median -- because the mean of 8 reps is destroyed by a single outlier, and
#'     the outliers here have an identified cause (see next bullet);
#'   * mean over reps whose Stage A grid cut was ADJACENT to the truth, i.e. whose
#'     bracket actually contained t*. When it does not, Stage B is clipped at a
#'     bracket edge and its error measures the handoff, not the jump mechanism.
#'     Separating those is the difference between "the mechanism does not work" and
#'     "the bracket was too narrow", which are opposite conclusions.
cell_summary <- function(res) {
  ns <- sort(unique(res$n))
  rows <- lapply(ns, function(nn) {
    d  <- res[res$n == nn, , drop = FALSE]
    dk <- d[d$topology_ok, , drop = FALSE]                       # topology recovered
    db <- dk[dk$bracket1_contains_truth, , drop = FALSE]         # ...and bracket valid
    data.frame(
      n = nn,
      reps = nrow(d),
      reps_topology_ok = nrow(dk),
      reps_bracket1_ok = nrow(db),
      n_bins = d$n_bins[[1L]],
      grid_mesh = 1 / d$n_bins[[1L]],
      mean_err1_stage_a = mean(dk$err1_stage_a),
      mean_err1_stage_b = mean(dk$err1_stage_b),
      median_err1_stage_b = stats::median(dk$err1_stage_b),
      mean_err1_stage_b_bracket_ok = if (nrow(db)) mean(db$err1_stage_b) else NA_real_,
      mean_err2_stage_a = mean(dk$err2_stage_a),
      mean_err2_stage_b = mean(dk$err2_stage_b),
      median_err2_stage_b = stats::median(dk$err2_stage_b),
      mean_leaves = mean(d$n_leaves),
      mean_secs = mean(d$secs_stage_a + d$secs_stage_b)
    )
  })
  do.call(rbind, rows)
}

#' log-log slope of a mean error against n
#'
#' Three sample sizes give one degree of freedom for the residual, so the
#' standard error is indicative only. This is a go/no-go read, not a
#' publication-quality rate curve.
loglog_slope <- function(n, err) {
  ok <- is.finite(err) & err > 0
  if (sum(ok) < 2L) return(list(slope = NA_real_, se = NA_real_, r2 = NA_real_, n_points = sum(ok)))
  m <- stats::lm(log(err[ok]) ~ log(n[ok]))
  s <- summary(m)
  list(
    slope = unname(stats::coef(m)[[2L]]),
    se = if (nrow(s$coefficients) >= 2L && ncol(s$coefficients) >= 2L) s$coefficients[2L, 2L] else NA_real_,
    r2 = s$r.squared,
    n_points = sum(ok),
    fit = m
  )
}

#' Verdict from the Stage B slope
verdict_for <- function(slope) {
  if (!is.finite(slope)) return("INCONCLUSIVE (slope not estimable)")
  d1  <- abs(slope - (-1))
  d13 <- abs(slope - (-1 / 3))
  if (d1 <= 0.25) {
    "SUPERCONSISTENCY CONFIRMED (slope close to -1)"
  } else if (d13 <= 0.20) {
    "NOT CONFIRMED -- error is grid-limited (slope close to -1/3)"
  } else if (slope > -0.25) {
    "NOT CONFIRMED -- error is barely shrinking in n"
  } else {
    "NOT CONFIRMED -- slope is not close to -1"
  }
}

analyze <- function() {
  cells <- load_cells()
  res <- cells$results
  summ <- cell_summary(res)

  cli::cli_h1("Per-cell threshold error (split 1: root, X1, t1* = {res$t1_star[[1L]]})")
  print(summ[, c("n", "reps", "reps_topology_ok", "reps_bracket1_ok", "n_bins",
                 "grid_mesh", "mean_err1_stage_a", "mean_err1_stage_b",
                 "median_err1_stage_b", "mean_err1_stage_b_bracket_ok",
                 "mean_leaves", "mean_secs")],
        row.names = FALSE, digits = 4)

  cli::cli_h1("Bracket diagnostics (split 1)")
  bad <- res[!res$bracket1_contains_truth, , drop = FALSE]
  cli::cli_inform("reps where the +/-1-grid-point bracket did NOT contain t1*: {nrow(bad)} of {nrow(res)}")
  if (nrow(bad)) {
    print(bad[, c("n", "rep", "n_leaves", "t1_grid", "t1_lo", "t1_hi",
                  "t1_refined", "err1_stage_a", "err1_stage_b")],
          row.names = FALSE, digits = 5)
    cli::cli_inform(c(
      "!" = "In these reps Stage B is CLIPPED: the scan cannot reach t1* because it lies outside the bracket.",
      "i" = "max |refined_cut - nearest bracket edge| = {signif(max(pmin(abs(bad$t1_refined - bad$t1_lo), abs(bad$t1_refined - bad$t1_hi))), 3)} -- i.e. the estimate sits ON the edge.",
      "i" = "This is a limitation of the Stage A -> Stage B handoff (bracket width), NOT evidence about the jump mechanism."
    ))
  }

  # Headline slope, plus the two robustness reads that isolate the bracket artifact.
  sb    <- loglog_slope(summ$n, summ$mean_err1_stage_b)
  sb_md <- loglog_slope(summ$n, summ$median_err1_stage_b)
  sb_bo <- loglog_slope(summ$n, summ$mean_err1_stage_b_bracket_ok)
  sa    <- loglog_slope(summ$n, summ$mean_err1_stage_a)
  s2    <- loglog_slope(summ$n, summ$mean_err2_stage_b)
  s2_md <- loglog_slope(summ$n, summ$median_err2_stage_b)

  cli::cli_h1("log-log slopes")
  cli::cli_inform(c(
    "*" = "split 1 STAGE B, mean, all topology-ok reps : slope = {signif(sb$slope, 4)} (se {signif(sb$se, 3)}, R2 {signif(sb$r2, 3)})",
    "*" = "split 1 STAGE B, median                     : slope = {signif(sb_md$slope, 4)} (se {signif(sb_md$se, 3)}, R2 {signif(sb_md$r2, 3)})",
    "*" = "split 1 STAGE B, mean, bracket-valid reps   : slope = {signif(sb_bo$slope, 4)} (se {signif(sb_bo$se, 3)}, R2 {signif(sb_bo$r2, 3)})",
    "*" = "split 1 Stage A (control)                   : slope = {signif(sa$slope, 4)}",
    "*" = "split 2 STAGE B, mean                       : slope = {signif(s2$slope, 4)} (se {signif(s2$se, 3)}, R2 {signif(s2$r2, 3)})",
    "*" = "split 2 STAGE B, median                     : slope = {signif(s2_md$slope, 4)}"
  ))

  # The verdict is taken from the bracket-valid read, because a clipped scan does
  # not test the claim at all. The unrestricted mean is reported alongside so the
  # contamination is visible rather than hidden.
  v_primary <- verdict_for(sb_bo$slope)
  v_mean    <- verdict_for(sb$slope)
  v_median  <- verdict_for(sb_md$slope)
  v_split2  <- verdict_for(s2$slope)

  cli::cli_h1("VERDICT")
  cli::cli_inform(c(
    "*" = "PRIMARY (split 1, bracket-valid reps): {v_primary}",
    "*" = "split 1, median: {v_median}",
    "*" = "split 1, unrestricted mean: {v_mean}",
    "*" = "split 2, mean: {v_split2}"
  ))
  cli::cli_inform("reference slopes: -1 (jump / O_p(1/n)), -1/2, -1/3 (grid mesh n^(-1/3))")

  cli::cli_h1("Memory (per-replicate, after full gc)")
  print(cells$memory, row.names = FALSE, digits = 4)

  out <- list(
    cell_summary = summ,
    slope_split1_stage_b_mean = sb[c("slope", "se", "r2", "n_points")],
    slope_split1_stage_b_median = sb_md[c("slope", "se", "r2", "n_points")],
    slope_split1_stage_b_bracket_ok = sb_bo[c("slope", "se", "r2", "n_points")],
    slope_split1_stage_a = sa[c("slope", "se", "r2", "n_points")],
    slope_split2_stage_b_mean = s2[c("slope", "se", "r2", "n_points")],
    slope_split2_stage_b_median = s2_md[c("slope", "se", "r2", "n_points")],
    verdict_primary = v_primary,
    verdict_median = v_median,
    verdict_unrestricted_mean = v_mean,
    verdict_split2 = v_split2,
    bracket_misses = bad,
    results = res,
    memory = cells$memory,
    meta = cells$meta,
    analyzed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  path <- file.path(STUDY, "results", "summary_slope.rds")
  saveRDS(out, path)
  cli::cli_alert_success("wrote {path}")
  invisible(out)
}

if (!interactive() && sys.nframe() == 0L) analyze()
