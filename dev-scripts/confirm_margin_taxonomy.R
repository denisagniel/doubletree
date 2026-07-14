# Confirm the "assumption-boundary, not difficulty-continuum" framing:
# On a MARGIN-SATISFYING DGP (unique oracle structure), M-split should track
# cross-fitting. On a MARGIN-VIOLATING DGP (diffuse outcome structure), M-split
# should show excess variance while cross-fitting stays fine.
#
# Metrics per (DGP, method): bias, SD(theta), coverage. For M-split also report
# median outcome-structure modal frequency (the margin proxy: high = unique mode).
# Matched seeds across methods within each rep.

suppressPackageStartupMessages({
  library(optimaltrees)
  devtools::load_all(".", quiet = TRUE)
})
source(file.path("simulations", "six_approach_comparison", "code", "dgps.R"))
set.seed(20260713)

REPS <- 30L
M <- 10L; K <- 5L
n <- 1000L

run_dgp <- function(gen, label) {
  cat(sprintf("\n=== %s  n=%d  %d reps ===\n", label, n, REPS)); flush.console()
  cf_th <- cf_cov <- ms_th <- ms_cov <- ms_freq <- rep(NA_real_, REPS)
  for (r in seq_len(REPS)) {
    dat <- gen(n); sb <- 1000L * r
    cf <- tryCatch(estimate_att(dat$X, dat$A, dat$Y, K = K, seed = sb,
                    use_rashomon = FALSE, outcome_type = "binary"),
                   error = function(e) NULL)
    ms <- tryCatch(estimate_att_msplit(dat$X, dat$A, dat$Y, M = M, K = K,
                    seed_base = sb, structure_selection = "modal",
                    outcome_type = "binary"),
                   error = function(e) NULL)
    if (!is.null(cf)) { cf_th[r] <- cf$theta
      cf_cov[r] <- (0.15 >= cf$ci_95[1]) && (0.15 <= cf$ci_95[2]) }
    if (!is.null(ms)) { ms_th[r] <- ms$theta
      ms_cov[r] <- (0.15 >= ms$ci_95[1]) && (0.15 <= ms$ci_95[2])
      ms_freq[r] <- ms$diagnostics$structure_frequency_m0 }
    cat(sprintf("  rep %2d: cf=%.4f  ms=%.4f  freq_m0=%.2f\n",
                r, cf_th[r], ms_th[r], ms_freq[r])); flush.console()
  }
  cat(sprintf("  CROSSFIT  bias=%+.4f SD=%.4f cov=%.2f\n",
              mean(cf_th,na.rm=TRUE)-0.15, sd(cf_th,na.rm=TRUE), mean(cf_cov,na.rm=TRUE)))
  cat(sprintf("  M-SPLIT   bias=%+.4f SD=%.4f cov=%.2f  median freq_m0=%.2f\n",
              mean(ms_th,na.rm=TRUE)-0.15, sd(ms_th,na.rm=TRUE), mean(ms_cov,na.rm=TRUE),
              median(ms_freq,na.rm=TRUE)))
  cat(sprintf("  --> SD ratio M-split/crossfit = %.2f\n",
              sd(ms_th,na.rm=TRUE)/sd(cf_th,na.rm=TRUE)))
  flush.console()
}

run_dgp(generate_dgp_simple,  "DGP1-simple (margin-satisfying)")
run_dgp(generate_dgp_complex, "DGP3-complex (margin-violating)")
