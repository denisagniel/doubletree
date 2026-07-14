# Diagnose M-split failure on DGP3 (complex, tree-representable) vs DGP4 (continuous).
# Question: is the failure statistical (diffuse mode -> arbitrary structure) or
# computational (CV/solver error -> dropped rep)?
#
# Per replication (single estimate_att_msplit call; structure stats come from its
# own $diagnostics field so we do NOT re-run Stage 1):
#   - errored: did estimate_att_msplit throw? (this is the "rep collapse")
#   - freq_e/m0 : plurality of the winning partition across M splits
#   - ndist_e/m0: number of distinct partitions found across M splits
#   - theta, covered: point estimate and whether 95% CI covers true ATT
#
# Per-rep line flushed to stdout for live progress. Prototype scale.

suppressPackageStartupMessages({
  library(optimaltrees)
  library(doubletree)
})

source(file.path("simulations", "six_approach_comparison", "code", "dgps.R"))

set.seed(20260713)

REPS <- 25L
M    <- 10L
K    <- 5L

diagnose_one <- function(gen, n, seed_base) {
  dat <- gen(n)
  fit <- tryCatch(
    estimate_att_msplit(dat$X, dat$A, dat$Y, M = M, K = K,
                        seed_base = seed_base, verbose = FALSE,
                        outcome_type = "binary"),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "fit_err")
  )
  if (inherits(fit, "fit_err")) {
    return(data.frame(errored = TRUE, freq_e = NA, freq_m0 = NA,
                      ndist_e = NA, ndist_m0 = NA, theta = NA, covered = NA,
                      msg = fit$msg, stringsAsFactors = FALSE))
  }
  d <- fit$diagnostics
  covered <- (0.15 >= fit$ci_95[1]) && (0.15 <= fit$ci_95[2])
  data.frame(
    errored = FALSE,
    freq_e = d$structure_frequency_e, freq_m0 = d$structure_frequency_m0,
    ndist_e = length(d$structure_counts_e), ndist_m0 = length(d$structure_counts_m0),
    theta = fit$theta, covered = covered, msg = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_grid <- function(gen, label, n) {
  cat(sprintf("\n=== %s, n=%d, %d reps ===\n", label, n, REPS)); flush.console()
  rows <- vector("list", REPS)
  for (r in seq_len(REPS)) {
    res <- diagnose_one(gen, n, seed_base = 1000L * r)
    res$rep <- r
    rows[[r]] <- res
    cat(sprintf("  rep %2d: %s freq_e=%s freq_m0=%s ndist=(%s,%s) theta=%s\n",
                r,
                if (res$errored) "ERR " else "ok  ",
                ifelse(is.na(res$freq_e), "  NA", sprintf("%.2f", res$freq_e)),
                ifelse(is.na(res$freq_m0), "  NA", sprintf("%.2f", res$freq_m0)),
                ifelse(is.na(res$ndist_e), "NA", res$ndist_e),
                ifelse(is.na(res$ndist_m0), "NA", res$ndist_m0),
                ifelse(is.na(res$theta), "NA", sprintf("%.4f", res$theta))))
    flush.console()
  }
  df <- do.call(rbind, rows)
  ok <- df[!df$errored, ]
  cat(sprintf("  -- error/drop rate: %.0f%% (%d/%d)\n",
              100 * mean(df$errored), sum(df$errored), nrow(df)))
  if (any(df$errored)) {
    for (mm in unique(df$msg[df$errored])) {
      cat(sprintf("     ERR: %s\n", substr(mm, 1, 100)))
    }
  }
  if (nrow(ok) > 0) {
    cat(sprintf("  -- modal freq e : median %.2f [%.2f, %.2f]\n",
                median(ok$freq_e), min(ok$freq_e), max(ok$freq_e)))
    cat(sprintf("  -- modal freq m0: median %.2f [%.2f, %.2f]\n",
                median(ok$freq_m0), min(ok$freq_m0), max(ok$freq_m0)))
    cat(sprintf("  -- # distinct e : median %.0f max %d | m0: median %.0f max %d\n",
                median(ok$ndist_e), max(ok$ndist_e),
                median(ok$ndist_m0), max(ok$ndist_m0)))
    cat(sprintf("  -- bias: %.4f  coverage: %.2f (n_ok=%d)\n",
                mean(ok$theta) - 0.15, mean(ok$covered), nrow(ok)))
  }
  flush.console()
  df$label <- label; df$n <- n
  df
}

all_res <- rbind(
  run_grid(generate_dgp_complex,    "DGP3-complex",    500L),
  run_grid(generate_dgp_continuous, "DGP4-continuous", 500L)
)

out <- file.path("dev-scripts", "msplit_dgp3_diagnostic.csv")
readr::write_csv(all_res, out)
cat(sprintf("\nWrote %s\n", out))
