# Scoped current-ground-truth check for the six-approach comparison.
# Uses the sim's OWN estimator wrappers (code/estimators.R) and the SAME seed
# formula as run_single_replication.R, so results are directly comparable to the
# archived results/combined/summary_inference.csv.
#
# Scope (local, fast): approaches 2 (crossfit) and 5 (msplit) -- the pair central
# to "does M-split match crossfit" -- across all 4 DGPs at n=1000, REPS reps.
# The full 500-rep x 6 x 4 x 3 grid remains a cluster job (launch_adaptive_cv.sh).
#
# Run from the six_approach_comparison/ directory.

suppressPackageStartupMessages({
  library(optimaltrees)
  devtools::load_all("../..", quiet = TRUE)  # local doubletree R/ (repo is doubletree/)
})

code_dir <- "code"
source(file.path(code_dir, "dgps.R"))
source(file.path(code_dir, "estimators.R"))
source(file.path(code_dir, "metrics.R"))

REPS <- 50L
N <- 1000L
approaches <- list(`2` = estimate_att_crossfit, `5` = estimate_att_msplit)
approach_names <- c(`2` = "crossfit_separate", `5` = "msplit")
# Binary DGPs 1-3 only: DGP4-continuous is a separate documented hard case
# (tree-depth blow-up on continuous covariates in the sim wrappers) and is not
# central to the M-split-vs-crossfit question. Deferred to the cluster run.
dgp_map <- list(`1` = generate_dgp_simple, `2` = generate_dgp_moderate,
                `3` = generate_dgp_complex)
dgp_names <- c(`1` = "simple", `2` = "moderate", `3` = "complex")

rows <- list()
for (ap in names(approaches)) {
  for (dg in names(dgp_map)) {
    est <- approaches[[ap]]; gen <- dgp_map[[dg]]
    th <- se <- cov <- rep(NA_real_, REPS)
    for (rep in seq_len(REPS)) {
      # SAME seed formula as run_single_replication.R:114
      seed <- 1000000 + as.integer(ap) * 1000000 + as.integer(dg) * 100000 + N + rep
      set.seed(seed)
      dat <- gen(n = N)
      r <- tryCatch(est(X = dat$X, A = dat$A, Y = dat$Y),
                    error = function(e) list(theta = NA_real_, se = NA_real_))
      if (!is.na(r$theta) && !is.na(r$se)) {
        th[rep] <- r$theta; se[rep] <- r$se
        cov[rep] <- (dat$true_att >= r$theta - 1.96 * r$se) &&
                    (dat$true_att <= r$theta + 1.96 * r$se)
      }
    }
    rows[[length(rows) + 1]] <- data.frame(
      approach = as.integer(ap), approach_name = approach_names[[ap]],
      dgp = as.integer(dg), dgp_name = dgp_names[[dg]], n = N,
      reps_ok = sum(!is.na(th)),
      bias = mean(th, na.rm = TRUE) - 0.15,
      sd_theta = sd(th, na.rm = TRUE),
      coverage = mean(cov, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    cat(sprintf("A%s %-9s bias=%+.4f SD=%.4f cov=%.2f (ok %d/%d)\n",
                ap, dgp_names[[dg]], tail(rows,1)[[1]]$bias,
                tail(rows,1)[[1]]$sd_theta, tail(rows,1)[[1]]$coverage,
                tail(rows,1)[[1]]$reps_ok, REPS)); flush.console()
  }
}

cur <- do.call(rbind, rows)
readr::write_csv(cur, "ground_truth_current.csv")

# Diff vs archived summary (bias + coverage per approach/dgp at n=1000)
arch_path <- "results/combined/summary_inference.csv"
if (file.exists(arch_path)) {
  arch <- readr::read_csv(arch_path, show_col_types = FALSE)
  cat("\n=== DIFF vs archived (n=1000) ===\n")
  cat(sprintf("%-9s %-9s | %8s %8s | %8s %8s\n",
              "approach","dgp","cur_cov","arch_cov","cur_bias","arch_bias"))
  for (i in seq_len(nrow(cur))) {
    a <- cur[i,]
    m <- arch[arch$approach == a$approach & arch$dgp == a$dgp &
              arch$n == a$n, ]
    if (nrow(m) >= 1) {
      cat(sprintf("A%-8d %-9s | %8.2f %8.2f | %+8.4f %+8.4f\n",
                  a$approach, a$dgp_name, a$coverage,
                  m$coverage[1], a$bias, m$bias[1]))
    }
  }
} else {
  cat("\n(archived summary not found for diff)\n")
}
cat("\nWrote ground_truth_current.csv\n")
