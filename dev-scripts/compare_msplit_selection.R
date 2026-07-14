# Head-to-head: modal vs lowest_risk structure selection for M-split on DGP3.
# Matched seeds per rep so the ONLY difference is the selection rule. We expect
# lowest_risk to reduce theta variance (SD across reps) on the diffuse-mode
# outcome nuisance, while keeping a single interpretable tree.

suppressPackageStartupMessages({
  library(optimaltrees)
  devtools::load_all(".", quiet = TRUE)  # pick up local R/ changes without install
})
source(file.path("simulations", "six_approach_comparison", "code", "dgps.R"))
set.seed(20260713)

REPS <- 30L
M <- 10L; K <- 5L

run <- function(gen, label, n) {
  cat(sprintf("\n=== %s n=%d ===\n", label, n)); flush.console()
  th_modal <- th_risk <- cov_modal <- cov_risk <- rep(NA_real_, REPS)
  for (r in seq_len(REPS)) {
    dat <- gen(n)
    sb <- 1000L * r
    fm <- tryCatch(estimate_att_msplit(dat$X, dat$A, dat$Y, M = M, K = K,
                    seed_base = sb, structure_selection = "modal",
                    outcome_type = "binary"), error = function(e) NULL)
    fr <- tryCatch(estimate_att_msplit(dat$X, dat$A, dat$Y, M = M, K = K,
                    seed_base = sb, structure_selection = "lowest_risk",
                    outcome_type = "binary"), error = function(e) NULL)
    if (!is.null(fm)) {
      th_modal[r] <- fm$theta
      cov_modal[r] <- (0.15 >= fm$ci_95[1]) && (0.15 <= fm$ci_95[2])
    }
    if (!is.null(fr)) {
      th_risk[r] <- fr$theta
      cov_risk[r] <- (0.15 >= fr$ci_95[1]) && (0.15 <= fr$ci_95[2])
    }
    cat(sprintf("  rep %2d: modal=%.4f  risk=%.4f\n", r,
                th_modal[r], th_risk[r])); flush.console()
  }
  summ <- function(th, cov, nm) {
    cat(sprintf("  %-11s bias=%+.4f  SD(theta)=%.4f  RMSE=%.4f  cov=%.2f (n=%d)\n",
                nm, mean(th, na.rm=TRUE) - 0.15, sd(th, na.rm=TRUE),
                sqrt(mean((th - 0.15)^2, na.rm=TRUE)), mean(cov, na.rm=TRUE),
                sum(!is.na(th))))
  }
  summ(th_modal, cov_modal, "modal")
  summ(th_risk,  cov_risk,  "lowest_risk")
  flush.console()
}

run(generate_dgp_complex, "DGP3-complex", 500L)
run(generate_dgp_complex, "DGP3-complex", 1000L)
