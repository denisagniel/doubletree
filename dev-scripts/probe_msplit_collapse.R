# Probe the M-split "rep collapse" that archived sims show at n>=1000 on DGP3
# (445->295->125 reps as n goes 500->1000->2000). n=500 showed 0% errors in the
# first diagnostic, so the collapse is an n>=1000 phenomenon. Capture the FULL
# error message + which stage throws, so we know if it's CV, GOSDT/solver, refit,
# or prediction.

suppressPackageStartupMessages({
  library(optimaltrees)
  library(doubletree)
})
source(file.path("simulations", "six_approach_comparison", "code", "dgps.R"))
set.seed(20260713)

REPS <- 15L
M <- 10L; K <- 5L

probe <- function(n) {
  cat(sprintf("\n=== DGP3 n=%d, %d reps ===\n", n, REPS)); flush.console()
  for (r in seq_len(REPS)) {
    dat <- generate_dgp_complex(n)
    res <- tryCatch({
      fit <- estimate_att_msplit(dat$X, dat$A, dat$Y, M = M, K = K,
                                 seed_base = 1000L * r, verbose = FALSE,
                                 outcome_type = "binary")
      sprintf("ok theta=%.4f", fit$theta)
    }, error = function(e) paste0("ERR: ", conditionMessage(e)))
    cat(sprintf("  rep %2d: %s\n", r, substr(res, 1, 150))); flush.console()
  }
}

probe(1000L)
probe(2000L)
