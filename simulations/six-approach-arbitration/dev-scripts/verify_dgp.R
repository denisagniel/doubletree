#!/usr/bin/env Rscript
# verify_dgp.R -- confirm the strengthened DGPs (2026-07-08) satisfy:
#   (1) empirical true ATT ~= 0.15 (unchanged by strengthening propensity)
#   (2) overlap/positivity: e_true away from 0/1 (no near-violations)
#   (3) propensity STRUCTURE recoverable by CV at n=500 (non-stump)
# Run from the study dir; uses installed/loaded optimaltrees + the study's dgp.R.
# Runs each CV probe in a killable subprocess (depth-capped; explosion already fixed).
suppressPackageStartupMessages({ library(callr); library(ps) })
PKG_OT <- normalizePath(file.path("..", "..", "..", "..", "optimaltrees"), mustWork = FALSE)
if (!dir.exists(PKG_OT)) PKG_OT <- "/Users/dagniel/RAND/rprojects/global-scholars/optimaltrees"
STUDY  <- normalizePath(file.path(".."))            # dev-scripts/ -> study dir
if (!file.exists(file.path(STUDY, "R", "dgp.R")))
  STUDY <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/six-approach-arbitration"

source(file.path(STUDY, "R", "dgp.R"))

# --- (1) & (2): estimand + overlap, analytic over a large sample -------------
# True ATT = E[Y1 - Y0 | A=1]. We know potential outcomes in the DGP, so compute
# directly on a big draw (this is a property of the DGP, independent of estimation).
true_att_and_overlap <- function(dgp, n = 200000L) {
  set.seed(20260708)
  # Re-run the generator but retain Y1,Y0,e_true by re-implementing the reveal.
  # Easiest: call the private generator, then recover ATT via a second draw of
  # potential outcomes at the same covariates is not exposed -> instead draw large
  # and use the identity ATT = mean(mu0+0.15 - mu0 | A=1) = 0.15 by construction,
  # but we VERIFY empirically by regenerating Y1,Y0 from the same mechanism.
  gen <- get(paste0(".dgp_", dgp))
  # Instrument: temporarily capture e_true and potential outcomes by re-sourcing a
  # patched generator is heavy; instead we recompute e_true here from the same coefs
  # by drawing X the same way is fragile. Simpler + robust: Monte Carlo the ATT via
  # the observed-data definition is not identified without a model. So we verify the
  # CONSTRUCTION invariant instead: draw many datasets, for each compute the
  # oracle ATT using the KNOWN potential outcomes exposed by a thin wrapper.
  NA
}

# Thin oracle wrapper: replicate each DGP's potential-outcome mechanism to get the
# TRUE ATT and overlap without depending on estimation. Mirrors dgp.R exactly.
oracle <- function(dgp, n = 200000L) {
  set.seed(20260708)
  expit <- function(x) 1 / (1 + exp(-x))
  if (dgp == "simple") {
    X <- data.frame(x1 = rbinom(n,1,.5), x2 = rbinom(n,1,.5), x3 = rbinom(n,1,.5))
    e <- expit(-0.7 + 1.0*X$x1 + 0.8*X$x2); mu0 <- 0.2 + 0.15*X$x1 + 0.15*X$x3
  } else if (dgp == "moderate") {
    X <- data.frame(x1=rbinom(n,1,.5),x2=rbinom(n,1,.5),x3=rbinom(n,1,.5),x4=rbinom(n,1,.5))
    e <- expit(-0.7 + 0.9*X$x1 + 0.7*X$x2 + 0.9*X$x1*X$x2)
    mu0 <- 0.2 + 0.2*X$x3 + 0.15*X$x4 + 0.2*X$x3*X$x4
  } else if (dgp == "complex") {
    X <- data.frame(x1=rbinom(n,1,.5),x2=rbinom(n,1,.5),x3=rbinom(n,1,.5),x4=rbinom(n,1,.5),x5=rbinom(n,1,.5))
    e <- expit(-0.9 + 0.6*(X$x1+X$x2+X$x3) + 0.8*X$x1*X$x2 + 0.6*X$x2*X$x3)
    mu0 <- 0.05 + 0.15*(X$x3+X$x4+X$x5) + 0.2*X$x3*X$x4 + 0.15*X$x4*X$x5
  } else if (dgp == "continuous") {
    X <- data.frame(x1=rbinom(n,1,.5),x2=rbinom(n,1,.5),x3=runif(n,-1,1),x4=rnorm(n,0,1))
    e <- expit(-0.5 + 0.8*X$x1 + 1.0*X$x3 + 0.4*X$x4 + 0.5*X$x1*X$x3)
    mu0 <- pmax(.01, pmin(.99, 0.2 + 0.15*X$x2 + 0.2*X$x3 + 0.15*(X$x4^2/2) + 0.1*X$x2*X$x3))
  } else stop("unknown dgp")
  mu1 <- pmin(mu0 + 0.15, 1)
  A <- rbinom(n, 1, e)
  # True ATT = E[mu1 - mu0 | A = 1]
  att <- weighted.mean(mu1 - mu0, w = A)      # among treated
  list(att = att, e_min = min(e), e_max = max(e), treated = mean(A),
       frac_extreme = mean(e < 0.02 | e > 0.98))
}

# --- (3): CV recovers propensity structure at n=500 (killable subprocess) ----
cv_structure <- function(dgp, n = 500L, wall_s = 90) {
  proc <- callr::r_bg(function(dgp, n, pkg, study) {
    devtools::load_all(pkg, quiet = TRUE); source(file.path(study, "R", "dgp.R"))
    set.seed(20260708)
    d <- generate_data(list(dgp = dgp, n = n))
    cv <- optimaltrees::cv_regularization_adaptive(
      X = d$X, y = d$A, loss_function = "log_loss", K = 5,
      max_lambda = 20*log(n)/n, refit = TRUE, verbose = FALSE, max_depth = 4L)
    p <- predict(cv$model, d$X, type = "prob")[, 2]
    list(cv_c = cv$best_lambda/(log(n)/n), leaves = length(unique(round(p, 6))))
  }, args = list(dgp=dgp, n=n, pkg=PKG_OT, study=STUDY), supervise = TRUE)
  t0 <- proc.time()[["elapsed"]]
  repeat { if (!proc$is_alive()) break
    if (proc.time()[["elapsed"]] - t0 > wall_s) { proc$kill_tree(); return(list(cv_c=NA, leaves=NA)) }
    Sys.sleep(0.3) }
  proc$wait(timeout = 5000)
  r <- tryCatch(proc$get_result(), error = function(e) list(cv_c=NA, leaves=NA))
  r
}

cat(sprintf("%-11s %8s %7s %7s %8s %9s | %8s %7s\n",
            "dgp", "ATT", "e_min", "e_max", "treated", "extreme", "CV_c@500", "leaves"))
cat(strrep("-", 78), "\n")
for (dgp in c("simple","moderate","complex","continuous")) {
  o <- oracle(dgp)
  s <- cv_structure(dgp)
  cat(sprintf("%-11s %8.4f %7.3f %7.3f %8.0f%% %8.2f%% | %8.2f %7s\n",
              dgp, o$att, o$e_min, o$e_max, 100*o$treated, 100*o$frac_extreme,
              ifelse(is.na(s$cv_c), -9, s$cv_c),
              ifelse(is.na(s$leaves), "KILLED", as.character(s$leaves))))
}
cat("\nPASS if: ATT ~= 0.15 (all), e in ~[.05,.95] (overlap ok), extreme ~0%, leaves > 1 (non-stump).\n")
