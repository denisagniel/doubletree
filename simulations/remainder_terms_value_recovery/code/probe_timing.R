# ---------------------------------------------------------------------------
# Timing / feasibility probe.
#
# Purpose: measure the wall-clock cost of ONE budget-constrained optimal-tree
# fit at each (n, Lbar) cell of the planned design, and calibrate the
# regularization lambda that lands on the budget, so the main run can use a
# single fit_tree() call per nuisance instead of a 40-fit bisection.
#
# Run: Rscript code/probe_timing.R           (from the study directory)
# ---------------------------------------------------------------------------

library(optimaltrees)
library(readr)
library(fs)

set.seed(20260908)

STUDY_DIR <- if (interactive()) "." else dirname(dirname(normalizePath(sys.frame(1)$ofile %||% "code/x")))
source(fs::path("code", "dgp.R"))

N_GRID <- c(2000L, 5000L, 10000L, 20000L)
L_GRID <- c(5L, 10L, 20L, 30L)
FIT_TIME_LIMIT <- 120                              # seconds per internal fit
PROBE_CELLS <- expand.grid(n = N_GRID, L = L_GRID, KEEP.OUT.ATTRS = FALSE)

max_depth_for <- function(L) max(2L, ceiling(log2(L)))

probe_one <- function(n, L) {
  truth <- make_truth(L)
  dat <- simulate_data(n, truth)
  lambda_n <- log(n) / n

  t_e <- system.time({
    e_fit <- optimaltrees::bisect_lambda_to_budget(
      as.data.frame(dat$Xbin), dat$A,
      leaf_budget = L, lambda_n = lambda_n, m_n = 1L,
      loss_function = "log_loss",
      max_depth = max_depth_for(L), depth_restricted = TRUE,
      fit_time_limit = FIT_TIME_LIMIT, verbose = FALSE
    )
  })[["elapsed"]]

  ctrl <- dat$A == 0L
  t_m <- system.time({
    m_fit <- optimaltrees::bisect_lambda_to_budget(
      as.data.frame(dat$Xbin[ctrl, , drop = FALSE]), dat$Y[ctrl],
      leaf_budget = L, lambda_n = lambda_n, m_n = 1L,
      loss_function = "squared_error",
      max_depth = max_depth_for(L), depth_restricted = TRUE,
      fit_time_limit = FIT_TIME_LIMIT, verbose = FALSE
    )
  })[["elapsed"]]

  list(n = n, L = L,
       lambda_e = e_fit$lambda, leaves_e = e_fit$n_leaves,
       lambda_m = m_fit$lambda, leaves_m = m_fit$n_leaves,
       sec_e = t_e, sec_m = t_m,
       e_names = paste(names(e_fit), collapse = ","))
}

cat("cells:", nrow(PROBE_CELLS), "\n")
out <- vector("list", nrow(PROBE_CELLS))
for (i in seq_len(nrow(PROBE_CELLS))) {
  cell <- PROBE_CELLS[i, ]
  cat(sprintf("[%2d/%2d] n=%6d L=%2d ... ", i, nrow(PROBE_CELLS), cell$n, cell$L))
  res <- tryCatch(probe_one(cell$n, cell$L), error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(res)) {
    cat(sprintf("e: %.1fs (%d leaves, lam=%.2e)  m: %.1fs (%d leaves, lam=%.2e)\n",
                res$sec_e, res$leaves_e, res$lambda_e,
                res$sec_m, res$leaves_m, res$lambda_m))
  }
  out[[i]] <- res
}

probe_tbl <- do.call(rbind, lapply(Filter(Negate(is.null), out), as.data.frame))
print(probe_tbl)
fs::dir_create("results")
readr::write_rds(probe_tbl, fs::path("results", "probe_timing.rds"))
