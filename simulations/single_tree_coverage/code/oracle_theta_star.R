# ============================================================
# oracle_theta_star.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
# Purpose: Compute the working-model ATT theta*_n -- the coverage target of
#   Theorem 1 -- per (DGP, analysis-n) cell, via a large oracle Monte-Carlo pass.
#   theta*_n =/= theta_0 in general: theta_0 (= true_att = 0.15) is the TRUE ATT,
#   whereas theta*_n is the ATT of the GRID-OPTIMAL piecewise-constant (PC) tree
#   model at the adaptive grid G_n the estimator uses at sample size n. On the
#   binary DGPs (structural-PC, Cor. 2) theta*_n = theta_0 up to grid resolution;
#   on `continuous` (approx-PC) the x4^2 term is not tree-representable so
#   theta*_n =/= theta_0 by construction (the discretization gap, Prop. 1).
# Inputs : none (draws from the DGPs).
# Outputs: results/theta_star_oracle.rds (cache read by the runners + analyze.R).
# Run    : from repo root ->
#   Rscript doubletree/simulations/single_tree_coverage/code/oracle_theta_star.R
# Run time: ~1-2 min (four 2e6 draws + grouped cell means).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(fs)
})
suppressMessages(devtools::load_all("doubletree", quiet = TRUE))

STUDY_DIR   <- fs::path("doubletree", "simulations", "single_tree_coverage")
source(fs::path(STUDY_DIR, "code", "dgps.R"))

# --- Study constants (shared with the runners via source()) ------------------
ORACLE_SEED <- 20260804L        # fixed oracle seed (per-DGP offset below)
N_ORACLE    <- 2000000L         # oracle Monte-Carlo sample size (~2e6)
N_GRID      <- c(500L, 1000L, 2000L, 5000L, 10000L)   # analysis sample sizes
DGP_NAMES   <- c("simple", "moderate", "complex", "continuous")
THETA0      <- 0.15             # true ATT (== true_att in every generate_dgp_*)

# Propensity clip mirrors estimate_att_single_tree (predict-time clamp to [0.01,0.99]).
PROP_LO <- doubletree:::.PROPENSITY_LOWER_BOUND
PROP_HI <- doubletree:::.PROPENSITY_UPPER_BOUND

#' Grid-cell id per observation at adaptive resolution b_n
#'
#' Reproduces the estimator's grid G_n: binary covariates pass through (a tree
#' fits the subcube exactly); continuous covariates are split at the SAME
#' population-quantile thresholds the package uses (`compute_thresholds`,
#' method "quantiles", b_n bins). The full product of per-feature bins is the
#' partition a deep-enough grid-optimal tree realizes, so cell membership defines
#' the grid-optimal PC nuisance.
#'
#' @param X data.frame of covariates (from a DGP oracle draw).
#' @param b_n Integer bins per continuous feature at the analysis n.
#' @return Character vector of length nrow(X): the grid-cell id per observation.
compute_grid_cell <- function(X, b_n) {
  ids <- purrr::map(X, function(x) {
    if (optimaltrees:::is_binary(x)) {
      return(as.integer(x))
    }
    uv <- unique(x[!is.na(x)])
    if (length(uv) <= 2) {
      return(as.integer(x == max(uv)))            # 2-level covariate -> {0,1}
    }
    thr <- optimaltrees:::compute_thresholds(x, method = "quantiles", n_bins = b_n)
    findInterval(x, sort(thr))                     # grid bin at population quantiles
  })
  do.call(paste, c(ids, sep = "|"))
}

#' theta*_n for one (DGP, analysis-n) cell, from the KNOWN truth (no MC noise)
#'
#' Computes the population working-model ATT theta*_n as the estimator TARGETS it.
#' Because the estimator cross-fits, it targets the population grid projections
#' (not an in-sample overfit), so the oracle must too. We therefore build the grid
#' projections and evaluate the ATT moment from the DGP's SMOOTH truth functions
#' (\code{e_true}, \code{mu0_true}) rather than the noisy Bernoulli draws, which
#' (a) removes the own-observation overfitting bias that otherwise grows with the
#' cell count (finer grid) and would spuriously inflate the continuous gap, and
#' (b) integrates out the outcome noise so theta*_n is pinned to ~1e-4 at 2e6.
#'
#' Population identities used (so this equals the observable-data projection):
#'   e_taus(cell)  = E[A|cell]     = E[e_true | cell]
#'   m0_taus(cell) = E[Y|A=0,cell] = E[(1-e_true) mu0_true | cell] / E[(1-e_true) | cell]
#' Taking conditional expectations given X in the orthogonal ATT score psi
#' (E[A|X]=e_true, E[Y|A=1,X]=mu1_true, E[Y|A=0,X]=mu0_true) and solving E[psi]=0:
#'   theta*_n = { E[e(mu1 - m0_taus)] - E[(e_taus/(1-e_taus))(1-e)(mu0 - m0_taus)] } / E[e]
#'
#' @param X data.frame of oracle covariates.
#' @param e_true True propensity e(X) = P(A=1|X) (returned by the DGP).
#' @param mu0_true True control mean mu0(X) = E[Y|A=0,X] (returned by the DGP).
#' @param dgp_name,n_analysis Cell identity (for the returned row).
#' @return One-row tibble: dgp, n, b_n, n_cells, theta_star, theta0, theta0_realized,
#'   gap, gap_vs_realized.
oracle_theta_star_cell <- function(X, e_true, mu0_true, dgp_name, n_analysis) {
  b_n  <- optimaltrees:::compute_bin_count("adaptive", n_analysis)
  cell <- compute_grid_cell(X, b_n)

  # E[Y|A=1,X]: the DGPs draw Y1 ~ Bern(pmin(mu0 + ATT, 1)) with ATT = 0.15, so the
  # treated conditional mean is the same clipped shift (the clip only bites in the
  # continuous DGP's upper tail, making the realized ATT slightly below nominal).
  mu1_true <- pmin(mu0_true + THETA0, 1)

  # Population grid projections from the smooth truth (no Bernoulli noise).
  e_taus_cell <- tapply(e_true, cell, mean)                  # E[e_true | cell]
  w0          <- 1 - e_true                                  # control-density weight
  num_cell    <- tapply(w0 * mu0_true, cell, sum)
  den_cell    <- tapply(w0, cell, sum)
  m0_taus_cell <- num_cell / den_cell                        # E[(1-e)mu0|cell]/E[(1-e)|cell]

  e_taus  <- pmin(pmax(as.numeric(e_taus_cell[cell]), PROP_LO), PROP_HI)  # mirror clip
  m0_taus <- as.numeric(m0_taus_cell[cell])
  # Guard degenerate cells (all-treated: no control mass) -> global control mean.
  if (anyNA(m0_taus) || any(!is.finite(m0_taus))) {
    global_m0 <- sum(w0 * mu0_true) / sum(w0)
    m0_taus[is.na(m0_taus) | !is.finite(m0_taus)] <- global_m0
  }

  # Population ATT moment E[psi(theta)] = 0 with (e_taus, m0_taus) plugged in.
  pi_true <- mean(e_true)
  num <- mean(e_true * (mu1_true - m0_taus) -
              (e_taus / (1 - e_taus)) * (1 - e_true) * (mu0_true - m0_taus))
  theta_star <- num / pi_true

  # Realized true ATT from the truth = E[mu1 - mu0 | A=1]; equals THETA0 exactly on
  # the binary DGPs, slightly below on `continuous` (upper-tail clipping of mu1).
  theta0_realized <- mean(e_true * (mu1_true - mu0_true)) / pi_true

  tibble::tibble(
    dgp = dgp_name, n = n_analysis, b_n = b_n, n_cells = length(unique(cell)),
    theta_star = theta_star, theta0 = THETA0, theta0_realized = theta0_realized,
    gap = theta_star - THETA0, gap_vs_realized = theta_star - theta0_realized
  )
}

#' Compute theta*_n across all (DGP, n) cells (one oracle draw reused per DGP)
#'
#' One 2e6 draw per DGP (the population is fixed; only the grid b_n varies with
#' the analysis n), so cells within a DGP share the same oracle sample.
#'
#' @param n_grid Analysis sample sizes.
#' @param dgp_names DGP names (generate_dgp_<name>).
#' @param n_oracle Oracle sample size.
#' @param base_seed Base oracle seed (per-DGP offset added).
#' @return Tibble of theta*_n rows across all cells.
compute_oracle_theta_star <- function(n_grid = N_GRID, dgp_names = DGP_NAMES,
                                      n_oracle = N_ORACLE, base_seed = ORACLE_SEED) {
  purrr::imap(dgp_names, function(dgp_name, j) {
    gen <- get(paste0("generate_dgp_", dgp_name))
    set.seed(base_seed + j)                        # deterministic per-DGP oracle draw
    d <- gen(n_oracle)
    purrr::map(n_grid, function(nn) {
      oracle_theta_star_cell(d$X, d$e_true, d$mu0_true, dgp_name, nn)
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind()
}

#' Load the cached oracle theta*_n table (compute + cache on first call)
#'
#' @param study_dir Study directory.
#' @param force Recompute even if the cache exists (default FALSE).
#' @return Tibble of theta*_n rows.
load_oracle <- function(study_dir = STUDY_DIR, force = FALSE) {
  cache <- fs::path(study_dir, "results", "theta_star_oracle.rds")
  if (!force && fs::file_exists(cache)) {
    return(readr::read_rds(cache))
  }
  tbl <- compute_oracle_theta_star()
  fs::dir_create(fs::path(study_dir, "results"))
  readr::write_rds(tbl, cache)
  tbl
}

# --- main: compute, cache, print theta*_n vs theta_0 -------------------------
# Guard so source()-ing this file (from the runners) defines the helpers WITHOUT
# recomputing; running via Rscript executes the pass.
if (sys.nframe() == 0L) {
  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("Computing oracle theta*_n (n_oracle = %s) ...\n",
              format(N_ORACLE, big.mark = ",", scientific = FALSE)))
  oracle_tbl <- load_oracle(force = TRUE)
  cat(sprintf("Done in %.1f s. Cached to %s\n\n",
              proc.time()[["elapsed"]] - t0,
              fs::path(STUDY_DIR, "results", "theta_star_oracle.rds")))

  cat("theta*_n vs theta_0 (= 0.15) by (DGP, n):\n")
  oracle_tbl |>
    dplyr::mutate(dplyr::across(c(theta_star, theta0_realized, gap, gap_vs_realized),
                                \(x) round(x, 4))) |>
    dplyr::select(dgp, n, b_n, n_cells, theta_star, theta0, theta0_realized,
                  gap, gap_vs_realized) |>
    as.data.frame() |>
    print(row.names = FALSE)

  cat("\nMax |gap vs realized ATT| by DGP (pure discretization bias;\n",
      "expect ~0 on binary, > 0 on continuous):\n", sep = "")
  oracle_tbl |>
    dplyr::summarise(max_abs_gap_realized = max(abs(gap_vs_realized)), .by = dgp) |>
    dplyr::mutate(max_abs_gap_realized = round(max_abs_gap_realized, 5)) |>
    as.data.frame() |>
    print(row.names = FALSE)
}
