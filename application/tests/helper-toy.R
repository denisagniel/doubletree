## ============================================================================
## application/tests/helper-toy.R
##
## Hand-built, deterministic toy frames shaped like the real tables. Sourced by
## testthat before every test-*.R in this directory.
##
## SCOPE, deliberately narrow. These fixtures exist to prove that THIS
## DIRECTORY's code runs and that its assertions fire -- the covariate selection,
## the prior_cost discretization, the design-matrix guard, and one estimator call
## at the real SHAPE. They are not, and cannot be, a stand-in for the data:
## smi_fixture() is blocked upstream (no census -- see 01_declare_requirements.R),
## and inventing a distribution here would be exactly the retired synthetic-data
## strategy smidata was built to replace.
##
## What is NOT tested here, on purpose:
##   * smidata's promoted smi_*() functions. They are tested in smidata's own
##     suite (tests/testthat/test-msr_aap.R, test-cost_claims.R). Re-testing them
##     from this repo would re-create the duplication their promotion removed.
##   * doubletree's estimators themselves. tests/testthat/ (the package suite)
##     covers those. This directory tests only that the APPLICATION calls them at
##     a shape they accept.
##
## Every expected value in the test files is HAND-CALCULATED and stated in a
## comment. An expectation produced by running the code proves nothing.
## ============================================================================

## Resolve application/ from wherever testthat set the working directory.
app_dir_for_tests <- local({
  candidates <- c("application", ".", "..", file.path("..", ".."))
  hit <- candidates[vapply(
    candidates, function(d) file.exists(file.path(d, "_config.R")), logical(1)
  )]
  if (length(hit) == 0L) {
    ## test_dir() runs with wd = the test directory, so application/ is one up.
    hit <- if (file.exists(file.path("..", "_config.R"))) ".." else NA_character_
  }
  if (is.na(hit[1])) {
    stop("Cannot locate application/_config.R from ", getwd())
  }
  hit[1]
})

source(file.path(app_dir_for_tests, "_config.R"))
source(file.path(app_dir_for_tests, "helpers", "covariate_blocks.R"))
source(file.path(app_dir_for_tests, "helpers", "discretize_prior_cost.R"))
source(file.path(app_dir_for_tests, "helpers", "design_matrix.R"))
source(file.path(app_dir_for_tests, "helpers", "complete_case.R"))

#' Toy `larger_smi_covariates`-shaped frame: the 12 confirmed columns, 0/1 coded.
#'
#' Column values are a deterministic bit pattern: column j alternates with period
#' 2^(((j - 1) mod 6) + 1). Every column therefore has both levels present at
#' n >= 64, and the whole frame is reproducible by hand from `n` alone. The period
#' cycles every 6 columns, so some column PAIRS are identical -- harmless here,
#' since no test asserts anything about between-column structure.
#'
#' @param n Number of patients. Must be at least 64 (= 2^6): the longest-period
#'   column needs 64 rows to show both levels, and a constant column would
#'   (correctly) trip `assert_binary_design_matrix()`.
#' @return Tibble with `ID`, `INDEX_DT`, `prior_cost`, and the 12 `_YN` columns.
toy_covariates <- function(n = 64L) {
  if (n < 64L) {
    stop("toy_covariates() needs n >= 64 so every column has both levels; got ", n)
  }
  yn_cols <- unlist(covariate_blocks(), use.names = FALSE)
  bits <- lapply(seq_along(yn_cols), function(j) {
    ## Period 2^(((j - 1) mod 6) + 1): both levels present for every column at
    ## n >= 64, and no column all-0 or all-1.
    period <- 2^((j - 1L) %% 6L + 1L)
    as.integer((seq_len(n) - 1L) %% period < period / 2)
  })
  names(bits) <- yn_cols

  tibble::as_tibble(c(
    list(
      ID = sprintf("toy_%03d", seq_len(n)),
      INDEX_DT = as.Date("2019-01-15"),
      ## Deterministic, strictly increasing, zero-free: 0 is exercised
      ## separately in test-discretize-prior-cost.R's zero-inflation cases.
      prior_cost = as.numeric(seq_len(n)) * 100
    ),
    bits
  ))
}

#' Toy `msr_aap`-shaped frame: one measurement period per patient.
#'
#' Shape only -- `MSR_YR` is smidata's confirmed `RY<yyyy>-<mm>` format, and
#' whether the values are plausible is beside the point. Used to confirm this
#' directory composes with `smidata::smi_select_aap_row()`'s column contract, not
#' to test that function.
#'
#' @param ids Patient IDs.
#' @return Tibble with `ID`, `MSR`, `MSR_YR`, `MSR_DEN`, `MSR_NUM`.
toy_msr_aap <- function(ids) {
  tibble::tibble(
    ID = ids,
    MSR = "AAP",
    MSR_YR = "RY2020-01",
    MSR_DEN = 10,
    ## Alternating 6/10 and 3/10 straddles the 0.5 achievement threshold, so
    ## both exposure arms are represented.
    MSR_NUM = rep(c(6, 3), length.out = length(ids))
  )
}

#' Toy cost-claims-shaped frame.
#'
#' Shape only: the five columns `SMI_COLS$cost_claims` names. One claim per
#' patient, one excluded Managed Care Invoice row, so the MCI exclusion has
#' something to match on if a caller wants to exercise it.
#'
#' @param ids Patient IDs.
#' @return Tibble with `ID`, `SRV_DT`, `AMOUNT_PAID`, `SEQ_ID`, `SS_DESC`.
toy_cost_claims <- function(ids) {
  tibble::tibble(
    ID = c(ids, ids[1]),
    SRV_DT = as.Date("2020-03-01"),
    AMOUNT_PAID = c(rep(50, length(ids)), 9999),
    SEQ_ID = seq_len(length(ids) + 1L),
    SS_DESC = c(rep("Inpatient", length(ids)), "Managed Care Invoice")
  )
}

#' A known-truth toy `(X, A, Y)` at the real 15-column binary shape.
#'
#' @description
#' Both nuisances are EXACTLY representable within `leaf_budget` leaves on this
#' grid, by construction: each depends on a single binary covariate, so a 2-leaf
#' tree represents it exactly and a 4-leaf budget is more than sufficient.
#'
#' \itemize{
#'   \item propensity  \eqn{e(X) = 0.3 + 0.4 \cdot} `HOUSING_YN`  (2 leaves)
#'   \item control outcome \eqn{\mu_0(X) = 100 + 900 \cdot} `VAR5G_DM2_YN` (2 leaves)
#'   \item \eqn{Y = \mu_0(X) + \tau A + \epsilon}, constant treatment effect
#'         \eqn{\tau}, so the ATT is \eqn{\tau} exactly.
#' }
#'
#' \strong{What this establishes, and what it does not.} It shows grid-exact
#' sparsity is SATISFIABLE at this covariate shape in principle -- 15 binary
#' columns and a 4-leaf budget are not self-contradictory. It says NOTHING about
#' whether sparsity holds on the real joint distribution of the 15 confirmed
#' covariates, which is not checkable from data at all (see
#' `OPEN_DECISIONS$diagnostics` and `07_estimate_att.R`'s header). The 13
#' covariates the DGP ignores are there precisely so the estimator faces the real
#' grid dimension rather than a 2-column toy.
#'
#' @param n Number of observations.
#' @param tau True constant treatment effect (the ATT).
#' @param sigma Residual SD. Small relative to `tau` by default, so the smoke
#'   test is about mechanics rather than power.
#' @param seed RNG seed.
#' @return List with `X` (data.frame, 15 binary columns in declared order), `A`,
#'   `Y`, and `tau`.
toy_att_data <- function(n = 400L, tau = 250, sigma = 25, seed = 1L) {
  set.seed(seed)
  cov <- toy_covariates(n = n)
  dummies <- discretize_prior_cost(cov$prior_cost)
  x_with_id <- assemble_design_matrix(cov, dummies)
  X <- as.data.frame(x_with_id[, design_matrix_columns(), drop = FALSE])

  e <- 0.3 + 0.4 * X$HOUSING_YN
  A <- stats::rbinom(n, 1L, e)
  mu0 <- 100 + 900 * X$VAR5G_DM2_YN
  Y <- mu0 + tau * A + stats::rnorm(n, 0, sigma)

  list(X = X, A = as.integer(A), Y = as.numeric(Y), tau = tau)
}
