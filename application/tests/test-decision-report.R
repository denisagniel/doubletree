## ============================================================================
## application/tests/test-decision-report.R
##
## The review report (08_decision_report.R + helpers/decision_report.R) on the
## toy intermediates from toy_report_inputs(). Every expected value is
## HAND-CALCULATED in a comment; none is copied from a run.
##
## Toy design (see toy_report_inputs() for the construction):
##   02: 200 base patients p001-p200; p001-p100 INDEX_DT 2020-03-01, p101-p200
##       2020-03-15. SAA rows p001-p180; a period starting exactly at s:
##       p001-p160; MSR_DEN >= 1: p001-p150; cost-file coverage: 136 remain.
##       Achieved: p001-p060 and p101-p130.
##   03: 136 patients, outcome [2020-06-01, 2021-06-01), prior [2018-06-01,
##       2019-06-01).
##   05: 136 -> 100 fully enrolled, 30 insufficient, 6 absent.
##   06: final sample 100, design column j (1..12) is 1 iff i %% (j+1) == 0
##       for i = 0..99; prior_cost = (row - 20) * 100.
## ============================================================================

skip_if_not_installed("smidata")
skip_if_not("smi_report" %in% getNamespaceExports("smidata"))

source(file.path(app_dir_for_tests, "helpers", "decision_report.R"))

toy <- toy_report_inputs()
rep_toy <- build_decision_report(toy, title = "toy")
tb <- rep_toy$tables

test_that("the levers are read from the config objects", {
  lv <- rep_toy$levers
  expect_match(lv$achievement_threshold, ">= 0.5", fixed = TRUE)
  expect_match(lv$period_rule, "anchor = \"ceiling_month\"", fixed = TRUE)
  expect_match(lv$negative_amount_policy, "^net")
  expect_match(lv$mci_exclusion, "Managed Care Invoice", fixed = TRUE)
  ## Windows: outcome on s, prior_cost on INDEX_DT, complete case 24 months.
  expect_match(lv$windows, "outcome [s+12mo, s+24mo)", fixed = TRUE)
  expect_match(lv$windows, "prior_cost [INDEX_DT-12mo, INDEX_DT)", fixed = TRUE)
  expect_match(lv$windows, "complete case: [s, s+24mo)", fixed = TRUE)
  expect_match(lv$enrollment_rule, "type_policy = require_agreement", fixed = TRUE)
  expect_match(lv$enrollment_rule, "gap_months = 0, gap_mode = total", fixed = TRUE)
  expect_match(lv$cost_families, "18 files", fixed = TRUE)
  expect_match(lv$prior_cost_floor, "set to 0", fixed = TRUE)
  expect_match(lv$winsorizing, "^not applicable")
  expect_match(lv$provider_attribution, "^not applicable")
  ## 15 design columns listed.
  expect_equal(length(strsplit(lv$covariate_set, ", ")[[1]]), 15L)
})

## Property: no masked cell is recoverable. Each identity below links counts that
## are published together; the number of unknown terms must not be exactly one.
n_unknown <- function(...) sum(is.na(c(...)))

test_that("the attrition table merges a small step instead of masking its drop", {
  fl <- tb$flow
  ## 200 -> 180 (20 with no SAA rows) -> 160 (20 with rows but none at s).
  expect_equal(fl$n[1:3], c(200L, 180L, 160L))
  expect_equal(fl$n_dropped[1:3], c(0L, 20L, 20L))
  ## 160 -> 150 drops 10 (< 11). Showing both 160 and 150 would publish the 10, so
  ## the 150 row is merged into the next step: 160 -> 136 drops 24.
  expect_false(150L %in% fl$n)
  expect_equal(fl$n[4], 136L)
  expect_equal(fl$n_dropped[4], 24L)
  expect_match(fl$step[4], "MSR_DEN >= 1", fixed = TRUE)
  expect_match(fl$step[4], "cost-file coverage", fixed = TRUE)
  expect_match(fl$reason[4], "MSR_DEN < 1", fixed = TRUE)
  ## complete case: 136 -> 100, 36 dropped, 100/200 = 50% of the start.
  expect_equal(fl$n[5], 100L)
  expect_equal(fl$n_dropped[5], 36L)
  expect_equal(fl$pct_of_start[5], 50)
  expect_match(fl$step[5], "[s, s+24mo)", fixed = TRUE)
  ## the asserted-zero step is listed and drops nobody.
  expect_equal(fl$n_dropped[6], 0L)
  ## no two published counts of the chain differ by 1..10.
  gaps <- diff(fl$n[!is.na(fl$n)])
  expect_false(any(gaps < 0 & gaps > -11))
})

test_that("hide_small_gaps() hides the inner count of a small drop, repeatedly", {
  ## 100 -> 95 (5): 95 is inner and is hidden; then 100 -> 90 (10) hides 90; then
  ## 100 -> 60 is fine.
  expect_equal(hide_small_gaps(c(100, 95, 90, 60), 11L), c(FALSE, TRUE, TRUE, FALSE))
  ## last count kept; the previous (inner) count is hidden instead.
  expect_equal(hide_small_gaps(c(100, 50, 45), 11L), c(FALSE, TRUE, FALSE))
  ## two counts only: the last is hidden (nothing inner to hide).
  expect_equal(hide_small_gaps(c(100, 95), 11L), c(FALSE, TRUE))
  ## drops of 0 or >= 11 are fine.
  expect_equal(hide_small_gaps(c(100, 100, 89, 50), 11L), rep(FALSE, 4))
  fl <- data.frame(step = c("a", "b", "c", "d"), n = c(100L, 95L, 90L, 60L),
                   reason = c("", "rb", "rc", "rd"), stringsAsFactors = FALSE)
  co <- coarsen_flow(fl, 11L)
  expect_equal(co$n, c(100L, 60L))
  expect_equal(co$reason[2], "rb; rc; rd")
})

test_that("exposure-period counts are published once; a hidden eligible count is not recoverable", {
  mp <- tb$measure_period.summary
  v <- stats::setNames(mp$n, mp$quantity)
  ## patients / matched / unmatched / den_below_1 live only in the flow.
  expect_equal(mp$quantity, c("eligible", "achieved"))
  ## eligible = 150 is the hidden count above (160 - 150 = 10), so it is not
  ## published here, and neither table that decomposes it is.
  expect_true(is.na(v[["eligible"]]))
  expect_equal(v[["achieved"]], 90L)             # 60 + 30
  expect_named(tb$measure_period.by_period_start, "note")
  expect_named(tb$measure_period.index_day, "note")
})

test_that("with no small step, eligible and its decompositions are published", {
  big <- toy
  big$p02$selection$elig <- seq_len(200) <= 140L   # 160 -> 140: drop 20, >= 11
  big$p02$selection$MSR_DEN <- as.numeric(big$p02$selection$elig)
  big$p02$selection$achieved <- dplyr::if_else(big$p02$selection$elig,
                                               seq_len(200) %in% c(1:60, 101:130), NA)
  big$p02$flow$n[4] <- 140L
  rb <- build_decision_report(big, title = "toy-big")$tables
  mp <- stats::setNames(rb$measure_period.summary$n, rb$measure_period.summary$quantity)
  expect_equal(mp[["eligible"]], 140L)
  expect_equal(mp[["achieved"]], 90L)
  expect_equal(rb$measure_period.by_period_start$n_eligible, c(100L, 40L))
  expect_equal(rb$flow$n[4], 140L)
})

test_that("cost tables: negative lines, floors, quantiles, MCI and categories", {
  co <- tb$costs_cohort
  pr <- co[co$window == "prior_cost", ]
  ou <- co[co$window == "outcome", ]
  ## prior_cost: +100 each for 136 (gross 13,600); p001-p019 also -150
  ## (19 lines, -2,850). Net 13,600 - 2,850 = 10,750; mean 10,750 / 136.
  expect_equal(pr$n_patients, 136L)
  expect_equal(pr$n_claims, 136L + 19L)
  expect_equal(pr$n_negative_lines, 19L)
  expect_equal(pr$negative_dollars, -2850)
  expect_equal(pr$gross_positive_dollars, 13600)
  expect_equal(pr$total_cost, 10750)
  expect_equal(pr$mean_cost, 10750 / 136)
  ## 19 patients net 100 - 150 = -50: 19 cells with a negative total, and 19/136
  ## = 14% of cells are at -50, so p01 and p05 are -50 and p25 is 100.
  expect_equal(pr$n_cells_negative_total, 19L)
  expect_equal(c(pr$p01, pr$p05, pr$p25), c(-50, -50, 100))
  expect_equal(pr$negative_share_of_gross, round(2850 / 13600, 5))
  ## outcome: 200 + 100 each = 300 (gross 136 * 300 = 40,800); p001-p012 also
  ## -50 (12 lines, -600). Net 40,200. 12 / 136 = 8.8% of cells net 250.
  expect_equal(ou$n_claims, 2L * 136L + 12L)
  expect_equal(ou$negative_dollars, -600)
  expect_equal(ou$total_cost, 40200)
  expect_equal(c(ou$p01, ou$p05, ou$p25), c(250, 250, 300))
  expect_equal(ou$n_excluded_mci_lines, 12L)     # p001-p012, 12 >= 11: shown
  expect_equal(ou$n_cells_negative_total, 0L)
  expect_false(any(c("min", "max", "min_claim_dt", "max_claim_dt") %in% names(co)))

  ## Final sample = p001-p100 only: all 19 prior reversals and all 12 outcome
  ## reversals are inside it. prior: 100 * 100 = 10,000 gross, -2,850, net 7,150;
  ## outcome: 100 * 300 = 30,000 gross, -600, net 29,400.
  cf <- tb$costs_final_sample
  expect_equal(cf$n_patients, c(100L, 100L))
  expect_equal(cf$total_cost[cf$window == "prior_cost"], 7150)
  expect_equal(cf$total_cost[cf$window == "outcome"], 29400)
  expect_equal(cf$mean_cost[cf$window == "prior_cost"], 71.5)   # 7,150 / 100

  ## Categories (outcome): Inpatient 136 lines, 27,200; Outpatient 136 + 12
  ## lines, 13,600 - 600 = 13,000; window net 27,200 + 13,000 = 40,200 so the
  ## shares are 27,200 / 40,200 = 0.67662 and 13,000 / 40,200 = 0.32338. The
  ## Managed Care Invoice lines (12 x 999 = 11,988) are shown but in no share.
  ca <- tb$cost_categories
  o <- ca[ca$window == "outcome", ]
  expect_equal(o$net_dollars[o$setting == "Inpatient"], 27200)
  expect_equal(o$net_dollars[o$setting == "Outpatient"], 13000)
  expect_equal(o$n_claims[o$setting == "Outpatient"], 148L)
  expect_equal(o$share_of_window_net[o$setting == "Inpatient"], 0.67662)
  expect_equal(o$share_of_window_net[o$setting == "Outpatient"], 0.32338)
  expect_true(o$excluded[o$setting == "Managed Care Invoice"])
  expect_equal(o$net_dollars[o$setting == "Managed Care Invoice"], 11988)
  expect_true(is.na(o$share_of_window_net[o$setting == "Managed Care Invoice"]))

  expect_equal(tb$cost_known$n_known, c(136L, 136L))
})

test_that("enrollment: complete-case tally, one arm-level partition", {
  es <- tb$enrollment.summary
  v <- stats::setNames(es$n, es$quantity)
  expect_equal(v[["n_cohort"]], 136L)
  expect_equal(v[["n_covered"]], 100L)
  expect_true(is.na(v[["n_absent"]]))            # 6 < 11
  ## 136 - 100 = 36 = insufficient (30) + absent (6): the 6 is masked, so the 30
  ## must be too, or 36 - 30 would give it back.
  expect_true(is.na(v[["n_insufficient"]]))
  bm <- tb$enrollment.by_months_covered
  expect_equal(bm$n_patients[bm$n_covered == 24L], 100L)
  expect_equal(bm$pct[bm$n_covered == 24L], round(100 * 100 / 136, 2))  # 73.53
  expect_equal(bm$window_months[1], 24L)
  ## by arm: A=0 45 / 60 = 0.75 (15 not complete); A=1 55 / 76 = 0.723684 (21).
  ba <- tb$enrollment_by_arm
  expect_equal(ba$n_complete, c(45L, 55L))
  expect_equal(ba$n_not_complete, c(15L, 21L))
  expect_equal(ba$rate, c(45 / 60, 55 / 76))
  ## the lag-by-arm cells (which had masked 10, 9, 6 and 10) are not published.
  expect_false("enrollment_by_lag_arm" %in% names(tb))
  ## the arm table adds up to the published totals.
  expect_equal(sum(ba$n_complete), v[["n_covered"]])
  expect_equal(sum(ba$n_not_complete), v[["n_cohort"]] - v[["n_covered"]])
})

test_that("a small arm-level cell is masked together with its complement", {
  small <- toy
  small$p05$by_arm <- tibble::tibble(
    A = c(0L, 1L), n = c(60L, 76L), n_complete = c(50L, 66L), n_unresolved = c(0L, 0L)
  ) |> dplyr::mutate(rate = .data$n_complete / .data$n)    # not complete: 10 and 10
  ba <- build_decision_report(small, title = "toy-small")$tables$enrollment_by_arm
  expect_true(all(is.na(ba$n_not_complete)))
  expect_true(all(is.na(ba$rate)))      # the rate would give n_complete / n back
})

test_that("a cohort-minus-final difference below 11 withholds the final-sample costs", {
  near <- toy
  near$p06$final_ids <- sprintf("p%03d", 1:130)      # 136 - 130 = 6 excluded
  near$p06$n_complete_case <- 130L
  near$p06$n_final <- 130L
  near$p06$design <- summarise_design(
    cbind(tibble::tibble(ID = sprintf("p%03d", 1:130), A = 1L, prior_cost = 0), 
          as.data.frame(matrix(0L, 130, 15, dimnames = list(NULL, design_matrix_columns())))),
    n_floored = 0L, floor_value = 0, cutpoints = c(p25 = 1, p50 = 2, p75 = 3))
  rn <- build_decision_report(near, title = "toy-near")$tables
  expect_named(rn$costs_final_sample, "note")
})

test_that("prior_cost floor, bins and cut points (hand-calculated)", {
  ## prior_cost = (row - 20) * 100, rows 1..100: rows 1-19 are negative (19
  ## floored to 0), row 20 is 0, rows 21-100 are 100..8,000. After flooring
  ## 20 patients sit at 0.
  pf <- tb$prior_cost_floor
  v <- stats::setNames(pf$n, pf$quantity)
  expect_equal(v[["n_final"]], 100L)
  expect_equal(v[["n_negative_prior_cost_floored"]], 19L)
  ## 20 at the floor - 19 floored = 1 patient with raw prior_cost exactly 0: a
  ## small cell recoverable from the two counts, so the at-floor count is hidden.
  expect_true(is.na(v[["n_prior_cost_at_floor_after_flooring"]]))
  ## Sorted floored values: x[1..20] = 0 and x[20 + m] = 100 m. Type-7 quantile
  ## position h = 1 + 99 p: p25 h = 25.75 -> x[25] + .75 (x[26] - x[25]) =
  ## 500 + .75 * 100 = 575; p50 h = 50.5 -> 3,000 + .5 * 100 = 3,050; p75 h =
  ## 75.25 -> 5,500 + .25 * 100 = 5,525.
  expect_equal(tb$prior_cost_cutpoints$cutpoint_usd, c(575, 3050, 5525))
  ## Bins (-Inf, 575], (575, 3050], (3050, 5525], (5525, Inf): Q1 = 20 zeros +
  ## 100..500 = 25; Q2 = 600..3,000 = 25; Q3 = 3,100..5,500 = 25; Q4 = 25.
  expect_equal(tb$prior_cost_bins$n, c(25L, 25L, 25L, 25L))
  expect_equal(tb$prior_cost_bins$pct, rep(25, 4))
})

test_that("design degeneracy table: prevalence shown, small splits masked", {
  dp <- tb$design_prevalence
  get <- function(nm) dp[dp$covariate == nm, ]
  ## column j is 1 iff i %% (j + 1) == 0 for i = 0..99, so n_yes =
  ## ceiling(100 / (j + 1)): HOUSING (j1) 50, FOOD 34, POVERTY 25, MH_IP 20,
  ## MH_ER 17, MED_IP 15, VAR5D_SZ_YN (j7) ceiling(100 / 8) = 13,
  ## VAR9E_ANTIPSY_YN (j8) ceiling(100 / 9) = 12, VAR9G (j9) ceiling(100 / 10) =
  ## 10, DM2 (j10) ceiling(100 / 11) = 10, HTN (j11) ceiling(100 / 12) = 9,
  ## GLU_OR_LIPID (j12) ceiling(100 / 13) = 8.
  expect_equal(get("HOUSING_YN")$n_yes, 50L)
  expect_equal(get("FOOD_YN")$n_yes, 34L)
  expect_equal(get("VAR5D_SZ_YN")$n_yes, 13L)
  expect_equal(get("VAR5D_SZ_YN")$n_no, 87L)
  expect_equal(get("VAR5D_SZ_YN")$prevalence, 0.13)
  expect_equal(get("VAR9E_ANTIPSY_YN")$n_yes, 12L)
  expect_equal(get("VAR9E_ANTIPSY_YN")$prevalence, 0.12)
  ## The last four have 10, 10, 9 and 8 patients with a 1: masked with their
  ## rates (the 0 count is derivable from the published n = 100).
  for (nm in c("VAR9G_ANTIDEP_YN", "VAR5G_DM2_YN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN")) {
    expect_true(all(is.na(unlist(get(nm)[c("n_yes", "n_no", "prevalence")]))), info = nm)
  }
  ## The quartile dummies are the Q2-Q4 counts of prior_cost_bins: published there.
  expect_false(any(grepl("^prior_cost_q", dp$covariate)))
  ## Grid: 2^15 atoms; the number of realized patterns is not asserted by hand
  ## here (see design_patterns() below); it must be between 1 and 100.
  g <- tb$design_grid
  expect_equal(g$n[1], 32768)                    # 2^15
  expect_lte(g$n[2], 100)
  expect_gte(g$n[2], 1)
  expect_equal(tb$final_arms$n, c(45L, 55L))     # A = 1 for rows 1..55
})

test_that("design_patterns() counts distinct rows and small-pattern patients", {
  ## Rows (0,1) (0,1) (1,1) (1,0) (1,0): three patterns of sizes 2, 1, 2.
  X <- data.frame(a = c(0, 0, 1, 1, 1), b = c(1, 1, 1, 0, 0))
  expect_equal(design_patterns(X, min_cell = 3L)$n_patterns, 3L)
  ## min_cell 3: all three patterns are smaller than 3 -> 2 + 1 + 2 = 5 patients.
  expect_equal(design_patterns(X, min_cell = 3L)$n_in_small_patterns, 5L)
  ## min_cell 2: only the size-1 pattern is smaller than 2 -> 1 patient.
  expect_equal(design_patterns(X, min_cell = 2L)$n_in_small_patterns, 1L)
})

test_that("mask_binary_split() hides both sides of a small split, keeps zero", {
  tb2 <- tibble::tibble(
    n_yes = c(5L, 50L, 95L, 0L), n_no = c(95L, 50L, 5L, 100L),
    prevalence = c(.05, .5, .95, 0), minority_share = c(.05, .5, .05, 0)
  )
  m <- mask_binary_split(tb2, min_cell = 11L)
  ## 5/95 and 95/5 are masked (either side below 11); 50/50 and 0/100 are not.
  expect_equal(is.na(m$n_yes), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(is.na(m$prevalence), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("the report refuses intermediates from different runs", {
  bad <- toy
  bad$p06$n_pre <- 999L
  expect_error(build_decision_report(bad), "disagree")
  bad2 <- toy
  bad2$p06$.meta$env <- "server"
  expect_error(build_decision_report(bad2), "different environments")
})

test_that("report_record() saves a list with environment metadata", {
  d <- withr::local_tempdir()
  path <- report_record("zz_test", list(x = 1:3), d)
  back <- readRDS(path)
  expect_equal(back$x, 1:3)
  expect_true(!is.null(back$.meta$env))
  expect_equal(basename(path), "zz_test.rds")
})

test_that("08 runs end to end off-server on the toy inputs and exports nothing", {
  skip_if(
    length(list.files(file.path(app_dir_for_tests, "output", "report_inputs"),
                      pattern = "[.]rds$")) > 0L,
    "real/fixture report intermediates exist locally; 08 would use or reject them"
  )
  skip_if(smidata::smi_is_server(), "never run the toy path on the server")
  root <- dirname(normalizePath(app_dir_for_tests))
  skip_if_not(basename(normalizePath(app_dir_for_tests)) == "application")
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(sprintf(
      "setwd(%s); source('application/08_decision_report.R')", shQuote(root, type = "cmd")
    ))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  expect_true(is.null(status) || status == 0L, info = paste(utils::tail(out, 15), collapse = "\n"))
  expect_true(any(grepl("Not exporting: not on the server", out, fixed = TRUE)))
  md <- file.path(app_dir_for_tests, "output", "decision_report_toy", "decision_report.md")
  expect_true(file.exists(md))
  lines <- readLines(md)
  expect_true(any(grepl("TOY INPUTS", lines, fixed = TRUE)))
  expect_true(any(grepl("^## prior_cost_floor", lines)))
})

test_that("smi_report_export() is only reachable behind smi_is_server()", {
  src <- readLines(file.path(app_dir_for_tests, "08_decision_report.R"), warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  at <- grep("smi_report_export\\(", src)
  expect_length(at, 1L)
  expect_true(any(grepl("smi_is_server\\(\\)", src[max(1L, at - 3L):at])))
})
