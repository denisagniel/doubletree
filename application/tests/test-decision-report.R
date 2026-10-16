## ============================================================================
## application/tests/test-decision-report.R
##
## The review report (08_decision_report.R + helpers/decision_report.R) on the
## toy intermediates from toy_report_inputs(). Every expected value is
## HAND-CALCULATED in a comment; none is copied from a run.
##
## DISCLOSURE. build_decision_report() hands smidata::smi_report() the identities
## it declares between published counts; smi_report() masks further cells until
## none returns a masked value, then audits. These tests assert (1) the audit
## passes on the toy and on every perturbed fixture, (2) the small value in the
## fixture is masked or not published at all, and (3) it is not the difference of
## two published counts, which is checked here by explicit subtraction as well as
## by the audit. The audit covers only DECLARED identities; (3) is a second,
## independent look at the cases the fixtures were built around.
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

## ---- helpers for the disclosure assertions ------------------------------------

audit_rows <- function(rep) nrow(smidata::smi_report_audit(rep, abort = FALSE))

## Every published (non-NA) value in the count columns of every table.
count_values <- function(rep) {
  unlist(lapply(rep$tables, function(tb) {
    cols <- intersect(attr(tb, "smi_n_cols"), names(tb))
    unlist(lapply(cols, function(cl) tb[[cl]][!is.na(tb[[cl]])]), use.names = FALSE)
  }), use.names = FALSE)
}

## Every published copy of the analysis-cohort count (the 02 -> 03 cohort) outside
## the flow.
cohort_copies <- function(rep) {
  t <- rep$tables
  c(t$cost_known$n_cohort, t$cost_known$n_known, t$costs_cohort$n_patients,
    t$costs_cohort$n_cells, t$enrollment.summary$n[t$enrollment.summary$quantity == "n_cohort"])
}

## Is any difference of two of `x` a count in 1..10?
has_small_difference <- function(x, lo = 1L, hi = 10L) {
  x <- sort(unique(x))
  d <- outer(x, x, "-")
  any(d >= lo & d <= hi)
}

test_that("the toy report passes smidata's audit, with the repair step's masks hand-counted", {
  expect_true(rep_toy$audited)
  expect_equal(audit_rows(rep_toy), 0L)
  expect_gt(length(rep_toy$identities), 0L)
  ## Cells masked by the repair step (beyond the tables' own primary masks):
  ## design_prevalence: n_yes = 10, 10, 9, 8 (HTN, DM2, ANTIDEP, GLU) are small and
  ## masked first; each is declared against n_final = 100 (n_yes + n_no = 100), so
  ## its n_no (90, 90, 91, 92) is masked too: 4 cells. prior_cost_floor: 20 patients
  ## at the floor - 19 floored = 1 is a small gap; smi_gap hides the inner (smaller)
  ## end, the 19: 1 cell. enrollment.summary: absent = 6 is masked (primary);
  ## n_cohort = covered + insufficient + absent is declared, so the smallest other
  ## positive part (insufficient, 30; covered is 100) is masked: 1 cell, together with
  ## its declared copy, the insufficient months-covered bucket (one repair step). The
  ## absent bucket (6 patients at 24 months) is a primary mask. Total 4 + 1 + 1 = 6.
  expect_equal(rep_toy$resolved, 6L)
  md <- smidata::smi_report_md(rep_toy)
  expect_true(any(grepl("Disclosure audit: passed against", md, fixed = TRUE)))
})

test_that("the attrition table merges small steps instead of masking drops", {
  fl <- tb$flow
  ## 200 -> 180 (20 with no SAA rows) -> 160 (20 with rows but none at s) -> 150
  ## (10 with MSR_DEN < 1) -> 136 (14 outside the cost files) -> 100 (36 not
  ## complete-case) -> 100 (asserted, none dropped).
  ## 160 - 150 = 10 is a small step, so smi_coarsen_flow() merges the 150 row into
  ## the next. And eligible (150) is published in measure_period.summary, so no
  ## flow count may sit within 10 of it: 160 - 150 = 10 would be the MSR_DEN < 1
  ## count. The 160 row is merged into the next step as well: 180 -> 136 = 44.
  expect_equal(fl$step_no, 1:5)
  expect_equal(fl$n, c(200L, 180L, 136L, 100L, 100L))
  expect_equal(fl$n_dropped, c(0L, 20L, 44L, 36L, 0L))
  expect_equal(fl$pct_of_start, c(100, 90, 68, 50, 50))     # 180/200, 136/200, 100/200
  expect_false(any(c(150L, 160L) %in% fl$n))
  expect_match(fl$step[3], "a period starts exactly at s", fixed = TRUE)
  expect_match(fl$step[3], "MSR_DEN >= 1", fixed = TRUE)
  expect_match(fl$step[3], "cost-file coverage", fixed = TRUE)
  expect_match(fl$reason[3], "MSR_DEN < 1", fixed = TRUE)   # nothing leaves the account
  expect_match(fl$step[4], "[s, s+24mo)", fixed = TRUE)
  ## no published count differs from another (flow n or eligible) by 1..10.
  expect_false(has_small_difference(c(fl$n, tb$measure_period.summary$n[1])))   # eligible
})

test_that("exposure-period counts: only the eligible split is published; the flow carries the rest", {
  mp <- tb$measure_period.summary
  v <- stats::setNames(mp$n, mp$quantity)
  ## Eligible p001-p150 (MSR_DEN >= 1) = 150. Achieved p001-p060 and p101-p130 =
  ## 60 + 30 = 90, so not achieved = 150 - 90 = 60.
  expect_equal(mp$quantity, c("eligible", "eligible_not_achieved", "eligible_achieved"))
  expect_equal(v[["eligible"]], 150L)
  expect_equal(v[["eligible_not_achieved"]], 60L)
  expect_equal(v[["eligible_achieved"]], 90L)
  ## patients (200), not matched (40) and MSR_DEN < 1 (10) are NOT republished
  ## here: the flow has them, and 10 is small.
  expect_false(any(c("patients", "not_matched", "matched_den_below_1") %in% mp$quantity))
  ## By period-start month: 2020-03 is p001-p100 (all eligible, 60 achieved, 40
  ## not); 2020-04 is p101-p150 (50 eligible, 30 achieved: p101-p130; 20 not).
  bp <- tb$measure_period.by_period_start
  expect_equal(bp$period_start, c("2020-03", "2020-04"))
  expect_equal(bp$n_eligible, c(100L, 50L))
  expect_equal(bp$n_achieved, c(60L, 30L))
  expect_equal(bp$n_not_achieved, c(40L, 20L))
  expect_equal(bp$achieved_rate, c(0.6, 0.6))
  ## INDEX_DT on the 1st: p001-p100 (100 patients, all 100 eligible); later in the
  ## month: p101-p200 (100 patients, p101-p150 = 50 eligible).
  expect_equal(tb$measure_period.index_day$n_patients, c(100L, 100L))
  expect_equal(tb$measure_period.index_day$n_eligible, c(100L, 50L))
})

test_that("a DEN<1 count of 5 is neither published nor the difference of two published counts", {
  ## p001-p155 eligible: 160 matched - 155 = 5 with MSR_DEN < 1. The flow's
  ## 160 -> 155 step is 5 and the 155 -> 136 coverage step is 19.
  small <- toy
  small$p02$selection$elig <- seq_len(200) <= 155L
  small$p02$selection$MSR_DEN <- as.numeric(small$p02$selection$elig)
  small$p02$selection$achieved <- dplyr::if_else(small$p02$selection$elig,
                                                 seq_len(200) %in% c(1:60, 101:130), NA)
  small$p02$flow$n[4] <- 155L
  r <- build_decision_report(small, title = "toy-den5")
  expect_equal(audit_rows(r), 0L)
  fl <- r$tables$flow
  ## chain 200, 180, 160, 155, 136, 100, 100: 160 - 155 = 5 hides 155; eligible
  ## (155) is published, so 160 (within 10 of it) is merged away too.
  expect_equal(fl$n, c(200L, 180L, 136L, 100L, 100L))
  mp <- stats::setNames(r$tables$measure_period.summary$n, r$tables$measure_period.summary$quantity)
  expect_equal(mp[["eligible"]], 155L)
  expect_equal(mp[["eligible_not_achieved"]], 65L)          # 155 - 90
  expect_false(has_small_difference(c(fl$n, mp[["eligible"]])))
  expect_false(5L %in% count_values(r))
  expect_false(160L %in% count_values(r))
})

test_that("a flow drop of 5 between the first two steps is merged, not masked", {
  small <- toy
  small$p02$flow$n[2] <- 195L          # 200 -> 195: 5 patients without an SAA row
  r <- build_decision_report(small, title = "toy-flow5")
  expect_equal(audit_rows(r), 0L)
  ## chain 200, 195, 160, 150, 136, 100, 100: 200 - 195 = 5 hides 195; 160 - 150 =
  ## 10 hides 150; eligible = 150 is published, so 160 goes too: 200, 136, 100, 100.
  expect_equal(r$tables$flow$n, c(200L, 136L, 100L, 100L))
  expect_equal(r$tables$flow$n_dropped, c(0L, 64L, 36L, 0L))
  expect_false(has_small_difference(c(r$tables$flow$n, 150L)))
})

test_that("with a large DEN<1 count the flow keeps its steps", {
  big <- toy
  big$p02$selection$elig <- seq_len(200) <= 140L   # 160 - 140 = 20 >= 11
  big$p02$selection$MSR_DEN <- as.numeric(big$p02$selection$elig)
  big$p02$selection$achieved <- dplyr::if_else(big$p02$selection$elig,
                                               seq_len(200) %in% c(1:60, 101:130), NA)
  big$p02$flow$n[4] <- 140L
  rb <- build_decision_report(big, title = "toy-big")
  ## 200, 180, 160, 140, 136, 100, 100: the 140 -> 136 coverage step is 4 (small),
  ## so 136 is hidden; eligible (140) is published and 136 is 4 below it. The
  ## coverage step 136 is therefore merged into the complete case step (100):
  ## 200, 180, 160, 140, 100, 100.
  expect_equal(audit_rows(rb), 0L)
  expect_equal(rb$tables$flow$n, c(200L, 180L, 160L, 140L, 100L, 100L))
  expect_equal(rb$tables$measure_period.summary$n, c(140L, 50L, 90L))   # 140 - 90 = 50
  expect_equal(rb$tables$measure_period.by_period_start$n_eligible, c(100L, 40L))
  ## 136 is hidden by the flow, so no other table may publish it either.
  expect_false(136L %in% cohort_copies(rb))
  expect_false(136L %in% rb$tables$flow$n)
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
  ## 19 patients net 100 - 150 = -50: 19 cells with a negative total; the other
  ## 117 cells are 100. Percentiles are order statistics (type 1): p25 is the
  ## smallest x with F(x) >= .25, and F(-50) = 19/136 = .14 < .25, so p25 = 100 (and
  ## p50, p75 = 100). A tail percentile is published only when floor(n_cells *
  ## min(p, 1 - p)) >= 11: with 136 cells p01 gives floor(1.36) = 1, p05 floor(6.8)
  ## = 6, p95 6 and p99 1, all withheld; p25 floor(34) = 34 and p50 floor(68) = 68.
  expect_equal(pr$n_cells_negative_total, 19L)
  expect_true(all(is.na(c(pr$p01, pr$p05, pr$p95, pr$p99))))
  expect_equal(c(pr$p25, pr$p50, pr$p75), c(100, 100, 100))
  expect_equal(pr$negative_share_of_gross, round(2850 / 13600, 5))
  ## outcome: 200 + 100 each = 300 (gross 136 * 300 = 40,800); p001-p012 also
  ## -50 (12 lines, -600). Net 40,200. 12 / 136 = 8.8% of cells net 250 (the other
  ## 124 net 300), so p25 = p50 = p75 = 300 (F(250) = .088 < .25).
  expect_equal(ou$n_claims, 2L * 136L + 12L)
  expect_equal(ou$negative_dollars, -600)
  expect_equal(ou$total_cost, 40200)
  expect_true(all(is.na(c(ou$p01, ou$p05, ou$p95, ou$p99))))
  expect_equal(c(ou$p25, ou$p50, ou$p75), c(300, 300, 300))
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
  ## Months covered is split by reason: an absent patient has n_covered =
  ## window_months (24), so the 24-month bucket holds covered (100) AND absent (6)
  ## patients. Undivided, bucket - covered would return the absent count.
  bm <- tb$enrollment.by_months_covered
  expect_equal(names(bm)[1], "reason")
  expect_equal(bm$reason, c("absent_from_flag_table", "covered", "insufficient_coverage"))
  expect_equal(bm$n_covered, c(24L, 24L, 20L))
  expect_true(is.na(bm$n_patients[1]))                       # 6 < 11
  expect_equal(bm$n_patients[bm$reason == "covered"], 100L)
  expect_true(is.na(bm$n_patients[3]))        # 30, masked as the complement of the 6
  expect_equal(bm$pct[bm$reason == "covered"], round(100 * 100 / 136, 2))  # 73.53
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

test_that("a small arm-level cell is masked, and its complement is not recoverable", {
  ## Consistent fixture: arms A=0 (60) and A=1 (76) with 50 and 64 complete, so 10
  ## and 12 not complete; 114 covered, 16 insufficient, 6 absent (16 + 6 = 22 =
  ## 10 + 12); months covered: 114 covered at 24, 16 insufficient at 20, 6 absent at
  ## 24. (The masked sums are 22, not 20: two masked parts summing to exactly
  ## 2 * 10 are themselves flagged by the audit.)
  small <- toy
  small$p05$by_arm <- tibble::tibble(
    A = c(0L, 1L), n = c(60L, 76L), n_complete = c(50L, 64L), n_unresolved = c(0L, 0L)
  ) |> dplyr::mutate(rate = .data$n_complete / .data$n)
  small$p05$coverage <- toy_coverage(114L, 16L, 6L)
  small$p06$n_complete_case <- 114L
  r <- build_decision_report(small, title = "toy-small")
  expect_equal(audit_rows(r), 0L)
  ba <- r$tables$enrollment_by_arm
  ## A=0 n_not_complete = 10 < 11 is masked, and the rate (which would return
  ## n_complete / n) with it; the pooled complement masks A=1's 12 (else
  ## 22 - 12 = 10, or 76 - 64 = 12 minus a published total, would return it).
  ## Their sum, 22 = 136 - 114, is not a small count.
  expect_true(all(is.na(ba$n_not_complete)))
  expect_true(all(is.na(ba$rate)))
  expect_equal(ba$n_complete, c(50L, 64L))
  ## absent = 6: masked together with insufficient (16, the smallest other part)
  ## so that 136 - 114 - 16 does not return it.
  es <- stats::setNames(r$tables$enrollment.summary$n, r$tables$enrollment.summary$quantity)
  expect_true(is.na(es[["n_absent"]]))
  expect_true(is.na(es[["n_insufficient"]]))
  expect_equal(es[["n_covered"]], 114L)
  expect_false(6L %in% count_values(r))
  ## The 24-month bucket is split by reason: 114 covered, absent masked.
  bm <- r$tables$enrollment.by_months_covered
  expect_equal(bm$n_patients[bm$reason == "covered"], 114L)
  expect_true(is.na(bm$n_patients[bm$reason == "absent_from_flag_table"]))
})

test_that("a cohort-minus-final difference below 11 hides the cohort count everywhere", {
  ## 136 -> 130 (6 excluded by the complete-case restriction). Consistent fixture:
  ## 130 covered, 6 insufficient, none absent; arms 56 + 74 = 130 complete.
  near <- toy
  near$p06$final_ids <- sprintf("p%03d", 1:130)
  near$p06$n_complete_case <- 130L
  near$p06$n_final <- 130L
  near$p06$design <- summarise_design(
    cbind(tibble::tibble(ID = sprintf("p%03d", 1:130), A = 1L, prior_cost = 0),
          as.data.frame(matrix(0L, 130, 15, dimnames = list(NULL, design_matrix_columns())))),
    n_floored = 0L, floor_value = 0, cutpoints = c(p25 = 1, p50 = 2, p75 = 3))
  near$p05$coverage <- toy_coverage(130L, 6L, 0L)
  near$p05$by_arm <- tibble::tibble(
    A = c(0L, 1L), n = c(60L, 76L), n_complete = c(56L, 74L), n_unresolved = c(0L, 0L)
  ) |> dplyr::mutate(rate = .data$n_complete / .data$n)
  rn <- build_decision_report(near, title = "toy-near")
  expect_equal(audit_rows(rn), 0L)
  ## Chain 200, 180, 160, 150, 136, 130, 130: 160 - 150 = 10 and 136 - 130 = 6 are
  ## small steps, and eligible (150) is published, so the flow publishes only
  ## 200, 180 and 130 (the final sample).
  expect_equal(rn$tables$flow$n, c(200L, 180L, 130L))
  ## 136 is not in the flow, so no other table may carry it (cost_known, the cost
  ## tables, the coverage tally): 136 - 130 = 6 would return the excluded count.
  expect_true(all(is.na(cohort_copies(rn))))
  expect_true(130L %in% count_values(rn))
  ## The final-sample costs are published (the old report withheld them).
  expect_equal(rn$tables$costs_final_sample$n_patients, c(130L, 130L))
  expect_false(has_small_difference(c(rn$tables$flow$n, 150L)))
  ## claim lines: the six excluded patients have six lines (prior_cost window:
  ## cohort 136 + 19 = 155 vs final 130 + 19 = 149, a difference of 6). The
  ## declared gap hides the inner (final-sample) count; the cohort's 155 stays.
  co <- rn$tables$costs_cohort
  cf <- rn$tables$costs_final_sample
  expect_true(is.na(cf$n_claims[cf$window == "prior_cost"]))
  expect_equal(co$n_claims[co$window == "prior_cost"], 155L)
  ## outcome: 2 * 136 + 12 = 284 vs 2 * 130 + 12 = 272, a difference of 12: both shown.
  expect_equal(c(co$n_claims[co$window == "outcome"], cf$n_claims[cf$window == "outcome"]),
               c(284L, 272L))
})

test_that("a service category with 5 claim lines is masked and not recoverable from the window total", {
  ## 5 extra "Dental" lines in the outcome window (p001-p005, 40 each).
  extra <- tibble::tibble(ID = sprintf("p%03d", 1:5), SRV_DT = as.Date("2020-07-01"),
                          AMOUNT_PAID = 40, SS_DESC = "Dental")
  rd <- build_decision_report(toy_report_inputs(extra_outcome_claims = extra), title = "toy-dental")
  expect_equal(audit_rows(rd), 0L)
  ## Window outcome claim lines: 136 Inpatient + 148 Outpatient + 5 Dental = 289.
  ca <- rd$tables$cost_categories
  o <- ca[ca$window == "outcome" & !ca$excluded, ]
  expect_true(is.na(o$n_claims[o$setting == "Dental"]))
  ## 289 - 136 - 148 = 5 would return it, so another part of the identity is
  ## masked as well (the smallest other published count, Inpatient's 136). Five
  ## lines are also five PATIENTS (p001-p005), so the patient counts are masked too.
  expect_true(is.na(o$n_claims[o$setting == "Inpatient"]))
  expect_false(5L %in% count_values(rd))
  ## the masked parts sum to 141, not a small count.
  co <- rd$tables$costs_cohort
  expect_equal(co$n_claims[co$window == "outcome"], 289L)
})

test_that("a category with 1 patient and 40 claim lines is masked on patients, not lines", {
  ## 40 extra "Dental" lines, all p001's: 40 >= 11 lines but ONE patient. Category
  ## rows are suppressed on n_patients, so the patient count (1) is masked, with
  ## the pooled complement (Inpatient's 136, the smallest other positive), and the
  ## row's dollars (40 * 40 = 1,600) and share are withheld. The line count 40 is
  ## itself >= 11 and is published (smidata masks each count column on its own).
  extra <- tibble::tibble(ID = rep("p001", 40), SRV_DT = as.Date("2020-07-01"),
                          AMOUNT_PAID = 40, SS_DESC = "Dental")
  rd <- build_decision_report(toy_report_inputs(extra_outcome_claims = extra), title = "toy-dental1")
  expect_equal(audit_rows(rd), 0L)
  ca <- rd$tables$cost_categories
  d <- ca[ca$window == "outcome" & ca$setting == "Dental", ]
  expect_equal(nrow(d), 1L)
  expect_true(is.na(d$n_patients))
  expect_equal(d$n_claims, 40L)
  expect_true(all(is.na(d[c("gross_positive", "net_dollars", "share_of_window_net")])))
  expect_false(1600 %in% unlist(ca[c("gross_positive", "net_dollars")]))
  ## Inpatient's patients (136) is the pooled complement; its dollars are masked too.
  expect_true(is.na(ca$n_patients[ca$setting == "Inpatient" & ca$window == "outcome"]))
  expect_true(is.na(ca$net_dollars[ca$setting == "Inpatient" & ca$window == "outcome"]))
  ## Window outcome lines: 136 + 148 + 40 = 324, shown in costs_cohort.
  expect_equal(rd$tables$costs_cohort$n_claims[rd$tables$costs_cohort$window == "outcome"], 324L)
})

test_that("a published count below 11 in a declared count column is a small_cell violation", {
  bad <- rep_toy
  bad$tables$msr_filter$n[1] <- 7L        # the masked 7 (non-SAA rows) put back
  v <- smidata::smi_report_audit(bad, abort = FALSE)
  expect_true("small_cell" %in% v$kind)
  expect_true(any(v$table == "msr_filter" & v$col == "n" & v$kind == "small_cell"))
  expect_error(smidata::smi_report_audit(bad), "disclosure audit")
})

test_that("bad stage inputs abort rather than publish a difference", {
  ## NA INDEX_DT: the INDEX_DT-day split would drop the patient.
  na_dt <- toy
  na_dt$p02$index_dt$INDEX_DT[1] <- as.Date(NA)
  expect_error(build_decision_report(na_dt), "INDEX_DT")
  ## A patient absent from index_dt.
  gone <- toy
  gone$p02$index_dt <- gone$p02$index_dt[-1, ]
  expect_error(build_decision_report(gone), "INDEX_DT")
  ## Duplicate IDs in the selection.
  dup <- toy
  dup$p02$selection <- rbind(dup$p02$selection, dup$p02$selection[1, ])
  expect_error(build_decision_report(dup), "unique")
})

test_that("a filtered cost table is refused: its setting summary describes the original", {
  filt <- toy
  filt$p03$cost <- filt$p03$cost[filt$p03$cost$id %in% sprintf("p%03d", 1:100), ]
  expect_error(build_decision_report(filt), "filtered or modified")
})

test_that("blank_cell() blanks a derived column only when the blanked count is its source", {
  tb0 <- smidata::smi_report_tag(
    tibble::tibble(k = c("x", "y"), n_a = c(50L, 60L), n_b = c(70L, 80L), r = c(.5, .6)),
    n_cols = c("n_a", "n_b"), derived_cols = "r", derived_from = list(r = "n_a")
  )
  out <- blank_cell(list(t = tb0), smidata::smi_cell("t", "n_b", k = "y"))$t
  expect_equal(out$n_b, c(70L, NA))
  expect_equal(out$r, c(.5, .6))                 # r depends on n_a only
  out <- blank_cell(list(t = tb0), smidata::smi_cell("t", "n_a", k = "y"))$t
  expect_equal(out$r, c(.5, NA))
})

test_that("prior_cost floor, bins and cut points (hand-calculated)", {
  ## prior_cost = (row - 20) * 100, rows 1..100: rows 1-19 are negative (19
  ## floored to 0), row 20 is 0, rows 21-100 are 100..8,000. After flooring
  ## 20 patients sit at 0.
  pf <- tb$prior_cost_floor
  v <- stats::setNames(pf$n, pf$quantity)
  expect_equal(v[["n_final"]], 100L)
  expect_equal(v[["n_prior_cost_at_floor_after_flooring"]], 20L)
  ## 20 at the floor - 19 floored = 1 patient with raw prior_cost exactly 0: a
  ## small cell recoverable from the two counts. smi_gap hides the inner (smaller)
  ## end, the floored count.
  expect_true(is.na(v[["n_negative_prior_cost_floored"]]))
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

test_that("a small yes/no split is masked on both sides, zero is not, and the rest stay", {
  ## n_yes + n_no = n_final = 100 is published (prior_cost_floor). VAR9G (10),
  ## DM2 (10), HTN (9) and GLU (8) have a small yes count: masked, and the declared
  ## identity n_yes + n_no = n_final masks the no side too (else 100 - n_no).
  dp <- tb$design_prevalence
  small_rows <- dp$covariate %in% c("VAR9G_ANTIDEP_YN", "VAR5G_DM2_YN", "VAR5G_HTN_YN", "GLU_OR_LIPID_YN")
  expect_true(all(is.na(dp$n_yes[small_rows]) & is.na(dp$n_no[small_rows])))
  expect_false(any(is.na(dp$n_yes[!small_rows])))
  ## A column with ZERO yes (a degeneracy flag) is published, not masked.
  zero <- toy
  zero$p06$design$prevalence$n_yes[1] <- 0L
  zero$p06$design$prevalence$n_no[1] <- 100L
  rz <- build_decision_report(zero, title = "toy-zero")$tables$design_prevalence
  expect_equal(rz$n_yes[1], 0L)
  expect_equal(rz$n_no[1], 100L)
  ## A single small split of 5 (95 no) in an otherwise large table.
  five <- toy
  five$p06$design$prevalence$n_yes[1] <- 5L
  five$p06$design$prevalence$n_no[1] <- 95L
  rf <- build_decision_report(five, title = "toy-five")
  expect_equal(audit_rows(rf), 0L)
  expect_true(is.na(rf$tables$design_prevalence$n_yes[1]) && is.na(rf$tables$design_prevalence$n_no[1]))
  expect_true(is.na(rf$tables$design_prevalence$prevalence[1]))
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
  expect_true(any(grepl("Disclosure audit: passed against", lines, fixed = TRUE)))
})

test_that("the toy report exports through smi_report_export() into a temporary outbox", {
  skip_if(smidata::smi_is_server(), "never export the toy on the server")
  withr::local_envvar(c(SMI_OUTBOX = file.path(withr::local_tempdir(), "outbox")))
  recs <- suppressMessages(suppressWarnings(
    smidata::smi_report_export(rep_toy, "doubletree-toy")
  ))
  expect_equal(nrow(recs), length(rep_toy$tables))
  ## A report that fails its audit is refused at the door.
  bad <- rep_toy
  bad$tables$prior_cost_floor$n[2] <- 19L      # un-mask the floored count: 20 - 19 = 1
  expect_error(smidata::smi_report_export(bad, "doubletree-bad"), "disclosure audit")
})

test_that("smi_report_export() is only reachable behind smi_is_server()", {
  src <- readLines(file.path(app_dir_for_tests, "08_decision_report.R"), warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  at <- grep("smi_report_export\\(", src)
  expect_length(at, 1L)
  expect_true(any(grepl("smi_is_server\\(\\)", src[max(1L, at - 3L):at])))
})
