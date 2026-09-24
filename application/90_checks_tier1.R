## ============================================================================
## 90_checks_tier1.R
##
## STUB. No executable Tier-1 logic -- every check here requires the real files
## on the secure server, which this repo has no access to. The comment blocks
## specify each check precisely enough to be typed out on the server without
## re-deriving it.
##
## "Tier 1" is smidata's middle tier: run on the server against a real
## subsample, PLUS the whole-data checks a subsample cannot answer (row counts,
## key overlap, totals).
##
## TWO OF THE THREE CHECKS BELOW ARE dual-bounds', REFERENCED RATHER THAN COPIED.
## This paper inherits dual-bounds' outcome and exposure design verbatim, so it
## inherits the same two risks, on the same files, with the same resolving
## procedure. Re-deriving them here would produce two specifications of one check
## that can drift apart and be answered differently -- which is exactly the
## failure mode the shared OPEN_DECISIONS entries (with their `shared_with` and
## `registry_ref` fields) exist to prevent. Run them ONCE, on the server, and
## record the answer in smidata's registry where both papers read it.
## ============================================================================


## ---------------------------------------------------------------------------
## CHECK 1 -- resolve OPEN_DECISIONS$cost_family_scope    [SHARED, REFERENCED]
##
## SPECIFICATION: ~/RAND/rprojects/smi/dual-bounds/analysis/real_data/
##                90_checks_tier1.R, CHECK 1.
## REGISTRY:      smidata::inst/analyses/dual-bounds__application.yml#cost_family_scope
##
## Do not re-derive it. That file specifies the whole procedure: pick one shared
## suffix year, read both aim3_* and new_* with (ID, SRV_DT, AMOUNT_PAID, SEQ_ID),
## then compare (a) row counts, (b) distinct claim keys AND their overlap, which
## is the decisive number, (c) total dollars including the de-duplicated union as
## a consistency check on (b), and (d) whether SEQ_ID is a real cross-file claim
## identifier or merely a within-file row sequence -- if the latter, (b) is
## meaningless and the fallback match key is specified there too.
##
## WHY THIS PAPER CARES *MORE* THAN dual-bounds DOES, and the one thing to add
## when running it: double counting here corrupts the COVARIATE GRID, not just
## the outcome scale. prior_cost's quartile cutpoints are computed within the
## analytic sample (06_assemble_analytic_data.R), so inflated dollars move the
## boundaries, which changes which patients are in which quartile, which changes
## the 2^15 grid that estimate_att()'s grid-exact sparsity assumption is stated
## over. An outcome inflated by a constant factor rescales an ATT; a covariate
## whose bins are drawn on inflated values is a different covariate.
##
## SO: when running dual-bounds' CHECK 1, additionally report the QUARTILE
## CUTPOINTS of pre-index total cost under each candidate scope (all four
## families summed, versus aim3_* only, versus the de-duplicated union) on the
## AAP-eligible cohort. If the three cutpoint vectors agree to within rounding,
## this paper is insensitive to the resolution and can proceed while dual-bounds
## waits. If they disagree, they disagree, and 04/06 cannot run until the scope
## is settled.
##
## SUPPRESSION. Counts and totals over millions of records -- no small cells.
## Cutpoints are quantiles of a large sample, likewise. Still route anything
## leaving the server through smidata::smi_export().
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 2 -- msr_code ambiguity, and the AAP gap distribution
##                                                        [SHARED, REFERENCED]
##
## SPECIFICATION: ~/RAND/rprojects/smi/dual-bounds/analysis/real_data/
##                90_checks_tier1.R, CHECK 2.
## REGISTRY:      smidata::inst/analyses/dual-bounds__application.yml#msr_code
##
## Do not re-derive it. That file specifies: the distribution of month_gap among
## matched patients, the tie frequency, the unmatched count, and cohort size at
## candidate gap cutoffs (6, 12, 18 months) so the eligibility/precision
## trade-off is a PI decision made on numbers.
##
## Its output resolves OPEN_DECISIONS$msr_code for BOTH papers simultaneously,
## because exposure A here is `status: inherited` from dual-bounds -- the same
## smidata::smi_select_aap_row() call on the same table with the same threshold.
##
## ONE ADDITION FOR THIS PAPER. Report the resulting A = 1 / A = 0 counts at each
## candidate gap cutoff, not just total cohort size. estimate_att() fits its mu
## tree on CONTROLS ALONE and requires at least leaf_budget * m_n control units
## (R/estimate_att.R's n_control check); a cutoff that keeps a large cohort but
## thins one arm can fail the estimator where it would not fail marbounds.
##
## SUPPRESSION. The per-gap patient counts WILL have small cells in the tail.
## Suppress n < 11 (CMS/HIPAA Safe Harbor) via smidata::smi_suppress() before
## export. Do not collapse the tail by hand, and do not widen buckets to clear 11
## without saying so.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## CHECK 3 -- are the 12 confirmed _YN columns genuinely {0, 1}?
##                                                      [THIS PAPER'S OWN]
##
## QUESTION. The PI-confirmed covariate selection (helpers/covariate_blocks.R)
## rests on the claim that these 12 columns are ALREADY BINARY and need no
## construction -- which is what makes a 15-column grid-exact design feasible at
## all. The fingerprint confirms the NAMES and the `_YN` suffix. It does not
## confirm the CODING. A SAS `_YN` flag is, in the wild, any of: 0/1 numeric,
## 1/2 numeric, "Y"/"N" character, "1"/"0" character, or TRUE/FALSE.
##
## WHY THIS IS NOT A COSMETIC CHECK. helpers/design_matrix.R's
## assert_binary_design_matrix() will catch a mis-coded column -- but the package
## itself will not, on the path that matters: check_att_data() (R/utils.R) does
## not examine X's coding at all, and estimate_att_crossfit() adds no further
## check, so a 1/2-coded column runs to completion and returns a number,
## discretized internally by optimaltrees at cutpoints nobody declared.
## estimate_att() does reject it, but only after the analytic data has been built
## -- i.e. after the multi-hour cost aggregation in 03. Knowing the coding BEFORE
## that run is the difference between a five-minute fix and a repeated run.
##
## THE CHECK. Read larger_smi_covariates with ID plus the 12 columns
## (unlist(covariate_blocks())). For each column report:
##
##   (a) CLASS and the full set of DISTINCT VALUES, including NA as its own
##       level. With 12 flags on 381,018 rows this is a 12-row table and is
##       cheap. Do NOT summarise to "n distinct = 2": 1/2 has two distinct
##       values and is wrong.
##
##   (b) MISSINGNESS RATE per column. An NA in a _YN flag is not automatically
##       "no": for a Z-code social-determinant flag (HOUSING_YN, FOOD_YN,
##       POVERTY_YN) it plausibly means "not screened", which is a different
##       thing from "screened, negative" and is itself correlated with
##       utilization. assert_binary_design_matrix() rejects any NA rather than
##       imputing, so a non-trivial rate on any of the 12 is a NEW PI decision
##       (recode, or replace that covariate) and would need its own registry
##       entry.
##
##   (c) PREVALENCE (mean, once coding is confirmed 0/1) per column, on the
##       AAP-eligible cohort rather than the full table. Two reasons. A flag at
##       0.1% prevalence is very nearly a constant column: it will not be
##       rejected by the constant-column check, but it spends a grid coordinate
##       to distinguish a handful of patients, which is a poor use of a scarce
##       budget. And prevalence within the A = 0 arm specifically bounds what the
##       mu tree can split on, since that tree is fit on controls alone.
##
##   (d) EXPOSURE-LEAKAGE CHECK on VAR9E_ANTIPSY_YN. This column is
##       antipsychotic use and the exposure is antipsychotic ADHERENCE. It is
##       included as a PRE-index confounder, which is correct only if its
##       measurement window precedes INDEX_DT. Confirm against the table's own
##       PRE_START / INDEX_DT columns that it does. Then, as a sanity check,
##       cross-tabulate VAR9E_ANTIPSY_YN against A: near-perfect agreement would
##       indicate the column is measured over a window overlapping the exposure
##       window, in which case it is not a confounder and must be dropped.
##       (VAR5G_AP_YN, another AP-related flag, is NOT in this paper's selection
##       -- but if a replacement covariate is ever needed, it carries the same
##       risk and needs the same check.)
##
## SUPPRESSION. (a) is a value-set listing, no cells. (b) and (c) are rates over
## a large cohort, but (c) crossed with A gives 2x2 tables per column and (d) is
## explicitly a cross-tabulation -- these CAN contain cells under 11 for a rare
## flag. Route through smidata::smi_suppress() before export.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## CHECK 4 -- resolve whether larger_smi_medicaid_monthly_flag$TYPE can be
## safely ignored (OPEN_DECISIONS$enrollment_source's IP/OP sub-question)
##   [SHARED, REFERENCED -- identical check to dual-bounds' CHECK 3]
##
## QUESTION. PI (2026-09-24) deferred the IP-vs-OP choice and directed:
## proceed WITHOUT filtering on TYPE for now -- every
## larger_smi_medicaid_monthly_flag row counts as enrollment evidence
## (toward complete-case status) regardless of TYPE. That is only safe if
## TYPE does not fragment a single patient's own enrollment signal in a way
## that changes the answer: if a patient's IP and OP rows for the same
## (ID, YEAR_MONTH) never DISAGREE on MEDICAID_FLAG, ignoring TYPE is a
## no-op; if they can disagree, "ignore TYPE" silently picks whichever row
## happens to be read first/last, which is not a decision anyone made.
##
## THE CHECK. Identical to dual-bounds' CHECK 3 -- run once, applies to both
## papers (same table, same question). On the server, group
## larger_smi_medicaid_monthly_flag by (ID, YEAR_MONTH) and compute:
##
##   (a) HOW MANY PATIENTS carry BOTH TYPE values at all (any row with
##       TYPE=='IP' AND any row with TYPE=='OP' for the same ID, anywhere in
##       their history) vs. patients with only one TYPE ever.
##
##   (b) WITHIN an (ID, YEAR_MONTH) cell that has both TYPE=='IP' and
##       TYPE=='OP' rows, do MEDICAID_FLAG values ever DISAGREE (one row 1,
##       the other 0)? Tabulate agree/disagree counts. Disagreement means
##       "ignore TYPE" is not neutral, and the deferred IP/OP choice will
##       change results, not just documentation.
##
##   (c) If (b) finds disagreement, report its size (patient-months
##       affected) so the PI can see whether the deferred choice is a
##       rounding error or a real modeling decision before
##       OPEN_DECISIONS$enrollment_source's interim "ignore TYPE" default is
##       relied on for a real number.
##
## SUPPRESSION. (a) is a per-category patient count and (c) is a
## patient-month count -- suppress any cell with n < 11 (CMS/HIPAA Safe
## Harbor) via smidata::smi_suppress() before export.
## ---------------------------------------------------------------------------
