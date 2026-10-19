Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile

── Open decisions (unconfirmed -- see application/_config.R) ──

• enrollment_source [confirmed, SHARED WITH dual-bounds]: Table/columns
providing continuous-enrollment spans through s+24mo (s = first of the month
                                                      on or after INDEX_DT), needed to determine COMPLETE-CASE STATUS (whether Y is
                                                                                                                       observed for the whole outcome window). NO table in the confirmed design
carries enrollment spans reaching that far: larger_smi_covariates is a
covariate snapshot (72 columns, no span fields), msr_aap is measure periods
(not enrollment), and the sibling aim1_smi_charac -- not elected by this
paper -- encodes only ~12mo post-index enrollment.
• cost_family_scope [confirmed, SHARED WITH dual-bounds]: SHARED WITH
dual-bounds -- see registry_ref. Whether all four cost-claims families are
disjoint claim sources or aim3_*/new_* are overlapping extracts of the same
claims (which would double-count cost if all four are summed).
• msr_code [confirmed, SHARED WITH dual-bounds]: SHARED WITH dual-bounds -- see
registry_ref. Which msr_saa$MSR code value identifies the SAA measure
specifically (msr_saa is not filtered to one measure, so 'the single
  nearest-month row' is ambiguous without it).
• coverage_eligibility [confirmed, SHARED WITH dual-bounds]: SHARED WITH
dual-bounds. Whether a patient whose study period (either DOUBLETREE_WINDOWS
                                                   window) extends beyond COST_SOURCES' file-date coverage is excluded from the
  analysis cohort.
• diagnostics [confirmed]: Which diagnostics gate a reportable ATT.
  estimate_att() returns sparsity-PROXY diagnostics (certified_e, certified_m0,
  used_search_e, used_search_m0, n_leaves_e, n_leaves_m0, gap_e, gap_m0) and
  never acts on them; which of them, at which thresholds, licenses reporting
  the flagship estimate rather than the cross-fit fallback is unspecified.
• prior_cost_negative_floor [open_decision]: How to treat patients whose
  pre-index prior_cost (cost summed over [INDEX_DT-12mo, INDEX_DT), NETTING
  negative AMOUNT_PAID lines) is NEGATIVE -- claim reversals in the window
  exceeding payments. discretize_prior_cost() refuses to bin a negative value
  by design (its own note: 'a decision, not a binning detail') -- net it to
  zero, drop the patient, or allow the signed value into the quantile grid are
  three genuinely different choices with different consequences for the
  covariate grid.

── 07 -- tree-based ATT estimation ─────────────────────────────────────────────
ℹ n = 51243; A = 1: 30962; A = 0: 20281;
Y in [0, 9140156]

── estimate_att() -- flagship, requires grid-exact sparsity ──

ℹ theta = -5764  (SE 1731)
ℹ 95% CI: [-9157, -2371]
ℹ certified_e = TRUE, certified_m0 = FALSE;
leaves: e = 8, m0 = 15 (budget
16)
ℹ used_search_e = FALSE, used_search_m0 =
FALSE; gaps: 0,
NA
ℹ any_truncated_e = FALSE, any_truncated_m0 =
FALSE; clip_rate_e = 0
ℹ Flagship licensed under the confirmed diagnostics rule
("certified_e && certified_m0"): FALSE
! Licence FAILS -- certified_e and certified_m0 do not both hold, so this
proxy does not support grid-exact sparsity for the flagship estimate.
There is NO fallback to estimate_att_crossfit(), silent or explicit:
estimate_att() remains the flagship (per config_estimator below) and
estimate_att_crossfit() remains its companion, unconditionally. A failed
licence is a signal to scrutinize the discrepancy reported below, not an
instruction to report estimate_att_crossfit()'s number in its place.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.

── estimate_att_crossfit() -- companion, no sparsity requirement ──

ℹ theta = -5319  (SE 1926)
ℹ 95% CI: [-9093, -1545]

── Discrepancy ──

ℹ estimate_att() - estimate_att_crossfit() = -444.9
(-0.17 pooled SEs)
! This is a DISAGREEMENT MEASURE, not a hypothesis test. The two estimates share
the same data, so pooled_se overstates the standard error of their difference;
and the estimators differ in more than the sparsity assumption (in-sample
                                                                versus cross-fit nuisances, different regularization). A large value says
'investigate', never 'sparsity is rejected'.
ℹ Flagship per `config_estimator = "att"`: theta =
  -5764
✔ Wrote the 07 result files to application/output.
