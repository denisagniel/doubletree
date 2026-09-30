Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile

── Open decisions (unconfirmed -- see application/_config.R) ──

• enrollment_source [confirmed, SHARED WITH dual-bounds]: Table/columns
providing continuous-enrollment spans through INDEX_DT+24mo, needed to
determine COMPLETE-CASE STATUS (whether Y is observed for the whole outcome
                                window). NO table in the confirmed design carries enrollment spans reaching
that far: larger_smi_covariates is a covariate snapshot (72 columns, no span
                                                         fields), msr_aap is measure periods (not enrollment), and the sibling
aim1_smi_charac -- not elected by this paper -- encodes only ~12mo post-index
enrollment.
• cost_family_scope [confirmed, SHARED WITH dual-bounds]: SHARED WITH
dual-bounds -- see registry_ref. Whether all four cost-claims families are
disjoint claim sources or aim3_*/new_* are overlapping extracts of the same
claims (which would double-count cost if all four are summed).
• msr_code [confirmed, SHARED WITH dual-bounds]: SHARED WITH dual-bounds -- see
registry_ref. Which msr_aap$MSR code value identifies the AAP measure
specifically (msr_aap is not filtered to one measure, so 'the single
  nearest-month row' is ambiguous without it).
• diagnostics [confirmed]: Which diagnostics gate a reportable ATT.
estimate_att() returns sparsity-PROXY diagnostics (certified_e, certified_m0,
                                                   used_search_e, used_search_m0, n_leaves_e, n_leaves_m0, gap_e, gap_m0) and
never acts on them; which of them, at which thresholds, licenses reporting
the flagship estimate rather than the cross-fit fallback is unspecified.
• prior_cost_negative_floor [confirmed]: How to treat patients whose pre-index
prior_cost (sum of AMOUNT_PAID over [INDEX_DT-12mo, INDEX_DT)) is NEGATIVE --
  claim reversals in the window exceeding payments. discretize_prior_cost()
refuses to bin a negative value by design (its own note: 'a decision, not a
  binning detail') -- net it to zero, drop the patient, or allow the signed
value into the quantile grid are three genuinely different choices with
different consequences for the covariate grid.

── 07 -- tree-based ATT estimation ─────────────────────────────────────────────
ℹ n = 227309; A = 1: 213150; A = 0: 14159;
Y in [0, 24064008]
bisect_lambda_to_budget: leaf_budget = 15 requires searching to depth 14 for full theoretical compliance (the worst-case chain-shaped topology) -- defaulting max_depth to 14. With 15 raw covariate(s) (more after discretization), this bisection can run up to 81 fits, each searching an unbounded-feature-count tree at depth 14 -- can be slow (see quality_reports/plans/2026-09-01_two-stage-package-defaults-session.md §1.3 for measured costs at this scale). Pass a smaller `max_depth` together with `depth_restricted = TRUE` to trade the full-class guarantee for speed.
Warning: High-dimensional data detected (estimated 465 binary features after discretization). Setting model_limit=1000000 to prevent memory issues. Override by passing model_limit=0 (unlimited) or model_limit=<value> explicitly.
leaf_budget = config_leaf_budget,


==> Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source doubletree

* installing to library 'C:/Users/DAgniel/AppData/Local/R/win-library/4.4/_build'
* installing *source* package 'doubletree' ...
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile
** help
Warning: C:/smi/projects/doubletree/man/estimate_att.Rd:233: unknown macro '\S'
*** installing help indices
** building package indices
Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile
** testing if installed package can be loaded from temporary location
Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile
** testing if installed package can be loaded from final location
Note: optimaltrees not found in ../optimaltrees
If optimaltrees is elsewhere, load it manually or update this .Rprofile
** testing if installed package keeps a record of temporary installation path
* DONE (doubletree)
