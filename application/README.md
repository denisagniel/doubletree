# `application/` — real-data application pipeline

<!-- NOTE (2026-10-09 audit): AAP is Adults' Access to Preventive/Ambulatory Health Services (HEDIS), not antipsychotic adherence (that is msr_saa). Framing under PI review; see smidata QUALITY_MEASURES.md. -->
Equity-blind, tree-based ATT estimation of antipsychotic-medication adherence (AAP achievement,
months `[s, s+12mo)`) on total SERVICE cost (no pharmacy; negative claim lines netted; nominal dollars; months `[s+12mo, s+24mo)`), in the NYS OMH Medicaid SMI cohort, via
`doubletree::estimate_att()` / `doubletree::estimate_att_crossfit()`.

**Read this before running anything.**

## What this is, and what it is not

This directory is a **deliberate methodological contrast** with
`~/RAND/rprojects/smi/dual-bounds/analysis/real_data/`, not an independent study. It reuses that
paper's population, exposure, outcome and covariate-source design **verbatim** — smidata's
`inst/analyses/doubletree__application.yml` marks sections 2–6 `status: inherited` — and differs
in three places (the 2026-10-15 project-wide PI rules — `s`, the anchors, service-only costs — apply to both papers alike):

| | dual-bounds | this paper |
|---|---|---|
| estimator | `marbounds` partial-identification bounds | `doubletree::estimate_att()` / `estimate_att_crossfit()` |
| outcome missingness | modelled response indicator `R` | **complete-case restriction**: full enrollment in `[s, s+24mo)`, `type_policy = "require_agreement"`, gap 0, `gap_mode = "total"` recorded (dual-bounds keeps everyone and models `R`) |
| covariates in `X` | ~68, four blocks, non-binary allowed | **15, all binary** (12 confirmed `_YN` + 3 `prior_cost` quartile dummies) |

The covariate count is the substantive consequence of the estimator choice, not a shortcut.
`estimate_att()`'s identifying assumption is **grid-exact sparsity**: a tree of at most
`leaf_budget` leaves represents *both* nuisances exactly *on the analyst's pre-specified grid*.
15 binary columns give a 2^15-atom grid; 68 would give 2^68, on which no 4-leaf tree represents
anything exactly. Selecting 3 covariates per clinical block **is** what makes the flagship
estimator's assumption plausible.

Inheriting a design is not sharing code. Nothing here reads anything from dual-bounds' repo. The
one thing genuinely shared is the **implementation**, and it is shared by calling it: the
dataset-semantics functions live in **smidata ≥ 0.2.0** as exported `smi_*()` and are called
directly. No local copy of `smi_select_measure_row()`, `smi_patient_windows()`,
`smi_sum_cost_in_windows()`, `smi_stream_cost_files()`, `smi_classify_is_mci()`,
`smi_window_coverage()` or `smi_read()` exists in this directory, and adding one would be exactly
the drift their promotion eliminated.

## Status

**Tier 0 (laptop, fixture):** reached for stages `02`–`04`, and as far as `06`'s own assembly
logic can be verified. `application/tests/` runs today with no data and no server — **165
assertions, 1 failure** (the 1 failure and 3 skips below are the local `optimaltrees.so`
security-tool block, unrelated to data/server availability — see `session_notes/2026-09-24.md`) —
and includes a call to `estimate_att()` at the real covariate shape (15 binary columns,
`outcome_type = "continuous"`, `leaf_budget = config_leaf_budget`) on a hand-built known-truth DGP, plus
`helpers/msr_filter.R`'s full fixture suite (added 2026-09-24) and the thinned complete-case
translation suite (2026-10-01 migration -- the month-by-month mechanics moved to smidata's own
104-assertion `test-enrollment.R`).

`Rscript application/run_pipeline.R` (new 2026-09-24, convenience only — every script cascades to
its own upstream dependency via `require_stage()` and is independently runnable) sources `01`–`07`
in order and now runs `01`–`04` fully against Tier-0 fixtures, then correctly **stops inside `05`'s
own limitation**: `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` column has no real-data census
levels, so `smi_fixture()` emits unparseable placeholder text and `smi_enrollment_panel()` cannot
be exercised against it (it correctly aborts if asked to try — that is the desired behaviour
against a genuinely wrong format on real data). This is an honest Tier-0 boundary, not a defect.
`06`'s own assembly logic (join, three-reason attrition accounting, structural checks, the
PROVISIONAL write) was verified separately, with hand-built inputs standing in for `02`–`05`'s
outputs (bypassing the fixture cascade at exactly the point it cannot reach) — see
`session_notes/2026-09-24.md` for what that caught. **Do not hand-write a stand-in fixture for
the pipeline itself** — the hand-built inputs used to verify `06`'s own logic are a unit-test
technique for one script's internal joins/assertions, not a substitute for `smi_fixture()`.

**Tier 1 (server, subsample + whole-data checks):** written, NOT YET RUN on the server. `90_checks_tier1.R`
specifies six checks. CHECKS 1, 2 (partly) and 4 are dual-bounds' (referenced, not re-derived).
CHECKS 3 (`_YN` coding and prevalence), 5 (cost-file coverage, which elig_coverage assumes) and 6
(msr_saa exposure convention and exact-start sizes) are this paper's own and are EXECUTABLE:
`Rscript application/90_checks_tier1.R` on the server (it aborts locally; each function was
exercised against Tier-0 fixtures only). CHECK 4 (added 2026-09-24) verifies the `enrollment_source` "ignore TYPE"
default's prerequisite at whole-table scope — `smidata::smi_enrollment_panel()` self-checks the identical
question at cohort/window-restricted scope, but that is not a substitute for CHECK 4's run.
Also worth a Tier-1 pass once on the server: `03`'s MCI-exclusion count (unverifiable at Tier 0 —
the census carries no real `SS_DESC` levels) and the coverage-truncation rate for `Y`/`prior_cost`
(under-exercised at Tier 0 — the fixture's `INDEX_DT` distribution happens to fall inside the
2016–2024 file coverage more often than real data plausibly will).

**Tier 2 (server, full cohort):** no longer blocked by an open decision (`diagnostics` resolved
2026-09-30) — but **not attempted**. `03`–`06` still need a real-data run, and `06` must be
re-run before `07` can read the un-suffixed `analytic_cohort.parquet` it now requires; the only
cohort file on the server is the `_PROVISIONAL` one written before the resolution.

```bash
Rscript application/run_tests.R          # all tests/ files (testthat); the count changes as tests are added
```

## Open decisions

Registered in `_config.R`, echoed at the top of every numbered script. Confirmed decisions are
read via `confirmed_value()` (no warning); unresolved placeholders are read via `open_value()`
(warns on every access — see `_config.R`). **No decision blocks the final artifact as of
2026-09-30** — five of the six entries are `status = "confirmed"`; `prior_cost_negative_floor`
was reopened 2026-10-15 (`status = "open_decision"`, still `blocking_final = FALSE`, read via
`open_value()` so every read warns). `blocking_final = FALSE` throughout. `assert_no_blocking_open()` in `07` is now a live regression
guard, not a known stop, and `06` writes the un-suffixed `analytic_cohort.parquet`.

| id | shared with dual-bounds | status | one line |
|---|---|---|---|
| `enrollment_source` | same table, different use | **confirmed** (2026-09-24) | `MEDICAID_FLAG==1` in `larger_smi_medicaid_monthly_flag` for every one of the 24 months in `[s, s+24mo)` (`s` = `T0`; PI 2026-10-15), no `TYPE` filter, zero gap tolerance — a PI-given, revisable default. Implemented via `smidata::smi_enrollment_panel()` + `smi_enrollment_coverage()` since the 2026-10-01 migration, translated as `type_policy = "require_agreement"` (not `"any"` — see `05_complete_case.R`'s header) and `gap_months = 0` / `gap_mode = "total"`; provisional on CHECK 4 (below) not yet running. |
| `cost_family_scope` | **yes** — identical question | **confirmed** (`service_two_families`) | The two SERVICE families (`aim3_svc`, `new_svc`) are disjoint claim sources; sum both. The four-family value is superseded: pharmacy files are not read (PI, 2026-10-15). |
| `msr_code` | **yes** — identical question | **confirmed** (`"SAA (no filter needed)"`; exposure switched AAP -> SAA 2026-10-09, SAA literal not yet seen in the real table) | `msr_saa` is expected to be an SAA-only extract — no `MSR` filter needed at the *registry* level. **Fixed 2026-09-24**: `02_population_and_eligibility.R`'s gate now reads this value (was stale since 2026-09-17). But "no filter needed" is not what `02` actually does — `smidata::smi_select_measure_row()`'s abort-on-tie does NOT guard against a other-measure row silently winning when it merely lands closer to the target month than any true row of the measure (no tie forms), so `02` filters to the literal `SAA_MSR_CODE` explicitly and reports any exclusion. See this decision's `note` for the full correction — the registry's "abort-on-tie is the safety net" claim was wrong, and is not yet fixed upstream (smidata, dual-bounds). |
| `diagnostics` | no — this paper's own | **confirmed** (`"certified_e && certified_m0"`, 2026-09-30) | Reported diagnostic flag, **not** an estimator switch: `certified_e` AND `certified_m0` both TRUE means this proxy supports grid-exact sparsity for the flagship `estimate_att()` result. **No fallback to `estimate_att_crossfit()`, silent or explicit** — the flagship stays the flagship regardless; a failed licence signals that `07`'s already-reported discrepancy deserves scrutiny, not that the crossfit number should be substituted in. No separate `gap_*` threshold — `certified_* == TRUE` already implies `gap == 0`. **Necessary, not sufficient**: `certified_*` says the search solved `eq:select` at the configured `leaf_budget`, not that the budget is large enough for the truth; the 2026-09-24 general standard (overlap/balance + nuisance fit vs OLS/GLM) applies on top. `07` **prints** whether the licence holds and never auto-switches. |
| `prior_cost_negative_floor` | no — this paper's own | **open_decision** (placeholder `0`; was confirmed 2026-09-28, reopened 2026-10-15 to match smidata's registry `covariates.prior_cost_floor`) — TODO | Behaviour left as is pending that decision: costs are netted (negative lines kept) first, THEN a negative window total is floored at 0. 2 patients had negative `prior_cost` (claim reversals exceeding pre-index payments) on the first real server run. Floored to 0 before `discretize_prior_cost()`, via `helpers/discretize_prior_cost.R::floor_negative_prior_cost()` — the raw dollar value retained in `06`'s output is unaffected; only the value fed to binning is floored. `discretize_prior_cost()` itself still aborts unchanged on a negative input. |

The two shared entries carry machine-readable `shared_with` and `registry_ref` fields instead of a
re-worded question. Two projects independently wording one question is how the two answers end up
differing.

**Deliberately absent, relative to dual-bounds' registry** — each for a stated reason, not by
oversight:

- `admin_gap_days` — superseded. There is no modelled `R` here to absorb a gap tolerance; the
  complete-case restriction is a hard inclusion criterion, so permitting a gap would make the
  retained patients' `Y` a sum over a window with holes in it. A tolerance would be a *new* PI
  decision with its own registry entry.
- `msr_tie_break` — not applicable. This pipeline calls `smi_select_measure_row()` at smidata's own
  default (`"earlier"`), and smidata aborts on any other value. No choice is being made here.
- `sl_lib` — not applicable. SuperLearner is `marbounds`' nuisance machinery. doubletree fits its
  nuisances with optimal trees via `optimaltrees`.

## Files

| file | runs today | what it does |
|---|---|---|
| `_config.R` | yes | Gates, `OPEN_DECISIONS`, `SMI_KEYS`, `SMI_COLS`, the 18 service cost dataset keys, `DOUBLETREE_WINDOWS` (with its `anchor` column), `doubletree_patient_windows()`, `config_estimator`, `config_leaf_budget`. Plain top-level variables, edited in place — not `Sys.getenv()`. |
| `helpers/covariate_blocks.R` | yes (tested) | The **PI-confirmed** 15-column selection, 4 blocks × 3 plus 3 `prior_cost` dummies, and the declared column ORDER. Authoritative, *not* a candidate — contrast dual-bounds' file of the same name. |
| `helpers/discretize_prior_cost.R` | yes (tested) | Pure: dollars → 3 quartile dummies, Q1 implicit reference, cutpoints returned as an attribute. Aborts on a collapsed boundary or an empty bin. |
| `helpers/design_matrix.R` | yes (tested) | `assemble_design_matrix()`, and `assert_binary_design_matrix()` — the four checks the package does **not** give you (see below). |
| ~~`helpers/complete_case.R`~~ | **DELETED 2026-10-01** | Promoted into smidata as `smi_parse_year_month()` / `smi_enrollment_panel()` / `smi_enrollment_coverage()`. Deleted rather than shimmed: it had been copied into dual-bounds on 2026-09-28 and the two copies had already drifted 82 lines apart. |
| `helpers/msr_filter.R` | yes (tested) | `filter_to_msr_code()` — the real `msr_code` contamination guard `02` uses (see that decision's row, above), pure and hand-fixture-tested. |
| `01_declare_requirements.R` | yes | Declares 21 dataset requirements (3 + the 18 service cost files) (added `medicaid_monthly_flag` 2026-09-24) + the 12 `_YN` columns via `smidata::smi_require()`; reports the census gap. Skips gracefully without smidata. |
| `02_population_and_eligibility.R` | yes (fixture) | **Fixed 2026-09-24** (was stale — see `msr_code`'s row above). Runs end to end locally against a Tier-0 fixture (`smidata::smi_read()` dispatches to `smi_fixture()` off-server automatically): builds `base_cohort`, filters `msr_saa` to `SAA_MSR_CODE` (`helpers/msr_filter.R`, tested), calls `smi_select_measure_row(anchor = "ceiling_month")` (a period must start exactly at `s`; no match or `MSR_DEN < 1` -> dropped), reports eligibility/exposure-arm counts. On the server it reads the real tables instead — same code path, `smidata::smi_env()` dispatches. |
| `03_cost_windows.R` | yes (fixture) | **Real as of 2026-09-24.** Cascade-guarded on `02` (`require_stage()`). Builds `windows`/`coverage`, then `smi_sum_cost_in_windows()` for real with `negative_amount_policy = "net"` over the 18 service files only (~10s against Tier-0 fixtures, verified before wiring live). `outcome_cost`/`prior_cost` carry an explicit `coverage` column and are `NA` (not coalesced to 0) whenever file-date coverage is not `"full"` — a patient with a partial-window sum is marked unknowable, not passed through as if complete; see this file's header. |
| `04_covariates.R` | yes (fixture) | **Real as of 2026-09-24.** Cascade-guarded on `03`. Reads the 12 confirmed `_YN` columns + `INDEX_DT`, joins `03`'s `prior_cost` (`[INDEX_DT-12mo, INDEX_DT)`) (no `coalesce()` — see `03`'s row), runs a PREVIEW discretize/assemble/assert pass on the pre-filter cohort. Explicitly a preview: `06` re-discretizes on the real, final (complete-case-filtered) sample. |
| `05_complete_case.R` | yes (fixture, up to a real Tier-0 boundary) | Computes complete-case status for real via `smidata::smi_enrollment_panel()` + `smi_enrollment_coverage()`, cascade-guarded on `02`. Locally, `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` has no real-data census levels, so the fixture emits unparseable placeholder text and the real computation is skipped with an explicit message (not a crash) — the hand-built two-patient demo still runs regardless, proving the function itself. On the server this limitation does not apply. Deliberately does **not** proxy enrollment from claims presence — `MEDICAID_FLAG` is a real enrollment indicator, independent of utilization; the file's header spells out why neither directional failure mode of a claims-presence proxy applies to it. |
| `06_assemble_analytic_data.R` | yes (fixture chain stops at `05`'s boundary; own logic verified separately) | Joins `02`–`05`, reports THREE separate attrition reasons (not complete-case; unknowable `Y`; unknowable `prior_cost` — see `03`'s coverage column), re-discretizes `prior_cost` on the real final sample, `assert_binary_design_matrix()` (unconditional), then writes with a `_PROVISIONAL` suffix while anything blocks, or without one once nothing does. **The PROVISIONAL release valve is now real** (fixed 2026-09-24, design-reviewed by `oracle`) — it was pure documentation before: `assert_no_blocking_open()` always aborted first, so the write was never reachable. Moved that call to `07` only (`diagnostics` gates the *estimate*, not the assembled *data*). **Found and fixed while verifying this script with hand-built inputs**: the original complete-case filter used `isTRUE(.data$complete_case)`, which is not vectorized — it would have silently zeroed the analytic sample on every real run, for a reason unrelated to complete-case status itself. Fixed to `.data$complete_case %in% TRUE`. |
| `07_estimate_att.R` | yes (needs real analytic data) | Gate (`assert_no_blocking_open()`) is the first statement after config — the only place it still runs, as of 2026-09-24 (moved out of `06`); now a live regression guard, since all decisions are confirmed. Reads the un-suffixed `analytic_cohort.parquet` (aborts by name if only the `_PROVISIONAL` file exists). Runs both estimators, prints whether the confirmed `diagnostics` licence holds (`certified_e && certified_m0`) without auto-switching on it, and reports their discrepancy. Writes `output/att_estimates.parquet`, `att_diagnostics.parquet`, `att_provenance.csv` and `att_fits.rds`. |
| `run_pipeline.R` | yes | New 2026-09-24. Convenience only (sources `01`–`07` in order) — every script is independently runnable via its own `require_stage()` cascade; this just saves typing seven commands. |
| `helpers/decision_report.R` | yes (tested) | Review-report plumbing: `report_record()` (how `02`/`03`/`05`/`06` save their counts), `summarise_design()`, `design_patterns()`, `report_levers()` (levers read from `_config.R`'s own objects), `build_decision_report()`, and `toy_report_inputs()` (hand-calculated toy intermediates). |
| `helpers/disclosure.R` | yes (tested) | Masking helpers smidata lacks (`mask_counts_primary()`, `mask_counts_partition()`, `tag_report_table()` over `smidata::smi_report_tag()`) and `tier1_saa_report()`; the identity-aware disclosure logic itself is smidata's (`smi_report()`, `smi_identity()`, `smi_gap()`, `smi_coarsen_flow()`). The project's former `hide_small_gaps()` / `coarsen_flow()` / `mask_partition_cells()` / `mask_binary_split()` / `tier1_disclosable_*()` were deleted. |
| `08_decision_report.R` | yes (toy off-server) | New 2026-10-15. Builds the review report from the intermediates `02`/`03`/`05`/`06` save to `output/report_inputs/`; see "Review report" below. Exports through `smi_report_export()` only when `smidata::smi_is_server()`. |
| `90_checks_tier1.R` | no | Server-only. Checks 1, 2 (partly) and 4 reference dual-bounds' (comment spec); checks 3, 5 and 6 are this paper's own and executable (`_YN` coding/prevalence, cost-file coverage, msr_saa convention and exact-start sizes), exporting through `smi_export()`. Aborts locally. Check 4 (added 2026-09-24) verifies the `enrollment_source` TYPE-ignore default's prerequisite. |
| `91_leaf_budget_feasibility.R` | yes (needs real analytic data) | New 2026-09-28. STANDALONE, NOT gated by `assert_no_blocking_open()` and NOT part of the 01–07 cascade — times `estimate_att()` alone at `config_leaf_budget` on the real analytic data and reports elapsed time. Reads the reportable `analytic_cohort.parquet` first, falls back to `analytic_cohort_PROVISIONAL.parquet` with an explicit warning (order reversed 2026-09-30 when `diagnostics` resolved — preferring PROVISIONAL unconditionally would silently time a stale cohort). Informational only; never a reportable estimate. Added after benchmarking `leaf_budget = 15` (raised from 4, PI 2026-09-28) on the 400-row Tier-0 toy fixture: 4.8s at `leaf_budget = 4`, 94.5s at `6`, still running past 5+ minutes at `15` (killed) — the real cohort's `07` risks discovering an infeasible budget only after `03`'s multi-hour cost-streaming pass. Invoke wrapped in the shell's `timeout` (the real bound; this script's own `setTimeLimit()` is best-effort only against optimaltrees' compiled C++ search — see the file's own header): `timeout 900 Rscript application/91_leaf_budget_feasibility.R`. |
| `tests/` | yes | testthat files, no data, no server. |
| `run_tests.R` | yes | `testthat::test_dir("application/tests")`. |

## Review report (`08_decision_report.R`)

Lets the PI check every lever and quote the counts in the paper. Each stage saves the raw counts it
already computes as it runs (`02` the attrition steps and the period selection, `03` the cost-window
object with its `setting_summary` attribute, `05` the enrollment coverage, `06` the final-sample
design counts) to `application/output/report_inputs/*.rds` (gitignored, local, never exported).
`08` reads those files, takes the levers from `_config.R`'s own objects (so the report cannot drift
from the code), and builds smidata's decision-report tables plus this paper's own.

**What the audit guarantees.** Every cell below 11 is masked. `build_decision_report()` then
declares the arithmetic between the published counts (`smidata::smi_identity()`, `smi_gap()`):
every count published in more than one table (a flow step against the cost, enrollment and
final-sample tables; eligible against its split; the final sample against the arm, bin and
prevalence tables), every total beside its breakdown (cohort = covered + insufficient + absent,
`n_yes + n_no` = final sample, the cost table's claim, negative, missing-amount and excluded lines
against the category table), and every difference of two nested counts (cohort minus final sample,
cohort minus known cost, at the prior-cost floor minus floored). `smidata::smi_report()` masks
further cells until no declared identity returns a masked value, then audits. The audit guarantees
that **no masked value is determined by the published cells through the declared identities**, and
that no declared sum of masked parts or declared difference of published counts is a count of 1..10.
It does **not** cover identities nobody declared, dollar dominance (a cost total dominated by one
patient inside a large group), or a masked cell narrowed to a range inside 1..10. The report's
header line says how many identities were audited and how many further cells the repair step
masked. Read the rendered report before releasing it; "audited" is not "safe to export". No
minimum, maximum or single-patient value is reported.

**Rules the stricter audit enforces** (smidata 15e2888). Every table must declare its count columns
(`tag_report_table()`, a wrapper over the exported `smidata::smi_report_tag()`; a derived column can name
the counts it depends on with `derived_from`), and the audit flags any *published* count of 1..10 in
one (`small_cell`), which the repair masks. A declared identity that is false on the data, or
inconsistent with the published cells, **aborts** rather than being audited; a cell-set over no rows
aborts too. A `smi_gap(big, small)` must have `big >= small` on the published cells, so the Tier-1 report
takes the direction from the raw counts. What aborts the build:

- `measure_period`: duplicate IDs, a patient absent from `index_dt` or with `NA` `INDEX_DT`, `NA`
  `matched`/`elig`, eligible but unmatched, an eligible patient with `NA` `achieved`;
- `cost_categories` / `cost_identities`: a cost table filtered or modified after
  `smi_sum_cost_in_windows()` (checksum mismatch), or one without its `negative_amount_policy`;
- an identity this report declares that is false on the unmasked counts (stage intermediates from
  different runs).

**Cost tables.** Percentiles are order statistics (type 1, an observed cell, never an interpolation),
and a tail percentile is published only when `floor(n_cells * min(p, 1 - p)) >= 11`: `p01`/`p99` need
`n_cells >= 1100`, `p05`/`p95` need 220, `p25`/`p75` need 44, `p50` needs 22. With the toy's 136 cells only
`p25`/`p50`/`p75` print. A masked negative-line count also masks `negative_dollars` and
`negative_share_of_gross`. `cost_categories` is suppressed on **distinct patients** (`n_patients`) per
(window, excluded) block, not on claim lines: one patient with 40 lines in a rare setting has their
patient count, dollars and share withheld (the line count 40 is published, being >= 11), and the
pooled complement stays inside the same block. Rows are sorted by window, excluded, setting.

**Enrollment.** `enrollment.by_months_covered` has a leading `reason` column and is split by it: an
absent patient has `n_covered = window_months`, so the full-coverage bucket holds covered and absent
patients and an undivided bucket would give `bucket - covered` = absent. `05` therefore records
`reason` beside `n_covered`; identities are per reason.

Where a count would otherwise have to be masked in several tables, it is published once or not at
all: a small attrition step is merged into the next one (`smi_coarsen_flow()`), eligible is
published in `measure_period.summary` so no flow count may sit within 10 of it (that step is merged
too), a count the flow does not publish is blanked in every other table, and the `n_matched`
column of `coverage_by_window`, the lag-by-arm enrollment cells and the quartile-dummy rows of the
prevalence table are not published separately.

```
Rscript application/run_pipeline.R          # 01-07, then 08 (or run 08 alone after 06)
Rscript application/08_decision_report.R    # Markdown -> application/output/decision_report/decision_report.md
```

On the server (`smidata::smi_is_server()` and intermediates recorded there) the tables are also
exported through `smidata::smi_report_export()` (the only door off the server). Off the server, with no
intermediates, it runs on the toy inputs, titles the report "TOY INPUTS", writes to
`output/decision_report_toy/`, and exports nothing. Stage intermediates from different runs or
environments abort the build.

| Lever / question | Where it is evidenced (report table) |
|---|---|
| Period rule (`s`, exact start, anchor), `MSR_DEN < 1` drop, achievement threshold, measure code, windows, cost policy, cost families, MCI rule, enrollment rule, floor, bins, covariate set, winsorizing / provider attribution (not applicable) | "Levers in force" (read from `_config.R`) |
| Every exclusion step, base table to analytic sample, with reasons (small steps merged into the next, text kept) | `flow` (now with a `step_no` column) |
| SAA cohort: eligible, split achieved / not achieved; by period-start month; `INDEX_DT` on the 1st vs later. Patients, exact-start matches, unmatched and `MSR_DEN < 1` are steps of `flow` (published once, merged when small) | `measure_period.summary` (`eligible`, `eligible_not_achieved`, `eligible_achieved`), `.by_period_start`, `.index_day`, `flow` |
| Cost-file coverage exclusions, per window (patients not fully covered; the matched denominator is the `flow` step) | `coverage_by_window` |
| Non-SAA `msr_saa` rows excluded; rows for IDs outside the cohort | `msr_filter` |
| Negative lines and dollars, cells with negative or zero totals, missing amounts, MCI lines excluded, distribution (p01-p99), prior vs outcome window | `costs_cohort` (analysis cohort), `costs_final_sample` (final sample) |
| Service-setting breakdown with net dollars and shares | `cost_categories` (analysis cohort; the attribute cannot be restricted) |
| Patients with a known `Y` / `prior_cost` | `cost_known` |
| Complete case over `[s, s+24mo)` and its attrition: covered, insufficient, absent; months covered; by `T0` lag and arm; by arm | `enrollment.summary`, `enrollment.by_months_covered`, `enrollment_by_arm` (the by-lag-by-arm cells are not published), and the last two `flow` steps |
| `prior_cost` floor: cells floored at 0, cells at the floor | `prior_cost_floor` |
| Quartile cut points and the count in each bin | `prior_cost_cutpoints`, `prior_cost_bins` |
| 15-covariate binary design degeneracy: prevalence of each `_YN` covariate (incl. `VAR5D_SZ_YN`, `VAR9E_ANTIPSY_YN`) in the final sample, realized patterns of the 2^15 grid | `design_prevalence`, `design_grid` |
| Final arm sizes | `final_arms` |

Not covered, and why: the costs tables describe the analysis cohort (everyone before the complete-case
filter) and the final sample, not the complete-case-only cohort in `cost_categories`, because
smidata's category table is read from an attribute of the unmodified cost object. The final-sample
cost table has no category table, so its claim, negative and missing-amount lines are linked to the
cohort table by gaps, not by a breakdown. Tests: `tests/test-decision-report.R` (hand-calculated toy
values and perturbed fixtures: a `MSR_DEN < 1` count of 5, a flow step of 5, a 6-patient
cohort-minus-final difference, a small arm cell, a 5-line service category, a 40-line one-patient
category, a small `_YN` split, an `NA` `INDEX_DT`, a duplicate ID, a filtered cost table, a published
count below 11, withheld tail percentiles, the by-reason months-covered split) and
`tests/test-disclosure.R` (the Tier-1 report).

**Tier-1 exports (`90_checks_tier1.R`).** Every table now leaves through `smi_report()` ->
`smi_report_export()`, named `doubletree_tier1__<table>` (it was `doubletree_tier1_<table>`).
The exact-start sizes, the MSR_YR offsets, the gap-by-offset table and the `_YN` prevalence table
publish related counts (an arm size is the sum of a column's yes / no / NA cells; the chain with-rows
>= matched >= eligible; the nearest-period match counts contain the exact-start ones), so they are
**one** audited report, written once by whichever of CHECK 6 / CHECK 3 runs first (the other no
longer writes them again). `yn_coding` and `file_coverage` are independent counts: primary
suppression, no declared identities, so only the per-cell guarantee applies to them.

## What the package does not check for you

Read `R/utils.R`'s `check_att_data()` before editing `helpers/design_matrix.R`. It validates
non-`NULL`ness, matching lengths, no `NA` in `X`/`A`/`Y`, that `A` is 0/1, and `Y` against
`outcome_type`. **It does not examine `X`'s coding at all.** The two entry points then diverge:

- `estimate_att()` **does** reject non-binary `X` itself, after `check_att_data()` — but only at
  estimation time, i.e. after 03's multi-hour cost aggregation and 06's write.
- `estimate_att_crossfit()` **does not.** `check_att_data()` is its only `X` validation, so a
  1/2-coded or continuous column runs to completion and returns a number, discretized internally
  by `optimaltrees` at cutpoints the analyst never declared.

Neither checks, on either path: the column count, the column **order**, or whether a column is
constant. `assert_binary_design_matrix()` checks all four, at assembly time.
`tests/test-design-matrix.R` contains a negative test that *demonstrates* the crossfit hole
rather than asserting it: it feeds a 1/2-coded column to `estimate_att_crossfit()` and confirms a
finite estimate comes back.

This is also why `prior_cost` is pre-discretized here rather than left to the estimator. Both
entry points forward `discretize_method` / `discretize_bins` to `optimaltrees`
(`R/estimate_att.R:384,402`; `R/estimate_att_crossfit.R:219`), so an un-preprocessed continuous
column is binned *adaptively, per fit*, at cutpoints that do not exist until after fitting — and
grid-exact sparsity would then be an assumption about an object that did not exist when it was
stated.

## Why both estimators are always reported

`estimate_att()` does not fail when grid-exact sparsity fails; it returns a confidently-stated
wrong number, because the assumption that licensed the plain Wald interval was false. It reports
proxies and, by explicit design, never switches estimator on the basis of them (auto-switching
would be data-dependent post-selection inference the paper does not analyse). So `07` always runs
`estimate_att_crossfit()` too and reports the **discrepancy**, turning a silent bias into a
visible disagreement. That discrepancy is a disagreement *measure*, not a hypothesis test — the
two estimates share data, and the estimators differ in more than the sparsity assumption.

## Relationship to the package (`R/`, `DESCRIPTION`, `R CMD check`)

`^application/` is in `.Rbuildignore`, so **nothing in this directory is inside `R CMD check`'s
scope**, and `application/tests/` is separate from `tests/testthat/` (the package suite) for
exactly that reason: this pipeline needs `dplyr`, `tibble`, `cli`, `tidyr`, `arrow`, `fs` and
`smidata`, none of which the package depends on, and none of which enters `DESCRIPTION`.

`DESCRIPTION` was **not** modified. `smidata` was deliberately *not* added to `Suggests`: it is an
internal, non-CRAN package, so listing it would make `R CMD check` emit a new note about an
unavailable dependency — a regression in the package's check status in exchange for documentation
that this README already provides. `_config.R` guards it with `requireNamespace()` and every
script reports its absence explicitly.

`R/` and `tests/testthat/` are untouched.

## Conventions worth knowing before editing

- **Configuration** is plain top-level R variables in `_config.R`, edited before a run.
  Environment variables are reserved for smidata's own per-machine path overrides
  (`SMI_DATA_DIR`, `ANALYSIS_ENV`, …), which must resolve identically across the five projects
  sharing the server.
- **Windows are half-open**: start inclusive, end exclusive, so no claim is counted twice. The
  `prior_cost` window ends *at* `INDEX_DT`, so a claim on `INDEX_DT` itself (the index event) is
  not in it.
- **Two anchors** (PI project-wide rules, 2026-10-15). `s` (called `T0` in the code) is the first
  day of the month on or after `INDEX_DT` (`INDEX_DT` itself when it is the 1st; `exposure_start()`
  in `helpers/period_start.R`). The `msr_saa` period must **start exactly at `s`**
  (`smidata::smi_select_measure_row(anchor = "ceiling_month")`: no nearest period, no gap
  tolerance, no ties); a patient with no such row, or `MSR_DEN < 1`, is ineligible and **dropped**.
  Achievement is `MSR_NUM / MSR_DEN >= 0.5`. Exposure `[s, s+12mo)`, outcome `[s+12mo, s+24mo)` and
  complete case `[s, s+24mo)` anchor on `s`; `prior_cost` `[INDEX_DT-12mo, INDEX_DT)` and the 12 `_YN`
  covariates anchor on `INDEX_DT` itself. The prior-cost and exposure windows are therefore not
  contiguous when `INDEX_DT` is not the 1st (the days up to `s` belong to neither).
  `doubletree_patient_windows()` (`_config.R`) builds both anchors; do not pass
  `DOUBLETREE_WINDOWS` to `smi_patient_windows()` directly. The earlier rule (`start > INDEX_DT`,
  offset 1, a gap cap `config_saa_max_month_gap`, `prior_cost` on `[T0-12mo, T0)`) is superseded
  and its code deleted.
- **Costs are SERVICE only** (`aim3_svc_cost_*`, `new_svc_cost_*`; 18 files), no pharmacy. Negative
  `AMOUNT_PAID` lines are kept (`negative_amount_policy = "net"`); Managed Care Invoice lines are
  excluded; nominal dollars.
- **The quartile cutpoints are sample-dependent and must be reported.** The confirmed design says
  quartiles are computed "within the analytic sample", which is the *post*-restriction sample —
  hence `discretize_prior_cost()` is called in `06`, after the complete-case filter, not only in
  `04`. Discretizing before the filter and then subsetting gives different boundaries.
- **No silent fallbacks.** No `%||%` substituting for a missing data value, no
  `tryCatch(→ NULL)`, no `try(silent = TRUE)`, no `possibly()`/`safely()`. Where the preferred
  behaviour is not implementable, the code aborts and names the gap.
- **The release valve is `_PROVISIONAL`, never a relaxed gate.** If a number is wanted before the
  decisions land, `06` names a `_PROVISIONAL`-suffixed artifact. A provisional file is
  self-labelling; a relaxed gate is invisible in the output.
- **Run time.** The cost-claims family is 171.5 GB across 36 files (service + pharmacy) — 91.8% of
  every byte in the secure tree; this pipeline reads only the 18 service files. `smi_sum_cost_in_windows()` streams one file at a time for that reason. Budget
  hours.

## Provenance

The confirmed design is recorded in
`~/RAND/tools/smidata/inst/analyses/doubletree__application.yml`; that file, not this directory,
is authoritative. Column names and row counts come from smidata snapshot `2026-09-15_2824080e`.
The 15-column selection and `prior_cost`'s quartile discretization were confirmed by the PI on
2026-09-17.

**Shared dataset semantics** (MSR_YR format, `msr_aap` as an access-to-care measure (the earlier schizophrenia-subgroup scoping was retracted 2026-10-09; see `~/RAND/tools/smidata/inst/docs/dictionary/QUALITY_MEASURES.md`),
SMI cohort ICD-10 definitions, cost-claims family structure) now live in
`~/RAND/tools/smidata/inst/docs/DATA_DICTIONARY.md` and its `dictionary/` family files, added
2026-09-28 — distinct from this file's own registry-recorded *decisions*.
