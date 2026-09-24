# `application/` — real-data application pipeline

Equity-blind, tree-based ATT estimation of antipsychotic-medication adherence (AAP achievement,
months 0–12 post-index) on total cost (months 13–24), in the NYS OMH Medicaid SMI cohort, via
`doubletree::estimate_att()` / `doubletree::estimate_att_crossfit()`.

**Read this before running anything.**

## What this is, and what it is not

This directory is a **deliberate methodological contrast** with
`~/RAND/rprojects/smi/dual-bounds/analysis/real_data/`, not an independent study. It reuses that
paper's population, exposure, outcome and covariate-source design **verbatim** — smidata's
`inst/analyses/doubletree__application.yml` marks sections 2–6 `status: inherited` — and differs
in exactly two places:

| | dual-bounds | this paper |
|---|---|---|
| estimator | `marbounds` partial-identification bounds | `doubletree::estimate_att()` / `estimate_att_crossfit()` |
| outcome missingness | modelled response indicator `R` | **complete-case restriction** through `INDEX_DT+24mo` |
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
directly. No local copy of `smi_select_aap_row()`, `smi_patient_windows()`,
`smi_sum_cost_in_windows()`, `smi_stream_cost_files()`, `smi_classify_is_mci()`,
`smi_window_coverage()` or `smi_read()` exists in this directory, and adding one would be exactly
the drift their promotion eliminated.

## Status

**Tier 0 (laptop, fixture):** reached for stages `02`–`04`, and as far as `06`'s own assembly
logic can be verified. `application/tests/` runs today with no data and no server — **165
assertions, 1 failure** (the 1 failure and 3 skips below are the local `optimaltrees.so`
security-tool block, unrelated to data/server availability — see `session_notes/2026-09-24.md`) —
and includes a call to `estimate_att()` at the real covariate shape (15 binary columns,
`outcome_type = "continuous"`, `leaf_budget = 4`) on a hand-built known-truth DGP, plus
`helpers/complete_case.R`'s and `helpers/msr_filter.R`'s full fixture suites (added 2026-09-24).

`Rscript application/run_pipeline.R` (new 2026-09-24, convenience only — every script cascades to
its own upstream dependency via `require_stage()` and is independently runnable) sources `01`–`07`
in order and now runs `01`–`04` fully against Tier-0 fixtures, then correctly **stops inside `05`'s
own limitation**: `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` column has no real-data census
levels, so `smi_fixture()` emits unparseable placeholder text and `compute_complete_case()` cannot
be exercised against it (it correctly aborts if asked to try — that is the desired behaviour
against a genuinely wrong format on real data). This is an honest Tier-0 boundary, not a defect.
`06`'s own assembly logic (join, three-reason attrition accounting, structural checks, the
PROVISIONAL write) was verified separately, with hand-built inputs standing in for `02`–`05`'s
outputs (bypassing the fixture cascade at exactly the point it cannot reach) — see
`session_notes/2026-09-24.md` for what that caught. **Do not hand-write a stand-in fixture for
the pipeline itself** — the hand-built inputs used to verify `06`'s own logic are a unit-test
technique for one script's internal joins/assertions, not a substitute for `smi_fixture()`.

**Tier 1 (server, subsample + whole-data checks):** not started. `90_checks_tier1.R` specifies
four checks; two are dual-bounds' own (referenced, not re-derived), one is this paper's own
(`_YN` coding), and CHECK 4 (added 2026-09-24) verifies the `enrollment_source` "ignore TYPE"
default's prerequisite at whole-table scope — `helpers/complete_case.R` self-checks the identical
question at cohort/window-restricted scope, but that is not a substitute for CHECK 4's run.
Also worth a Tier-1 pass once on the server: `03`'s MCI-exclusion count (unverifiable at Tier 0 —
the census carries no real `SS_DESC` levels) and the coverage-truncation rate for `Y`/`prior_cost`
(under-exercised at Tier 0 — the fixture's `INDEX_DT` distribution happens to fall inside the
2016–2024 file coverage more often than real data plausibly will).

**Tier 2 (server, full cohort):** blocked behind `diagnostics` (the one open decision still
`blocking_final = TRUE`). `02`'s gate is fixed (see `msr_code`'s row, above); stages `03`–`06`
still need attempting against real data once `diagnostics` resolves.

```bash
Rscript application/run_tests.R          # 165 passed, 1 failed (optimaltrees.so block), 3 skipped
```

## Open decisions

Registered in `_config.R`, echoed at the top of every numbered script. Confirmed decisions are
read via `confirmed_value()` (no warning); unresolved placeholders are read via `open_value()`
(warns on every access — see `_config.R`). **Only `diagnostics` still blocks the final artifact**
as of 2026-09-24; the other three resolved on 2026-09-17/2026-09-24 (see each row).

| id | shared with dual-bounds | status | one line |
|---|---|---|---|
| `enrollment_source` | same table, different use | **confirmed** (2026-09-24) | `MEDICAID_FLAG==1` in `larger_smi_medicaid_monthly_flag` for every one of the 24 months in `[INDEX_DT, INDEX_DT+24mo)`, no `TYPE` filter, zero gap tolerance — a PI-given, revisable default. Implemented in `helpers/complete_case.R`; provisional on CHECK 4 (below) not yet running. |
| `cost_family_scope` | **yes** — identical question | **confirmed** (`all_four`) | All four cost-claims families are disjoint claim sources; sum all four. |
| `msr_code` | **yes** — identical question | **confirmed** (`"AAP (no filter needed)"`) | `msr_aap` is already an AAP-only extract — no `MSR` filter needed at the *registry* level. **Fixed 2026-09-24**: `02_population_and_eligibility.R`'s gate now reads this value (was stale since 2026-09-17). But "no filter needed" is not what `02` actually does — `smidata::smi_select_aap_row()`'s abort-on-tie does NOT guard against a non-AAP row silently winning when it merely lands closer to the target month than any true AAP row (no tie forms), so `02` filters to the literal `AAP_MSR_CODE` explicitly and reports any exclusion. See this decision's `note` for the full correction — the registry's "abort-on-tie is the safety net" claim was wrong, and is not yet fixed upstream (smidata, dual-bounds). |
| `diagnostics` | no — this paper's own | **open, BLOCKS FINAL ARTIFACT** | Which sparsity **proxies** (`certified_*`, `n_leaves_*`, `gap_*`) license reporting the flagship estimate. Note sparsity itself is *not* checkable from data. |

The two shared entries carry machine-readable `shared_with` and `registry_ref` fields instead of a
re-worded question. Two projects independently wording one question is how the two answers end up
differing.

**Deliberately absent, relative to dual-bounds' registry** — each for a stated reason, not by
oversight:

- `admin_gap_days` — superseded. There is no modelled `R` here to absorb a gap tolerance; the
  complete-case restriction is a hard inclusion criterion, so permitting a gap would make the
  retained patients' `Y` a sum over a window with holes in it. A tolerance would be a *new* PI
  decision with its own registry entry.
- `msr_tie_break` — not applicable. This pipeline calls `smi_select_aap_row()` at smidata's own
  default (`"earlier"`), and smidata aborts on any other value. No choice is being made here.
- `sl_lib` — not applicable. SuperLearner is `marbounds`' nuisance machinery. doubletree fits its
  nuisances with optimal trees via `optimaltrees`.

## Files

| file | runs today | what it does |
|---|---|---|
| `_config.R` | yes | Gates, `OPEN_DECISIONS`, `SMI_KEYS`, `SMI_COLS`, the 36 cost dataset keys, `DOUBLETREE_WINDOWS`, `config_estimator`, `config_leaf_budget`. Plain top-level variables, edited in place — not `Sys.getenv()`. |
| `helpers/covariate_blocks.R` | yes (tested) | The **PI-confirmed** 15-column selection, 4 blocks × 3 plus 3 `prior_cost` dummies, and the declared column ORDER. Authoritative, *not* a candidate — contrast dual-bounds' file of the same name. |
| `helpers/discretize_prior_cost.R` | yes (tested) | Pure: dollars → 3 quartile dummies, Q1 implicit reference, cutpoints returned as an attribute. Aborts on a collapsed boundary or an empty bin. |
| `helpers/design_matrix.R` | yes (tested) | `assemble_design_matrix()`, and `assert_binary_design_matrix()` — the four checks the package does **not** give you (see below). |
| `helpers/complete_case.R` | yes (tested) | `parse_year_month()`, `ym_month_index()`, `compute_complete_case()` — the real `enrollment_source` logic, pure and hand-fixture-tested (no server, no `smidata`). See its file header for the split rationale. |
| `helpers/msr_filter.R` | yes (tested) | `filter_to_msr_code()` — the real `msr_code` contamination guard `02` uses (see that decision's row, above), pure and hand-fixture-tested. |
| `01_declare_requirements.R` | yes | Declares 39 dataset requirements (added `medicaid_monthly_flag` 2026-09-24) + the 12 `_YN` columns via `smidata::smi_require()`; reports the census gap. Skips gracefully without smidata. |
| `02_population_and_eligibility.R` | yes (fixture) | **Fixed 2026-09-24** (was stale — see `msr_code`'s row above). Runs end to end locally against a Tier-0 fixture (`smidata::smi_read()` dispatches to `smi_fixture()` off-server automatically): builds `base_cohort`, filters `msr_aap` to `AAP_MSR_CODE` (`helpers/msr_filter.R`, tested), calls `smi_select_aap_row()`, reports eligibility/exposure-arm counts. On the server it reads the real tables instead — same code path, `smidata::smi_env()` dispatches. |
| `03_cost_windows.R` | yes (fixture) | **Real as of 2026-09-24.** Cascade-guarded on `02` (`require_stage()`). Builds `windows`/`coverage`, then `smi_sum_cost_in_windows()` for real (~10s against 36 Tier-0 fixtures, verified before wiring live). `outcome_cost`/`prior_cost` carry an explicit `coverage` column and are `NA` (not coalesced to 0) whenever file-date coverage is not `"full"` — a patient with a partial-window sum is marked unknowable, not passed through as if complete; see this file's header. |
| `04_covariates.R` | yes (fixture) | **Real as of 2026-09-24.** Cascade-guarded on `03`. Reads the 12 confirmed `_YN` columns + `INDEX_DT`, joins `03`'s `prior_cost` (no `coalesce()` — see `03`'s row), runs a PREVIEW discretize/assemble/assert pass on the pre-filter cohort. Explicitly a preview: `06` re-discretizes on the real, final (complete-case-filtered) sample. |
| `05_complete_case.R` | yes (fixture, up to a real Tier-0 boundary) | Computes complete-case status for real via `helpers/complete_case.R::compute_complete_case()`, cascade-guarded on `02`. Locally, `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` has no real-data census levels, so the fixture emits unparseable placeholder text and the real computation is skipped with an explicit message (not a crash) — the hand-built two-patient demo still runs regardless, proving the function itself. On the server this limitation does not apply. Deliberately does **not** proxy enrollment from claims presence — `MEDICAID_FLAG` is a real enrollment indicator, independent of utilization; the file's header spells out why neither directional failure mode of a claims-presence proxy applies to it. |
| `06_assemble_analytic_data.R` | yes (fixture chain stops at `05`'s boundary; own logic verified separately) | Joins `02`–`05`, reports THREE separate attrition reasons (not complete-case; unknowable `Y`; unknowable `prior_cost` — see `03`'s coverage column), re-discretizes `prior_cost` on the real final sample, `assert_binary_design_matrix()` (unconditional), then writes with a `_PROVISIONAL` suffix while anything blocks, or without one once nothing does. **The PROVISIONAL release valve is now real** (fixed 2026-09-24, design-reviewed by `oracle`) — it was pure documentation before: `assert_no_blocking_open()` always aborted first, so the write was never reachable. Moved that call to `07` only (`diagnostics` gates the *estimate*, not the assembled *data*). **Found and fixed while verifying this script with hand-built inputs**: the original complete-case filter used `isTRUE(.data$complete_case)`, which is not vectorized — it would have silently zeroed the analytic sample on every real run, for a reason unrelated to complete-case status itself. Fixed to `.data$complete_case %in% TRUE`. |
| `07_estimate_att.R` | **aborts by design** | Gate (`assert_no_blocking_open()`) is the first statement after config — the only place it still runs, as of 2026-09-24 (moved out of `06`). Then both estimators, and their discrepancy. |
| `run_pipeline.R` | yes | New 2026-09-24. Convenience only (sources `01`–`07` in order) — every script is independently runnable via its own `require_stage()` cascade; this just saves typing seven commands. |
| `90_checks_tier1.R` | no | Comment-only spec. Checks 1–2 reference dual-bounds'; check 3 (`_YN` coding) is this paper's own. |
| `tests/` | yes | 165 assertions, no data, no server. |
| `run_tests.R` | yes | `testthat::test_dir("application/tests")`. |

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
- **Windows are half-open**: start inclusive, end exclusive. `prior_cost` ends *at* `INDEX_DT`, so
  the index day belongs to neither window and no claim is counted twice.
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
- **Run time.** The cost-claims family is 171.5 GB across 36 files — 91.8% of every byte in the
  secure tree. `smi_sum_cost_in_windows()` streams one file at a time for that reason. Budget
  hours.

## Provenance

The confirmed design is recorded in
`~/RAND/tools/smidata/inst/analyses/doubletree__application.yml`; that file, not this directory,
is authoritative. Column names and row counts come from smidata snapshot `2026-09-15_2824080e`.
The 15-column selection and `prior_cost`'s quartile discretization were confirmed by the PI on
2026-09-17.
