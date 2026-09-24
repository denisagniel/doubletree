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

**Tier 0 (laptop, fixture):** partially reached. `application/tests/` runs today with no data and
no server — **135 assertions, 1 failure** (the 1 failure and 3 skips below are the local
`optimaltrees.so` security-tool block, unrelated to data/server availability — see
`session_notes/2026-09-24.md`) — and includes a call to `estimate_att()` at the real covariate
shape (15 binary columns, `outcome_type = "continuous"`, `leaf_budget = 4`) on a hand-built
known-truth DGP, plus `helpers/complete_case.R`'s full fixture suite (added 2026-09-24). What is
*not* reached is a fixture-backed run of stages 02–06: `smi_fixture()` cannot build one, because
`smi_census()` has not been run on `larger_smi_covariates`, `msr_aap`, or the cost-claims family.
The unblocker is a server-side `smi_census()` run, not anything in this repo. **Do not hand-write
a stand-in fixture.**

**Tier 1 (server, subsample + whole-data checks):** not started. `90_checks_tier1.R` specifies
four checks; two are dual-bounds' own (referenced, not re-derived), one is this paper's own
(`_YN` coding), and CHECK 4 (added 2026-09-24) verifies the `enrollment_source` "ignore TYPE"
default's prerequisite at whole-table scope — `helpers/complete_case.R` self-checks the identical
question at cohort/window-restricted scope, but that is not a substitute for CHECK 4's run.

**Tier 2 (server, full cohort):** blocked behind `diagnostics` (the one open decision still
`blocking_final = TRUE`) and, separately, `02_population_and_eligibility.R`'s own stale gate (see
Files, above — `msr_code` resolved 2026-09-17 but `02`'s code was not updated to match).

```bash
Rscript application/run_tests.R          # 135 passed, 1 failed (optimaltrees.so block), 3 skipped
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
| `msr_code` | **yes** — identical question | **confirmed** (`"AAP (no filter needed)"`) | `msr_aap` is already an AAP-only extract — no `MSR` filter needed. **Not yet mirrored into `02_population_and_eligibility.R`'s own gate**, which still aborts pending this value; see that file's status below. |
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
| `01_declare_requirements.R` | yes | Declares 39 dataset requirements (added `medicaid_monthly_flag` 2026-09-24) + the 12 `_YN` columns via `smidata::smi_require()`; reports the census gap. Skips gracefully without smidata. |
| `02_population_and_eligibility.R` | **aborts by design** | Aborts on the unknown AAP `MSR` code. **Stale as of 2026-09-24**: `msr_code` is now `confirmed`, but this script's own gate has not been updated to read that value — a separate, pre-existing gap, not fixed by this pass. Intended pipeline in comments. |
| `03_cost_windows.R` | partly | Builds the 36-row `cost_sources`; smoke-checks `smi_patient_windows()`/`smi_window_coverage()` on two hand-written patients. The outcome window and the `prior_cost` window are two slices of **one** pass. |
| `04_covariates.R` | yes | Asserts the confirmed 15-column selection against the source column list. No fallback — there is nothing to fall back from. |
| `05_complete_case.R` | yes (real logic; server read still blocked on `02`, above) | Computes complete-case status for real via `helpers/complete_case.R::compute_complete_case()`. `enrollment_source` resolved 2026-09-24; the demo on this file's own tail runs on a hand-built fixture, not real data — the real orchestration is written out (commented) pending `analysis_cohort` from a working `02`. Deliberately does **not** proxy enrollment from claims presence; the two directional failure modes are spelled out in the file, and why `MEDICAID_FLAG` avoids both. |
| `06_assemble_analytic_data.R` | **aborts by design** | Joins 02–05, re-discretizes `prior_cost` **on the final sample**, `assert_binary_design_matrix()`, then `assert_no_blocking_open()`. |
| `07_estimate_att.R` | **aborts by design** | Gate is the first statement after config. Then both estimators, and their discrepancy. |
| `90_checks_tier1.R` | no | Comment-only spec. Checks 1–2 reference dual-bounds'; check 3 (`_YN` coding) is this paper's own. |
| `tests/` | yes | 135 assertions, no data, no server. |
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
