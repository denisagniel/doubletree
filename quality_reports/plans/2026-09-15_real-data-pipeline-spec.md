# Plan: Real-Data Application Pipeline Spec (smidata-shaped, uniform template)

**Date:** 2026-09-16 (supersedes the 2026-09-15 spec-shell — see "History" at the bottom).
**Status:** Design decided by the PI (2026-09-16); reuses dual-bounds' design wholesale for
population/exposure/outcome/covariate-source, differing only in estimator. Not yet built —
doubletree has zero real-data application code.
**Scope:** doubletree's real-data application only (the package's methods/simulation content,
mature and out of scope here — see `README.md`).
**Why this doc exists:** Fifth of five pipeline specs written this session against one uniform
8-point template. Unlike the shell it replaces, this version has a real design to record: the
PI decided (2026-09-16) that doubletree reuses dual-bounds' entire real-data design
(`smi/dual-bounds/quality_reports/plans/2026-09-15_real-data-pipeline-spec.md`) rather than
building its own from scratch, applying doubletree's own tree-based ATT estimator to the same
applied dataset instead of dual-bounds' `marbounds` partial-identification method.

---

## 1. Population / cohort

**Reused from dual-bounds' spec (§1), PI-confirmed 2026-09-16.** Same base table, same anchor,
same eligibility filter:

- **Base table: `larger_smi_covariates`, ~381,191 rows.** Same table as the finalized paper,
  telehealth-analysis, and dual-bounds. NOT `aim1_smi_charac` (that's the applied
  paper/health-affairs-case-mix-paper's smaller cohort — see §6).
- **Anchor: `INDEX_DT`**, taken directly from `larger_smi_covariates`. Pre-computed in the
  extract; no bespoke diagnosis-date derivation needed.
- **Eligibility filter: AAP eligibility, REPLACING age/enrollment filters entirely.**
  `elig_aap = 1` iff `MSR_DEN >= 1` in the `msr_aap` file. PI direction (2026-09-16): AAP
  eligibility supersedes the original pseudocode's age/enrollment filters — this was an open
  question in dual-bounds' own spec (§7 item 1) at the time this doc was written; apply the
  resolution here regardless of whether dual-bounds' file has been updated to say so yet.

Doubletree adds no eligibility criteria of its own beyond this. Whatever open questions remain
about the AAP-eligibility mechanics (see dual-bounds' spec §7) apply identically here.

## 2. Exposure/treatment

**Reused from dual-bounds' spec (§2), PI-confirmed 2026-09-16.** `A` = AAP achievement in
months 0–12 post-`INDEX_DT`: whether the patient achieved the `msr_aap` measure
(`MSR_NUM/MSR_DEN >= 0.5` among those eligible) during `[INDEX_DT, INDEX_DT+12mo]`. Binary, no
change needed for doubletree's estimator — `estimate_att()` and `estimate_att_crossfit()` both
take a binary treatment vector `A` as-is.

No raw `pharmacy_claims` needed for exposure construction, for the same reason dual-bounds
doesn't need it: `A` is sourced from the existing `msr_aap` measure infrastructure, not built
from fill-level claims.

## 3. Outcome

**Reused from dual-bounds' spec (§3), PI-confirmed 2026-09-16.** `Y` = total cost, measured in
months 13–24 post-`INDEX_DT` (`[INDEX_DT+12mo, INDEX_DT+24mo]`), same disjoint exposure-then-
outcome window structure as the finalized paper and applied paper.

**Managed Care Invoice exclusion applies identically here**: cost computed from
`cost_claims`-family files (see §5) must exclude Managed Care Invoice records, matching every
other cost computation in this cross-project effort.

`Y` is continuous. doubletree's estimator supports this directly via `outcome_type =
"continuous"` on either `estimate_att()` or `estimate_att_crossfit()` (confirmed in
`README.md`, "Continuous outcome" note) — no outcome-type mismatch with dual-bounds' choice.

**`R`/death handling: RESOLVED (PI, 2026-09-16) — complete-case restriction.** Dual-bounds' `R`
(response indicator) exists because `marbounds`' joint-sensitivity model explicitly models
informative missingness/censoring as part of its partial-identification approach.
`estimate_att()`/`estimate_att_crossfit()` take `(X, A, Y)` with no censoring argument — they
assume `Y` is observed for the estimation sample. Confirmed: doubletree's real-data application
**restricts the estimation sample to patients with `Y` observed through `INDEX_DT+24mo`** — a
complete-case restriction, not a censoring-aware estimator. This is a genuine difference from
dual-bounds, not a copy-forward gap: it REPLACES dual-bounds' `R`/death-ascertainment claims
rather than inheriting them, since those have no place to attach in doubletree's estimator
signature.

## 4. Covariates

**Reused from dual-bounds' spec (§4), PI-confirmed 2026-09-16, WITH ONE NEW CONSTRAINT dual-
bounds does not have.** Covariates come entirely from `larger_smi_covariates`'s own columns —
no separate raw-claims-derived covariate blocks. Exact column list pending a real server
fingerprint (per every other spec in this series); dual-bounds' spec notes ~79 candidate
columns, unverified.

**New open question specific to doubletree, not present in dual-bounds' spec: grid-exact
sparsity / covariate dimensionality.** doubletree's flagship no-sample-splitting estimator,
`estimate_att()`, requires:
- Covariates in **binary (0/1) form** (confirmed in `README.md`: "Covariates must be binary
  (0/1); see `?estimate_att` for why").
- **Structural/grid-exact sparsity**: a tree of at most `leaf_budget` leaves must represent both
  true nuisance functions (propensity and outcome model) exactly on the analyst's chosen grid.

Neither constraint applies to dual-bounds' `marbounds` estimator, which is a genuine difference
in what doubletree needs from the same underlying data. Whether `larger_smi_covariates`'s ~79
columns, once confirmed, will need discretization (continuous → binary indicators) and/or
dimensionality reduction to make grid-exact sparsity plausible is **not yet resolved** — see §7.
If sparsity is implausible for the real covariate set, `estimate_att_crossfit()` (the cross-
fitting fallback, which does not require structural sparsity but still needs binary covariates
per the same README note) is the fallback path, not a change to the covariate list itself.

## 5. Dataset/measure inventory

**Same inventory as dual-bounds' spec (§5), PI-confirmed 2026-09-16:**

| dataset_key | Role | Status |
|---|---|---|
| `larger_smi_covariates` | Cohort base, `INDEX_DT`, baseline covariates | Confirmed elsewhere this session (finalized paper, telehealth-analysis, dual-bounds) |
| `msr_aap` | AAP eligibility (`MSR_DEN`) + exposure (`MSR_NUM`/`MSR_DEN`) | Confirmed elsewhere this session (columns: `ID`, `MSR_YR`, `MSR_NUM`, `MSR_DEN`) |
| `cost_claims` family (`aim3_svc_cost_*`, `new_svc_cost_*`, possibly `aim3_pharm_cost_*`/`new_pharm_cost_*`) | Outcome (cost) | Same naming discrepancy flagged in vector-incremental-effects's applied-paper spec and dual-bounds' spec — PI-confirmed ground truth is `aim3_svc_cost_*`/`new_svc_cost_*`, not the `service_claims`/`svc-cost-data_*` names some project code documents |

No dataset_key is unique to doubletree — this application reuses dual-bounds' entire data
footprint. Doubletree's own contribution is the estimator applied to it, not new data.

## 6. Relationship to other papers

**To dual-bounds: same applied dataset, different estimator — this is the point of doubletree's
own contribution.** Same population (`larger_smi_covariates`, AAP-eligible), same exposure (AAP
achievement, months 0–12), same outcome (total cost, months 13–24). The two papers give the
reader a shared applied example analyzed two ways: dual-bounds' partial-identification/joint-
sensitivity bounds via `marbounds`, doubletree's point-identified ATT via interpretable trees
(`estimate_att()`/`estimate_att_crossfit()`). This is deliberate — a genuine methodological
contrast on the same data, not a coincidence of convenience.

**To the finalized paper and telehealth-analysis:** same base-cohort table
(`larger_smi_covariates`), same window convention (exposure 0–12, outcome 13–24 post-index).

**To the applied paper (vector-incremental-effects) and health-affairs-case-mix-paper:**
DIFFERENT base cohort. Those two use `aim1_smi_charac` (~281,941 rows), the smaller cohort — the
still-open D001 cohort discrepancy (see the smidata plan and dual-bounds' spec §6) applies to
that pairing, not to doubletree/dual-bounds/finalized-paper/telehealth-analysis, which all agree
on `larger_smi_covariates`.

## 7. Open questions

Nearly everything in §1–3 that was open in dual-bounds' spec is either resolved by PI direction
(the AAP-replaces-filters question) or simply inherited unchanged (AAP-eligibility join
mechanics, `cost_claims` naming). Carried forward, plus doubletree-specific additions:

1. **AAP-eligibility join mechanics** (inherited from dual-bounds' spec §7 item 1, now resolved
   as "replaces" per PI direction 2026-09-16 — the mechanical join to `larger_smi_covariates`
   still isn't built).
2. ~~**Framing**: does AAP-eligibility + AAP-achievement-as-exposure need reconciling with
   doubletree's own manuscript framing?~~ **RESOLVED (PI, 2026-09-16), addressed by the same
   manuscript edit as item 7 below** — `manuscript.tex`'s Discussion section now names the
   AAP-adherence exposure and cost outcome explicitly, so the framing is already reconciled as
   drafted, not left implicit.
3. **Covariate source**: can `larger_smi_covariates`'s ~79 columns populate a usable covariate
   set at all? (Inherited from dual-bounds' spec §7 item 3.)
4. **NEW, doubletree-specific — grid-exact sparsity / covariate dimensionality (§4).** Does the
   real covariate set, once confirmed, support a small-leaf-budget tree that represents both
   nuisance functions exactly? If not plausible, does discretization of continuous columns into
   binary indicators get there, or is `estimate_att_crossfit()` the required fallback instead?
   This question has **no counterpart in dual-bounds' spec** — `marbounds` has no sparsity or
   binary-covariate requirement, so this is a genuine additional burden doubletree's real-data
   application carries that dual-bounds' does not. **Still open — genuinely blocked on a real
   `smi_fingerprint()` run; nothing to decide until the covariate column list exists.**
5. ~~**Outcome missingness/attrition through month 24 (§3).**~~ **RESOLVED (PI, 2026-09-16):
   complete-case restriction — the estimation sample is restricted to patients with `Y` observed
   through `INDEX_DT+24mo`.** This REPLACES dual-bounds' `R`/death-ascertainment items rather
   than inheriting them, since those have no place to attach in doubletree's estimator signature.
6. **`cost_claims` naming discrepancy** — shared with vector-incremental-effects's applied paper
   and dual-bounds; needs a real `smi_fingerprint()` run.
7. ~~**Manuscript reconciliation.**~~ **RESOLVED (PI, 2026-09-16): updated now, not deferred.**
   `manuscript.tex`'s `\section{Discussion}` (line ~701) no longer says the real-data application
   is generically "left to future work" — it now names the specific planned application (same
   population/exposure/outcome as dual-bounds, contrasting identification strategy), and states
   the grid-exact-sparsity/binary-covariate question as an open, checkable question the
   application must answer rather than a precondition assumed to hold. Verified: `xelatex`
   compiles the manuscript cleanly (21 pages, 0 errors) after this edit.

## 8. Current code state

**Pure methods package. No real-data application code exists yet — this design is the starting
point for building it, not a description of anything runnable.** Confirmed:

- R package structure (DESCRIPTION/NAMESPACE/R/tests), `optimaltrees (>= 0.4.1)` as a hard
  `Imports` dependency.
- `estimate_att()`, `estimate_att_crossfit()`, `estimate_att_rashomon()` (superseded) are
  implemented and exercised by the mature simulation study (3 DGPs × 4 methods × 3 sample sizes
  × 500 replications) — but none of that harness has ever been pointed at real data.
- `data/` remains an empty placeholder directory (`.gitkeep` only, per the shell doc this spec
  replaces) — no fixture, no extract, no schema contract exists for this application.
- `inst/paper/manuscript.tex`'s Discussion section was updated 2026-09-16 (see §7 item 7) to
  name this application specifically — the manuscript text is ahead of the code here: the
  design is now described in prose before any pipeline script exists.

**Contrast with dual-bounds' own §8**: dual-bounds' spec describes a design sketch not yet
runnable, but with an existing 437-line pipeline pseudocode file (`analysis/real_data/
pipeline_pseudocode.R`) implementing an earlier, now-superseded version of a real-data pipeline.
Doubletree has no equivalent file at all — not even a superseded one. This spec is the first
real-data design artifact doubletree has ever had.

**Next step**: build a `smi_require_gate()`-declared skeleton (matching vector-incremental-
effects's pattern, per dual-bounds' spec's own stated next step) once §7 item 4 (grid-exact
sparsity feasibility) is checked against the real covariate fingerprint — that question gates
whether `estimate_att()` or `estimate_att_crossfit()` is the right entry point to build around.

---

## History

**2026-09-15**: original spec-shell (`2026-09-15_real-data-pipeline-spec-shell.md`) recorded
that doubletree's real-data application was genuinely undecided — no population, exposure,
outcome, or covariate set had ever been specified, and no file in the repo mentioned SMI, OMH,
Medicaid, or any real dataset. That file's own §7 posed the design questions this spec now
answers (whether doubletree reuses the SMI population, what treatment/outcome, and what the
grid-exact-sparsity constraint implies for covariate choice).

**2026-09-16**: PI decided doubletree's real-data application reuses dual-bounds' entire design
(population, exposure, outcome, covariate source) rather than building its own, applying
doubletree's own tree-based ATT estimator instead of dual-bounds' `marbounds` method. This file
replaces the shell as the authoritative spec. The shell is retained unmodified for the historical
record of the pre-decision blank state; a one-line pointer has been added to its header (see
that file) directing readers here.

**2026-09-16, second revision (same day, via smidata's cross-project analysis registry
walkthrough).** PI resolved three of §7's remaining open items: outcome missingness (complete-
case restriction through `INDEX_DT+24mo`, replacing dual-bounds' `R`/death-ascertainment claims
rather than inheriting them), manuscript reconciliation (updated `manuscript.tex` now, see §7
item 7 and §8), and, as a consequence of that same edit, the framing-reconciliation question
(§7 item 2). **Deliberately left open**: diagnostics. **Still open, genuinely blocked, not a PI
decision**: grid-exact sparsity/covariate dimensionality (needs the real column list) and
`cost_claims` naming (needs a real fingerprint).
