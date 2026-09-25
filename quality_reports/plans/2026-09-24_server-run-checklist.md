# Plan: Server-Run Checklist for `application/` (real-data pipeline)

**Date:** 2026-09-24
**Status:** **NOT YET ACTIONABLE.** Blocked on data ingest, not on this pipeline's own
readiness — see §0. This doc is written now so nothing is guessed at once ingest lands;
it is a readiness checklist, not a green light.
**Scope:** doubletree's `application/` real-data pipeline (`01`–`07` + `90_checks_tier1.R`).
Does not cover dual-bounds'/eqwt's own pipelines, though several items below are shared
with dual-bounds (marked `[SHARED]`) and resolving them resolves both at once.
**Why this doc exists:** Today's session got the pipeline to run end-to-end against
Tier-0 fixtures, up to a real, honest boundary (`05`'s `YEAR_MONTH` gap — see §2). The
next real progress needs the actual server, which requires (a) the data actually being
there in a form `smi_read()` can use quickly, and (b) someone with server/Citrix access
to run it. Neither is available to the agent that did today's work. This doc is the
handoff: exactly what to check, run, and report back, in the right order, so the person
running it doesn't have to re-derive any of today's findings.

---

## 0. Read this before doing anything

**As of 2026-09-24 13:15 (smidata's own session notes, `~/RAND/tools/smidata/session_notes/2026-09-24.md`),
the frozen-data ingest (`smi_ingest()` / `06_ingest.R`) is in progress at roughly 16 of 64
datasets per 24 hours, with an estimated ~3 more days remaining as of that timestamp** — i.e.
plan on this being ready no earlier than **2026-09-27**. That same entry states directly:
**"no real run is possible today regardless"** of any of the three sibling pipelines'
(dual-bounds, doubletree, eqwt) own code readiness. This is not a doubletree-specific
blocker — it applies identically to everything below.

**Why this matters mechanically, not just as a "wait" instruction.** `smidata::smi_read()`
resolves the frozen ingest cache (fast, cached Parquet) first, and falls back to the RAW
`.sas7bdat` file only when no ingest holds that dataset yet. Per
`~/RAND/tools/smidata/inst/docs/SERVER_SESSION_RUNBOOK.md`: **opening a raw `.sas7bdat` file
costs HOURS, regardless of how much is subsequently read from it** — this is a per-file-open
cost, not a per-byte cost, and it scales with file size. **Before ingest completes, every
`smi_read()` call in this pipeline that touches a not-yet-ingested dataset is an hours-scale
operation**, and the runbook's own guidance is to budget **one** such operation per server
session, not several.

**CORRECTION (2026-09-25, `session_notes/2026-09-25.md`): checking "per-dataset progress" is
not actually useful, and the paragraph below originally implied otherwise.** Read
`smi_ingest()`'s actual implementation (`~/RAND/tools/smidata/R/ingest.R`) and `smi_read()`'s
server-side lookup (`.smi_newest_ingest_path()`, `R/read.R`): `smi_ingest()` stages every
requested dataset into `.staging/<ingest_id>/` and publishes the WHOLE batch via a single
atomic rename plus a `_COMPLETE` sentinel, written only once **every** dataset in that run has
finished. The code's own comment states the design intent directly: *"Never cherry-picks
survivors. A partial `<ingest_id>` invites someone to pin an ingest that is silently missing a
dataset."* `smi_read()`'s lookup filters strictly to `complete == TRUE` ingests, so **an
in-progress run is invisible to `smi_read()` regardless of how many individual files it has
already converted internally** — a dataset finishing early inside the batch (e.g.
`aim3_svc_cost_16` done on day one) does not make it readable via the fast path one hour
sooner than the last dataset in that same run. `06_ingest.R` calls `smi_ingest()` exactly once
for the whole run, not once per file, so "16/64 in 24h" is one batch's internal progress, not
16 independently-published datasets.

**What this means for the table below:** it still tells you WHICH four datasets matter and
roughly how expensive each is, but "check whether dataset X specifically is done" is not a
useful question to ask while the run is in progress — the only question that matters is
whether the WHOLE run (all datasets requested in this `smi_ingest()` call, not just these four)
has published. Ask whoever is running it for the run's overall completion status, not a
per-file breakdown.

**Action before running anything below:** confirm the batch has fully published (the ETA below,
or ask directly) for these four datasets this pipeline touches:

| Dataset | Needed by | Size (2026-09-15 fingerprint) |
|---|---|---|
| `larger_smi_covariates` | `02`, `04` | moderate (72 cols × 381,018 rows) |
| `msr_aap` | `02` | large (19,012,819 rows) |
| `larger_smi_medicaid_monthly_flag` | `05`, CHECK 4 | 5.8 GB / ~24M rows |
| the 36 cost-claims files (`aim3_svc_cost_*`, `new_svc_cost_*`, `aim3_pharm_cost_*`, `new_pharm_cost_*`) | `03` | **171.5 GB — 91.8% of every byte in the whole tree** |

If the batch hasn't fully published yet, you have two real options, not just "wait":

- **Wait for full publish** (~3 more days per the ETA above) — safest, and every `smi_read()`
  call across the whole pipeline becomes fast afterward.
- **Use `smi_read()`'s raw-`.sas7bdat` fallback today** for ONE specific dataset if there's real
  urgency — it works independent of the in-progress ingest (reads its own copy, no conflict),
  at the hours-scale-per-file cost the runbook warns about. Only worth it if starting `02`
  (say) today matters more than waiting; still budget one such operation per session.

Either way, see §1 for why attempting several stages via the raw-file path across separate
sessions is riskier than it looks.

---

## 1. The one thing that must be decided BEFORE the first real run

**`03_cost_windows.R`'s expensive step has no cross-session checkpoint.** Each numbered script
now cascades to whatever upstream object it needs via `require_stage()` (`_config.R`) — this
is exactly what makes `Rscript application/06_....R` work standalone against fixtures. But on
the server, if `03` is run in one session and `04`/`05`/`06` are run in a **separate** process
later (a new `Rscript` invocation, or a fresh Citrix session), `require_stage()` will find
`outcome_cost`/`prior_cost` missing and **re-run all of `03` from scratch** — re-streaming all
36 cost-claims files (171.5 GB) a second time. Combined with §0's "one hours-scale operation
per session" budget, this is the single most expensive mistake available in this checklist.

**Two ways to avoid it — pick one before starting, don't discover it mid-run:**

- **(a) Run `01`–`07` in ONE continuous R session** (`Rscript application/run_pipeline.R`, or
  `source()` the seven scripts in order at an interactive console) that stays alive through
  the whole chain, so `03`'s in-memory `outcome_cost`/`prior_cost` are still there when `04`
  needs them. Simplest, no code change, but means accepting one long sitting (budget for the
  cost-claims pass alone being multi-hour, plus everything else) and a Citrix-disconnect risk —
  see the runbook's own Step 8 ("Citrix logoff survival test") for whether a long job survives
  a disconnect on this VDI.
- **(b) Implement a checkpoint first**: have `03` write `outcome_cost`/`prior_cost` to
  `application/output/` (parquet, matching `06`'s own convention), and extend `04`'s
  `require_stage()` call to check that file before falling back to re-running `03`. Not done as
  part of today's session — flagged, not fixed (see `session_notes/2026-09-24.md`'s 15:19
  entry, "Next"). Worth 30–60 minutes of implementation if the server run is going to be split
  across sessions regardless (e.g. because `03` alone doesn't fit in a single sitting).

**Recommendation: do (b) before attempting the real run**, since a Citrix session surviving a
171.5GB streaming pass uninterrupted is not something to bet a "no real run possible" multi-day
wait on. Ask for this to be implemented (it's a small, self-contained change) rather than
discovering the gap mid-run.

---

## 2. What already works, and where it honestly stops (today's state)

Confirmed against Tier-0 fixtures, 2026-09-24 (commits `3874761`, `26fb071`):

- `01`–`04` run end-to-end for real (fixture-backed locally; same code path reads real tables
  on the server via `smidata::smi_env()` dispatch — no code difference between fixture and
  real runs).
- `05` computes complete-case status for real via `helpers/complete_case.R::compute_complete_case()`,
  **except** it cannot be exercised against the Tier-0 fixture specifically because
  `larger_smi_medicaid_monthly_flag`'s `YEAR_MONTH` column has no real-data census levels (the
  fixture emits placeholder text `FIXTURE_TEXT_<n>`, which the real parser correctly refuses to
  parse — that refusal is the desired behavior against a genuinely wrong format on real data).
  **This limitation is fixture-only; it does not apply on the server.**
- `06` assembles the analytic dataset, with three separately-reported attrition reasons
  (not complete-case; unknowable `Y`; unknowable `prior_cost`, both from file-date coverage
  truncation — see `03`'s row below) and a working `_PROVISIONAL` release valve. Verified with
  hand-built inputs (the fixture chain can't reach `06` locally, per `05`'s limitation above) —
  this caught a real bug (§3, item 4) that predates today's session.
- `07` still correctly refuses to run (`assert_no_blocking_open()`): `diagnostics` is the one
  `OPEN_DECISIONS` entry still `blocking_final = TRUE`.

**Net effect:** every line of code that will run on the server has already been exercised as
far as fixtures allow, and one real, previously-latent bug was found and fixed in the process
(§3, item 4) — the server run's job is to exercise the parts fixtures structurally cannot
(§0's four datasets' real content), not to debug the pipeline's own logic from scratch.

---

## 3. Bugs found and fixed today that the server run should specifically confirm

1. **`msr_code`'s stale gate** (`02`) — was calling `open_value()` on an already-confirmed
   decision, silently never reaching the real assembly. Fixed; also found and fixed a real gap
   in `smi_select_aap_row()`'s documented safety net (it does not protect against a non-AAP row
   silently winning when no tie forms — see `OPEN_DECISIONS$msr_code$note`). **Confirm on the
   server:** `02`'s MSR-filter warning should fire ZERO exclusions if the "always AAP" finding
   holds; if it fires with a nonzero count, that is a real finding, not a fixture artifact —
   report the exact count and the other MSR values seen.
2. **`cost_family_scope`'s identical stale-gate bug** (`03`) — same fix, same pattern.
3. **Y/`prior_cost` coalesce bug** (`03`) — a patient with PARTIAL file-date coverage but a
   real (truncated) partial cost sum was passing that sum through unmarked. Fixed to force `NA`
   whenever `coverage != "full"`, regardless of whether `smi_sum_cost_in_windows()` returned a
   sum. **Confirm on the server:** the real coverage-truncation rate — Tier-0 fixtures
   under-exercise this (fixture `INDEX_DT` happens to fall inside 2016–2024 more than real data
   plausibly will). Report the `outcome_cost`/`prior_cost` known-vs-unknown counts `03` already
   prints.
4. **`06`'s `isTRUE()` bug — the most important one.** `dplyr::filter(analytic_pre,
   isTRUE(.data$complete_case))` is inherited near-verbatim from the ORIGINAL pre-existing
   pseudocode (predates today). `isTRUE()` is not vectorized — it collapses a whole column to
   one logical (almost always a single `FALSE`), and `dplyr::filter()` then applies that ONE
   value to every row, silently dropping ALL of them. **This would have zeroed the analytic
   sample on the very first real run**, for a reason unrelated to whether `complete_case` was
   actually resolved — i.e., it would have looked exactly like "complete-case status isn't
   computable yet," which is a plausible-sounding wrong explanation that could easily have gone
   unquestioned. Fixed to `.data$complete_case %in% TRUE`. **Confirm on the server:** `06`'s
   attrition report ("N eligible -> M complete-case -> ...") should show a plausible,
   nonzero `M`, not zero.
5. **MCI-exclusion count is unverifiable at Tier 0** (`03`) — the census has no real `SS_DESC`
   levels, so the fixture can never exercise the "Managed Care Invoice" exclusion regex.
   **Confirm on the server:** `n_excluded_mci` across the cost-window sum should be nonzero
   (03 already warns if it's exactly 0 when `!is_local`). If it IS zero on real data, that is a
   real finding — check `config_cost_setting_col`/the exclude regex before trusting any cost
   number.

---

## 4. Step-by-step, once §0/§1 are cleared

Run everything below in ONE continuous session (§1) unless the checkpoint fix has been
implemented.

```bash
cd application
# Confirm environment resolves to "server", not "local":
Rscript -e 'cat(smidata::smi_env(), "\n")'
```

If that prints `local`, something is wrong with the environment (see the runbook's own
troubleshooting table for `smi_paths()` / `SMI_DATA_DIR` — drive-letter/casing issues are the
usual cause on a fresh VDI) — fix that before proceeding; every subsequent step silently
becomes a fixture run otherwise, which is a hard failure mode to notice ("it worked" but on
the wrong data).

```bash
Rscript run_pipeline.R
```

Expect, in order: `01` reports 39/39 declared against the real contract (seconds). `02` reads
`larger_smi_covariates` + `msr_aap` for real (hours if not yet ingested, per §0). `03` runs the
36-file cost-claims streaming pass — THE expensive step (hours; 171.5 GB). `04` is fast (joins
already-computed tables). `05` reads `larger_smi_medicaid_monthly_flag` for real (hours if not
yet ingested) and, unlike the fixture run, should NOT hit the Tier-0-only limitations in §2 —
if it does, that's a new finding, not expected. `06` writes `analytic_cohort_PROVISIONAL.parquet`
(still provisional — `diagnostics` remains open). `07` still correctly aborts.

---

## 5. The four Tier-1 checks (`90_checks_tier1.R`) — separate from the pipeline run above

These are comment-only specs today (server-only, never executed) — typed out precisely enough
to run without re-deriving anything. Each needs its own dataset read, so each is its own
budget-one-operation-per-session item under §0, independent of whether the pipeline run in §4
has happened yet.

| Check | Needs | Resolves | Status |
|---|---|---|---|
| CHECK 1 | cost-claims (one shared suffix year) | `cost_family_scope` (already confirmed, PI 2026-09-17 — this is corroboration, not a blocker) | Not run |
| CHECK 2 | `msr_aap` | `msr_code` (already confirmed — corroboration) | Not run |
| CHECK 3 | `larger_smi_covariates` | Whether the 12 confirmed `_YN` columns are genuinely `{0,1}` (this paper's own; not shared) | Not run |
| **CHECK 4** | `larger_smi_medicaid_monthly_flag` | The prerequisite for `enrollment_source`'s "ignore TYPE" default — **the highest-priority of the four**, since `05`'s real output depends on it and it's currently an unverified assumption, not a confirmed one | Not run |

**CHECK 4 is the one to prioritize** if only one Tier-1 check fits in a session before the full
pipeline run: it's the prerequisite for trusting `05`'s output at all, per
`OPEN_DECISIONS$enrollment_source`'s own note. `compute_complete_case()` already self-checks
the identical question at cohort/window-restricted scope and **aborts** if it ever finds
disagreement — so if CHECK 4 finds TYPE disagreement is rare/localized but the pipeline run in
§4 didn't abort, that's not a contradiction (different scope, see `helpers/complete_case.R`'s
own Details); if CHECK 4 finds disagreement is common, expect `05`'s real run to abort with a
message naming it, not to silently produce a wrong number.

Full specification for each check is written out in `90_checks_tier1.R` itself — read it
directly on the server rather than retyping it from this table; this table is a priority/status
summary, not the spec.

---

## 6. What to bring back / report format

For each script run, capture (copy-paste is fine, these are `cli`-formatted console messages
already designed to be self-explanatory):

- **`01`**: the "Declared N of N" line and the Tier-0/environment lines (confirms real vs. fixture).
- **`02`**: the MSR-exclusion warning (or its absence), `Eligible: X of Y`, `A=1`/`A=0` counts.
- **`03`**: the coverage-truncation warning's exact numbers, the MCI-exclusion warning (or its
  absence), `outcome_cost`/`prior_cost` known-vs-unknown counts.
- **`04`**: the preview cutpoints line (informational only — `06`'s are the real ones).
- **`05`**: whether it reaches the real computation at all (§2's limitation should NOT apply on
  the server) — the `Complete-case: N of M` line and the by-arm attrition table.
- **`06`**: the full attrition chain (`N eligible -> M complete-case -> ... -> final`), the
  realized quartile cutpoints (methods-section material), arm counts, and confirm the written
  file is `analytic_cohort_PROVISIONAL.parquet` (expected — `diagnostics` still blocks the
  reportable one) plus its manifest's contents.
- **CHECK 3/CHECK 4**: the exact tables specified in `90_checks_tier1.R` (coding/prevalence
  per column for CHECK 3; the (a)/(b)/(c) disagreement tabulation for CHECK 4), suppressed per
  each check's own stated `smi_suppress()` rule before anything leaves the server.

Bring back console output/exported tables, not just a verbal summary — several of today's
fixture-only findings (the truncation-rate numbers, the MCI-exclusion count) are exactly the
kind of thing that's easy to mis-remember as "looked fine" without the actual number.

---

## 7. Safety reminders (mostly already enforced by the code, restated here for the human)

- **Never** report a number computed from a Tier-0 fixture as a real finding — the code already
  warns loudly every time this could happen, but the discipline extends to what gets written in
  a report/email afterward too.
- **Suppression**: every count/table named in §5's checks has an explicit `n < 11`
  (CMS/HIPAA Safe Harbor) suppression rule stated in `90_checks_tier1.R` — route through
  `smidata::smi_suppress()` before anything leaves the server, per each check's own spec.
- **Do not hand-write a stand-in fixture** for anything the real run can't reach — if `05`'s
  Tier-0 limitation (§2) turns out to somehow also apply on the server (it should not), stop and
  report that as a finding rather than working around it.
- **If `03` needs to be re-run** because of §1's checkpoint gap, say so explicitly in whatever
  gets reported back — a second 171.5 GB pass is a cost worth knowing about, not absorbing
  silently into "it took a while."
