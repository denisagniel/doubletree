# S5 (`head_to_head_comparison`) claims audit — 2026-09-24

**Trigger.** The S5 cluster sweep finished 2026-09-23 (37/37 cells, per
`session_notes/2026-09-23.md`). Before drafting the manuscript's planned §4
"Simulation evidence" section against it, this pass re-derived every regime's
numbers directly from the on-disk `.rds` checkpoints (via the study's own
`analyze.R::summarise_cell()`/`rmse_ratio_table()`), rather than trusting the
prose summary in the session note, and cross-checked the result against
`claims.md`'s existing verdicts and `manuscript.tex`'s exact claim statements.
An `oracle` consultation (session `ses_f2c6c78b1ffelW5AfBYWwwKQMQ`) pressure-tested
the four findings below before any fix was applied; two of the four write-ups
here differ from the auditor's first draft because of that consultation
(marked inline).

**Scope check.** `manuscript.tex` currently has **zero** occurrences of
`head_to_head` / "Simulation evidence" (confirmed by grep, 2026-09-24) — no
manuscript prose yet cites any S5 number, so nothing here is a correction to
already-published text. This audit is timely precisely because it precedes
that drafting.

---

## Finding 1 — stale checkpoint contamination in `results/` (fixed)

19 of the 37 result cells still had a 2026-09-09 pilot/smoke-test checkpoint
(`_r2.rds` or `_r30.rds`) sitting alongside the 2026-09-24 full-scale one
(`_r1000.rds`/`_r2000.rds`). `analyze.R::load_results()` did an unconditional
`list.files()` + `rbind()` with no dedup, so any call to it — including the
one behind the 2026-09-23 "first substantive look" — silently pooled both
vintages for those 19 cells.

Quantified directly: **1,942** duplicate `(regime, dgp, n, arm, rep)` rows
across the 19 cells, of which **70 are NOT byte-identical** between the two
vintages (i.e. genuinely different `optimaltrees` code states pooled as if
they were exchangeable draws), concentrated almost entirely in R4/DGP-C
(64 of 70 — see Finding 4's mechanism).

**Numerical impact, checked directly (not asserted):** recomputing every
regime's `summarise_cell()` output from only the correct full-scale file per
cell moved no coverage/bias/CI-width number by more than ~0.1 percentage
point anywhere, and `rmse_ratio_table()`'s paired-bootstrap output for R1/R2
(the headline C1/C3 ratio statistics) was byte-identical naive-vs-clean.
Because `rep_seed()` is a pure function of `(regime, dgp, n, rep)`, the
over-weighted subset (reps 1–30) is an *arbitrary but exchangeable* slice of
the replication distribution, which is why over-weighting it can't bias a
point estimate — this was confirmed empirically, not just asserted.

**Fix applied 2026-09-24:**
- `analyze.R::load_results()` now dedups per cell (`dedup_checkpoints()`),
  keeping only the checkpoint with the largest `_r<nsim>` vintage suffix, and
  cross-checks every stale vintage against it: any `(arm, rep)` key shared
  between a stale and the kept checkpoint must match on `theta` to
  floating-point tolerance, or the function `cli::cli_abort()`s rather than
  silently pooling. Verified: running it against the (still-contaminated at
  the time) directory correctly aborted on the known R3/`B_eps15` mismatch;
  after archiving the stale files it loads cleanly and reproduces the
  hand-verified "clean" numbers exactly.
- Stale `_r2.rds`/`_r30.rds` files moved (not deleted) to
  `simulations/head_to_head_comparison/results/archive/pilot_2026-09-09_superseded/`.
  `results/` now contains exactly 37 `cell_*.rds` files, one per cell.

---

## Finding 2 — `rem:two-trees-why`'s "Covered" verdict oversteps the study's own licensing text

`claims.md` cited S5's R2 (C3: flagship ≈ its own cross-fit fallback at
matched `n`) as "Covered" evidence for *"no-splitting works **because** the
candidate class is fixed/finite, not from the bilinear product-rate condition
alone."* `head_to_head_comparison/README.md` §0 pre-commits: *"No comparative
sentence beyond 'matched X' or 'retained coverage where Y did not' is
licensed by this design."* C3 is a precision-matching **consequence**; the
claim cited for it is a proof-**mechanism** attribution — R2 cannot
distinguish "no-splitting works because of the fixed-finite-class argument"
from any other explanation producing the same observable consequence, and
structurally cannot even in principle, since R2 also confounds splitting
(yes/no) with nuisance-class depth (`max_depth=4` vs. the flagship's own
`Lbar`) at the same time.

`rem:two-trees-why` (`manuscript.tex:232-234`) actually bundles four
separable claims. Per the auditor's `oracle` consultation, the remark is
*not* purely a mechanism claim — its last sentence makes one falsifiable,
currently-untested empirical prediction that deserves its own GAP row rather
than being lumped in with the untestable pieces:

1. **(a)** The bias decomposition's product-rate condition,
   $D_wD_\mu\to0$, is indifferent to whether $\hat w$ is itself tree-based —
   an analytic consequence of `lem:biasbound`, same footing as `prop:bilinear`.
   **Not sim-testable as stated** (a population-identity claim, same class as
   `cor:single-tree`).
2. **(b)** No-splitting's own asymptotic linearity additionally needs the
   own-observation empirical-process term $o_p(n^{-1/2})$, and *that* is what
   the fixed/finite candidate class + exact structural sparsity actually
   controls — the "why no-splitting works" mechanism attribution.
   **Not sim-testable as stated**; R2/C3 shows only the consequence
   consistent with it.
3. **(b)′** *"An arbitrary non-tree $\hat w$ fit on the same $n$ observations
   and satisfying only the product-rate condition is not in general covered
   by this argument"* — falsifiable and observable: a no-splitting, non-tree
   propensity arm satisfying only the product-rate condition should show
   inference degradation the tree-based flagship does not. No such arm
   exists in any current study. **GAP** — new, not previously registered.
4. The `hahn1998` propensity-informativeness point (second paragraph of the
   remark) is a citation, not an empirical claim. **N/A — theoretical.**

**Fix applied 2026-09-24:** `claims.md` Table 1's single `rem:two-trees-why`
row split into these four; see the diff for exact wording.

---

## Finding 3 — R4 (`ass:causal` stress): real undercoverage, wrong initial framing, unresolved pre-registered open item

**What's true.** Plain `doubletree`'s coverage under R4 (DGP-C, near-boundary
propensity, worst-cell $e_0=0.982$): **89.7%** at $n=500$ (empirical MC SE
0.96pp $\to$ ~5.5 SEs below the 95% nominal target), 91.3% at $n=1000$, 93.7%
at $n=2000$ (~1.7 SEs below) — a real, only partially-resolving-with-$n$
shortfall, not Monte Carlo noise at the smaller $n$'s.

**What the auditor got wrong on the first pass, and `oracle` corrected.**
`ass:causal` is **not violated** in R4 — DGP-C's worst $e_0=0.982<1$, so the
uniform bound $c>0$ still holds throughout; this is a **small-$c$ regime**,
not a positivity failure. Calling it an "`ass:causal` violation" (the
auditor's first draft) overstates what the DGP does.

The reason R4 has no anchor arm is *not* (as first drafted) "the anchor's
proof might not go through when `ass:causal` is stressed" — `thm:anchor`'s
hypothesis is a primitive on the anchor pair itself
($\sqrt n(\tilde\theta-\theta_0)\rightsquigarrow N(0,\tilde V)$, etc.), not a
chain resting on `ass:causal`, and the anchor CI provably contains its own
Wald interval (triangle inequality; the harness asserts this every
replication), so it can only fail to help if the *underlying* failure is one
its widening mechanism targets. The anchor widens by
$|\hat\theta_{\text{full}}-\hat\theta_{\text{anchor}}|$, which targets the
**bilinear bias remainder** (`prop:bilinear`/`lem:biasbound`). R4's own
numbers argue this isn't primarily a bias problem: relative bias shrinks fast
(pilot-scale 44.3%$\to$10.8%, $n=500\to1000$) while undercoverage persists —
a dissociation pointing at variance/tail behaviour near a thin propensity
margin ($\tilde V$ contains an $e_0^2/(1-e_0)$ term, inflated ~55$\times$ at
$e_0=0.982$ — `manuscript.tex:171`), which the anchor's bias-targeted
widening is not aimed at.

This also means R4's shortfall is most naturally read as **empirical
confirmation of a limitation the manuscript already discloses**
(`manuscript.tex:786`: guarantees are asymptotic and pointwise in the DGP,
not a promise that a given leaf budget/grid/$n$ falls inside the regime where
the asymptotics are a good guide) — not an undisclosed gap. What *is*
unresolved: `head_to_head_comparison/README.md` §10's own pilot-stage open
item 2 pre-registered *"if 0.900 persists it is a real finding"* — it
persisted at the full-scale nsim=1000 — and nothing in `claims.md` or the
README discharges that open item as of this audit.

**Compounding factor (not incidental).** Of Finding 1's 70 genuinely
mismatched stale-vs-fresh rows, **64 are in R4/DGP-C**. Mechanism: DGP-C's
near-empty control cells (README §5: $P(\text{zero controls})=7.5\%$ at
$n=500$) make the tree-splitting objective tie, and the full sweep ran under
an emergency `optimaltrees` tie-breaking revert (pre-`140d110`, restoring the
old tied-tree-enumeration path, because the intended deterministic
`single_model()` extraction has a live, unresolved glibc heap-corruption bug —
`session_notes/2026-09-23.md`, 11:00 entry). So **R4's coverage number is the
single least re-run-stable number in the study, and it is also the source of
the study's headline undercoverage finding.** The revert caveat currently
lives only in a session-notes file; it did not previously propagate to
`claims.md`.

**Fixes applied 2026-09-24:**
- `claims.md`'s `ass:causal` row reworded: "consequence tested" replaced with
  an explicit small-$c$/undercoverage framing, the persisted README open item
  noted, and the revert-caveat/R4-sensitivity cross-reference added.
- `head_to_head_comparison/code/common.R`'s R4 regime definition: added
  `doubletree_crossfit` and `doubletree_anchor` arms, with a **pre-registered,
  falsifiable prediction** recorded in-line (before any new number exists):
  the anchor will widen but still fall short of nominal, because it targets
  bias and R4's shortfall is variance/tail-driven. A confirmed
  predicted-negative is stronger evidence than the current silence; if the
  anchor instead *does* restore coverage here, that is itself a finding worth
  surfacing (evidence the anchor's reach extends beyond `ass:rate`
  violations), not something to quietly absorb.
- **Re-run status:** code change only. Actually executing R4 at full scale
  (nsim=1000, $n\in\{500,1000,2000\}$, now with the crossfit-heavy anchor arm)
  requires the O2 cluster — confirmed 2026-09-24 that `optimaltrees` cannot
  currently be loaded locally by either path: `pkgload::load_all()` on the
  dev source tree fails to compile ("Could not find tools necessary to
  compile a package"), and the already-installed binary
  (`~/Library/R/arm64/4.5/library/optimaltrees`) fails at `dyn.load()` with
  `mmap(...) errno=1` — the same security-tool-blocks-local-.so signature
  logged earlier today (07:57 entry, `session_notes/2026-09-24.md`) for a
  different binary. **Not yet run.**
  - **Gotcha for whoever runs this on O2:** `run_sweep.R`'s resumability
    check is keyed on `file.exists()` for each `(regime, dgp, n)` cell's
    output path only — it has no way to know the *arm set* inside an
    existing checkpoint changed. Resubmitting R4 as-is will silently SKIP
    all three cells, because `cell_R4_C_n{500,1000,2000}_r1000.rds` already
    exist from the pre-change run. The three existing R4 checkpoint files
    must be removed or moved aside on O2 **before** resubmitting, or nothing
    will happen. Exact procedure:
    ```sh
    # on O2, from the doubletree package root
    git pull <github-remote> main   # never bare `git pull` -- see 2026-09-23 note
    mkdir -p simulations/head_to_head_comparison/results/archive/pre-anchor-arm_2026-09-24
    mv simulations/head_to_head_comparison/results/cell_R4_C_n{500,1000,2000}_r1000.rds \
       simulations/head_to_head_comparison/results/archive/pre-anchor-arm_2026-09-24/
    bash simulations/head_to_head_comparison/slurm/launch_subset.sh R4
    # then, once complete: check_progress.sh, combine_results.R, rsync back
    ```
  - Ideally this happens *after* the upstream `single_model()` heap-corruption
    bug is fixed (not under the emergency revert), given the concentration
    finding above — otherwise the new anchor/crossfit arm's numbers inherit
    the same tie-breaking exposure the flagship arm already has in this cell.

---

## Finding 4 — a session-notes number doesn't match the data (fixed by correction, not code)

`session_notes/2026-09-23.md` states doubletree is "130-1000x faster" than
its cross-fit fallback (R2/C3). Direct computation of `secs_mean` ratios
across $n\in\{200,500,1000,2000,4000,8000\}$: approximately 136x, 149x, 163x,
163x, 156x, 147x — a narrow band, no trend toward 1000x anywhere. Likely a
transcription/misremembering error, not a simulation defect. No code change;
flagged so the correct range (~136x-163x) is what gets carried into any
future manuscript §4 draft, not "130-1000x."

---

## Priority ranking (per `oracle`)

1. Finding 3 (R4) — the only one that could put a wrong number in the paper.
2. Finding 2 (`rem:two-trees-why`) — cheap fix, would be an easy referee catch.
3. Finding 1 (stale checkpoints) — real hygiene defect, verified immaterial numerically.
4. Finding 4 (typo) — cosmetic.

## Not yet done (tracked, not silently dropped)

- The R4 anchor+crossfit re-run itself (needs O2; see Finding 3).
- `glm_main`'s R1 relative bias plateaus at 25-30% through $n=4000$ and never
  approaches the closed-form asymptotic prediction of 45% — the pilot's
  "not a discrepancy to chase, full sweep will resolve it" dismissal does not
  look right against the full-scale data. Not investigated further this pass;
  flagged for a follow-up.
- R2's two-factor confound (splitting **and** `max_depth`) limits what C3 can
  ever license about *why* not-splitting works, independent of Finding 2's
  registry fix; a matched-nuisance-class R2 variant would be needed to close
  this properly, and is not implemented.
