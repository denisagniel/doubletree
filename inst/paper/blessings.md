# Blessings — `doubletree`

**Instance path:** `inst/paper/blessings.md`
**Governed by:** `.claude/rules/paper-sequencing-gate.md` HARD GATE 4 (author blessing),
`.claude/rules/paper-protocol.md` "Never block on style" (this file is not a style gate — see
that section's "What this rule does not cover")
**Status:** live, seeded 2026-09-25 (registry created; no author blessing recorded yet)

> **Why this file exists.** Resuming collaborative drafting, the author set two standing
> requirements: (1) no section, subsection, or paragraph is final/stable/settled/done until the
> author has said so explicitly; (2) every drafted paragraph is audited for AI-speak/tone
> (`audit-drafted-prose`) before the author ever sees it. This file is the persistent record for
> requirement (1) — it survives across sessions, so a future agent cannot report something as
> "done" just because a prior agent drafted it, ran a reviewer, or the author said "looks good,
> keep going" about something else.

**Default is `unblessed`.** A region absent from this table is unblessed. Only the author's
explicit statement creates or upgrades a row. An agent records a blessing; it **never** infers
one. Silence is not blessing. A merged commit is not blessing. A clean `domain-reviewer` /
`proofreader` report is not blessing. Prior `domain-reviewer`/`proofreader` audit passes
recorded in `quality_reports/reviews/` (e.g. the 2026-09-21 and 2026-09-24 reports) verified
correctness and fidelity — they are **not** author sign-off, and are not treated as such here.

**Three states only:**

| State | Meaning |
|---|---|
| `unblessed` | Default. No author sign-off recorded, or none ever given. |
| `blessed` | Author explicitly approved this region as final, pinned to a commit sha. |
| `blessed-stale` | Was `blessed`, but the manuscript has changed inside that region's line range since the pinned sha — treated as `unblessed` until re-blessed. |

**Mechanical check** (run by `edit-paper-with-context` before editing, and by any agent before
reporting a region as "done"): resolve the target passage to its enclosing `\label` (fall back
to the verbatim `\section`/`\subsection` title when unlabelled) → look up the row → if `blessed`,
run `git diff --unified=0 <sha>..HEAD -- inst/paper/manuscript.tex` and check whether any changed
hunk falls inside the region's current line range → if yes, downgrade to `blessed-stale` in this
same turn and treat as unblessed.

---

## Regions

Seeded from `manuscript.tex`'s current section/subsection structure as of commit `796dd26`
(2026-09-21, last commit to touch the file). Rows below are `unblessed` by default; any
exception is recorded explicitly with the author's own sign-off, per the contract above.

| Region (`\label` or heading) | § | Status | Date | Sha | Note |
|---|---|---|---|---|---|
| "Abstract" (no `\label`, precedes §1) | — | **blessed** | 2026-09-25 | `e03a868` | Author-approved after one revision round (fixed "audible"→"auditable" typo; restored a short point-estimate-stays-auditable clause; knowingly dropped the old pointwise-not-uniform/no-ML-ensemble-flexibility scope sentence). Committed in `e03a868`. |
| "Introduction" (no `\label`) | 1 | **blessed** | 2026-09-25 | `e03a868` | Author-approved after a full-section (not paragraph-by-paragraph) review: duplicate `\citep{chernozhukov2018}` removed (`:38`); AI-speak/correctness audit found no blocking issues, two low-severity notes accepted as-is (borderline "Furthermore" transition; `:42`'s "no sample splitting" omitting "or cross-fitting" vs. the abstract's fuller phrasing). Committed in `e03a868`. |
| "Data, estimand, and notation" (no own `\label`; falls under `\label{sec:setup}`) | 2.1 | **blessed** | 2026-09-25 | `e03a868` | Author-blessed directly (author's own review, not a full agent audit pass this session — recorded per "author records, never infers" §2.4 of this file's contract: the author's explicit statement is sufficient on its own). Committed in `e03a868`. |
| `sec:trees` | 2.2 | **blessed** | 2026-09-25 | `e03a868` | Author-blessed after an extensive edit history this session: full audit (AI-speak clean; found 3 correctness issues — bold-vector inconsistency, unstated `\Pi^\nu_\tau` reference measure, `rem:saturated-general`'s binary-only leaf-cost formula); bold-vectors fix (document-wide); `p_{\bx}` removed, `\Pi^\nu_\tau` realigned to `theory.tex`'s measure convention (author-confirmed intent); `rem:saturated-general` broadened beyond GLM then independently `domain-reviewer`-audited (`quality_reports/reviews/2026-09-25_sec2.2-saturated-general-broadening.md`, 6 blocking findings, all fixed same day); real citations added for the discretization-motivation sentence (`pope2004hcc`, `wilson1998framingham`, `breyer1988sentencing`, verified via `librarian`); "axis-aligned splits" explained at first mention. `claims.md`, `notation.md`, `story.md` (S8 + new reconciliation-log entry) all kept consistent. Committed in `e03a868`. |
| "Identification and influence function" (no own `\label`; falls under `sec:setup`) | 2.3 | **blessed** | 2026-09-25 | `e03a868` | Author-blessed directly (author's own review, no additional agent audit pass requested for this region). Content unchanged this session; the only edit was the document-wide bold-vector pass touching the positivity-bound display and its surrounding prose — mechanical only, no content change. Committed in `e03a868`. |
| "The doubletree algorithm" (no own `\label`; falls under `sec:setup`) | 2.4 | unblessed | — | — | **Audited and fixed 2026-09-25, this session, uncommitted — awaiting your review, not yet blessed.** Two-pass audit (`domain-reviewer` for correctness, a paired proofreader pass for tone) found 3 blocking correctness defects, all fixed same day: (i) `eq:select` defined the selected partition as `\that_j`, a macro rendering as $\hat t_j$ (hat over Latin *t*), disconnected from every downstream `\tauhat_j` use — fixed by replacing `\that_j`→`\tauhat_j` throughout; (ii) the Clipping assumption (`ass:construct`) was unscoped, literally clipping the unclipped-by-design leaf refit used in structure selection and leaving $\thetahat$ technically undefined on an all-treated leaf — fixed by adopting `theory.tex:235`'s own overloading convention ($\hat e_{\tauhat_e}$ denotes the clipped value only at the selected partition); (iii) $R^{(\mu)}_n$ was on the control-conditional scale ($\div n_0$), which `theory.tex` explicitly rules out — fixed to the unnormalised control-weighted form matching §2.1's $\|\cdot\|_{2,w}$. Two MAJOR overclaims also narrowed: the unlabeled remark's unsupported "more efficient use of...data" claim (replaced with what the construction actually licenses), and `rem:two-trees-why`'s closing sentence (narrowed to name exactly `rem:single-tree-se`/`cor:single-tree-saturated`, with the matching overclaim at §3.5's own opening, line ~684, corrected to match). Full audit: `quality_reports/reviews/2026-09-25_sec2.4-algorithm-audit.md`. Recompiled clean (4-pass XeLaTeX, 0 undefined refs/citations). **Known open item, not fixed this pass:** `ass:construct`'s label covers only Clipping, not `theory.tex`'s full four-part bundle its downstream citations imply — flagged in `notation.md`, deferred. `claims.md`, `notation.md`, `story.md` all updated. |
| `sec:partition-recovery` | 3.1 | unblessed | — | — | |
| `sec:instantiation1` | 3.2 | unblessed | — | — | |
| `sec:honest-manuscript` | 3.3 | unblessed | — | — | |
| `sec:finite-discussion` | 3.4 | unblessed | — | — | |
| `sec:single-tree` | 3.5 | unblessed | — | — | |
| `sec:lasso-comparison` | 3.6 | unblessed | — | — | Confirmed complete per `outline.md`, 2026-09-17 — that note is a drafting-completeness claim, not a blessing. |
| Simulation evidence *(not yet drafted)* | 4 | — | — | — | No row: nothing exists yet to bless or leave unblessed. Add a row when §4 is first drafted. |
| `sec:discussion` | 5 | unblessed | — | — | Items (i)–(iii) audited 2026-09-21 (`quality_reports/reviews/2026-09-21_discussion-completion.md`). **New draft addition, 2026-09-25 (this session, uncommitted, NOT yet author-reviewed):** a sentence restoring the "no claim of matching an unrestricted ML ensemble's flexibility" disclosure (`:788`, in the "What this paper has not done" paragraph) — this claim existed only in the pre-edit abstract and had disappeared from the entire manuscript once the abstract was trimmed. Self-audited for AI-speak before this handoff (categories 1–6, none flagged); not yet blessed — still `unblessed` until you review this specific addition. |
| `app:ate` | Appendix | unblessed | — | — | Audited 2026-09-21 (`quality_reports/reviews/2026-09-21_ate-appendix.md`). |

---

## Maintenance contract

- **A drafting or editing skill creates or confirms a row for the region it touches, in the same
  turn** — same discipline as `notation.md`/`claims.md`'s registry-delta requirement.
- **Only the author's own words upgrade a row to `blessed`.** "Looks good," "keep going," or the
  absence of objection is not sufficient — the author must say something equivalent to "this
  section is final" or "I bless this" for the specific region named.
- **Unblessed prose is labelled `DRAFT` in every handoff.** Never report an unblessed region as
  done, final, settled, or stable, no matter how many reviewer passes it has been through.
- **Blessed prose is never silently rewritten.** An edit to a `blessed` region requires naming
  the intended change and asking first, before editing.
- **Do not narrate blessing changes in the manuscript** (`paper-protocol.md` §1) — this table is
  the record; the prose is the argument.
