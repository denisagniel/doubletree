# Review log — Discussion (§5) completion, 2026-09-21

**Scope:** `inst/paper/manuscript.tex`, `\section{Discussion}\label{sec:discussion}`. Three items
decided in the 2026-09-18 session (per `outline.md`'s §5 row, `structure-reviewer`'s original
findings) were drafted this session: (i) `prop:spectest`'s unvalidated size/power caveat, (ii)
narrowing the audit-trail-availability claim, and (iii) a new ATT-vs-ATE scope paragraph.

## What changed and why

**Paragraph 1 (audit-trail paragraph).** `structure-reviewer` (2026-09-18 session) found this
paragraph presented `prop:spectest` as an unqualified operational decision rule and claimed the
audit trail is "available whether or not structural sparsity actually holds" without
acknowledging that the interval's width, when sparsity fails, depends on a possibly
non-interpretable anchor. Drafted a revision adding the caveat and narrowing the claim, then
commissioned a `domain-reviewer` audit (session `ses_f3bafbc32ffeJ3KT6f13vDo4u3`) and a
`proofreader` pass (session `ses_f3bafbbb6ffex4obMlBarDmKzD`) against the draft.

`domain-reviewer` found the first draft's fixes did not go far enough and introduced new
problems:
- The caveat named the finite-sample size/power gap but omitted `prop:spectest`'s dependence on
  `ass:rate` (the anchor's fold-wise nuisance rate condition), which `notation.md` classifies as
  not empirically verifiable — the larger of the two limitations.
- Framing the diagnostic as a basis for "choosing which interval to report" implied the selected
  interval carries a coverage guarantee; `prop:spectest` proves only size→0 and consistency under
  an alternative, not coverage for a diagnostic-selected interval — the same post-selection-inference
  hazard `sec:lasso-comparison` cites `leeb2005` against.
- The closing sentence's first branch ("complete for the interval... when structural sparsity
  holds") was self-undercutting: the diagnostic's test statistic $\widehat\delta$ is itself built
  from the anchor, so running the diagnostic pulls the anchor into the audit trail even when
  sparsity does in fact hold — only assuming sparsity outright, without running the diagnostic,
  keeps the interval anchor-free.
- "The anchor... may use $K$-fold cross-fitting... inherits exactly the cost cross-fitting
  imposes" was a modal mismatch (the anchor is only "permitted," not required, to cross-fit) that
  also named the wrong cost when the anchor is a non-cross-fit black-box model instead.

`proofreader` separately found: inconsistent terminology ("bound" vs. the manuscript's
established "width"), a redundant "wider width," an overly dense sentence better split at an
existing colon, and an ambiguous closing clause ("is not in question").

All of the above were fixed in the final version: the caveat now names `ass:rate` explicitly, the
diagnostic is reframed as evidence informing confidence rather than a licensed selection rule
with an explicit non-coverage disclaimer and a stated safe-harbor (report the anchor throughout),
the closing sentence's logic is corrected to state the interval is anchor-free only when sparsity
is assumed outright without running the diagnostic, and the cross-fitting-cost sentence now
states the cost depends on which kind of anchor the analyst actually builds. Terminology and
sentence-length issues fixed to match.

**New paragraph (ATT-vs-ATE scope).** Drafted per the 2026-09-18 decision (ATT stays the primary
estimand; main machinery plausibly generalizes to the ATE but this is unverified; the single-tree
family and `rem:two-trees-why` stay ATT-only per `hahn1998`'s asymmetry). `domain-reviewer` found
a citation-fidelity problem in the first draft: the single-tree corollary family was said to
"rest on" `hahn1998`, but `cor:single-tree` itself is a pure algebraic identity with no
dependence on any efficiency-bound asymmetry — Hahn's result explains why the family is of
*interest* for the ATT, not a premise it uses. Separately, the draft attributed
`rem:two-trees-why`'s *entire* case for a propensity tree to Hahn, when the remark explicitly
separates two points and only the second (propensity-nuisance informativeness) rests on Hahn; the
first (propensity tree as an interpretability, not validity, choice) is not ATT-specific. Both
misattributions fixed. Also added, per `domain-reviewer`'s minor finding, that the ATE requires
both outcome regressions (three trees, not two) — a real structural consequence the scope
decision should not leave implicit. `proofreader` additionally fixed an ambiguous adverb
placement ("exactly lowers" → matched the manuscript's own existing phrasing at `:195`),
standardized "single-outcome-tree corollary family" to the manuscript's established "single-tree"
short form, added missing commas around a parenthetical, and trimmed one redundant restatement.

**Provenance-leakage check (both paragraphs, both audits):** clean — no reference to the
manuscript's own editing/review history, no auditor vocabulary in the prose.

## Registry deltas (same turn)

- `claims.md`, `prop:spectest` row: added a note that the Discussion now discloses this GAP to the
  reader (prose-level disclosure; simulation-evidence status itself unchanged).
- `claims.md`, ATT-vs-ATE scope row: status updated `Decided, not yet drafted` → `Drafted`
  (`manuscript.tex:739`); softened "is asserted to generalize" to match the prose's explicit
  non-claim; corrected the Hahn attribution to match the fixed prose; dropped "the CLT" from the
  machinery list per `domain-reviewer`'s recommendation to reconcile the registry down to the more
  defensible prose rather than the reverse.
- `outline.md` §5 (Discussion row): marked items (i)-(iii) drafted and audited, with a pointer to
  this log. Also corrected a separate, pre-existing error found while updating this row: item
  (iv) (the trees-vs-linear representational trade-off paragraph) was listed as "not yet drafted"
  but was in fact already drafted in the 2026-09-18 commit alongside `rem:saturated-general` —
  fixed the listing, not the manuscript (nothing there needed a change).

## Verification

Label/`\ref`/`\eqref` integrity re-checked after every edit in this batch (Python regex sweep of
`manuscript.tex`): 44 unique labels, 35 unique refs, zero missing. No LaTeX compile run this
session.

## Not done this session

- LaTeX compilation was not verified (per `edit-paper-with-context`'s own checklist, this needs a
  `compile-latex`/`verifier` pass before treating the edit as fully final).
- Item (iv) needed no new drafting (see correction above); no further action pending on it.
- The manuscript's Section 2 gap found in parallel this session (`ass:sparsity`, decided
  2026-09-17/18, never actually drafted into `sec:trees` despite `outline.md` claiming the whole
  restructuring pass was "EXECUTED") is unrelated to this log's scope and is tracked separately.
