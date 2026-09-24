# Cleanup Manifest: 2026-09-24 12:45:48

**Project:** doubletree
**Archived:** 21 files, 3 category directories

## Summary

- Old versions: 1 file
- Temporary/generated output: 20 files
- Misc (dead scaffolding): 1 file

## Detailed Archive

### Old Versions (1 file)
Moved to: `old_versions/`

- `inst/paper/manuscript-bak.tex` → Stale manual backup (Aug 18, 926 lines) of
  the tracked `inst/paper/manuscript.tex` (currently 866 lines, substantially
  restructured since — see commits `939daef`, `896a405`, `796dd26`,
  `b9142bb`, `c478733`). Superseded by git history; kept here only in case
  anyone wants a pre-restructuring diff without walking commits.

### Temporary/Generated Output (20 files)
Moved to: `temp_files/`

All files timestamped `20260909-*` — the same 2026-09-09 pilot/smoke-test
vintage that `quality_reports/reviews/2026-09-24_s5-claims-audit.md`
identified as stale and archived on the `.rds`-checkpoint side (see that
audit's item 2, "Stale-checkpoint contamination"). These are the pilot's
derived summary/verification tables, regenerable from `analyze.R` against
the now-current (deduplicated) `results/` directory if ever needed again.

**`head_to_head_comparison_tables/`** (16 files):
- `verify_probe_C_20260909-065906.csv`
- `verify_cells_A_20260909-065906.csv`
- `verify_remainder_20260909-065906.csv`
- `verify_confounding_20260909-065906.csv`
- `verify_dgps_20260909-065906.csv`
- `pilot_summary_20260909-070000.csv`, `_070017.csv`, `_070128.csv`, `_071859.csv`
- `pilot_projection_20260909-070000.csv`, `_070017.csv`, `_070128.csv`, `_071859.csv`
- `pilot_rmse_ratios_20260909-070000.csv`, `_070017.csv`, `_070128.csv`

**`single_tree_corollaries_tables/`** (4 files):
- `pilot_gates_20260909-145931.csv`
- `pilot_regimeE_20260909-145931.csv`
- `population_20260909-145931.csv`
- `pilot_summary_20260909-145931.csv`

### Misc — Dead Scaffolding (1 file)
Moved to: `misc/session_logs/`

- `inst/paper/quality_reports/session_logs/README.md` → The directory held
  only this README. Its own text: "historical archive, HISTORICAL, FROZEN
  2026-09-01 ... do not add to them" and points its "full diagnosis" link at
  a *different* project (`missing-data-did/...`). No actual session-log
  files ever lived here in doubletree — this repo's real session record is
  `session_notes/YYYY-MM-DD.md`, which is untouched. Archiving the empty
  scaffold, not any content.

## Not Archived (left in place, on purpose)

- `.omo/` — opencode/omo tooling run-continuation state (session JSON files),
  not project content. Added to `.gitignore` instead of archiving so it stops
  showing up in `git status` going forward.
- `quality_reports/`, `session_notes/` — never archived per this skill's
  rules; two stray uncommitted files/edits from these directories were
  instead **committed** (see commits `6af4dfd`), not moved here, because they
  were real, unlanded work product rather than clutter:
  `quality_reports/plans/2026-09-15_real-data-pipeline-spec.md` (the actual
  PI-confirmed spec, previously untracked) and a restored entry in
  `session_notes/2026-09-17.md` that appears to have been dropped in an
  earlier merge.

## Recovery

To recover any file: copy it from its category directory above back to the
original path shown in the "Moved to" line.

```bash
# Example: restore the pilot tables
cp -r archive/cleanup-2026-09-24-124548/temp_files/head_to_head_comparison_tables/* \
   simulations/head_to_head_comparison/tables/
```

## Recommendations

- Add `simulations/*/tables/` to `.gitignore` (regenerable output, same
  pattern already applied to `simulations/*/results/`) to prevent
  re-accumulation of untracked pilot-era tables.
- `.omo/` added to `.gitignore` this pass — see diff.
