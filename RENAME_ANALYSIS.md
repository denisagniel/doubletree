# Rename Analysis: dmltree → doubletree

**Created:** 2026-03-18
**Status:** Scoping analysis

## Summary

Renaming `dmltree` to `doubletree` would require changes in **62 files** across the package, with impacts ranging from critical package infrastructure to documentation and simulations. Based on your experience with the `optimaltrees` rename (formerly `treefarmr`), this analysis categorizes changes by criticality and identifies potential bug sources.

---

## Criticality Tiers

### 🔴 **TIER 1: Critical (Package Will Break Without These)**

These must be perfect or the package won't load/function:

1. **`DESCRIPTION`** (line 1)
   - Package name: `dmltree` → `doubletree`

2. **`tests/testthat.R`** (lines 10, 12)
   - `library(dmltree)` → `library(doubletree)`
   - `test_check("dmltree")` → `test_check("doubletree")`

3. **Git remote URL** (`.git/config` line 7)
   - Current: `git@code.rand.org:agniel-projects/methods_projects/dmltree.git`
   - Need to update to: `doubletree.git` (or keep old URL if not renaming repo)

4. **Bibliography filename** (`inst/paper/dmltree-refs.bib`)
   - Rename to: `doubletree-refs.bib`
   - Update references in manuscript: `\bibliography{dmltree-refs}` → `\bibliography{doubletree-refs}`

### 🟡 **TIER 2: High Priority (User-Facing)**

These affect how users interact with the package:

5. **README.md** (6 occurrences)
   - Title, installation examples, library calls
   - Lines: 1, 3, 10, 59, 72

6. **R documentation** (`man/*.Rd` files)
   - Package references in descriptions (auto-generated, will update with `devtools::document()`)

7. **Manuscript** (`inst/paper/manuscript.tex`, 3+ occurrences)
   - Line 43: Software package mention
   - Line 584: `\texttt{dmltree}` code reference
   - Line 1364: Bibliography reference
   - Also in backup files: `manuscript-before-merge.tex`, `manuscript.tex.backup`

### 🟢 **TIER 3: Medium Priority (Development & Simulations)**

Won't break the package but affects development workflow:

8. **Simulation scripts** (~30 files in `simulations/`)
   - Function calls: `dmltree::dml_att()` → `doubletree::dml_att()`
   - Comments referencing the package
   - Example files:
     - `simulations/exploration/run_simulations.R` (3 occurrences)
     - `simulations/production/methods/*.R` (comparison methods)
     - `simulations/diagnostics/*.R`

9. **Simulation documentation** (`simulations/docs/*.md`, `simulations/production/README.md`)
   - Package references in quickstart guides
   - Method comparison descriptions

10. **Session notes** (`session_notes/*.md`)
    - Historical references (5 files)
    - Can update for consistency but not critical

---

## File-by-File Breakdown

### Core Package Files (7 files)
```
✓ DESCRIPTION                    [Line 1: Package name]
✓ tests/testthat.R              [Lines 10, 12: library + test_check]
✓ README.md                     [6 references throughout]
✓ inst/paper/dmltree-refs.bib        [Filename + manuscript \bibliography{}]
✓ inst/paper/manuscript.tex          [3+ textual/code references]
✓ inst/paper/manuscript-before-merge.tex [3+ references]
✓ inst/paper/manuscript.tex.backup   [3+ references]
```

### Simulation Files (30+ files)
```
simulations/production/launch_production_sims_v2.sh
simulations/production/methods/method_linear_dml.R
simulations/production/methods/method_forest_dml.R
simulations/production/launch_production_sims.sh
simulations/production/BATCH_SIMULATION_README.md
simulations/production/run_primary.R
simulations/production/test_*.R (15 files)
simulations/production/profile_*.R (5 files)
simulations/production/debug_*.R (3 files)
simulations/production/README.md
simulations/docs/*.md (4 files)
simulations/diagnostics/*.R (2 files)
simulations/exploration/*.R (3 files)
```

### Documentation & Session Notes (10 files)
```
session_notes/2026-03-13.md
session_notes/2026-03-06.md
session_notes/2026-03-03.md
session_notes/2025-02-12.md
session_notes/2025-02-05.md
simulations/docs/LOG_BLOAT_PREVENTION_SUMMARY.md
simulations/docs/LOGGING_PROTOCOL.md
simulations/docs/BATCH_SIMULATION_README.md
simulations/docs/MEMORY_SAFE_SIMULATIONS.md
simulations/docs/QUICKSTART.md
```

---

## Potential Bug Sources (Based on optimaltrees Experience)

### 1. **Function Calls in Scripts**
**Risk:** Calls like `dmltree::dml_att()` scattered across 30+ simulation files.
**Impact:** Scripts will fail at runtime when they can't find the package.
**Detection:** Not caught until you run the script.
**Mitigation:** Systematic grep + replace; test key simulation scripts.

### 2. **Package Loading in .Rprofile**
**Risk:** `.Rprofile` doesn't directly reference `dmltree`, but simulation scripts may expect it loaded.
**Impact:** Scripts fail with "object not found" errors.
**Detection:** Only when running simulations in fresh R sessions.
**Mitigation:** Test with `R --vanilla` after rename.

### 3. **Git Remote URL**
**Risk:** If you rename the git repository on the remote server, local `.git/config` becomes stale.
**Impact:** Cannot push/pull until fixed.
**Detection:** First git push attempt.
**Mitigation:**
  - Option A: Keep remote repo as `dmltree.git` (package name != repo name is fine)
  - Option B: Rename remote repo + update `.git/config` + inform collaborators

### 4. **Bibliography References**
**Risk:** Manuscript has `\bibliography{dmltree-refs}` hardcoded.
**Impact:** LaTeX compilation fails (can't find .bib file).
**Detection:** `/compile-latex` will error.
**Mitigation:** Rename file + update all `\bibliography{}` commands.

### 5. **Roxygen-Generated Documentation**
**Risk:** `man/*.Rd` files may contain package name in descriptions.
**Impact:** Documentation mentions wrong package name.
**Detection:** `R CMD check` warnings; user confusion.
**Mitigation:** Run `devtools::document()` after rename (auto-regenerates from roxygen comments).

### 6. **NAMESPACE Auto-Generation**
**Risk:** `NAMESPACE` is hand-written, but typically auto-generated. After rename, need to verify exports.
**Impact:** Functions may not export correctly.
**Detection:** `devtools::check()` or loading package.
**Mitigation:** Run `devtools::document()` after rename.

### 7. **Cross-References in Comments**
**Risk:** Comments in R code or simulations may say "see dmltree::function()".
**Impact:** Confusing but non-breaking.
**Detection:** Manual review or systematic grep.
**Mitigation:** Search and replace in comments.

---

## Recommended Workflow

### Phase 1: Critical Infrastructure (MUST BE PERFECT)
1. ✅ **Backup everything** (commit current state)
2. ✅ **DESCRIPTION**: Change package name
3. ✅ **tests/testthat.R**: Update library() and test_check()
4. ✅ **Rename**: `inst/paper/dmltree-refs.bib` → `inst/paper/doubletree-refs.bib`
5. ✅ **Manuscript**: Update all `\bibliography{dmltree-refs}` references
6. ✅ **Run**: `devtools::document()` (regenerate man pages + NAMESPACE)
7. ✅ **Test**: `devtools::load_all()` - package should load
8. ✅ **Test**: `devtools::test()` - tests should pass

### Phase 2: User-Facing Documentation
9. ✅ **README.md**: Update all package references (6 locations)
10. ✅ **Manuscript**: Update textual references to package name (3+ locations)
11. ✅ **Verify**: `/compile-latex` compiles successfully

### Phase 3: Simulations & Scripts
12. ✅ **Systematic search/replace**: `dmltree::` → `doubletree::` in all `.R` files
13. ✅ **Simulation docs**: Update references in `simulations/docs/*.md`
14. ✅ **Test**: Run 2-3 key simulation scripts to verify they work
15. ✅ **Session notes**: Update for consistency (optional but clean)

### Phase 4: Git & Infrastructure
16. **Decision point**: Rename git remote repository or keep as `dmltree.git`?
    - If keeping: No changes needed
    - If renaming: Update `.git/config` + coordinate with collaborators
17. ✅ **Commit**: All changes with clear message
18. ✅ **Test fresh install**: `devtools::install()` + `library(doubletree)` in fresh session

---

## Estimated Effort

| Phase | Files | Effort | Risk |
|-------|-------|--------|------|
| Phase 1: Critical | 7 files | 15-20 min | **HIGH** - breaks if wrong |
| Phase 2: Documentation | 5 files | 10-15 min | Medium - user-facing |
| Phase 3: Simulations | 35+ files | 30-45 min | Low - runtime errors only |
| Phase 4: Git | 1 file | 5-10 min | Medium - coordination |
| **Total** | **~50 files** | **60-90 min** | - |

**Quality check time:** +30-60 min (running tests, compiling LaTeX, testing simulations)

**Total time estimate:** **2-2.5 hours** with systematic verification at each phase.

---

## Comparison to optimaltrees Rename

Based on your mention of "keep discovering bugs" from the `treefarmr` → `optimaltrees` rename:

### What's Similar:
- Package name in DESCRIPTION, tests, README
- Function calls scattered across multiple files
- Git infrastructure considerations

### What's Different (BETTER for doubletree):
- **Smaller codebase**: 62 files vs potentially more in optimaltrees
- **No C++ code**: No compiled code that might cache the old name
- **Clearer structure**: Simulations are in dedicated directory (easier to batch-update)
- **You've learned**: You now know to check for ALL namespace-qualified calls (`pkg::func`)

### What's Different (MORE COMPLEX for doubletree):
- **Manuscript integration**: Paper directly references package name in 3+ places
- **Bibliography file**: Filename must match `\bibliography{}` command exactly
- **More simulation scripts**: 30+ simulation files vs potentially fewer test scripts

---

## Recommended Decision Framework

### ✅ **Rename NOW if:**
- You plan to publish/share the package soon (better to rename before external users)
- "doubletree" better communicates the method (double ML + trees)
- You're okay spending 2-2.5 hours on systematic rename + verification

### ⏸️ **Wait to rename if:**
- Active collaboration where rename would disrupt others' workflows
- Manuscript under review (don't introduce package name changes mid-review)
- Other high-priority bugs/features more urgent
- You want to batch this with other breaking changes

### ❌ **Don't rename if:**
- "dmltree" is already established in papers/talks/GitHub stars
- Backward compatibility with existing user code is critical
- The name is already good enough and change doesn't add value

---

## Next Steps

If you decide to proceed:

1. **Approval**: Confirm you want to rename
2. **Plan Mode**: I can enter plan mode and create detailed rename plan
3. **Execute**: Systematic rename with verification at each phase
4. **Quality Gate**: Run full test suite + compile LaTeX + test key simulations
5. **Commit**: Single atomic commit with all changes

**Questions to answer before proceeding:**
- Do you want to rename the git remote repository too, or just the package?
- Are there active collaborators I should be aware of?
- Any other related changes you want to bundle with this rename?
