# Quick Test Instructions: Validating Code Review Fixes

**Date:** 2026-03-19
**Goal:** Confirm code review fixes (commit 1a16da0) restored proper 95% coverage

---

## Background

**Quick test results (n=800, DGP1, 100 sims):**
- Fixed λ: 96% coverage ✓
- CV λ: 96% coverage ✓

**Conclusion:** Code review fixes likely solved the problem. Testing with fixed λ first (simpler, faster).

---

## Running the Quick Test

### All DGPs in Parallel (~2-4 hours per DGP)

```bash
cd /Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/production

# Launch all 4 DGPs
./run_all_test100.sh
```

**What it does:**
- Runs 100 reps (vs 500 production) across 3 sample sizes (n=400/800/1600)
- Tests 4 methods (tree, rashomon, forest, linear)
- Total per DGP: 100 × 3 × 4 = 1,200 simulations
- Expected time: 2-4 hours per DGP in parallel

### Monitor Progress

```bash
# Watch logs
tail -f logs/dgp1_test100_*.log

# Check all convergence rates
grep "Convergence rate" logs/dgp*_test100_*.log

# Check all coverage results (after completion)
grep -A 20 "Summary Statistics" logs/dgp*_test100_*.log
```

### Stop if Needed

```bash
# Get PIDs from launch output, then:
kill PID1 PID2 PID3 PID4
```

---

## Expected Results

### Success Criteria
- **Convergence:** >95% of simulations converge
- **Coverage:**
  - Fold-specific (tree): 93-95%
  - Rashomon: 90-93% (DGPs 1-3), 70-80% (DGP4)
  - Oracle: ~95% (validates implementation)
  - Forest/Linear: comparison baselines

### If Coverage is Good (≥93%)
✓ **Code review fixes worked!**
- Original 84.7% was due to bugs (silent corruption, bad fallbacks)
- Fixed λ = log(n)/n is appropriate
- Proceed to full production runs (500 reps)

### If Coverage is Still Low (<90%)
✗ **Need further investigation:**
- Check which DGPs/sample sizes have problems
- May need CV after all (or different tree implementations)
- Review nuisance model diagnostics

---

## After Test Completion

### Analyze Results

```bash
cd simulations/production

# Each test creates results in results/dgpX_batch_YYYY-MM-DD/
# Look for files like dgpX_summary.csv

# Quick summary
for f in results/dgp*_batch_*/dgp*_summary.csv; do
  echo "=== $f ==="
  cat "$f"
  echo ""
done
```

### Expected Output Format
```
method,n,n_valid,bias,rmse,coverage
tree,400,100,0.005,0.045,0.94
tree,800,100,0.003,0.035,0.95
tree,1600,100,0.002,0.025,0.95
rashomon,400,100,0.006,0.048,0.92
...
```

---

## Decision Tree

```
Test complete?
├─ Yes
│  ├─ Coverage ≥93% across DGPs 1-3?
│  │  ├─ Yes → ✓ Run full production (500 reps)
│  │  └─ No → Investigate which DGPs/methods fail
│  └─ Convergence <95%?
│     └─ Yes → Check error logs, may have bugs
└─ No → Wait for completion (check logs)
```

---

## Full Production Run (After Successful Test)

If test shows good coverage:

```bash
cd simulations/production

# Launch full batch (500 reps each)
nohup Rscript run_dgp1_batch.R > logs/dgp1_prod_$(date +%Y%m%d).log 2>&1 &
nohup Rscript run_dgp2_batch.R > logs/dgp2_prod_$(date +%Y%m%d).log 2>&1 &
nohup Rscript run_dgp3_batch.R > logs/dgp3_prod_$(date +%Y%m%d).log 2>&1 &
nohup Rscript run_dgp4_batch.R > logs/dgp4_prod_$(date +%Y%m%d).log 2>&1 &

# Expected time: 6-12 hours per DGP
```

---

## Key Differences: Test vs Production

| Aspect | Quick Test | Production |
|--------|------------|------------|
| N_REPS | 100 | 500 |
| Runtime per DGP | 2-4 hours | 6-12 hours |
| Total sims per DGP | 1,200 | 6,000 |
| Purpose | Validate coverage | Final results |

---

## Files Created

- `run_dgp1_test100.R` - Quick test version of DGP1
- `run_dgp2_test100.R` - Quick test version of DGP2
- `run_dgp3_test100.R` - Quick test version of DGP3
- `run_dgp4_test100.R` - Quick test version of DGP4
- `run_all_test100.sh` - Launch all tests in parallel

---

## Commits

1. **1a16da0** - Fix 15 critical code review issues (CRITICAL)
2. **33a84d7** - Update simulations to use CV (exploratory, later reverted)
3. **[current]** - Revert to fixed λ after successful quick test

---

## Contact

Questions or unexpected results? Check:
1. Logs in `simulations/production/logs/`
2. Git history: `git log --oneline`
3. Code review fixes: Review commit 1a16da0 changes
