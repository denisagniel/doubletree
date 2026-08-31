# Approach 4 Auto-Tuning Fix (2026-05-29)

## Problem

Approach 4 (doubletree_averaged) was failing with 499/500 replications at n=500:

```
Rashomon intersection empty for propensity.
Rashomon intersection empty for outcome.
```

**Root cause:** Default `rashomon_bound_multiplier = 0.05` is too tight for small n (500) with complex DGPs. The K=5 fold-specific Rashomon sets don't overlap when epsilon is this small.

**This is not a bug** - the code correctly detects and reports empty intersection. But it means approach 4 doesn't work at these settings.

## Solution

Enable `auto_tune_intersecting = TRUE` in the simulation wrapper.

**What it does:** Automatically increases epsilon (0.05 → 0.10 → 0.15 → 0.20) until a non-empty intersection is found.

**Trade-off:**
- ✓ Approach 4 will succeed at n=500
- ✓ Still gets interpretable single tree
- ⚠️ May use epsilon > 0.05 (includes more suboptimal trees)
- ⚠️ Slightly slower (tries multiple epsilon values)

## Changes Made

**File:** `code/estimators.R`

```r
# Before:
result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y, K = K,
  regularization = regularization,
  outcome_type = "binary",
  verbose = FALSE
)

# After:
result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y, K = K,
  regularization = regularization,
  outcome_type = "binary",
  auto_tune_intersecting = TRUE,  # NEW
  verbose = FALSE
)
```

## Testing

**Local test:** `test_approach4_autotune.R`
- Tests DGP 3 (complex), n=500 - the configuration that was failing
- Verifies auto-tuning finds an intersection
- Checks that result is reasonable

**Run:**
```bash
cd simulations/six_approach_comparison
Rscript test_approach4_autotune.R
```

## Expected Behavior After Fix

**Before (with epsilon=0.05 fixed):**
- Approach 4 succeeds at n=1000, n=2000
- Approach 4 fails at n=500 (empty intersection)
- Success rate at n=500: ~0.2%

**After (with auto-tuning):**
- Approach 4 succeeds at all sample sizes
- May use epsilon=0.10 or 0.15 at n=500
- Success rate at n=500: ~95%+
- Result object will contain final epsilon_n used

## Relaunch Plan

After local test passes:

1. **Commit changes:**
   ```bash
   git add code/estimators.R test_approach4_autotune.R APPROACH4_AUTOTUNE_FIX.md
   git commit -m "Enable auto-tuning for approach 4 to handle empty intersection"
   git push
   ```

2. **Update cluster:**
   ```bash
   cd ~/doubletree
   git pull
   R CMD INSTALL .
   ```

3. **Relaunch approach 4:**
   ```bash
   cd simulations/six_approach_comparison
   ./slurm/cleanup_approach4.sh
   sbatch slurm/relaunch_approach4.sh
   ```

## Alternative Considered

**Option A:** Increase epsilon to 0.10 manually (fixed)
- Simpler, but loses ability to use tighter epsilon at larger n
- Auto-tuning is better: uses 0.05 when possible, increases only if needed

**Option B:** Accept failures at n=500
- Would leave approach 4 incomplete in results
- Missing important comparison at smallest sample size

**Option C (chosen):** Enable auto-tuning
- Best of both worlds: tight epsilon when possible, looser when needed
- Approach 4 works at all sample sizes
- Maintains interpretability (single tree output)

## Documentation for Results

When reporting results, note:
- Approach 4 uses adaptive epsilon selection
- Final epsilon_n may vary by configuration (larger for small n, complex DGP)
- This is theory-compliant: larger epsilon includes more near-optimal trees but maintains Rashomon validity
