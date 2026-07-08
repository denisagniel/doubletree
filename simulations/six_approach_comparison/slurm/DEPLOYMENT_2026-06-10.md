# Deployment Instructions — 2026-06-10

## What Changed

**doubletree** commit `34e5d68`:
- `estimate_att_msplit.R`: Added `discretize_bins="adaptive", discretize_method="quantiles"`
  to Stage 1 CV calls (was defaulting to `discretize_bins=2` = median split)
- `estimate_att_averaged.R`: Same fix; plus Stage 2 now stores cross-fitted predictions
  in `predictions_all_splits` return field

**optimaltrees** commit `a2b6321`:
- Refactored `cross_fitted_rashomon.R` (simplified, ~317 line reduction)
- Consolidated auto-tune regularization binary search
- All 1043 tests pass

**Root cause of fixed failures:**
- DGP4 approach 6 "argument is of length zero" (0/500 success) → stump modal structures
  from `discretize_bins=2` on continuous features
- DGP3 approach 6 propensity bound failures (70% failure) → same cause

**Validation:** 10/10 DGP4 success after fix; DGP1 100-rep bias unchanged (+0.003).

---

## Step 1: Push Changes to GitHub

On your local machine:
```bash
cd ~/RAND/rprojects/global-scholars/optimaltrees
git push

cd ~/RAND/rprojects/global-scholars/doubletree
git push
```

---

## Step 2: Update Packages on Cluster

SSH to O2, then:
```bash
# Pull both packages
cd ~/optimaltrees
git pull

cd ~/doubletree
git pull

# Reinstall both (order matters: optimaltrees first, doubletree depends on it)
module load gcc/14.2.0 R/4.4.2
R CMD INSTALL ~/optimaltrees
R CMD INSTALL ~/doubletree
```

Verify install:
```bash
Rscript -e "library(optimaltrees); library(doubletree); cat('OK\n')"
```

---

## Step 3: Quick Smoke Test on Cluster

Before submitting 600 jobs, run a single rep to confirm packages work:
```bash
cd ~/doubletree/simulations/six_approach_comparison

Rscript code/run_single_replication.R \
  --approach 6 --dgp 4 --n 500 \
  --rep_start 1 --rep_end 2 \
  --output /tmp/test_a6_dgp4.rds

# Should complete without errors. Check output:
Rscript -e "r <- readRDS('/tmp/test_a6_dgp4.rds'); cat('n_success:', sum(!is.na(r\$theta)), '\n')"
```

---

## Step 4: Clean Up Old Approach 6 Results

The previous fast-array run (jobs 25-36 of run_fast_approaches.sh) produced
partial/failed results for approach 6. Clean before resubmitting:
```bash
cd ~/doubletree/simulations/six_approach_comparison
rm -f results/raw/approach6_*.rds
rm -f logs/approach6_*.out logs/approach6_*.err
```

---

## Step 5: Submit Approach 6 Full Array

```bash
cd ~/doubletree/simulations/six_approach_comparison
sbatch slurm/run_approach6.sh
```

**Job details:**
- Array: 1-600 (600 jobs)
- Config: 12 configs (4 DGPs × 3 n values) × 50 batches
- Reps: 20 per batch = 1000 reps per config
- Time: 6h limit, 16G memory
- Output: `results/raw/approach6_<job>.rds`

---

## Step 6: Resubmit Approach 4 DGP3 Large-n (Optional)

Approach 4 DGP3 n=1000 and n=2000 previously timed out. These are slots
config_idx=8 (DGP3, n=1000) and config_idx=9 (DGP3, n=2000) in
`run_approach4.sh`. Per-rep time on complex DGP is ~12-15 min.

Consider increasing time limit before resubmitting:
- Current: 6h for approach 4
- Each rep: ~15 min × 20 reps = 300 min = 5h → tight at 6h
- Recommended: 8h for DGP3 large-n

---

## Step 7: Monitor

```bash
squeue -u $USER
# or
bash slurm/check_progress.sh
```

---

## Step 8: After Completion

```bash
cd ~/doubletree/simulations/six_approach_comparison
Rscript code/combine_results.R
```

---

## Notes on DGP4 Approaches 3/4 Performance

Local smoke test showed approach 3 on DGP4 continuous takes ~2600s/rep and
approach 4 takes ~3500s/rep. This is because adaptive discretization creates
more thresholds → larger Rashomon search space on continuous features.

The existing `run_approach3_dgp4.sh` and `run_approach4_dgp4.sh` have separate
time limits — check those if DGP4 Rashomon jobs are timing out again.
