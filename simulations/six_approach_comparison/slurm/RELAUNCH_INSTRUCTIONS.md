# Relaunch Instructions (2026-05-29)

## Context

Fixed S7 type error in `estimate_att_doubletree_averaged()` (approach 4). Need to rerun approach 4 jobs with updated package.

## Step 1: Update Package on Cluster

```bash
cd ~/doubletree
git pull
R CMD INSTALL .
```

## Step 2: Choose Relaunch Strategy

### Option A: Relaunch Only Approach 4 (Recommended)

Fastest option - only reruns the 12 failed jobs (approach 4 = doubletree_averaged).

```bash
cd simulations/six_approach_comparison

# Clean up old approach 4 results/logs
./slurm/cleanup_approach4.sh

# Relaunch approach 4 jobs (13-24)
sbatch slurm/relaunch_approach4.sh
```

**Jobs:** 12 jobs, ~3 hours max runtime
**Output:** Overwrites `results/raw/fast_approach_{13..24}.rds`

### Option B: Relaunch Everything

If you want to ensure all results use the latest package version.

```bash
cd simulations/six_approach_comparison

# Clean up all results/logs (interactive prompt)
./slurm/cleanup_all.sh

# Relaunch all three job arrays
sbatch slurm/run_fast_approaches.sh      # 36 jobs, 3h max
sbatch slurm/run_medium_approaches.sh    # 24 jobs, 4h max
sbatch slurm/run_msplit_approach.sh      # 60 jobs, 4h max
```

**Jobs:** 120 jobs total, ~4 hours max runtime
**Output:** Fresh `results/raw/*.rds` files

## Step 3: Monitor Progress

```bash
# Check queue status
squeue -u $USER

# Check specific job
squeue -j <job_id>

# Monitor logs in real-time
tail -f logs/approach4_rerun_13.out  # for approach 4 rerun
tail -f logs/fast_13.out             # for full relaunch
```

## Step 4: Check Results

After jobs complete:

```bash
cd simulations/six_approach_comparison

# Check for failures
grep -l "Exit code: [^0]" logs/*.out

# Count successful results
ls results/raw/*.rds | wc -l  # Should be 120

# Combine results
Rscript code/combine_results.R
```

## Job Mappings

### Fast Approaches (run_fast_approaches.sh)
- Jobs 1-12: Approach 1 (full_sample)
- **Jobs 13-24: Approach 4 (doubletree_averaged)** ← Fixed bug here
- Jobs 25-36: Approach 6 (msplit_singlefit)

### Medium Approaches (run_medium_approaches.sh)
- Jobs 1-12: Approach 2 (crossfit_separate)
- Jobs 13-24: Approach 3 (doubletree)

### M-Split Approach (run_msplit_approach.sh)
- Jobs 1-60: Approach 5 (msplit) with 5 batches per configuration

## What Was Fixed

**Bug:** Lines 758-759 in `R/estimate_att_averaged.R` called `extract_tree_structure(cf_e@intersecting_trees[[1]])`, but `@intersecting_trees` contains nested list trees, not OptimalTreesModel objects.

**Fix:** Removed those lines, set `structures = NULL` (averaged trees already contain structure info).

**Commit:** ca1070f
