# Functional Consistency Simulation - O2 Cluster

**Purpose:** Test whether averaged tree approach achieves perfect functional consistency while maintaining valid asymptotic inference.

**Key Research Questions:**
1. Does averaged tree achieve perfect FC (max_diff ≈ 0) for all n, K?
2. Does coverage degrade as n increases (testing if bias = O(1/n_ℓ) breaks √n-consistency)?
3. Does larger K (more cross-fitting) improve coverage?
4. Does standardized bias √n × bias grow for averaged tree but stay constant for standard M-split?

---

## Simulation Design

**Parameter Grid:**
- **n**: {200, 400, 800, 1600, 3200} - test asymptotic behavior
- **K**: {2, 3, 5} - proportion of cross-fit predictions
- **DGP**: {"simple", "complex", "sparse"} - robustness across data patterns
- **method**: {"standard_msplit", "averaged_tree", "pattern_aggregation"} - three-way comparison
- **M**: 10 (fixed) - number of sample splits

**Total Configurations:** 5 × 3 × 3 × 2 = 90
**Replications per Config:** 500
**Total Replications:** 67,500

**Expected Runtime:** ~8-12 hours on O2 short partition (30 min per job, ~135 jobs)

---

## Prerequisites

### Local Machine

1. **Packages installed:**
   ```r
   devtools::install("optimaltrees")
   devtools::install("doubletree")
   ```

2. **Test locally:**
   ```bash
   cd doubletree/simulations/functional_consistency/slurm
   bash quick_test.sh
   ```

3. **Commit and push:**
   ```bash
   git add -A
   git commit -m "Add functional consistency simulation infrastructure"
   git push
   ```

### O2 Cluster

1. **SSH access:** `ssh username@o2.hms.harvard.edu`
2. **R modules available:** `module avail R` (need R/4.2.1 or later)
3. **Git access:** Repository accessible from O2

---

## Deployment Steps

### 1. Connect to O2

```bash
ssh username@o2.hms.harvard.edu
```

### 2. Pull Repository

```bash
cd ~/global-scholars
git pull
```

### 3. Install Packages (if first time)

```bash
module load gcc/9.2.0 R/4.2.1

# Install from local source
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree
```

Verify installation:
```bash
Rscript -e "library(optimaltrees); library(doubletree); cat('Packages loaded successfully\n')"
```

### 4. Test on O2

```bash
cd doubletree/simulations/functional_consistency/slurm
bash quick_test.sh
```

Expected output: 3 successful replications with summary statistics

### 5. Launch Full Simulation

```bash
bash launch_all_simulations.sh
```

This submits 90 SLURM jobs (one per configuration).

---

## Monitoring

### Check Progress

```bash
bash check_progress.sh
```

Shows:
- Running/pending jobs
- Completed replications
- Progress percentage

Run periodically (e.g., every 30 minutes).

### Check Specific Job

```bash
# View job status
squeue -u $USER -n fc_sim

# View job output
cat logs/fc_sim_JOBID_TASKID.out

# View job errors
cat logs/fc_sim_JOBID_TASKID.err
```

### Cancel Jobs (if needed)

```bash
# Cancel all fc_sim jobs
scancel -u $USER -n fc_sim

# Cancel specific job
scancel JOBID
```

---

## After Completion

### 1. Combine Results

```bash
cd doubletree/simulations/functional_consistency/slurm
Rscript combine_results.R
```

Creates:
- `results/combined_fc_simulations.rds` (R format)
- `results/combined_fc_simulations.csv` (CSV format)

### 2. Transfer Results Back

**Option A: Git (if results are small enough)**
```bash
git add results/combined_fc_simulations.csv
git commit -m "Add functional consistency simulation results"
git push
```

**Option B: SCP (for larger files)**
```bash
# On local machine
scp username@o2.hms.harvard.edu:~/global-scholars/doubletree/simulations/functional_consistency/results/combined_fc_simulations.rds ./
```

### 3. Analyze Results

On local machine:
```r
results <- readRDS("combined_fc_simulations.rds")

# Key analyses:
# 1. Coverage vs n by method
# 2. Functional consistency by method
# 3. Standardized bias √n × bias vs n
# 4. Effect of K on coverage
```

---

## File Structure

```
functional_consistency/
├── run_fc_simulation.R           # Main simulation function
├── slurm/
│   ├── run_single_replication.R  # CLI for single rep
│   ├── run_simulations.slurm     # SLURM batch script
│   ├── launch_all_simulations.sh # Submit all configs
│   ├── quick_test.sh             # Local test (3 reps)
│   ├── check_progress.sh         # Monitor progress
│   ├── combine_results.R         # Aggregate results
│   ├── README_O2.md              # This file
│   └── logs/                     # SLURM output logs
└── results/                      # Output directory
    ├── [individual .rds files]
    ├── combined_fc_simulations.rds
    └── combined_fc_simulations.csv
```

---

## SLURM Configuration

**Resource Allocation:**
- Memory: 4G per task
- Time: 30 minutes per task
- Partition: short
- Array: 1-50 (50 tasks × 10 reps = 500 reps per config)

**Modules:**
- gcc/9.2.0
- R/4.2.1

**Scratch Directory:** `/n/scratch/users/${USER:0:1}/${USER}/fc_sim`

Results are copied from scratch to permanent storage after each task completes.

---

## Troubleshooting

### Problem: Jobs fail with "package not found"

**Solution:** Install packages on O2:
```bash
module load gcc/9.2.0 R/4.2.1
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree
```

### Problem: Jobs fail with "object not found" errors

**Solution:** Check that:
1. Repository is up to date: `git pull`
2. All simulation files are present
3. File paths in SLURM script are correct

### Problem: Jobs stuck in pending

**Solution:**
- Check partition limits: `squeue -u $USER`
- Consider switching to medium or long partition for large jobs
- Check O2 status: https://rc.hms.harvard.edu/

### Problem: Results files missing

**Solution:**
- Check SLURM logs: `cat slurm/logs/*.err`
- Verify scratch directory: `ls /n/scratch/users/${USER:0:1}/${USER}/fc_sim`
- Re-run failed configurations

---

## Expected Results

### Hypotheses

**H1:** Averaged tree achieves max_diff ≈ 0 (machine precision) for all n, K
**H2:** Coverage for averaged tree degrades as n increases (if β ≤ 1/2)
**H3:** Larger K improves coverage for averaged tree
**H4:** Standardized bias √n × bias grows with n for averaged tree, constant for standard M-split

### Key Metrics

1. **Functional consistency:** `max(max_diff_e, max_diff_m0)`
2. **Coverage:** Proportion of 95% CIs containing truth
3. **Standardized bias:** `sqrt(n) * (att_est - att_true)`
4. **Bias/SE ratio:** `|bias| / se` (should be << 1.96 for valid CIs)

---

## Contact

For issues or questions:
- Check O2 documentation: https://rc.hms.harvard.edu/
- Submit SLURM support ticket
- Review research constitution: `meta-spec/RESEARCH_CONSTITUTION.md`

---

## Version History

- **2026-04-13:** Initial simulation infrastructure created
