# DML-ATT Simulations on O2 Cluster

Complete infrastructure for running the primary DML-ATT simulations on Harvard's O2 cluster using SLURM array jobs.

## Current Status (2026-03-20)

**Working methods:** ✅ forest, linear (fully tested)
**Broken methods:** ❌ tree, rashomon (S7 predict issues from yesterday's migration)

**Recommendation:** Run forest + linear only (9,000 sims, 15-30 min), fix S7 issues later.

```bash
# Quick test (30 seconds)
bash slurm/quick_test.sh

# Launch forest + linear only (9,000 simulations)
bash slurm/launch_forest_linear.sh

# Or launch all methods (18,000 sims, includes broken tree/rashomon)
bash slurm/launch_all_simulations.sh
```

## Overview

**Goal:** Run 18,000 DML-ATT simulations for manuscript Table 1

**Design:**
- 3 DGPs (binary, continuous, moderate)
- 4 methods (tree-DML, rashomon-DML, forest-DML, linear-DML)
- 3 sample sizes (400, 800, 1600)
- 500 replications per configuration
- Total: 3 × 4 × 3 × 500 = 18,000 simulations

**Parallelization:**
- 36 SLURM array jobs (one per configuration)
- Each array job: 500 tasks (one per replication)
- Estimated time: 30-60 minutes (vs. 15 hours sequential)

---

## Setup (One-Time)

### 1. Install R packages on O2

```bash
# SSH to O2
ssh username@o2.hms.harvard.edu

# Load R
module load gcc/14.2.0
module load R/4.4.2

# Set library path
export R_LIBS_USER="${HOME}/R/library"
mkdir -p $R_LIBS_USER

# Install packages
R
```

In R:
```r
install.packages(c("dplyr", "optparse", "ranger"), lib = Sys.getenv("R_LIBS_USER"))
```

### 2. Build optimaltrees package

```bash
# On O2, navigate to package directory
cd /path/to/global-scholars/optimaltrees

# Load R
module load gcc/14.2.0
module load R/4.4.2

# Build and install
R CMD INSTALL .
```

Or use devtools (if available):
```r
devtools::install("/path/to/global-scholars/optimaltrees")
```

### 3. Transfer code to O2

**Using GitHub (recommended):**
```bash
# On local machine: commit and push
cd ~/RAND/rprojects/global-scholars
git add doubletree/simulations/
git commit -m "Add O2 simulation infrastructure"
git push

# On O2: pull latest code
ssh username@o2.hms.harvard.edu
cd /path/to/global-scholars
git pull
```

**Alternative - rsync (if not using git):**
```bash
# Transfer entire doubletree directory
rsync -avz --exclude 'results/' --exclude '.git/' \
  ~/RAND/rprojects/global-scholars/doubletree/ \
  username@o2.hms.harvard.edu:/path/to/global-scholars/doubletree/
```

---

## Running Simulations

### Quick Start

```bash
# SSH to O2
ssh username@o2.hms.harvard.edu

# Navigate to simulation directory
cd /path/to/global-scholars/doubletree/simulations/production

# Create logs directory
mkdir -p logs

# Launch all 36 configurations
bash slurm/launch_all_simulations.sh
```

This will submit 36 SLURM array jobs (18,000 total tasks).

### Launch Single Configuration (Testing)

To test with one configuration before launching all:

```bash
# Test with dgp1, n=400, tree method (500 replications)
sbatch --export=DGP=dgp1,N=400,METHOD=tree slurm/run_dml_simulations.slurm
```

### Test with Small Array (Dry Run)

To test with just 10 replications:

```bash
# Edit slurm script temporarily
sed -i 's/#SBATCH --array=1-500/#SBATCH --array=1-10/' slurm/run_dml_simulations.slurm

# Launch one config
sbatch --export=DGP=dgp1,N=400,METHOD=tree slurm/run_dml_simulations.slurm

# Check results after ~5 minutes
bash slurm/check_progress.sh

# Restore original array size
sed -i 's/#SBATCH --array=1-10/#SBATCH --array=1-500/' slurm/run_dml_simulations.slurm
```

---

## Monitoring

### Check Progress

```bash
# Overall progress (counts completed files)
bash slurm/check_progress.sh
```

Output shows:
- Completed replications by configuration
- Overall progress percentage
- Running SLURM jobs
- Failed tasks (if any)

### Check SLURM Queue

```bash
# View your running jobs
squeue -u $USER

# Watch queue in real-time
watch -n 10 'squeue -u $USER | grep dml_sim'

# Count running tasks
squeue -u $USER -n dml_sim | wc -l
```

### Check Individual Logs

```bash
# View recent output log
tail logs/dml_JOBID_TASKID.out

# View recent error log
tail logs/dml_JOBID_TASKID.err

# Check for failures
grep -l "exit code 1" logs/dml_*.err

# View failed task
cat logs/dml_JOBID_TASKID.err
```

---

## Combining Results

After simulations complete:

```bash
# Combine all individual .rds files into single dataset
Rscript slurm/combine_results.R
```

This creates:
- `results/o2_primary_combined/combined_results_YYYY-MM-DD.rds` (full data)
- `results/o2_primary_combined/combined_results_YYYY-MM-DD.csv` (for inspection)
- `results/o2_primary_combined/summary_stats_YYYY-MM-DD.csv` (performance metrics)
- `results/o2_primary_combined/metadata_YYYY-MM-DD.rds` (run information)

---

## Troubleshooting

### Job Fails Immediately

**Symptom:** Tasks fail within seconds

**Check:**
1. Module loading: `module load gcc/14.2.0 && module load R/4.4.2`
2. R library path: `echo $R_LIBS_USER`
3. Package installation: `R -e "library(optimaltrees)"`
4. File paths in `run_single_replication.R`

**Common fixes:**
- Incorrect relative paths → Use absolute paths or verify working directory
- Missing packages → Install to `$R_LIBS_USER`
- Wrong R version → Try `module load R` (loads default)

### Job Runs But No Output

**Symptom:** Task completes but no .rds file in output directory

**Check:**
1. Output directory permissions: `ls -ld $OUTPUT_DIR`
2. Scratch space quota: `du -sh /n/scratch/users/${USER:0:1}/${USER}`
3. Error logs: `tail logs/dml_*.err`

**Common fixes:**
- Permission denied → Create directory manually: `mkdir -p $OUTPUT_DIR`
- Quota exceeded → Clean up old files: `rm -rf /n/scratch/...`
- R error → Check logs for error messages

### Job Times Out

**Symptom:** Task hits 1-hour time limit

**Solutions:**
1. Increase time limit in SLURM script: `#SBATCH --time=02:00:00`
2. Reduce complexity:
   - Reduce K folds: `--k-folds 3` (default is 5)
   - Use fewer trees for forest: edit `method_forest_dml.R`

### Memory Error

**Symptom:** "Cannot allocate memory" or job killed

**Solutions:**
1. Increase memory: `#SBATCH --mem=8G` (default is 6G)
2. Monitor memory usage: `sacct -j JOBID --format=JobID,MaxRSS`

### Duplicate Replications

**Symptom:** Same replication appears twice

**Cause:** Job restarted or resubmitted

**Solution:**
```bash
# combine_results.R automatically deduplicates
# Or manually remove duplicates:
cd $OUTPUT_DIR
# Keep only first occurrence of each replication
```

---

## Managing Jobs

### Cancel All Jobs

```bash
# Get job IDs
squeue -u $USER -n dml_sim -o "%A"

# Cancel all dml_sim jobs
scancel -u $USER -n dml_sim

# Or cancel specific job
scancel JOBID
```

### Resubmit Failed Tasks

If some replications fail:

```bash
# Identify missing replications
bash slurm/check_progress.sh

# Resubmit specific array indices
sbatch --export=DGP=dgp1,N=400,METHOD=tree --array=45,67,89 \
  slurm/run_dml_simulations.slurm
```

### Hold/Release Jobs

```bash
# Hold job (prevent from running)
scontrol hold JOBID

# Release job
scontrol release JOBID
```

---

## File Locations

### On O2

**Scratch directory (temporary, fast):**
```
/n/scratch/users/${USER:0:1}/${USER}/global-scholars/results/o2_primary/
```
- Individual replication .rds files
- ~10 MB per replication → ~180 GB total
- **Cleaned periodically (30 days)** → Combine results ASAP!

**Logs directory:**
```
doubletree/simulations/production/logs/
```
- SLURM stdout: `dml_JOBID_TASKID.out`
- SLURM stderr: `dml_JOBID_TASKID.err`
- Check these for debugging

### Combined Results

**Final output:**
```
doubletree/simulations/production/results/o2_primary_combined/
```
- Permanent storage
- Combined dataset + summaries
- Safe to transfer back to local machine

---

## Expected Resource Usage

### Per Task
- **CPU:** 1 core
- **Memory:** 2-4 GB (allocated 6 GB for safety)
- **Time:** 10-30 minutes (allocated 1 hour)
- **Disk:** ~10 MB output per replication

### Total Job
- **Tasks:** 18,000
- **Parallelization:** ~500-1000 concurrent (depending on O2 availability)
- **Total time:** 30-60 minutes
- **Total disk:** ~180 GB (scratch), ~2-5 GB (combined)

---

## Comparison: O2 vs Local

| Metric | Local (Sequential) | O2 (Parallel) |
|--------|-------------------|---------------|
| Time | ~15 hours | ~30-60 minutes |
| Memory | 12-16 GB peak | 6 GB per task (distributed) |
| Risk | Single point of failure | Fault tolerant |
| Cost | Laptop unavailable | Free (O2 allocation) |

**Recommendation:** Use O2 for production runs (500+ reps), local for testing (10-50 reps)

---

## References

- O2 documentation: https://harvardmed.atlassian.net/wiki/spaces/O2
- SLURM commands: https://slurm.schedmd.com/quickstart.html
- R on O2: `module spider R`
- Simulation design: `doubletree/simulations/production/run_primary.R`

---

## Quick Reference

```bash
# Launch all simulations
bash slurm/launch_all_simulations.sh

# Check progress
bash slurm/check_progress.sh

# Monitor queue
watch -n 10 'squeue -u $USER | grep dml_sim'

# Combine results
Rscript slurm/combine_results.R

# Cancel all
scancel -u $USER -n dml_sim
```

---

## Support

- O2 helpdesk: rchelp@hms.harvard.edu
- Project issues: See `session_notes/` and `CLAUDE.md`
