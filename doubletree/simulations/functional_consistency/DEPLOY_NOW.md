# O2 Deployment Instructions - READY TO RUN

**Status:** Code committed (c6659e4) and pushed to GitHub
**Ready:** Yes, execute commands below

---

## Step 1: SSH to O2

```bash
ssh yourusername@o2.hms.harvard.edu
```

---

## Step 2: Pull Latest Code

```bash
cd ~/global-scholars/doubletree
git pull
```

Expected output: "Updating... Fast-forward... 17 files changed, 2803 insertions(+)"

---

## Step 3: Verify Packages Installed

```bash
module load gcc/9.2.0 R/4.2.1

Rscript -e "library(optimaltrees); library(doubletree); cat('Packages OK\n')"
```

If packages not found, install:
```bash
cd ~/global-scholars
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree
```

---

## Step 4: Quick Test (30 seconds)

```bash
cd ~/global-scholars/doubletree/simulations/functional_consistency/slurm
bash quick_test.sh
```

Expected: 3 successful replications with summary statistics showing:
- Standard M-split: max_diff > 0
- Averaged tree: max_diff = 0
- Pattern aggregation: max_diff = 0

---

## Step 5: Deploy Full Simulation

```bash
bash launch_all_simulations.sh
```

This submits **135 SLURM jobs** (one per configuration).

Expected output:
```
==========================================
Launching Functional Consistency Simulations
==========================================
Total configurations: 135
Replications per config: 500
Total replications: 67,500
==========================================

[1/135] Submitting: n=200, dgp=simple, method=standard_msplit, K=2, M=10
  Job ID: 12345678
[2/135] Submitting: n=200, dgp=simple, method=standard_msplit, K=3, M=10
  Job ID: 12345679
...
```

---

## Step 6: Monitor Progress

Run periodically (every 30-60 minutes):

```bash
bash check_progress.sh
```

Shows:
- Running/pending jobs
- Completed replications / 67,500
- Progress percentage

Or check directly:
```bash
squeue -u $USER -n fc_sim
```

---

## Step 7: After Completion (~8-12 hours)

When all jobs complete:

```bash
Rscript combine_results.R
```

Creates:
- `results/combined_fc_simulations.rds`
- `results/combined_fc_simulations.csv`

---

## Step 8: Transfer Results Back

**Option A: Git (if small enough)**
```bash
git add results/combined_fc_simulations.csv
git commit -m "FC simulation results: 67,500 replications"
git push
```

**Option B: SCP to local machine**
```bash
# On your local machine
scp username@o2.hms.harvard.edu:~/global-scholars/doubletree/simulations/functional_consistency/results/combined_fc_simulations.rds ~/Desktop/
```

---

## What to Expect

**Runtime:** ~8-12 hours
- 135 configurations
- Each: 50 array tasks × 10 reps = 500 reps
- Tasks run in parallel (depends on O2 availability)

**Storage:**
- Intermediate: ~2GB in scratch
- Final: ~200MB in results/

**Critical test:** Does coverage stay at 95% as n increases to 3200?
- If yes: Perfect FC is "free"!
- If no: Trade-off between FC and valid inference

---

## Troubleshooting

**Problem: Packages not found**
```bash
module load gcc/9.2.0 R/4.2.1
cd ~/global-scholars
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree
```

**Problem: Jobs fail immediately**
```bash
# Check error logs
cat slurm/logs/fc_sim_*.err
```

**Problem: Jobs stuck in queue**
```bash
# Check partition status
sinfo
# Consider switching partition in run_simulations.slurm
```

---

## After Analysis

Key plots to create (on local machine after transferring results):

1. **Coverage vs n by method** (critical!)
2. Functional consistency verification
3. Standardized bias trends
4. Method comparison at n=3200
5. Effect of K

See `doubletree/simulations/functional_consistency/DEPLOYMENT_CHECKLIST.md` for detailed analysis plan.

---

## Questions?

- Simulation design: See `SIMULATION_SUMMARY.md`
- Methods comparison: See `METHODS_COMPARISON.md`
- Full documentation: See `slurm/README_O2.md`
- Session log: See `quality_reports/session_logs/2026-04-13_functional-consistency-simulation.md`

---

**Ready to deploy!** Execute Step 1 above to begin.
