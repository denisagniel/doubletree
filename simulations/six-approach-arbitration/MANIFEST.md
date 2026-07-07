# Run Manifest -- six-approach-arbitration

_No run submitted yet._ This file is overwritten by `slurm/submit.sh` on each
submission with the run-id, git SHA, code hash, sizing, scratch/log paths, and
SLURM job ids for the most recent run.

Workflow:
```bash
# 1. Locally: size the jobs from a quick profile
Rscript slurm/profile_timing.R --target-hours 2

# 2. On O2: submit (mints run-id, writes this manifest)
bash slurm/submit.sh

# 3. On O2: monitor, then combine into home results/
bash slurm/monitor.sh
Rscript slurm/combine.R --run-id <run-id> --scratch-dir <scratch> --study-dir .

# 4. On O2: clean superseded scratch runs
bash slurm/clean.sh
```
