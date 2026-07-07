# O2 Simulation Infrastructure -- global-scholars

Repeatable, memory-efficient infrastructure for running R simulation studies on
the Harvard O2 SLURM cluster. Every study uses the **same layout** and the same
lifecycle, so once you learn one you know them all.

See `TRANSFER.md` for exact file-transfer commands (SSH-key auth, no passwords).

## Layout

```
simulations/
├── TRANSFER.md            # shared: O2 transfer commands
├── README_O2.md           # shared: this file
└── <study-name>/          # one self-contained subtree per study
    ├── config/
    │   ├── grid.R         # SINGLE SOURCE OF TRUTH: grid, reps, seeds, unit table
    │   ├── sizing.env     # written by profile_timing.R (bash-sourceable)
    │   └── sizing.json    # written by profile_timing.R (human-readable)
    ├── R/
    │   ├── dgp.R          # data-generating processes (you edit)
    │   ├── estimators.R   # methods under comparison (you edit)
    │   └── run_one.R      # run one (config, rep) unit -> one-row data frame
    ├── slurm/
    │   ├── run_replication.R  # runs one array task's block of units -> scratch
    │   ├── profile_timing.R   # LOCAL profiling -> job sizing (1-3 hr target)
    │   ├── array.slurm        # SLURM array job (module R, scratch I/O)
    │   ├── submit.sh          # chunk <=1000/array, throttle <=10000, run-id, waves
    │   ├── combine.R          # scratch task files -> ONE result in home results/
    │   ├── monitor.sh         # progress + failed tasks + log discovery
    │   └── clean.sh           # remove stale scratch runs (keeps home results)
    ├── logs/latest        # -> current run's scratch log dir (symlink)
    ├── results/           # FINAL combined results ONLY (home, backed up)
    └── MANIFEST.md        # current run-id, git SHA, sizing, job ids
```

## Storage discipline (scratch vs home)

- **Scratch** `/n/scratch/users/d/dma12/global-scholars/<study>/<run-id>/`
  holds all per-task result files and logs. Fast, large, **not backed up**,
  purged periodically. New code never mixes with old output because each run-id
  gets its own subtree.
- **Home** `results/<run-id>.rds` holds only the single combined final result
  per run (100 GiB quota, backed up). `combine.R` is the only thing that writes
  there.

## Lifecycle (per study)

```bash
cd simulations/<study-name>

# 1. Edit the science: config/grid.R, R/dgp.R, R/estimators.R, R/run_one.R.
#    Include at least one STRESS regime (Constitution Section 9).

# 2. Profile to size the array for a 1-3 hr/task target.
#    (A) Uniform sizing (single global reps/job) -- OK when methods have similar cost:
Rscript slurm/profile_timing.R --target-hours 2      # -> config/sizing.env + sizing.json
#    (B) Per-method sizing -- USE THIS for six-approach-arbitration, where the
#    methods differ ~50x in cost (M-split does M*K=50 fits/unit). Profiles each
#    method at its worst cell (n=2000, continuous) and sizes one array per method.
#    That worst cell is the documented Rashomon blow-up, so profile_per_method.R
#    runs each probe unit in a KILLABLE subprocess with a per-unit memory cap
#    (--mem-cap-gb, default 24): a runaway unit is killed and recorded (--mem sized
#    to a floor + loud warning) instead of OOM-killing the whole session. Set the
#    cap BELOW your salloc --mem. Run inside a generously sized interactive alloc:
salloc -p interactive -c 4 --mem=48G -t 0-04:00
Rscript slurm/profile_per_method.R --target-hours 2 --mem-cap-gb 24  # -> config/sizing_<method>.env (x7)

# 3. Commit + push; on O2 `git pull`. (See TRANSFER.md.)
#    IMPORTANT: profile on O2 (matches run hardware) OR re-profile there; per-unit
#    time is machine-dependent and drives --time. Run steps 2+ from an O2 session.

# 4. ON O2: submit. Mints a run-id, sets up scratch, writes MANIFEST.md.
bash slurm/submit.sh              # uniform sizing (A)
# OR, for per-method sizing (B):
bash slurm/submit_per_method.sh   # one array per method, per-method scratch subdirs

# 5. ON O2: monitor. Shows progress, newest logs, and failed tasks.
bash slurm/monitor.sh                # add --tail-failures to debug
#    Re-running submit.sh resumes: tasks whose scratch file exists are skipped.

# 6. ON O2: combine into the single home result (errors loudly if code changed
#    since submission, or if tasks are missing).
Rscript slurm/combine.R --run-id <run-id> \
  --scratch-dir /n/scratch/users/d/dma12/global-scholars/<study>/<run-id> \
  --study-dir .

# 7. Pull results/<run-id>.rds back to your machine (see TRANSFER.md).

# 8. ON O2: clean superseded scratch runs (home results kept).
bash slurm/clean.sh
```

## Adding another study

Re-invoke the `/setup-cluster-simulations` skill with a new study name. It
scaffolds a new `simulations/<new-study>/` subtree with the identical structure
and leaves existing studies untouched. Studies never collide: scratch, logs, and
results are all namespaced by study name.

## How this addresses common O2 pain points

| Pain point | Mechanism |
|------------|-----------|
| Different folder structures per project | One enforced layout scaffolded every time |
| Intermediate files clutter home | Scratch for per-task files/logs; home for final only |
| Hard to find recent logs | Run-id-namespaced log dir + `logs/latest` + `monitor.sh` newest-first listing |
| Stale results/logs after code changes | Per-run-id scratch subtree + source-code-hash guard in `combine.R` (loud fail) |
| Array > 1000 tasks | `submit.sh` chunks into <=1000-task arrays |
| > 10000 jobs at once | `submit.sh` throttles with `%cap` and submits dependent waves |
| Guessing wall time | `profile_timing.R` sizes reps/job from measured per-unit time (1-3 hr) |
| Out-of-memory | `profile_timing.R` sizes `--mem` from measured peak memory |
| Transfer command format | `TRANSFER.md` with exact, SSH-key-based O2 commands |

## Known-good practices baked in

- `library()` (not `devtools::load_all()`) -- module R has no devtools.
- Explicit `dest=` on every optparse option -- avoids hyphenated-arg parse bugs.
- Zero-padded task ids (`task_000123.rds`) -- consistent sort/glob.
- One-row **data frame** results (never lists) -- clean `rbind` in `combine.R`.
- Deterministic per-unit seeds (`BASE_SEED + unit`) -- reproducible & independent.
- Idempotent tasks -- resubmit to resume; existing task files are skipped.
- Atomic writes (`.tmp` then rename) -- no half-written result files.
