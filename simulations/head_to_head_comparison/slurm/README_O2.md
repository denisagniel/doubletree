# head_to_head_comparison on Harvard O2 — deployment guide

**Study:** `doubletree/simulations/head_to_head_comparison` (S5)
**Spec:** `quality_reports/specs/2026-09-08_head-to-head-comparison.md`
**Study README (design, arms, regimes, pilot results):** `../README.md`
**Created:** 2026-09-09

Read `../README.md` §0 and spec §0 before quoting any number this infrastructure
produces. **This study is not a leaderboard.** It checks three narrow claims
(C1/C3 "matched X", C2 "retained coverage where Y did not") and states in advance
where a competitor is expected to do as well or better.

---

## 1. What this deploys, and why it is batched the way it is

The sequential sweep (`code/run_sweep.R`) is **54.2 h** on one core. The pilot
(2026-09-09, `../README.md` §10) measured where that time goes:

| regime | nsim | cells | h (sequential) | share | dominant cost |
|---|---|---|---|---|---|
| R1 | 1000 | 4 | 0.94 | 2 % | flagship + GLM + lasso |
| R2 | 2000 | 6 | 22.84 | 42 % | `estimate_att_crossfit()` |
| R3 | 1000 | 8 | 11.64 | 21 % | `estimate_att_crossfit()` (anchor) |
| R4 | 1000 | 3 | 0.45 | 1 % | flagship + GLM + RF |
| R5 | 1000 | 4 | 18.14 | 33 % | **two** crossfit anchors per rep |
| R6 | 2000 | 12 | 0.23 | 0.4 % | flagship only |
| **total** | | **37** | **54.2** | | 77 % in R2 + R5 |

**The unit of work is a cost-sized block of replications of one
(regime, dgp, n) cell** — not one replication, and not one whole cell.

* Not one replication: R6 costs ~7–14 ms/rep, so per-replication tasks would
  spend more wall time starting R than estimating, and R6 alone would submit
  24,000 tasks for 14 minutes of arithmetic.
* Not one cell: R2's `n = 8000` cell is **11.6 h by itself**. Cell-level
  parallelism floors the wall clock at that one cell and spreads 37 tasks over
  four orders of magnitude of cost.

Block size is therefore `target_secs / (measured secs per rep at this n)`, with
`target_secs = 1800` (30 min) by default. `slurm/units.R` holds the map, the
measured cost anchors, and the full reasoning; `slurm/print_sizing.R` prints it.

### Job array size and budget

```
$ Rscript slurm/print_sizing.R

target seconds per array task: 1800  (override with H2H_TARGET_SECS)
walltime requested = 3 x projected + 900 s, rounded up to the minute

 regime n_tasks max_task_secs total_hours walltime mem_gb partition
     R1       5         904.8        0.94 01:01:00      4     short
     R2      48        1759.4       22.84 01:43:00      8     short
     R3      28        1675.9       11.64 01:39:00      6     short
     R4       3         928.7        0.45 01:02:00      4     short
     R5      38        1749.9       18.14 01:43:00      6     short
     R6      12         148.4        0.23 00:23:00      4     short

total array tasks: 134    total compute: 54.2 h    longest task: 29.3 min
wall clock if every task runs concurrently: ~29 min (the longest task)
```

* **134 tasks in 6 arrays** (one array per regime, sized independently). The
  regimes differ ~350× in per-replication cost, so a single `--time`/`--mem`
  would either time out R2/R5 or reserve 1.7 h and 8 GB for R6 tasks that finish
  in seconds.
* **Wall clock:** ~30 min if O2 schedules everything at once; realistically
  45–90 min on `short` with normal queueing. Against 54 h sequential.
* **Peak allocation** if all 134 run at once: 134 CPUs, ~740 GB. Cap it with
  `CONCURRENCY_CAP` (default 200 = uncapped here); `CONCURRENCY_CAP=40` costs
  ~4× the wall clock and is a good neighbour on a busy cluster.
* **Partition:** `short` throughout (O2 `short` allows 12 h; the longest request
  is 1 h 43 m).

### Two honesty notes on the sizing

1. **The cost model under-states two things, on purpose rather than silently.**
   It extrapolates *linearly in n* from the pilot's measured anchor, but
   (a) `estimate_att_crossfit()`'s nested CV is worse than linear in `n`, so R2's
   `n = 8000` and R5's `n = 2000` blocks are the likeliest to over-run; and
   (b) there is a fixed per-fit cost that does not scale with `n`, so *small*-`n`
   cheap cells run **longer** than projected — measured here, R6's `n = 200`
   block took 24.0 s against 14.8 s projected (1.6×). The `3×` safety factor on
   `--time` absorbs both, and partial flushes make an over-run recoverable rather
   than lost.
2. **Memory requests are not profiled.** The pilot ran on the dev box without RSS
   instrumentation. 4–8 GB is a conservative ~100× the data (every arm runs with
   `worker_limit = 1L`, no Rashomon enumeration, single-threaded BLAS). **Check
   the first wave and revise `MEM_GB` in `slurm/units.R`:**
   ```bash
   sacct -j <JOBID> --format=JobID,JobName%22,State,Elapsed,MaxRSS,ReqMem
   ```

---

## 2. Prerequisites

### 2.1 The one thing that is different from local development

`code/common.R`'s **default** package loading is
`pkgload::load_all()` on this source tree *and* on the sibling `../optimaltrees`
source tree. That is a dev-box pattern and it does not survive a SLURM job:
module R has no working `pkgload`/`devtools` dev-load, and `../optimaltrees` is a
relative path that is only correct when both repositories happen to be checked
out side by side.

`slurm/run_single_replication.R` therefore sets **`H2H_USE_INSTALLED=1`** before
sourcing `common.R`, which switches its section-0 block to

```r
library(optimaltrees)
library(doubletree)
```

against **installed** packages. Nothing else in `common.R` changes: every package
call in it is either `doubletree::`-qualified or reached through `asNamespace()`.
Which path ran is stamped into every cell's metadata as `meta$pkg_load_path`
(`"installed"` on the cluster, `"source"` locally), and `quick_test.sh` **fails**
if a shard reports anything but `"installed"`.

**Consequence: both packages MUST be `R CMD INSTALL`ed on O2, and re-installed
after every `git pull` that touches either package.** A forgotten re-install
means the whole sweep silently runs stale estimator code. `launch_subset.sh`
preflights that both packages are *present*, but it cannot know whether they are
*current* — re-install unconditionally.

### 2.2 Local, before pushing

```bash
cd /path/to/global-scholars
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree

cd doubletree
bash simulations/head_to_head_comparison/slurm/quick_test.sh
```

To test against a throwaway library instead of your user library:

```bash
mkdir -p /tmp/h2h_lib
R CMD INSTALL --library=/tmp/h2h_lib optimaltrees
R CMD INSTALL --library=/tmp/h2h_lib doubletree
R_LIBS=/tmp/h2h_lib bash simulations/head_to_head_comparison/slurm/quick_test.sh
```

Then commit and push (the study code and this `slurm/` directory live in the
`doubletree` repository):

```bash
cd doubletree
git add -A simulations/head_to_head_comparison
git commit -m "head_to_head_comparison: O2 deployment infrastructure"
git push
```

### 2.3 On O2

```bash
ssh dma12@o2.hms.harvard.edu
cd ~/global-scholars/doubletree
git pull

module load gcc/14.2.0
module load R/4.4.2

# Both packages, every time either repo changed. optimaltrees FIRST (doubletree
# links against it).
R CMD INSTALL ../optimaltrees
R CMD INSTALL .

Rscript -e 'library(optimaltrees); library(doubletree);
            cat("optimaltrees", as.character(packageVersion("optimaltrees")), "\n");
            cat("doubletree  ", as.character(packageVersion("doubletree")), "\n")'
```

CRAN dependencies the arms need, if the module R library does not already carry
them (`att_lasso_dml()` uses **glmnet**, `att_forest()` uses **ranger**, the
harness uses **optparse**, **digest**, **cli**, **knitr**; `data.table` is
optional):

```bash
Rscript -e 'pkgs <- c("optparse","digest","cli","knitr","glmnet","ranger","data.table")
            miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
            if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
            else cat("all present\n")'
```

Module versions above are **this project's current cluster pairing**, copied from
`../../six-approach-arbitration/slurm/{array.slurm,submit_per_method.sh}` (2026-07).

> **Stale-reference warning.** `../../functional_consistency/slurm/README_O2.md`
> still documents `gcc/9.2.0` + `R/4.2.1`. That is an April-2026 pairing and is
> **not** what this study uses. It was not modelled on and should not be copied.

---

## 3. Deploy

### 3.1 Smoke test on O2 first

```bash
cd ~/global-scholars/doubletree
bash simulations/head_to_head_comparison/slurm/quick_test.sh
```

Checks, in order: packages installed → sizing table resolves → one real block of
the cheapest regime runs end to end → the shard's contract holds (rows, arms,
`pkg_load_path == "installed"`, finite `theta`/`sigma`, no error rows). Takes
seconds. **Exit 0 is the gate; do not launch without it.**

Then inspect what would be submitted, without submitting:

```bash
DRY_RUN=1 bash simulations/head_to_head_comparison/slurm/launch_all_simulations.sh
```

### 3.2 Launch

```bash
cd ~/global-scholars/doubletree/simulations/head_to_head_comparison

# everything (134 tasks, 6 arrays)
bash slurm/launch_all_simulations.sh

# or a subset, e.g. the three cheap regimes first (20 tasks, ~1.6 h compute)
bash slurm/launch_subset.sh R1 R4 R6

# the two expensive ones on their own (86 tasks, ~41 h compute)
bash slurm/launch_subset.sh R2 R5
```

The launcher preflights (package presence, DGP-A's bilinear remainder), mints a
**run id**, creates the run's scratch tree, submits one array per regime, and
writes `../MANIFEST.md` with the run id, git SHA, sizing, job ids, and the exact
monitor and combine commands for *this* run.

**Environment gates** (all optional, all echoed at start):

| gate | default | effect |
|---|---|---|
| `RUN_ID` | fresh timestamp + git SHA | **resume** an existing run instead of starting one |
| `H2H_TARGET_SECS` | `1800` | target seconds/task; re-batches the whole array |
| `CONCURRENCY_CAP` | `200` | max simultaneously-running tasks per array |
| `SCRATCH_ROOT` | O2 scratch (below) | relocate scratch (used by local pipeline tests) |
| `DRY_RUN` | `0` | `1` prints the `sbatch` lines and submits nothing |

### 3.3 Storage layout

```
/n/scratch/users/d/dma12/global-scholars/head_to_head_comparison/<RUN_ID>/
├── sizing.txt                       # the sizing this run was submitted with
├── logs/<REGIME>_<jobid>_<task>.{out,err}
├── R1/  batch_R1_A_n500_r1000_b0001.rds ...   # one shard per array task
├── R2/  ...
└── R6/  ...
```

Scratch is **temporary** (O2 purges it). `combine_results.R` copies the assembled
cells into the study's own `results/` directory in `$HOME`, which is backed up.
`../logs/latest` is symlinked to the run's log dir for convenience.

---

## 4. Monitor

```bash
cd ~/global-scholars/doubletree/simulations/head_to_head_comparison
bash slurm/check_progress.sh                     # newest run
RUN_ID=<id> bash slurm/check_progress.sh         # a specific run
WATCH=1 bash slurm/check_progress.sh             # re-print every 60 s
```

Reports `squeue` state, shards-on-disk against the count `slurm/units.R` says
each array should produce, partials in flight, tasks killed at the wall limit
(the `received SIGTERM` sentinel), and distinct R errors across the logs.

Ad-hoc:

```bash
squeue -u $USER
scancel -u $USER --name=head_to_head_comparison-R2   # cancel one regime
tail -n 40 ../logs/latest/R2_*_7.out
```

### If tasks hit their wall limit

Nothing is lost. Each task flushes a partial every `ceiling(reps/20)`
replications, so a timeout costs at most 5 % of one block. **Resubmit with the
same run id**: finished tasks exit immediately, interrupted ones resume from
their last flush.

```bash
RUN_ID=<same id> bash slurm/launch_subset.sh R2 R5
```

If a regime times out repeatedly, the cost anchor is wrong for it: lower
`H2H_TARGET_SECS` (smaller blocks) or raise `SAFETY_FACTOR` in `slurm/units.R`.
Do **not** simply raise `--time` by hand — the sizing must stay derived from one
table, or the launcher and the workers can build different array maps.

---

## 5. Combine and retrieve

```bash
cd ~/global-scholars/doubletree/simulations/head_to_head_comparison

# what is complete, and which task ids are missing — writes nothing
Rscript slurm/combine_results.R --scratch-dir <SCRATCH_DIR> --pkg-root <PKG_ROOT> --dry-run

# assemble (the exact command with paths filled in is in ../MANIFEST.md)
Rscript slurm/combine_results.R \
  --scratch-dir /n/scratch/users/d/dma12/global-scholars/head_to_head_comparison/<RUN_ID> \
  --pkg-root ~/global-scholars/doubletree \
  --run-id <RUN_ID>
```

**Output shape is deliberately not a new format.** Each cell is reassembled into
the *same* payload `run_cell()` writes — `list(results, dgp_summary, meta)` — at
the *same* path `cell_path()` chooses, e.g.
`results/cell_R2_A2_n8000_r2000.rds`. `code/analyze.R`'s `load_results()`,
`summarise_cell()` and `rmse_ratio_table()` then work on cluster output with no
cluster-specific branch. Cluster provenance (run id, hosts, shard list, SLURM job
ids) is added under `meta$cluster`.

A cell is written **only** when its shards cover replications `1..nsim` with no
gaps; short cells are reported with the task ids to resubmit and skipped.
`--allow-incomplete` writes them anyway, named for the replications actually
present so the filename cannot overstate the study.

### The pilot-checkpoint collision — read this once

`load_results()` globs **every** `cell_*.rds` in `results/` and rbinds them. The
pilot left `cell_R1_A_n500_r30.rds` and friends there. Those do not collide by
*name* with `cell_R1_A_n500_r1000.rds`, but they do collide in
`load_results()`'s output, silently double-counting the pilot's 30 replications
alongside the sweep's 1000. `combine_results.R` **refuses to write** into a
`results/` directory that still holds pilot-scale checkpoints for cells this run
produced, and prints the `mv` command. `--archive-pilot` does the move for you
into `results/archive_pilot_<stamp>/`.

### Tables and transfer back

```bash
# regenerate the sweep's summary/ratio tables from the combined cells
Rscript simulations/head_to_head_comparison/code/run_sweep.R
# (every cell checkpoint now exists, so it computes no replications and only
#  writes tables/sweep_summary_*.csv and tables/sweep_rmse_ratios_*.csv)
```

```bash
# small: CSV tables through git
git add simulations/head_to_head_comparison/tables
git commit -m "head_to_head_comparison: full-sweep tables (run <RUN_ID>)"
git push

# large: the cell checkpoints by scp, from the LOCAL machine
scp -r dma12@o2.hms.harvard.edu:'~/global-scholars/doubletree/simulations/head_to_head_comparison/results/cell_*.rds' \
  doubletree/simulations/head_to_head_comparison/results/
```

---

## 6. File structure

```
head_to_head_comparison/
├── README.md                     # design, arms, regimes, pilot results
├── MANIFEST.md                   # written by the launcher: run id, sizing, commands
├── code/
│   ├── common.R                  # harness; section 0 gates dev vs installed loading
│   ├── dgps.R  analyze.R  verify_dgps.R  run_pilot.R  run_sweep.R
├── slurm/
│   ├── units.R                   # THE array map + pilot cost anchors + sizing
│   ├── print_sizing.R            # print the array size / walltime / memory
│   ├── run_single_replication.R  # one array task = one block of one cell
│   ├── run_simulations.slurm     # array script (modules, SIGTERM trap)
│   ├── launch_subset.sh          # the workhorse launcher (preflight, run id, sbatch)
│   ├── launch_all_simulations.sh # thin wrapper: all six regimes
│   ├── quick_test.sh             # local/cluster smoke test of the installed path
│   ├── check_progress.sh         # monitor a run
│   ├── combine_results.R         # shards -> cell_*.rds payloads analyze.R reads
│   └── README_O2.md              # this file
├── results/                      # cell_*.rds checkpoints (pilot + sweep)
├── tables/  figures/
```

---

## 7. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `there is no package called 'doubletree'` | packages not installed on O2 | §2.3 — `R CMD INSTALL ../optimaltrees` then `R CMD INSTALL .` |
| `could not find function "estimate_att"` / arg missing | stale installed package | re-install both; `common.R` asserts the expected formals and will say so |
| shard says `pkg_load_path: source` | `H2H_USE_INSTALLED` did not take | you are not running through `run_single_replication.R`; do not call `run_sweep.R` on the cluster |
| `No DESCRIPTION in ...` | wrong working directory | `--pkg-root` must be the `doubletree` package root |
| `Reused source ... not found` | partial checkout | `common.R` reuses S1/S2's DGP code from `simulations/`; `git pull` the whole repo |
| tasks TIMEOUT | crossfit cost model under-states this cell | resubmit same `RUN_ID`; then lower `H2H_TARGET_SECS` |
| tasks OUT_OF_MEMORY | `MEM_GB` guess too low | raise it in `slurm/units.R`, resubmit same `RUN_ID`; check `sacct MaxRSS` |
| combine reports cells short | tasks failed or never ran | `--dry-run` lists the task ids; resubmit same `RUN_ID` |
| combine aborts on pilot checkpoints | §5 collision guard | `--archive-pilot` |
| jobs pending forever | `short` partition busy | lower `CONCURRENCY_CAP`, or launch cheap regimes first |

**Reproducibility.** `rep_seed()` is a pure function of
`(SEED_MASTER, regime, dgp, n, rep)`, so replication *r* of a cell draws
byte-identical data whether it runs as the *r*-th replication of a sequential
`run_cell()` or as the 7th replication of batch 3 on a compute node. Batching is
invisible in the results, and the entire 54 h sweep remains a function of the
single integer `SEED_MASTER = 20260908`.

---

## 8. Version history

- **2026-09-09:** Created. Modelled on `six-approach-arbitration/slurm/` (this
  project's most current cluster pattern: `library()`-based loading, run-id
  scratch lifecycle, per-unit flushes, SIGTERM trap, code-freshness preflight),
  with the per-regime unit table replacing that study's global unit index plus
  `--unit-offset`/`--max-unit` reconciliation.
