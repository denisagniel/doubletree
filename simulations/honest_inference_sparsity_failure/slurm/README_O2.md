# O2 deployment — `honest_inference_sparsity_failure` (S2)

**Study:** `doubletree/simulations/honest_inference_sparsity_failure/`
**Spec:** `quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md`
(the 2026-09-01 changelog at the top of that file is authoritative)
**Scaffolded:** 2026-09-09
**Cluster:** Harvard O2 (SLURM), `gcc/14.2.0` + `R/4.4.2`

---

## 1. What this deployment is for

The spec's grid is **5 DGP variants × n ∈ {500, 2000, 8000}**:

| variant | `eps` | role |
|---|---|---|
| `eps0` | 0.00 | favourable / exact sparsity (sanity check) |
| `eps_small` | 0.03 | small violation |
| `eps_mod` | 0.09 | moderate violation |
| `eps_sev` | 0.15 | severe violation |
| `orth_blind` | 0.15 | `prop:spectest`'s **blind spot** — orthogonal misspecification, true bias exactly 0 |

Replications per cell (`code/common.R`, `reps_for_n()`):

- **R = 300** at n = 500 and n = 2000
- **R = 150** at n = 8000 — the spec's **own sanctioned exception** ("consider
  dropping the `n = 8000` cell to a smaller `R` if timing is prohibitive, stating
  the reduced `R` explicitly rather than silently"). It is stated at
  `REPS_DEFAULT`/`REPS_AT_MAX_N` in `code/common.R`, in `code/run_sweep.R`'s
  header, in every checkpoint's `meta$reps`, and in `code/analyze.R`'s output.

**Status as of 2026-09-09.** The five **n = 500** cells are complete on disk
(`results/cell_*_n500_r300.rds`, written sequentially on 2026-09-01). The study was
then blocked from 2026-09-04 on `optimaltrees >= 0.4.1` (0.4.0 was installed); that
blocker is resolved.

**The ten remaining cells — 5 variants × n ∈ {2000, 8000} = 2,250 replications
(1,500 + 750) — are what this directory exists to run.** n = 500 stays a valid
argument to `launch_subset.sh` so a backfill goes through the same machinery, but it
is never a default.

Every replication is **two** estimator calls on the same draw:
`estimate_att()` (the method under test) and `estimate_att_crossfit()` (the anchor).

---

## 2. Prerequisites, in order

### 2.1 Both packages must be **R CMD INSTALLed** on the cluster

Module R carries no working `devtools`/`pkgload` dev-load, and
`pkgload::load_all("../optimaltrees")` is a **relative sibling path** that is only
correct when both repositories happen to be checked out side by side. Neither
assumption is safe inside a SLURM job.

```bash
module load gcc/14.2.0
module load R/4.4.2

# From the parent directory holding both repos:
R CMD INSTALL optimaltrees     # must end up >= 0.4.1
R CMD INSTALL doubletree
```

`optimaltrees` has compiled code (`src/`, Rcpp), so this needs the `gcc` module
loaded first. `doubletree` is pure R.

Verify:

```bash
Rscript --no-init-file -e 'for (p in c("optimaltrees","doubletree"))
  cat(p, as.character(packageVersion(p)), dirname(find.package(p)), "\n")'
```

`launch_subset.sh` re-checks this as a preflight and **refuses to submit** if either
package is missing or if `optimaltrees < 0.4.1`. A forgotten `R CMD INSTALL` after a
`git pull` means every task silently runs stale estimator code and the whole wave has
to be thrown away; the preflight is the cheapest place to catch it.

### 2.2 `code/common.R`'s loading gate

`code/common.R` section 0 has **two** loading paths behind one environment variable
(added 2026-09-09; the same gate `head_to_head_comparison/code/common.R` carries):

| `HIS_USE_INSTALLED` | path |
|---|---|
| unset or `"0"` (default) | **DEV**: `pkgload::load_all()` on this source tree and on `../optimaltrees`. Unchanged from before this deployment existed. |
| `"1"` | **CLUSTER**: `library(optimaltrees); library(doubletree)` against installed packages. |

Every script in this directory sets `HIS_USE_INSTALLED=1` before sourcing
`common.R`. Nothing else in `common.R` differs between the paths: every package call
in it is either `doubletree::`-qualified or reached through `asNamespace()`. The path
that actually ran is stamped into each shard as `runtime$pkg_load_path` and into each
combined cell as `meta$cluster$pkg_load_paths`.

### 2.3 ⚠ The tracked `.Rprofile` hazard — read this one

**`doubletree/.Rprofile` is checked into git** and runs

```r
devtools::load_all("../optimaltrees")
```

at startup. R reads `./.Rprofile` as the user profile whenever `R_PROFILE_USER` is
unset. So a plain `Rscript` launched from (or near) the package root — which is what
SLURM does — **dev-loads `optimaltrees` from source before any of our code runs**.
`common.R`'s `isNamespaceLoaded("optimaltrees")` guard then finds it already present
and `library(optimaltrees)` is a no-op on the dev namespace.

The consequence is not an error. It is a completed sweep whose numbers came from
**uninstalled source**, reported as `pkg_load_path = "installed"`. Two defences:

1. **`Rscript --no-init-file`** everywhere in this directory
   (`run_simulations.slurm`, `launch_subset.sh`, `quick_test.sh`,
   `check_progress.sh`, `combine_results.R`). `.Renviron` is *still* read, so a
   per-machine `R_LIBS_USER` keeps working.
2. **`assert_installed_load()`** in `run_single_replication.R` compares each loaded
   namespace's directory against `.libPaths()` and **aborts** if either was
   dev-loaded. Verified to trip on 2026-09-09 by deliberately omitting
   `--no-init-file`.

`check_progress.sh` reports a "dev-load guard trips" count; anything other than 0
means tasks were launched without `--no-init-file`.

If you prefer a belt-and-braces shell-level fix: `export R_PROFILE_USER=/dev/null`.

### 2.4 Data dependencies

None. The DGP family is generated analytically by `code/dgps_delta_dial.R` with
exact enumeration over the `2^p` covariate grid, reusing S1's
`simulations/partition_recovery_clt/code/enumerate_sufficient_class.R`. No binary
data files need to reach the cluster — only the two git repositories.

---

## 3. Files

| file | role |
|---|---|
| `units.R` | **The array-index → unit-of-work map**, and the only copy of it. Cost anchors, batch boundaries, walltime/memory sizing, shard filenames. |
| `run_single_replication.R` | One array task: one contiguous **block** of replications of one `(dgp, n)` cell. Resume-from-partial, atomic flush, installed-load assertion. |
| `run_simulations.slurm` | The sbatch script. Module loads, SIGTERM trap + sentinel, `--no-init-file`. Not run directly. |
| `launch_subset.sh` | **The workhorse launcher.** Preflight, run-id, one array per wave, manifest. |
| `launch_all_simulations.sh` | One-line wrapper: `launch_subset.sh n2000 n8000`. |
| `print_sizing.R` | Prints the sizing in `table` / `env` / `units` formats. The launcher and the workers read the *same* table through this. |
| `quick_test.sh` | ~1-minute wiring smoke test at **both** missing sample sizes. No SLURM. |
| `check_progress.sh` | Monitor: squeue, shards-vs-expected, partials, timeouts, errors, spec-grid completeness. |
| `combine_results.R` | Shards → the `cell_*.rds` checkpoints `code/analyze.R` already reads. |

### What one unit of work is

**One unit = one `(dgp, n, contiguous block of replications)`.** Not one
replication (a SLURM task pays seconds of module-load and R-startup overhead against
a 3.8 s replication, and one-rep tasks would submit 2,250 tasks for 4 hours of
arithmetic), and not one whole cell (that floors the wall clock at the slowest cell
and wastes the array on 10 tasks).

**Batching is exactly reproducible against the sequential sweep.** `rep_seed()` is a
pure function of `(SEED_MASTER, dgp, n, rep)`, so replication 137 of a cell draws
byte-identical data whether it is the 137th replication of a sequential `run_cell()`
or the 62nd of batch 2 here. The cells this produces are directly comparable to the
five n = 500 cells produced sequentially on 2026-09-01. **Do not change the seed
scheme to anything sequential** — a `set.seed(1); for (r in ...)` stream would make
batch boundaries visible in the numbers.

### One array per **wave** (= per sample size)

`n` is the only axis the cost and the replication count vary along: all five variants
are the same size at the same `n` (same `p`, same `leaf_budget`, same estimator pair
— they differ only in `eps` and in which coordinate the residual lives on, neither of
which changes fit cost), while `n` changes per-replication cost 3× and `R` from 300 to
150. A single `--time`/`--mem` would either time out n = 8000 or reserve an hour for
n = 2000 tasks that finish in ten minutes.

Per-wave arrays also make the task-id → block map **offset-free**: each array's task
ids restart at 1 and index directly into that wave's own table.
(`six-approach-arbitration` carries a *global* unit index plus
`--unit-offset`/`--max-unit` to reconcile it, and its own comments record the
2026-07-10 contamination bug the reconciliation existed to fix. Per-wave tables
remove the class of bug rather than guarding against it.)

---

## 4. Job array size and budget

### Cost model: measured per `n`, **not** extrapolated

`head_to_head_comparison/slurm/units.R` anchors one measured cost per regime and
scales it linearly in `n`. That is deliberately **not** done here: this study's cost
is measurably **sublinear** over the relevant range, so a linear extrapolation from
the n = 2000 anchor would predict 15.3 s/rep at n = 8000 and over-reserve by 36 %.
Both needed sizes were measured directly, so the model is a lookup.

| n | measured s/rep | provenance |
|---|---|---|
| 500 | 1.90 | study's own 2026-09-01 pilot (recorded at `REPS_DEFAULT` in `code/common.R`) |
| 2000 | 3.83 | measured 2026-09-09, 3 reps of `eps_sev` via the installed-library path (3.75, 4.16, 3.57) |
| 8000 | 11.25 | measured 2026-09-09, same protocol (10.79, 11.36, 11.59) |

The 2026-09-09 figures reproduce the 2026-09-01 pilot's 3.8 / 11.1 to within 1–2 %,
on a *different* package-loading path — which is why three replications are enough to
size an array from.

**98 % of the cost is the ANCHOR, not the method under test.** At n = 8000,
`secs_full ≈ 0.07` against `secs_cf ≈ 11.2`: `estimate_att_crossfit()` runs nested CV
(K = 5 outer folds × an adaptive lambda search × `cv_K = 5` inner folds × 2
nuisances) against a single pair of tree fits for `estimate_att()`. Any future
re-timing should be read as a measurement of `estimate_att_crossfit()`, and a change
in its CV defaults is the one package change that invalidates the table.

### Sizing (`HIS_TARGET_SECS = 900`, the default)

```
 wave    n reps_per_cell n_tasks max_task_secs total_hours walltime mem_gb partition
n2000 2000           300      10         574.5        1.60 00:44:00      6     short
n8000 8000           150      10         843.8        2.34 00:58:00      8     short

cells: 10   total array tasks: 20   total compute: 3.94 h   longest task: 14.1 min
wall clock if every task runs concurrently: ~14 min
```

Reproduce with `Rscript --no-init-file slurm/print_sizing.R` (add
`--format units` for the full per-task map).

**Why 900 s and not `head_to_head_comparison`'s 1800.** The whole remaining sweep is
only ~3.9 h, so an 1800 s target maps each of the ten cells onto exactly one task and
floors the wall clock at the slowest single cell (28 min). At 900 s every cell splits
in two, the array doubles to 20 tasks, and the longest task drops to ~14 min. The
cost is one extra R startup per cell (~10 s, under 2 % of a task).

**Walltime** = `3 × projected + 900 s`, rounded up to the minute. The safety factor
is 3 even though the cost model is a direct measurement rather than an extrapolation,
because the measurement was taken on the **dev box** and an O2 compute node's
per-core throughput is not the dev box's. Over-reserving inside `short` (12 h limit)
is free; a timeout costs a resubmit.

**Memory is NOT profiled on O2** — stated so nobody reads `MEM_GB_BY_WAVE` as
measured. It is conservative: both estimator calls run `worker_limit = 1L`, no
Rashomon enumeration anywhere, `parallel = FALSE` on the cross-fit call, BLAS/OpenMP
pinned to one thread. `p = 5` binary covariates at n ≤ 8000 is a trivially small
design matrix; n = 8000 gets more only because K = 5 nested CV holds several fold fits
at once. **Revise these from `sacct -o MaxRSS` after the first wave** rather than
leaving them as guesses (§7).

---

## 5. Deploy

### 5.1 Locally, before pushing

```bash
cd doubletree
R CMD INSTALL ../optimaltrees      # only if not already installed
R CMD INSTALL .
bash simulations/honest_inference_sparsity_failure/slurm/quick_test.sh
```

Exit 0 means the `library()`-based path works end to end at **both** n = 2000 and
n = 8000. Verified 2026-09-09: 4.10 s/rep and 11.08 s/rep, all three intervals
finite, all coverage indicators non-NA, shard schema **identical** to the production
n = 500 cells (67 columns).

Then dry-run the launcher (works on a dev box; the module loads are guarded on
`command -v module`):

```bash
DRY_RUN=1 SCRATCH_ROOT=/tmp/his_test \
  bash simulations/honest_inference_sparsity_failure/slurm/launch_all_simulations.sh
```

### 5.2 On O2

```bash
ssh dma12@o2.hms.harvard.edu
cd ~/global-scholars/doubletree
git pull

module load gcc/14.2.0 R/4.4.2
R CMD INSTALL ../optimaltrees && R CMD INSTALL .    # §2.1

cd simulations/honest_inference_sparsity_failure

bash slurm/quick_test.sh                            # ~1 min, must exit 0
DRY_RUN=1 bash slurm/launch_all_simulations.sh       # inspect the plan
bash slurm/launch_all_simulations.sh                 # submit 2 arrays, 20 tasks
```

**Staged rollout is recommended for the first submit:** run the cheaper wave alone,
confirm shards land and `sacct` memory looks sane, then submit the expensive one.

```bash
bash slurm/launch_subset.sh n2000
# ... check, then:
bash slurm/launch_subset.sh n8000
```

`launch_subset.sh` writes `MANIFEST.md` in the study directory recording the run id,
git SHA, gates, sizing, job ids, and the exact monitor/combine commands with paths
already substituted.

### 5.3 Monitor

```bash
bash slurm/check_progress.sh          # newest run
WATCH=1 bash slurm/check_progress.sh  # re-print every 60 s
RUN_ID=<id> bash slurm/check_progress.sh
```

It reports shards-on-disk against the number `units.R` says the array should produce
(so "13/20" cannot drift from what was submitted), partials in flight, **wall-limit
kills via the `received SIGTERM` sentinel**, R-level errors, dev-load guard trips,
and the spec-grid view of `results/`.

### 5.4 Combine

```bash
Rscript --no-init-file slurm/combine_results.R \
  --scratch-dir /n/scratch/users/d/dma12/global-scholars/honest_inference_sparsity_failure/<RUN_ID> \
  --pkg-root $PWD/../.. --run-id <RUN_ID>
```

Add `--dry-run` first: it lists exactly which cells are short and which task ids to
resubmit. This writes `results/cell_<dgp>_n<n>_r<reps>.rds` — **the same payload
`run_cell()` writes, at the same path `cell_path()` chooses** — so `code/analyze.R`
works on cluster output with no cluster-specific branch.

### 5.5 Analyse

```bash
Rscript --no-init-file simulations/honest_inference_sparsity_failure/code/analyze.R
```

`analyze.R` **executes the spec's Decision rule**: if `ci_95_anchor` coverage falls
below `1-alpha` at any `eps` it prints a STOP banner, because that contradicts
`thm:anchor`'s unconditional proof and must be escalated to Oracle before anything
else in the study is interpreted. It also carries the Arm 0 flag (does
`estimate_att_crossfit()` itself satisfy `ass:rate` at this cell?) into every later
table.

---

## 6. Crash safety and resume

- **Idempotent tasks.** If a task's shard already exists it exits 0 immediately, so
  resubmitting a run-id costs nothing for finished tasks.
- **Atomic partial flushes.** Every `FLUSH_EVERY` (25) replications the block's
  completed rows are written to `part_<shard>.rds` via write-to-tmp-then-rename, so a
  kill mid-write cannot leave a truncated `.rds` a later resume would read as valid.
  A wall-time timeout loses at most 25 replications.
- **Resume is unit-scoped.** A partial is trusted only if it belongs to the exact
  same `(wave, dgp, n, rep_from, rep_to)`. A stale partial from a differently-sized
  array (different `--target-secs`, hence different batch boundaries) is reported and
  ignored, never silently absorbed.
- **To resume, resubmit with the same run id:**
  ```bash
  RUN_ID=<id> bash slurm/launch_subset.sh n8000
  ```
- **No silent NA.** `on_error = "record"`: a failed replication persists **one full
  row** carrying the condition message in `error_message`, per
  `common.R::rep_row_template()`. Never a plausible-looking NA estimate with no
  explanation, and never a short row that drops out of an `rbind`.

### What `combine_results.R` refuses to do

- **Write a short cell under a full name.** A file called `..._r300.rds` holding 210
  replications is *worse* than a missing one: every MC standard error is right for 210
  and reported as if 300, and `analyze.R` looks the file up by that exact name and
  would accept it silently. Short cells are skipped and reported with the task ids to
  resubmit. `--allow-incomplete` writes them named for the replications actually
  present, so the filename cannot lie (and `analyze.R` will then correctly not find
  them).
- **Combine a drifted map.** If any replication appears in more than one shard it
  aborts — two shards claiming the same block means the array map changed between
  launch and combine.
- **Combine under a different `R`.** If the shards' `reps_total` disagrees with
  `reps_for_n(n)` in the combining session (i.e. `HIS_REPS`/`HIS_REPS_MAX_N` changed)
  it aborts, because the cell would be written under a filename describing a
  replication count it does not have. Use the values recorded in the run's
  `MANIFEST.md`.
- **Combine dev-loaded results.** If any shard reports `pkg_load_path != "installed"`
  it aborts.
- **Overwrite an existing cell** unless `--overwrite` is passed.

---

## 7. After the first wave: revise the memory request

```bash
sacct -j <ARRAY_JOB_ID> --format=JobID,JobName%30,State,Elapsed,MaxRSS,ReqMem
```

Then edit `MEM_GB_BY_WAVE` in `slurm/units.R` from what was actually used. Leaving
guesses in place is how a later, larger wave gets OOM-killed.

Also compare `runtime$elapsed_secs` against `unit$est_secs` in a few shards. If O2
nodes are systematically slower than the dev box, **re-measure one cost anchor** in
`COST_ANCHOR_BY_N` — that is the fix, not raising `SAFETY_FACTOR` forever.

---

## 8. Environment gates (all defaulted, all echoed)

| gate | default | effect |
|---|---|---|
| `HIS_USE_INSTALLED` | `0` | `1` = `library()` cluster path in `code/common.R` (every script here sets it) |
| `HIS_TARGET_SECS` | `900` | target wall cost per array task; **must** match between launcher, worker and combine |
| `HIS_REPS` | `300` | replications per cell for n < 8000 |
| `HIS_REPS_MAX_N` | `150` | replications per cell at n = 8000 |
| `RUN_ID` | fresh timestamp + git SHA | set to **resume** an existing run |
| `SCRATCH_ROOT` | `/n/scratch/users/d/dma12/global-scholars/honest_inference_sparsity_failure` | scratch root |
| `CONCURRENCY_CAP` | `200` | max simultaneously-running tasks per array |
| `DRY_RUN` | `0` | `1` = print the sbatch lines, submit nothing |
| `WATCH` | `0` | `1` = `check_progress.sh` re-prints every 60 s |
| `WAVES`, `TASK_ID`, `REPS` | `n2000 n8000`, `1`, `2` | `quick_test.sh` scope |

`HIS_TARGET_SECS`, `HIS_REPS` and `HIS_REPS_MAX_N` are **exported** by the launcher
into the array's environment, not merely documented, so the map the worker builds
cannot differ from the map the array size was computed from.

---

## 9. Storage

| tier | path | contents |
|---|---|---|
| scratch (fast, temporary) | `/n/scratch/users/d/dma12/global-scholars/honest_inference_sparsity_failure/<RUN_ID>/<wave>/` | `batch_*.rds` shards, `part_*.rds` partials, `logs/` |
| repo (permanent) | `doubletree/simulations/honest_inference_sparsity_failure/results/` | `cell_*.rds` — what `analyze.R` reads |

Shards stay in scratch; only the combined cells land in the repository. Total scratch
footprint for this deployment is small (each shard is a ≤150-row data.frame plus
population summary — tens of KB), so no cleanup is urgent, but O2 scratch is purged on
a retention policy: **combine before the retention window closes.**

Bringing results back:

```bash
# from the dev box
scp -i ~/.ssh/id_rsa_o2 \
  'dma12@transfer.rc.hms.harvard.edu:~/global-scholars/doubletree/simulations/honest_inference_sparsity_failure/results/cell_*_n2000_r300.rds' \
  doubletree/simulations/honest_inference_sparsity_failure/results/
scp -i ~/.ssh/id_rsa_o2 \
  'dma12@transfer.rc.hms.harvard.edu:~/global-scholars/doubletree/simulations/honest_inference_sparsity_failure/results/cell_*_n8000_r150_.rds' \
  doubletree/simulations/honest_inference_sparsity_failure/results/
```

(Use `transfer.rc.hms.harvard.edu`, not the login node, for bulk copies.)

---

## 10. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `package(s) optimaltrees were DEV-LOADED from source` | tracked `.Rprofile` ran `devtools::load_all("../optimaltrees")` | launch with `Rscript --no-init-file` (`run_simulations.slurm` does) or `export R_PROFILE_USER=/dev/null`. §2.3 |
| `package doubletree is NOT installed on this node` | forgot `R CMD INSTALL` after `git pull` | §2.1 |
| `optimaltrees 0.4.0 is installed but this study needs >= 0.4.1` | stale sibling install | reinstall `optimaltrees` |
| `Loaded doubletree::estimate_att has no <arg> argument` | a stale installed copy shadows the source tree | reinstall both packages |
| `cli` error `Invalid cli literal: {.arg} starts with a dot` | pre-existing latent bug in `common.R`'s fail-fast formals check (fires only when the check *fails*, i.e. when a package is stale) | it is telling you the packages are stale — reinstall. Same latent bug exists in `head_to_head_comparison/code/common.R`. |
| tasks time out | dev-box cost anchor too optimistic for an O2 core | partials are safe; resubmit the same `RUN_ID`. Then re-measure the anchor (§7) |
| `replication(s) ... appear in more than one shard` | `HIS_TARGET_SECS` differed between launch and combine | re-launch with a fresh run id |
| `shards ... were produced with reps_total = X, but reps_for_n(n) = Y` | `HIS_REPS`/`HIS_REPS_MAX_N` differ from launch | use the values in the run's `MANIFEST.md` |
| `analyze.R` aborts naming missing checkpoints | some cells not yet combined | `check_progress.sh`, then `combine_results.R` |
| `ci_95_anchor failed to cover while the anchor's own Wald interval did` | deterministic impossibility → `anchor_ci_harness()` is misconstructed | this is a real bug, not a statistical outcome. See `common.R`'s note on the spec's spurious `/sqrt(n)` |

`squeue -u $USER`, `scancel <jobid>`, `scontrol show job <jobid>` as usual.

---

## 11. Provenance recorded per run

Each combined cell's `meta` carries `run_cell()`'s own stamp (seed scheme,
`leaf_budget`, `lambda_n`, anchor settings, exact `(delta_e, delta_mu)`,
`bias_pseudo`, `D_w`/`D_mu`, `<g,h>_nu`, clip bounds, git SHAs, package versions)
**plus** `meta$cluster`:

`run_id`, `wave`, `scratch_dir`, `target_secs`, `n_shards`, `shard_files`, `hosts`,
`slurm_job_ids`, `pkg_load_paths`, `elapsed_secs_sum`, `complete`, `reps_target`.

A results file that is argued over can name the run, the hosts and the shards it came
from without anyone consulting a shell history.

---

## 12. What this deployment must **not** be used for

- **Do not** re-run the five n = 500 cells as part of a normal launch. They are
  complete and comparable; `launch_all_simulations.sh` excludes them deliberately.
- **Do not** change `rep_seed()` (§3). Batch reproducibility against the existing
  n = 500 cells depends on its purity.
- **Do not** raise `leaf_budget` at larger `eps`. It is fixed at 2 across every
  variant on purpose (spec, "Estimand and method"): letting it grow would let the
  tree "solve around" the manufactured sparsity violation and defeat the sweep.
  `dgps_delta_dial.R` design note (1) records why 2 specifically.
- **Do not** enable Rashomon enumeration or any parallelism inside a task. That was
  the documented cause of every previous memory-exhaustion incident in this codebase
  (`simulations/docs/MEMORY_SAFE_SIMULATIONS.md`, `MEMORY.md`
  `[LEARN:rashomon-memory]`).
- **Do not** quote a number from an incomplete cell. `combine_results.R` refuses to
  write one under a full name for exactly this reason.
