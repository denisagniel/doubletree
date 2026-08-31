# Archive Note: six_approach_comparison

**Archived:** 2026-08-03
**Canonical successor:** `simulations/six-approach-arbitration/`

## What this was

The `six_approach_comparison` study compared six approaches to doubletree-based
ATT estimation across 4 DGPs and 3 sample sizes (n = 500, 1000, 2000), 1000 reps
per config. It was the primary simulation study driving the approach-selection
decision for the doubletree manuscript.

## When it ran

Initial runs: May–Jun 2026. Final deployment (approach 6 fix for DGP4 failures
caused by `discretize_bins=2` stump structures): 2026-06-10. Last result file
activity: Jul 2026.

## Results

Final combined results in `study/results/combined/all_results.rds`.
Summary CSVs: `summary_inference.csv`, `summary_timing.csv`.
Deployment notes: `study/slurm/DEPLOYMENT_2026-06-10.md`.

## Why superseded

`six-approach-arbitration` replaces this study structurally:
- Uses the current skill-template scaffold (per-unit checkpoint, stale-code hash,
  per-stratum sizing, `combine.R` coverage assertion)
- Cleaner stratum-based approach decomposition
- Fixes the per-method contamination bug present in the old per-approach `.sh` pattern

## Do not delete

The `all_results.rds` file and deployment notes provide context for manuscript
decisions. Keep until the paper is submitted.
