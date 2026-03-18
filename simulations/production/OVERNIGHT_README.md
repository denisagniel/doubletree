# Running Simulations Overnight

## Quick Start

```bash
cd /Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/production

# Start simulations in background
nohup ./run_overnight.sh > logs/nohup.out 2>&1 &

# Check status
./check_status.sh

# Monitor in real-time
tail -f logs/primary_*.log
```

## What Gets Run

**Primary Simulations (`run_primary.R`):**
- 3 DGPs × 4 methods × 3 sample sizes × 500 reps = **18,000 runs**
- Methods: tree-DML, rashomon-DML, forest-DML, linear-DML
- Sample sizes: n ∈ {400, 800, 1600}
- **Estimated runtime:** ~4 hours with 4 cores

## Output

Results saved to:
```
results/primary_YYYY-MM-DD/
├── simulation_results.rds    # Full replication data
└── summary_stats.csv          # Aggregated metrics (coverage, bias, RMSE)
```

Logs saved to:
```
logs/
├── primary_YYYYMMDD_HHMMSS.log   # Full output
└── nohup.out                     # Background process output
```

## Monitoring

### Check if running
```bash
./check_status.sh
```

### Watch progress in real-time
```bash
tail -f logs/primary_*.log | grep -E "(Simulations:|Memory:|Completed:)"
```

### Check results so far
```bash
Rscript -e 'readRDS("results/primary_2026-03-12/simulation_results.rds") |> nrow()'
```

## Stopping

### Graceful stop
```bash
pkill -f "Rscript run_primary.R"
```

### Force stop (if graceful doesn't work)
```bash
pkill -9 -f "Rscript run_primary.R"
rm logs/primary.pid
```

## Troubleshooting

### Out of memory
If simulations crash with memory errors:
1. Reduce N_CORES in `run_primary.R` (line 72): change from 4 to 2
2. Reduce N_REPS for testing (line 58): change from 500 to 50

### Check for errors
```bash
grep -i "error\|failed\|warning" logs/primary_*.log
```

### Resume from checkpoint
If simulations crash partway through, they can be resumed (results are appended incrementally).

## After Completion

1. Check results:
   ```bash
   Rscript analyze_manuscript.R
   ```

2. Generate tables for manuscript:
   ```bash
   # Results will be in results/primary_*/summary_stats.csv
   ```

3. Commit results:
   ```bash
   git add results/primary_*/
   git commit -m "Add primary simulation results (YYYY-MM-DD)"
   ```

## Status Indicators

- ✓ **Running**: PID file exists, process active
- ✗ **Not running**: No PID file or stale PID
- **Progress**: Look for "Simulations: X/18000" lines in log
- **Memory**: Look for "Memory:" lines in log
- **Errors**: Look for "FAILED" or "Error:" in log

## Expected Timeline

- **Start:** When you run `nohup ./run_overnight.sh &`
- **First checkpoint:** ~15 minutes (first 100 sims complete)
- **25% complete:** ~1 hour
- **50% complete:** ~2 hours
- **75% complete:** ~3 hours
- **Completion:** ~4 hours

Check status periodically with `./check_status.sh` or let it run fully overnight.
