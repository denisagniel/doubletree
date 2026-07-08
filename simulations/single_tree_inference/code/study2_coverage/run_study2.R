# Study 2: Full Coverage Test
# Created: 2026-04-29
#
# Purpose: Full simulation testing bias-adjusted CI coverage
#
# Settings:
# - 3 DGPs: simple, moderate, complex
# - 3 n: 500, 1000, 2000
# - 500 reps per setting
# - 3 CI methods per rep
# Total: 4500 reps (~4-5 hours)

# Load functions
source("code/dgps.R")
source("code/estimators.R")
source("code/metrics.R")
source("code/utils.R")

# TODO: Implement full study script
# 1. Define full grid
# 2. Run replications (parallel)
# 3. Compute coverage for all CI methods
# 4. Save timestamped results
# 5. Print summary
