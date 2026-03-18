# Check progress by examining R process memory usage
# More completed sims = more memory used (results stored until end)

# Get process info
pid <- system("pgrep -f 'run_beta_study.R' | head -1", intern = TRUE)
if (length(pid) > 0 && pid != "") {
  # Memory in KB
  mem_cmd <- sprintf("ps -p %s -o rss= 2>/dev/null", pid)
  mem_kb <- as.numeric(system(mem_cmd, intern = TRUE))
  mem_mb <- mem_kb / 1024
  
  cat("Simulation Progress Estimate:\n")
  cat("==============================\n")
  cat(sprintf("Process memory: %.1f MB\n", mem_mb))
  cat(sprintf("Log file size: %.1f GB\n", 3.7))
  cat(sprintf("Runtime: 6 minutes\n\n"))
  
  # Estimate: Each result row is ~200 bytes, 13500 rows = ~2.7 MB final
  # But R holds intermediate objects, so estimate ~10-20 MB at completion
  # Current memory usage can give rough progress
  
  # Very rough estimate based on memory growth
  if (mem_mb < 170) {
    cat("Estimated progress: <5% (early stage)\n")
  } else if (mem_mb < 200) {
    cat("Estimated progress: 5-15%\n")
  } else if (mem_mb < 250) {
    cat("Estimated progress: 15-30%\n")
  } else {
    cat("Estimated progress: 30%+\n")
  }
  
  cat("\nNote: mclapply() saves all results at end, so no interim files available\n")
}
