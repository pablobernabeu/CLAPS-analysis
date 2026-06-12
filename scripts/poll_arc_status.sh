#!/usr/bin/env bash
# scripts/poll_arc_status.sh  
# Quick HPC status poll that writes output for automated parsing.
# Usage from design_analysis/: bash scripts/poll_arc_status.sh

set -euo pipefail

SSH="/c/Windows/System32/OpenSSH/ssh.exe"
OUTPUT_FILE="outputs/arc_queue_status.txt"
mkdir -p "$(dirname "$OUTPUT_FILE")"

{
  echo "=== ARC Queue Status ==="
  echo "Timestamp: $(date -Iseconds 2>/dev/null || date +%Y-%m-%dT%H:%M:%S%z)"
  echo ""
  
  echo "--- Active Queue ---"
  $SSH arc 'squeue -u "$(whoami)" --format="%.18i %.10j %.8T %.10M %.10l %.6D %R"' || echo "ERROR: squeue failed"
  
  echo ""
  echo "--- Recent Job Accounting (last 24h) ---"
  $SSH arc 'sacct -u "$(whoami)" --format=JobID%-20,JobName%-30,State,ExitCode,Elapsed,MaxRSS --starttime=now-24hours | head -100' || echo "ERROR: sacct failed"
  
  echo ""
  echo "--- Completed Design Analysis Outputs ---"
  find outputs/design_analysis -name "*.qs" -type f -printf "%f\n" 2>/dev/null | wc -l || echo "0"
  echo " outputs found locally"
  
} > "$OUTPUT_FILE" 2>&1

echo "[poll] Status written to: $OUTPUT_FILE"
cat "$OUTPUT_FILE"
