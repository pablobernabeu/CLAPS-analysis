#!/usr/bin/env bash
# scripts/analyze_logs.sh
# Analyze local SLURM logs for errors, failures, and diagnostics.
# Usage from design_analysis/: bash scripts/analyze_logs.sh

set -euo pipefail

LOG_DIR="outputs/logs"
SUMMARY_FILE="outputs/log_analysis.txt"

if [[ ! -d "$LOG_DIR" ]] || [[ -z "$(ls -A "$LOG_DIR" 2>/dev/null)" ]]; then
  echo "ERROR: No logs found in $LOG_DIR"
  echo "Run: bash scripts/fetch_arc_logs.sh"
  exit 1
fi

echo "========================================"
echo "  Analyzing SLURM Logs"
echo "========================================"
echo ""

{
  echo "=== SLURM Log Analysis ==="
  echo "Timestamp: $(date -Iseconds 2>/dev/null || date +%Y-%m-%dT%H:%M:%S%z)"
  echo ""
  
  # Count logs
  OUT_COUNT=$(find "$LOG_DIR" -name '*.out' -type f 2>/dev/null | wc -l | tr -d ' ')
  ERR_COUNT=$(find "$LOG_DIR" -name '*.err' -type f 2>/dev/null | wc -l | tr -d ' ')
  echo "Total .out files: $OUT_COUNT"
  echo "Total .err files: $ERR_COUNT"
  
  echo ""
  echo "--- OOM (Out of Memory) Failures ---"
  grep -l -E '(Out of memory|oom_kill|exit code 137|Killed|OOM)' "$LOG_DIR"/*.{out,err} 2>/dev/null | head -20 || echo "None found"
  
  echo ""
  echo "--- Timeout Failures ---"
  grep -l -E '(TIMEOUT|TIME LIMIT|DUE TO TIME LIMIT)' "$LOG_DIR"/*.{out,err} 2>/dev/null | head -20 || echo "None found"
  
  echo ""
  echo "--- R Errors ---"
  grep -l -E '(Error in|Error:|Execution halted|cannot open connection)' "$LOG_DIR"/*.out 2>/dev/null | head -20 || echo "None found"
  
  echo ""
  echo "--- Stan Compilation Errors ---"
  grep -l -E '(COMPILATION ERROR|Could not compile model|Stan compilation failed)' "$LOG_DIR"/*.out 2>/dev/null | head -10 || echo "None found"
  
  echo ""
  echo "--- Divergences (Posterior Pathologies) ---"
  grep -E 'divergent transition|max_treedepth|Rhat' "$LOG_DIR"/*.out 2>/dev/null | head -20 || echo "None found"
  
  echo ""
  echo "--- Successful Completions ---"
  grep -l 'COMPLETED\|Success\|Model fitted successfully' "$LOG_DIR"/*.out 2>/dev/null | wc -l | tr -d ' ' || echo "0"
  echo " jobs completed successfully"
  
  echo ""
  echo "--- Most Recent Error Lines ---"
  for errfile in $(find "$LOG_DIR" -name '*.err' -type f -size +0 2>/dev/null | sort -r | head -5); do
    echo ""
    echo "==> $(basename "$errfile") <=="
    tail -5 "$errfile"
  done
  
  echo ""
  echo "--- Largest Log Files (potential issues) ---"
  find "$LOG_DIR" -type f \( -name '*.out' -o -name '*.err' \) -exec ls -lh {} \; 2>/dev/null | sort -k5 -hr | head -10 || echo "None"
  
} > "$SUMMARY_FILE"

echo "Analysis complete!"
echo "Full report: $SUMMARY_FILE"
echo ""
cat "$SUMMARY_FILE"
