#!/usr/bin/env bash
# scripts/parse_arc_status.sh
# Parse the HPC status file and provide a summary.
# Usage from design_analysis/: bash scripts/parse_arc_status.sh

set -euo pipefail

STATUS_FILE="outputs/arc_queue_status.txt"

if [[ ! -f "$STATUS_FILE" ]]; then
  echo "ERROR: No status file found at $STATUS_FILE"
  echo "Run: bash scripts/poll_arc_status.sh"
  exit 1
fi

echo "========================================"
echo "  ARC HPC Status Summary"
echo "========================================"
echo ""

# Extract timestamp
TIMESTAMP=$(grep "^Timestamp:" "$STATUS_FILE" | cut -d' ' -f2- || echo "unknown")
echo "Poll time: $TIMESTAMP"
echo ""

# Count active jobs
ACTIVE_COUNT=$(awk '/^--- Active Queue ---$/,/^$/ {if (NF>0 && !/^---|^Timestamp|^===|^JOBID|ERROR/) print}' "$STATUS_FILE" | wc -l | tr -d ' ')
echo "Active queue jobs: $ACTIVE_COUNT"

# Show running jobs
echo ""
echo "--- Running Jobs ---"
awk '/^--- Active Queue ---$/,/^$/ {if ($3 ~ /^R/ || $3 ~ /RUNNING/) print}' "$STATUS_FILE" | head -20

# Show pending jobs 
echo ""
echo "--- Pending Jobs ---"
awk '/^--- Active Queue ---$/,/^$/ {if ($3 ~ /^PD/ || $3 ~ /PENDING/) print}' "$STATUS_FILE" | head -10

# Show recent failures
echo ""
echo "--- Recent Failures (last 24h) ---"
awk '/^--- Recent Job Accounting/,/^$/ {if ($3 ~ /FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY/) print}' "$STATUS_FILE" | head -15

# Local output count
echo ""
echo "--- Local Outputs ---"
LOCAL_COUNT=$(find outputs/design_analysis -name "*.qs" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "Completed analysis outputs: $LOCAL_COUNT"

echo ""
echo "========================================"
echo "Full details: $STATUS_FILE"
echo "========================================"
