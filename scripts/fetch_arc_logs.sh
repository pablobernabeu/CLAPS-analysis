#!/usr/bin/env bash
# scripts/fetch_arc_logs.sh
# Sync SLURM logs from ARC to local outputs/logs/ for analysis.
# Usage from design_analysis/: bash scripts/fetch_arc_logs.sh [--recent N]

set -euo pipefail

SSH="/c/Windows/System32/OpenSSH/ssh.exe"
RSYNC_SSH="/c/Windows/System32/OpenSSH/ssh.exe"

# Parse arguments
RECENT_LIMIT=50
if [[ "${1:-}" == "--recent" ]] && [[ -n "${2:-}" ]]; then
  RECENT_LIMIT="$2"
fi

LOCAL_LOG_DIR="outputs/logs"
SUMMARY_FILE="outputs/log_fetch_summary.txt"

mkdir -p "$LOCAL_LOG_DIR"

echo "========================================"
echo "  Fetching ARC Logs"
echo "========================================"
echo ""

# Get list of recent logs on ARC - search multiple possible locations
echo "[1/4] Detecting remote working directory..."
REMOTE_WD=$($SSH arc 'echo $HOME')
echo "Remote home: $REMOTE_WD"

echo ""
echo "[2/4] Searching for log directories..."
# Search common locations: ~/design_analysis/outputs/logs, ~/outputs/logs
SEARCH_PATHS=(
  "design_analysis/outputs/logs"
  "outputs/logs"
)

FOUND_PATH=""
for path in "${SEARCH_PATHS[@]}"; do
  if $SSH arc "test -d ~/${path}" 2>/dev/null; then
    FOUND_PATH="$path"
    echo "✓ Found logs at: ${REMOTE_WD}/${path}"
    break
  fi
done

if [[ -z "$FOUND_PATH" ]]; then
  echo "ERROR: Could not find outputs/logs directory on ARC"
  echo "Searched:"
  for path in "${SEARCH_PATHS[@]}"; do
    echo "  - ~/$path"
  done
  echo ""
  echo "Check where you submitted jobs from and ensure logs were created."
  exit 1
fi

echo ""
echo "[3/4] Scanning logs in ${FOUND_PATH}..."
REMOTE_FILES=$($SSH arc "cd ~/ && find ${FOUND_PATH} -name '*.out' -o -name '*.err' | sort -r | head -${RECENT_LIMIT}" 2>/dev/null || echo "")

if [[ -z "$REMOTE_FILES" ]]; then
  echo "No log files found in arc:~/${REMOTE_LOG_DIR}"
  echo "Have you submitted any jobs yet?"
  exit 0
fi

LOG_COUNT=$(echo "$REMOTE_FILES" | wc -l | tr -d ' ')
echo "Found $LOG_COUNT recent log files on ARC"

# Sync logs to local
echo ""
echo "[4/4] Syncing logs..."
rsync -avz -e "$RSYNC_SSH" \
  "arc:~/${FOUND_PATH}/" \
  "${LOCAL_LOG_DIR}/" \
  --include='*.out' \
  --include='*.err' \
  --exclude='*' \
  2>&1 | grep -E '(sending|receiving|sent|received|speedup)' | tail -5 || true

echo ""
echo "Generating summary..."

# Create summary
{
  echo "=== ARC Log Fetch Summary ==="
  echo "Timestamp: $(date -Iseconds 2>/dev/null || date +%Y-%m-%dT%H:%M:%S%z)"
  echo "Logs fetched: $LOG_COUNT"
  echo ""
  echo "--- Most Recent Logs ---"
  find "$LOCAL_LOG_DIR" -type f \( -name '*.out' -o -name '*.err' \) -printf '%T+ %p\n' 2>/dev/null | sort -r | head -20 || \
    ls -lt "$LOCAL_LOG_DIR"/*.{out,err} 2>/dev/null | head -20 || echo "No local logs yet"
  
  echo ""
  echo "--- Error Log Preview (last 10 lines of .err files) ---"
  for errfile in $(find "$LOCAL_LOG_DIR" -name '*.err' -type f 2>/dev/null | sort -r | head -5); do
    if [[ -s "$errfile" ]]; then
      echo ""
      echo "==> $errfile <=="
      tail -10 "$errfile"
    fi
  done
  
  echo ""
  echo "--- Job Status from Logs ---"
  grep -h -E '(COMPLETED|FAILED|TIMEOUT|ERROR|Exit code|slurmstepd)' "$LOCAL_LOG_DIR"/*.out 2>/dev/null | tail -20 || echo "No status lines found"
  
} > "$SUMMARY_FILE"

echo ""
echo "========================================"
echo "✓ Logs synced to: $LOCAL_LOG_DIR"
echo "✓ Summary written to: $SUMMARY_FILE"
echo "========================================"
echo ""
cat "$SUMMARY_FILE"
