#!/usr/bin/env bash
# scripts/poll_arc_status.sh
# Quick HPC status poll that writes outputs/arc_queue_status.txt for automated parsing.
#
# Usage from design_analysis/:
#   bash scripts/poll_arc_status.sh
#   ARC_SSH_HOST=arc bash scripts/poll_arc_status.sh      # override the host
#
# The host defaults to <arc-host>, which is the alias carrying the agent key. The
# script hardcoded "arc" until 2026-09-03; that alias closes the connection here,
# so the script truncated the file it was meant to refresh and left an ERROR line
# where the queue should be. See docs/hpc_polling_workflow.md.
#
# The queue is read from both clusters and filtered to CLAPS jobs. The previous
# version polled one cluster and wrote different section headings from the ones in
# the committed snapshot, so it could not reproduce its own output file.
#
# The poll is assembled in a temporary file and only replaces the committed
# snapshot once it is known to carry queue data. Writing straight to the output,
# as this script used to, means a refused connection overwrites a good record with
# an ERROR line, which is how a working snapshot was lost twice on 2026-09-03. The
# login node also rate-limits repeated connections, so a failed poll is often
# transient and the right response is to keep the previous snapshot and retry.

set -euo pipefail

SSH="${ARC_SSH:-/c/Windows/System32/OpenSSH/ssh.exe}"
ARC_SSH_HOST="${ARC_SSH_HOST:-<arc-host>}"
SACCT_SINCE="${SACCT_SINCE:-2026-08-01}"
OUTPUT_FILE="outputs/arc_queue_status.txt"
mkdir -p "$(dirname "$OUTPUT_FILE")"
TMP_FILE="$(mktemp "${OUTPUT_FILE}.XXXXXX")"
trap 'rm -f "$TMP_FILE"' EXIT

# Keep the cluster banners and the column header, then CLAPS jobs only.
claps_filter() { grep -E '^CLUSTER:|JOBID|claps' || true; }

{
  echo "=== ARC Queue Status ==="
  echo "Timestamp: $(date -Iseconds 2>/dev/null || date +%Y-%m-%dT%H:%M:%S%z)"
  echo ""

  echo "--- Active Queue (CLAPS jobs only, both clusters) ---"
  "$SSH" "$ARC_SSH_HOST" \
    'squeue -M all -u "$(whoami)" --format="%.18i %.16j %.8T %.10M %.10l %.6D %R"' \
    2>/dev/null | claps_filter || echo "ERROR: squeue failed"

  echo ""
  echo "--- Other users' jobs on this account are omitted; non-CLAPS jobs of this user are omitted ---"
  echo ""
  echo "--- Recent CLAPS Job Accounting (since ${SACCT_SINCE}) ---"
  "$SSH" "$ARC_SSH_HOST" \
    "sacct -M all -u \"\$(whoami)\" --starttime=${SACCT_SINCE} -X \
       --format=JobID%22,JobName%20,State%20,End%17,Elapsed%10,ExitCode%8" \
    2>/dev/null | grep -E 'JobID|^-|claps' || echo "ERROR: sacct failed"

  echo ""
  echo "--- Completed Design Analysis Outputs (local) ---"
  find outputs/design_analysis -name "*.rds" -type f 2>/dev/null | wc -l || echo "0"
  echo " outputs found locally"
} > "$TMP_FILE" 2>&1

# A poll counts as usable only if the queue section came back with real content.
if grep -qE '^CLUSTER:|claps' "$TMP_FILE" && ! grep -q "ERROR: squeue failed" "$TMP_FILE"; then
  mv "$TMP_FILE" "$OUTPUT_FILE"
  echo "[poll] Status written to: $OUTPUT_FILE"
else
  echo "[poll] Poll failed; $OUTPUT_FILE left unchanged. Collected output was:" >&2
  sed -n '1,12p' "$TMP_FILE" >&2
  rm -f "$TMP_FILE"
  exit 1
fi
