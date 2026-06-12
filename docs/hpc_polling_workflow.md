# HPC Polling Workflow

## Root cause
Windows process isolation + VS Code terminal context prevents agent-run terminals from accessing your unlocked SSH agent socket, even with correct environment variables.

## Solution
Two-script workflow:

### 1. Poll (you run in your authenticated terminal)
```bash
cd design_analysis
bash scripts/poll_arc_status.sh
```
Writes to `outputs/arc_queue_status.txt`

### 2. Parse (agent runs automatically)
```bash
bash scripts/parse_arc_status.sh
```
Reads the status file and shows summary.

## Setup (one-time)
```bash
# Add your key to agent (enter passphrase once)
ssh-add ~/.ssh/id_ed25519_arc

# Verify it worked
ssh-add -l

# Test connection
ssh arc "echo Connected"
```

## Fetch and Analyze Logs

### Fetch logs from ARC
```bash
cd design_analysis
bash scripts/fetch_arc_logs.sh          # Get 50 most recent logs
bash scripts/fetch_arc_logs.sh --recent 100  # Get 100 most recent
```

### Analyze local logs
```bash
bash scripts/analyze_logs.sh
```
Checks for: OOM failures, timeouts, R errors, Stan issues, divergences, completions.

## Scripts Created
- `scripts/poll_arc_status.sh` — Fetch queue status (you run)
- `scripts/parse_arc_status.sh` — Parse status file (agent runs)
- `scripts/fetch_arc_logs.sh` — Sync logs from ARC (you run)
- `scripts/analyze_logs.sh` — Analyze local logs (agent runs)

## Changes made
- `.bashrc`: Disabled automatic `ssh-add` on every shell startup
- `.ssh/config`: Changed `AddKeysToAgent` to `ask` (interactive sessions only)
