# ARC Submission Guide — CLAPS Bayesian Ordinal Workflow

## Overview

This guide covers:
1. Interactive node setup (software builds only)
2. Long-partition production runs
3. HTC CPU high-memory rescue
4. Monitoring and restart logic
5. Memory troubleshooting

---

## 1. Interactive Node Builds

**Never run analysis or model fitting on login nodes.**

Request an interactive session on the `devel` partition for software setup:

```bash
srun --partition=devel \
     --nodes=1 --ntasks=1 --cpus-per-task=4 \
     --mem=16G --time=01:00:00 --pty bash
```

Once on the interactive node:

```bash
# Load R module (verify name with: module spider R)
module purge
module load R   # or: module load R/4.3.1-foss-2023a

# Restore renv library
Rscript -e "renv::restore(prompt = FALSE)"

# Install CmdStan if not available as a module
Rscript -e "cmdstanr::install_cmdstan(cores = 4)"

# Verify CmdStan
Rscript -e "cmdstanr::check_cmdstan_toolchain(); cat(cmdstanr::cmdstan_version(), '\n')"

# Run smoke test
bash hpc/submit_devel_smoke_test.sh
```

---

## 2. Long-Partition Production Runs

All production model fits use `--partition=long` with explicit long wall times.
**The ARC default 1-day wall time is not adequate for ordinal mixed models.**

Example: submit the design analysis array:

```bash
sbatch hpc/submit_design_analysis_array.sh
```

This requests:
- Partition: `long`
- Wall time: `7-00:00:00` (7 days)
- CPUs: 16 per task
- Memory: 64 GB

### Adjusting array size

Before submitting, check `config/design_grid.csv`:

```bash
N=$(tail -n +2 config/design_grid.csv | wc -l)
echo "Design grid has $N rows"
```

Update `#SBATCH --array=1-${N}%50` in `hpc/submit_design_analysis_array.sh` to match.

### Chained submission (recommended)

```bash
# Step 1: pilot models
PILOT_JOB=$(sbatch hpc/submit_pilot_models_array.sh | awk '{print $4}')
echo "Pilot job: $PILOT_JOB"

# Step 2: prior sensitivity (after pilot)
SENS_JOB=$(sbatch \
  --dependency=afterok:$PILOT_JOB \
  hpc/submit_prior_sensitivity_array.sh | awk '{print $4}')

# Step 3: design analysis (after sensitivity)
DESIGN_JOB=$(sbatch \
  --dependency=afterok:$SENS_JOB \
  hpc/submit_design_analysis_array.sh | awk '{print $4}')

# Step 4: calibration (after design)
CALIB_JOB=$(sbatch \
  --dependency=afterok:$DESIGN_JOB \
  hpc/submit_bf_calibration_array.sh | awk '{print $4}')

# Step 5: aggregate (after design + calibration)
AGG_JOB=$(sbatch \
  --dependency=afterok:$DESIGN_JOB:$CALIB_JOB \
  hpc/submit_aggregate_afterok.sh | awk '{print $4}')

# Step 6: render (after aggregate)
sbatch --dependency=afterok:$AGG_JOB hpc/submit_render_reports.sh

echo "Full pipeline submitted."
```

---

## 3. HTC CPU High-Memory Rescue

Use `htc` **only** when arc jobs fail with OOM (exit code 137 or "killed" status).
HTC is CPU-only in this project. No GPU resources are requested.

### Identify OOM failures

```bash
# Check sacct for killed jobs
sacct -j $DESIGN_JOB \
  --format=JobID,ArrayTaskID,State,ExitCode,MaxRSS \
  --noheader | grep -E "FAILED|OUT_OF_ME|137"
```

### Populate failed cells list

```bash
# Extract failed array task IDs and append to failed_cells.txt
sacct -j $DESIGN_JOB \
  --format=ArrayTaskID,State --noheader \
  | awk '$2 ~ /FAILED|CANCELLED/ {print $1}' \
  >> outputs/failed_cells.txt
sort -u outputs/failed_cells.txt -o outputs/failed_cells.txt
```

### Submit rescue

```bash
N_FAILED=$(wc -l < outputs/failed_cells.txt)
sbatch --clusters=htc \
       --array=1-${N_FAILED}%5 \
       hpc/submit_highmem_rescue_htc_cpu.sh
```

---

## 4. Monitoring

```bash
# Job status
squeue -u $USER --format="%.18i %.10j %.8T %.10M %.10l %.6D %R"

# Array task breakdown
squeue -j $DESIGN_JOB --array --format="%.18i %.10T %.10M"

# Check progress against design grid
Rscript scripts/08_submit_status_report.R
# Output: outputs/job_status_report.csv

# Tail a specific task log
tail -f outputs/logs/design_analysis_${DESIGN_JOB}_42.out
```

---

## 5. Memory Troubleshooting

| Exit code | Likely cause | Action |
|-----------|-------------|--------|
| 137       | OOM kill    | Add to failed_cells.txt; rescue on HTC |
| 1         | R error     | Check `.err` log; fix R code; resubmit |
| Timeout   | Slow model  | Resubmit on `long` with larger `--time` |
| Divergences | Sampler issue | Increase `adapt_delta` in config (max 0.995); document |
| Rhat ≥ 1.05 | Convergence | Model falls back to next ladder level automatically |

### Increasing adapt_delta

Edit `config/analysis_config.yaml`:

```yaml
model:
  adapt_delta: 0.995   # increased from 0.99
```

Then resubmit with `OVERWRITE=1`:

```bash
OVERWRITE=1 sbatch hpc/submit_design_analysis_array.sh
```

---

## 6. Checking Output Completeness

```bash
# Count completed design cells
ls outputs/design_analysis/*.qs | wc -l

# Compare against expected grid size
tail -n +2 config/design_grid.csv | wc -l

# Full status report
Rscript scripts/08_submit_status_report.R
```

---

## Notes

- Module names (`R`, `cmdstan`, `quarto`) must be verified on your ARC system with `module spider`.
- See `config/arc_modules.yaml` for the documented known-unknowns.
- The `TMPDIR` variable on ARC points to fast local scratch; `CMDSTANR_OUTPUT_DIR` uses it by default.
