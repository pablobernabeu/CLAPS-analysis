# ARC Submission Guide — CLAPS Bayesian Ordinal Workflow

## Overview

This guide covers:
1. Interactive node setup (software builds only)
2. Long-partition production runs
3. HTC CPU high-memory rescue
4. Monitoring and restart logic
5. Memory troubleshooting

Read the section below first if you are trying to understand what an existing
submission script does, rather than to submit one.

---

## 0. Anatomy of a submission script

Every script in `hpc/` shares the same preamble. It is repeated rather than
factored into a common file because SLURM copies the submitted script to a spool
directory before running it, so a `source`d helper alongside it would not
reliably be found. What each part is for:

**`#SBATCH` block.** The resource request. `%A` and `%a` in the log paths expand
to the job ID and the array task ID, giving one log pair per task; without `%a`
every task in an array would write over the same file. `--mail-type=FAIL,END`
reports the array as a whole, not each task.

**`set -euo pipefail`.** Abort on an unset variable, on any command failure, and
on a failure anywhere in a pipeline. Without it a missing module or a failed
`module load` would be reported only as a warning and the job would proceed to
run R against the wrong environment, producing results that look valid. The CI
check in `.github/workflows/static-checks.yaml` enforces its presence.

**`cd "$HOME/design_analysis"`.** All scripts anchor to an absolute path rather
than deriving one from `$BASH_SOURCE`, because at run time `$BASH_SOURCE` points
into the non-writable SLURM spool copy, not the repository.

**`module purge` then `module load`.** Purge first so the environment does not
depend on whatever the submitting shell happened to have loaded, which is the
usual cause of a job that runs for one person and not another. The R module is
read from `$ARC_R_MODULE` with a pinned default, so a cluster-wide module change
does not silently alter the R version behind published results.

**`$DATA` guard.** Fails immediately with a clear message if the project storage
variable is unset. Everything below writes under `$DATA`, so continuing without it
would scatter outputs into the home quota and fill it.

**`R_LIBS_USER` and `RENV_PATHS_CACHE`.** Point R at the project library and the
renv cache on project storage. Home directories are too small for a Stan
toolchain, and a shared cache means parallel array tasks do not each rebuild
packages.

**`OMP_NUM_THREADS` and `STAN_NUM_THREADS`.** Both are set from
`$SLURM_CPUS_PER_TASK` so the process uses exactly the cores it was allocated.
This is the single most common cause of an HPC job being slower than expected:
without it, the threading libraries default to the *physical core count of the
node*, and every task on a shared node oversubscribes it. `production_sampling()`
in `R/04_model_formulas.R` reads `STAN_NUM_THREADS` for its `cores` argument, which
is how the R code inherits the allocation.

**`CMDSTANR_OUTPUT_DIR`.** Sends Stan's intermediate CSV output to `$TMPDIR`, on
node-local disk. Stan writes one CSV per chain per fit; on shared storage this is
slow and, across a large array, enough to disturb the filesystem for other users.

**`CMDSTAN` discovery.** `ls -d … | sort -V | tail -1` picks the highest installed
CmdStan version (`sort -V` compares version strings, so 2.36 sorts after 2.9).
The `stanc` executable check that follows turns a missing or half-installed
toolchain into an immediate, named failure instead of a compilation error inside R
several minutes later.

**The `$ROW -gt $N_ROWS` clean exit.** An array is often submitted with a range
wider than the grid, or split across accounts. A task beyond the last grid row
exits 0 with a message rather than failing, so an over-wide `--array` does not
fill the mail queue with spurious failures. `tail -n +2` skips the CSV header when
counting rows, which is why grid rows are 1-based and `--array` must start at 1.

**The `awk` echo of the grid row.** Prints the cell's parameters into the log, so
the log alone identifies what the task ran. Note that these `awk` field indices
are positional, so a change to a grid's column order will silently mislabel the
log line without affecting the analysis, which reads the CSV by column name.

**`${OVERWRITE:+--overwrite}`.** Passes the flag only when `$OVERWRITE` is set and
non-empty, so a resubmission can force recomputation with `OVERWRITE=1 sbatch …`
without a separate script.

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
ls outputs/design_analysis/*.rds | wc -l

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
