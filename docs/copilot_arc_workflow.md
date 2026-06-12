# CLAPS Bayesian Workflow — Copilot ARC Workflow Guide

This document describes the implemented ARC workflow for GitHub Copilot and future contributors.

## Repository overview

```
R/                  Helper modules (00–10)
scripts/            Runnable pipeline scripts (00–08)
hpc/                SLURM batch and array scripts
config/             YAML and CSV configuration files
reports/            Quarto (.qmd) report sources
docs/               Documentation
references.bib      BibTeX bibliography (audit-verified only)
targets.R           {targets} pipeline definition
renv.lock           Package lockfile
```

## Execution order

1. **Interactive node**: build CmdStan, restore renv, run smoke test
   ```
   srun --partition=devel --cpus-per-task=4 --mem=16G --time=01:00:00 --pty bash
   bash hpc/submit_devel_smoke_test.sh
   ```

2. **Prior predictive checks** (arc, short partition):
   ```
   sbatch hpc/submit_pilot_models_array.sh  # or run prior_predictive_checks first
   ```

3. **Pilot model fitting** (arc, long partition, SLURM array):
   ```
   PILOT_JOB=$(sbatch hpc/submit_pilot_models_array.sh | awk '{print $4}')
   ```

4. **Prior sensitivity** (arc, long partition, afterok pilot):
   ```
   SENS_JOB=$(sbatch --dependency=afterok:$PILOT_JOB \
     hpc/submit_prior_sensitivity_array.sh | awk '{print $4}')
   ```

5. **Design analysis** (arc, long partition, large array):
   ```
   DESIGN_JOB=$(sbatch --dependency=afterok:$SENS_JOB \
     hpc/submit_design_analysis_array.sh | awk '{print $4}')
   ```

6. **BF calibration** (arc, long partition, small subset):
   ```
   CALIB_JOB=$(sbatch --dependency=afterok:$DESIGN_JOB \
     hpc/submit_bf_calibration_array.sh | awk '{print $4}')
   ```

7. **Aggregate** (arc, afterok all):
   ```
   AGG_JOB=$(sbatch --dependency=afterok:$DESIGN_JOB:$CALIB_JOB \
     hpc/submit_aggregate_afterok.sh | awk '{print $4}')
   ```

8. **Render reports** (arc, afterok aggregate):
   ```
   sbatch --dependency=afterok:$AGG_JOB hpc/submit_render_reports.sh
   ```

9. **HTC rescue** (htc, CPU only, 256 GB, for OOM failures):
   ```
   # After design jobs finish, check for OOM failures:
   Rscript scripts/08_submit_status_report.R
   # Then rescue failed cells:
   sbatch --clusters=htc hpc/submit_highmem_rescue_htc_cpu.sh
   ```

## GPU policy

**There are no GPU resources in this project.**
Do not add `--gres=gpu` to any SLURM script.
Stan/CmdStan does not benefit from GPU for the ordinal mixed models used here.

## Module names (unresolved)

The exact R module name on ARC must be verified. Common patterns:
- `R/4.3.1-foss-2023a`
- `R/4.4.0`
- `R`

Set `ARC_R_MODULE` in your environment or update `config/arc_modules.yaml`.

## Memory troubleshooting

| Symptom | Action |
|---------|--------|
| Job killed (OOM, exit 137) | Add row index to `outputs/failed_cells.txt`; resubmit via HTC rescue |
| Job timed out | Resubmit with same cell and longer `--time` on `long` partition |
| Divergences > 50 | Run with increased `adapt_delta` (up to 0.995) via config; document escalation |
| Rhat ≥ 1.05 | Model falls back to next ladder level; document in `_ladder_log.csv` |

## Restart logic

All output scripts:
- Check for existing `.qs` / `.rds` output and skip unless `OVERWRITE=1`
- Write to `.tmp` first, rename on success (atomic write)
- Log job ID, array task ID, Git SHA, seed, cell metadata, and runtime

## Key constraints

- Never run analysis on login nodes
- Never request GPU resources
- Never use `htc` as the primary cluster (only OOM rescue)
- Never choose a simpler model opportunistically
- Never include pilot data in confirmatory inference
- Never invent references or DOIs
