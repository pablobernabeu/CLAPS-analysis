# Design / Power Analysis Pipeline

How the Bayes-factor **power** analysis is generated, submitted, aggregated, and
rendered on ARC. See `preregistration_decisions.md` §8–§9 for the locked
methodological decisions this implements.

## Two-stage design

1. **Feasibility / convergence study** (one-off, heavy sampler — 16 chains ×
   5000 iter). Establishes the *maximal feasible model* and per-fit runtimes per
   language. Run once; results archived under
   `$DATA/PROJECT_GROUP/outputs/design_analysis_feasibility_v1` (+ matching
   `design_summary_feasibility_v1`).
2. **Powered Monte-Carlo study** (this pipeline). Each design *condition* is
   replicated across `n_simulations_per_cell` independently seeded simulated data
   sets so that `P(BF₁₀ > threshold)` is a genuine power estimate, not a 0/1
   indicator. Replicates use the lighter `replication_sampler` (4 chains × 3000
   iter, `adapt_delta = 0.99`); every replicate is still convergence-checked.

## Config knobs (`config/analysis_config.yaml`)

| Key | Value | Meaning |
|-----|-------|---------|
| `design_analysis.n_simulations_per_cell` | 200 | replicates / single-language & gender point (MC SE ≈ 0.028 at p = 0.8) |
| `design_analysis.n_simulations_per_cell_cross` | 50 | replicates / cross-language point (8–12 h/fit) |
| `replication_sampler` | 4 × 3000 (1000 warmup) | per-replicate sampler |
| `model` (iter 5000, chains 16) | — | heavy sampler, feasibility study only |

## Step 1 — generate grids

```bash
module load R/4.4.2-gfbf-2024a
export R_LIBS_USER=$DATA/PROJECT_GROUP/R/library_4.4
Rscript scripts/generate_design_grid.R --out_dir config
```

Produces three grids (each < cluster `MaxArraySize` = 5000), one SLURM array task
per row:

| Grid | Rows | Conditions |
|------|------|-----------|
| `config/design_grid_single.csv` | 4200 | 3 langs × (N∈{30,40,50,60} primary + {weak,literature_centred,heavy_tailed}@N50) × 200 |
| `config/design_grid_gender.csv` | 2400 | 3 langs × N∈{30,40,50,60} × gender variation × 200 |
| `config/design_grid_cross.csv`  | 200  | AllLanguages L4 cross-uncorrelated × N∈{30,40,50,60} × 50 |

## Step 2 — submit the chained pipeline

Three design arrays (resources overridden at submit time; `GRID` selects the
grid), then aggregation (`afterany` — runs even if some replicates fail), then
report rendering (`afterok`). Run from `$HOME/design_analysis`, with `$DATA` set:

```bash
JS=$(sbatch --parsable --partition=short --time=12:00:00 --cpus-per-task=4 --mem=16G \
     --array=1-4200%40 --export=ALL,GRID=config/design_grid_single.csv \
     hpc/submit_design_analysis_array.sh)
JG=$(sbatch --parsable --partition=short --time=12:00:00 --cpus-per-task=4 --mem=16G \
     --array=1-2400%40 --export=ALL,GRID=config/design_grid_gender.csv \
     hpc/submit_design_analysis_array.sh)
JC=$(sbatch --parsable --partition=short --time=12:00:00 --cpus-per-task=4 --mem=64G \
     --array=1-200%12  --export=ALL,GRID=config/design_grid_cross.csv \
     hpc/submit_design_analysis_array.sh)
JA=$(sbatch --parsable --dependency=afterany:$JS:$JG:$JC hpc/submit_aggregate_afterok.sh)
JR=$(sbatch --parsable --dependency=afterok:$JA            hpc/submit_render_reports.sh)
echo "single=$JS gender=$JG cross=$JC aggregate=$JA render=$JR"
```

## Outputs (project storage, not the repo)

`$DATA/PROJECT_GROUP = /data/PROJECT_GROUP/PROJECT_ID/PROJECT_GROUP`

- per-cell `.rds` → `$DATA/PROJECT_GROUP/outputs/design_analysis/`
- summary CSVs → `$DATA/PROJECT_GROUP/outputs/design_summary/`
  (`bf_exceedance.csv` = power curves; `recommended_sample_size.csv`,
  `maximal_feasible_model.csv`, `runtime_summary.csv`, `failure_summary.csv`)
- rendered report → `$DATA/PROJECT_GROUP/outputs/reports/` (qmd reads CSVs via
  `CLAPS_OUTPUTS_ROOT`, exported by `submit_render_reports.sh`).

## Monitoring

```bash
ssh <arc-host> "sacct -j <JOBID> --format=JobID,State,Elapsed --noheader | grep -vE 'batch|extern'"
```

A SLURM array reports `COMPLETED` even when a cell catches a fit error and exits 0:
always confirm cell health from the `.rds` `$summary$status` / `$diagnostics$convergence_ok`,
or from `failure_summary.csv` after aggregation. `.out` logs are block-buffered
(no Stan iteration lines until the fit ends); use rising SLURM elapsed + clean
`.err` as the liveness signal.

## Gotchas baked into the scripts

- All sbatch scripts anchor to `$HOME/design_analysis` (not `$(dirname BASH_SOURCE)`,
  which resolves to a non-writable SLURM spool copy) and export
  `R_LIBS_USER`/`RENV_PATHS_CACHE` to the `$DATA` tree.
- Cross-language fits need more memory (`--mem=64G`); single-language run fine at 16G.
- **Quarto** is not an ARC module. A standalone build is installed at
  `$DATA/PROJECT_GROUP/quarto-1.9.38/` and linked to `~/bin/quarto` (on PATH). Its
  bundled `deno` **SIGTRAPs on the login node** (hardened environment) but renders
  fine on **compute nodes**, so report rendering only works inside a SLURM job
  (or `srun`), never via a bare login-node `quarto render`. `submit_render_reports.sh`
  puts `~/bin` on PATH and sets `QUARTO_PATH`/`TMPDIR` accordingly. To render
  ad hoc: `srun --partition=devel --time=00:09:00 --mem=8G bash -lc '… quarto render …'`.
