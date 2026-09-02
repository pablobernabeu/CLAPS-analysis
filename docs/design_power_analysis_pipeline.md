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
| `design_analysis.n_simulations_per_cell_cross` | 50 | replicates / cross-language point (29–96 h/fit, measured) |
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

## Seed allocation across grids

Every replicate is independently seeded, and the seed is the only thing that
distinguishes one replicate of a design point from another. Each grid generator
therefore allocates from its own numeric base, so that two grids do not hand the
same RNG stream to two different cells. A shared stream would mean the two cells
draw identical data whenever their design parameters also coincide, which would
present one simulated data set as two independent replicates.

Every base below is exclusive to one generator. The ranges are pairwise disjoint,
and `tests/testthat/test-seed-disjointness.R` fails if that ever stops being true.

| Base | Range | Generator / grid |
|------|-------|------------------|
| 1e5 | 100000–102399 | `generate_design_grid.R` power curves → `design_grid_single.csv` |
| 2e5 | 200000–201799 | `generate_design_grid.R` prior sensitivity → `design_grid_single.csv` |
| 3e5 | 300000–302399 | `generate_design_grid.R` gender → `design_grid_gender.csv` |
| 4e5 | 400000–402399 | `generate_design_grid.R` extended N → `design_grid_extend.csv` |
| 5e5 | from 500000 | `generate_feasibility_grid.R` |
| 6e5 | 600000–… | `generate_corrected_power_grid.R` |
| 7e5 | 700000–701049 | `generate_floor50_power_grid.R` |
| 8e5 | 800000–… | `generate_databased_grid.R` |
| 9e5 | 900000–900199 | `generate_design_grid.R` cross → `design_grid_cross.csv` |
| 9.1e5 | 910000–… | `generate_databased_v2_extendedN_grid.R` |
| 9.29e5 | 929000 | `generate_pooled_v2_grid.R` smoke test → `design_grid_pooled_v2_TEST.csv` |
| 9.3e5 | 930000–931490 | `generate_pooled_v2_grid.R` |
| 9.32e5 | 932000–932290 | `generate_pooled_v2_N80_grid.R` |
| 9.4e5 | 940000–941079 | `generate_databased_v2_decision_grid.R` |
| 9.49e5 | 949000 | `generate_pooled_model_comparison_grid.R --smoke` → `design_grid_pooled_modelcmp_SMOKE.csv` |
| 9.5e5 | 950000–951190 | `generate_pooled_model_comparison_grid.R` |
| 9.59e5 | 959000 | `generate_pooled_model_comparison_grid.R --smoke --verbblock` → `design_grid_pooled_verbblock_SMOKE.csv` (not committed) |
| 9.6e5 | 960000–960490 | `generate_pooled_model_comparison_grid.R 50 --verbblock`; at the default 120 groups it would reach 961190 |
| 1.0e6 | 1000000–1002099 | `generate_corrected_scale_grid.R` — superseded arm, see below |
| 1.1e6 | 1100000–1100439 | `generate_safeguard_grid.R` — superseded arm, see below |
| 1.2e6 | 1200000–1200719 | `generate_databased_v2_grid.R` |
| 1.3e6 | 1300000–1303790 | `generate_pooled_v2_N80_precision_grid.R` |

Bases from 1.4e6 upward are free. Leave a gap of at least 1e5 above a new base, so
that raising a grid's replicate count later cannot run it into its neighbour.

Four of the rows above name files kept only in the private repository, so a reader of
the public mirror will not find them there. They are `generate_pooled_model_comparison_grid.R`
with its two committed grids, and the two smoke grids at 9.29e5 and 9.49e5. The
model-comparison arm has never run and no figure in the report comes from it. Its ranges are
recorded here because the registry has to reserve every base that is already spoken for,
whether or not the grid holding it is published.

Smoke grids are the one exception, sitting 1000 below the production base they
rehearse (9.29e5 under 9.3e5, 9.49e5 under 9.5e5, 9.59e5 under 9.6e5), so that a
smoke cell and the real cell it rehearses are never confused. A smoke grid is one
group by construction and cannot grow into its neighbour. The grids that can grow
are the production ones, upward into whichever base sits next above them, and two
of those have less headroom than the 1e5 rule implies. At spacing 10 the 9.5e5
block fits 900 groups before it reaches the 9.59e5 smoke base, which leaves 780
beyond its current default of 120; the 9.3e5 block fits 200 before it reaches
9.32e5, and 150 are already committed. Raising either replicate count near those
limits means allocating a fresh base rather than extending in place.

Two early grids predate the scheme and allocate below 1e5. `design_grid.csv` holds
59 seeds scattered between 1001 and 9012, and `design_grid_fast.csv` 39 between
11001 and 13005, in small per-condition blocks rather than the dense runs above.
Both were written by hand rather than by a generator. They are recorded here
because the test compares every committed grid, not only those with a base in the
table.

### Two superseded arms

`generate_safeguard_grid.R` and `generate_corrected_scale_grid.R` keep their seed
bases in the registry above, but neither is part of the reported analysis. Commit
3c5a872 (2026-07-24) abandoned them as "literature-anchored refinement arms" that had
been "promised as forthcoming, but their jobs have been held (admin/user) and idle for
a month and are superseded", and reworded the report accordingly.

Their held cluster jobs were cancelled on 2026-07-30: `claps_safeguard` 12409496
(tasks 123–440, 122 outputs already on disk) and `claps_corrscale` 12419137 (which had
produced nothing). Nothing reads `outputs/design_safeguard` or
`outputs/design_corrected_scale` — the four aggregators read `design_databased_v2`,
`design_pooled_v2`, `design_corrected` and `design_analysis` respectively — and the
report cites neither directory.

Beware a name collision that has misled more than once: the report's "safeguard"
figures come from the **safeguard mode** of the pilot-grounded engine
(`R/10_simulate_from_pilot_v2.R`, reaching the report via
`outputs/design_summary_pilot/joint_power_pilot.csv`), not from the safeguard array.
The generators are retained so the record of what was tried survives, and so either
arm can be revived; if one is, its seeds will no longer match anything on disk.

### The 2026-07-30 re-basing

Two collisions existed until 2026-07-30 and are now resolved. `corrected_scale`,
`floor50` and `safeguard` all began at 700000, so the latter two ranges sat inside
the first; and `databased_v2` began at 900000, covering the 200 seeds already used
by the committed `design_grid_cross.csv`. Comments in those generators had described
each base as collision-free.

No published result was affected. The grids committed under `config/` were checked
pairwise and share no seeds, and the computed outputs on ARC were checked for
identical cell IDs across directories, of which there were none: the overlapping
seeds always went with a different verb count, so no two cells ever produced the
same simulated data.

Which generator moved was decided by how many cells each already had computed,
because the seed forms part of every output filename and re-basing therefore
orphans existing outputs:

| Generator | Old base | New base | Computed cells | Orphaned |
|---|---|---|---|---|
| `generate_floor50_power_grid.R` | 7e5 | 7e5 (kept) | 1028, in `outputs/design_corrected` | 0 |
| `generate_corrected_scale_grid.R` | 7e5 | 1.0e6 | 0, never run | 0 |
| `generate_safeguard_grid.R` | 7e5 | 1.1e6 | 122, in `outputs/design_safeguard` | 122 |
| `generate_databased_v2_grid.R` | 9e5 | 1.2e6 | 719, in `outputs/design_databased_v2` | 719 |

The 841 orphaned cells are the unavoidable minimum: one of the three 7e5 grids had
to move, and `databased_v2` had to move because the cross grid it collided with is
committed. Orphaned outputs are not deleted, only no longer matched by the
resume-by-existing-output check, so regenerating either grid and resubmitting
recomputes those cells under the new seeds. The 374 extended-N cells sharing the
`design_databased_v2` directory use base 9.1e5 and are unaffected.

### The 2026-09-02 re-basing

The test that the 2026-07-30 work introduced had been failing on two further
collisions since 2026-08-19. `main` was red throughout, and the failure went unread.

`generate_pooled_v2_N80_precision_grid.R`, added on 2026-08-19, took base 9.4e5,
which the table above had already assigned to `generate_databased_v2_decision_grid.R`.
The decision grid runs consecutively over 940000–941079, so 108 of the precision
grid's 380 rows fell inside it. Each row hands seed, seed+1 and seed+2 to the three
languages, so 324 distinct seed values collided. The precision generator's comment
described its base as clear of the pooled grids, which it was; what it did not do
was consider the decision grid.

`generate_pooled_model_comparison_grid.R` applied its `--smoke` branch to the group
count and the sampler budget but not to the seed base, so the single smoke group took
950000, which is group 1 of the committed full grid. That mattered more than the
count of one shared seed suggests. `R/11_simulate_pooled_v2.R` builds the output name
from the prior source, mode, sample size, model level and seed, so a smoke fit at 500
iterations and the real 3000-iteration cell resolve to the same file, and the
resume-by-existing-output check would have accepted whichever was written first.

No published result was affected, and nothing was orphaned: neither arm had been run.

| Generator | Old base | New base | Computed cells | Orphaned |
|---|---|---|---|---|
| `generate_databased_v2_decision_grid.R` | 9.4e5 | 9.4e5 (kept) | 1080, on ARC under `$DATA/outputs/design_databased_v2_decision`; none held locally | 0 |
| `generate_pooled_v2_N80_precision_grid.R` | 9.4e5 | 1.3e6 | 0, arrays 8607044 and 8607045 cancelled before any task started | 0 |
| `generate_pooled_model_comparison_grid.R --smoke` | 9.5e5 | 9.49e5 | 0, never run; the one matching three-task array was cancelled the same day, before any task started | 0 |

The direction was decided by the rule applied on 2026-07-30, that the grid with fewer
computed cells moves. The decision arm supplies all of the 80-participant cell the
recommendation rests on and 120 of the 144 replicates at 70 and 100, while the
precision arm has produced nothing.

Both generators previously restated other grids' ranges in a hand-copied list. Every
range in those lists was correct; what neither covered was the whole registry, and the
precision generator's omission of the decision grid is what let the collision through.
Those two comments now point at the table above instead of restating it.
`generate_databased_v2_decision_grid.R`, `generate_pooled_v2_grid.R` and
`generate_pooled_v2_N80_grid.R` still carry lists of their own, and are worth the same
treatment when they are next touched.

The same commit gave `--smoke --verbblock` its own output path. It had shared the
`--smoke` path, so running both flags overwrote the committed model-comparison smoke
grid with verbblock arms.

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
- Quarto is not an ARC module. A standalone build is installed at
  `$DATA/PROJECT_GROUP/quarto-1.9.38/` and linked to `~/bin/quarto` (on PATH). Its
  bundled `deno` hits a SIGTRAP on the hardened login node but renders fine on
  compute nodes, so report rendering only works inside a SLURM job
  (or `srun`) and never from a bare login-node `quarto render`. `submit_render_reports.sh`
  puts `~/bin` on PATH and sets `QUARTO_PATH`/`TMPDIR` accordingly. To render
  ad hoc, use `srun --partition=devel --time=00:09:00 --mem=8G bash -lc '… quarto render …'`.
