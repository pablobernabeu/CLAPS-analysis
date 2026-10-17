# CLAPS: analysis materials for the Registered Report

CLAPS (A Crosslinguistic Approach to Passive Semantics) asks whether the acceptability of
passive sentences tracks how strongly the event affects the person it happens to, more than
the acceptability of actives and pseudo-passives does, and whether this holds across languages.
Teams in at least 20 languages will each collect acceptability judgments from 80 participants
and affectedness ratings from 20 more. This repository holds the analysis registered in the
Stage 1 Registered Report, the design analysis behind its sampling plan, and the report of
that design analysis.

> **Data.** No participant data are included here. The design analysis rests on simulated data,
> and the summaries it reads are under `outputs/`. Scripts that process the pilot data are
> included for transparency but cannot be run without them.

## The hypotheses

| | Coefficient | Prediction | Tested on |
|---|---|---|---|
| H1 | Active by affectedness interaction | Negative: the affectedness slope is steeper for passives than for actives | The mean across all languages |
| H2 | Pseudo-passive by affectedness interaction | Negative: the slope is steeper for passives than for pseudo-passives | The mean across the languages with a pseudo-passive |

Each language is analysed with a Bayesian cumulative-logit mixed-effects model, and the
per-language estimates are combined in a Bayesian random-effects meta-analysis. Each hypothesis
is decided by a directional Bayes factor: 10 or more is strong evidence for it, 1/10 or less
strong evidence that the interaction runs the other way, and anything in between is
inconclusive. In the code and output files, H1 and H2 carry the labels `h1b` and `h2b`.

## The registered analysis

`stage2/` holds the confirmatory analysis as six numbered scripts, run in this order from the
repository root:

| Order | Script | What it does |
|---|---|---|
| 1 | `stage2/01_check_and_exclude.R` | Checks the teams' data and applies the registered exclusions |
| 2 | `stage2/02_affectedness.R` | Builds the affectedness predictor and runs the outcome-neutral check |
| 3 | `stage2/03_fit_language_models.R` | Fits each language's model by the registered ladder |
| 4 | `stage2/04_meta_analysis.R` | Tests H1 and H2 |
| 5 | `stage2/05_robustness_single_model.R` | Robustness check: a single model across languages |
| 6 | `stage2/06_robustness_families.R` | Robustness check: languages nested within families |

`stage2/README.md` describes the input format, the decisions each step makes and the cluster
scripts, and `stage2/dry_run/` runs every step end to end on simulated data. The functions are
in `R/14_stage2.R`, and their tests in `tests/testthat/test-stage2.R`. The registered settings
are in the `stage2` block of `config/analysis_config.yaml`.

## The design analysis

With 20 languages, 15 of them with a pseudo-passive, and 80 participants per language, the
registered test reaches a directional Bayes factor of 10 in 88% of simulated studies for H1 and
89% for H2 (`outputs/design_summary_pilot/klanguage_registered.csv`, from
`scripts/run_klanguage_registered.R`). These rates take the spread of the effect between
languages from the seven languages measured so far: the three CLAPS pilots (English, Norwegian
and Turkish) and the four further languages of Ambridge, Arnon and Bekman (2023). The spread is
what limits the test, and the report gives the rates under other plausible spreads.

## The report

`reports/CLAPS_preliminary_sample_size_analysis_<date>.pdf` is the current report, and
`reports/preliminary_sample_size_analysis.qmd` its source. The public source reads only the
summaries committed under `outputs/`, a disclosure-reviewed raster of its first figure and
aggregate sample counts, so it rebuilds without participant data. With R, Quarto and a LaTeX
distribution installed, restore the locked R environment and render from the repository root:

```sh
Rscript -e "renv::restore(prompt = FALSE)"
quarto render reports/preliminary_sample_size_analysis.qmd
```

The rebuild should reproduce the content, tables and pagination of the committed PDF, though not
byte for byte: the source carries the render date, TeX records creation metadata, and the public
build embeds the raster of Figure 1 where the private build draws it from the pilot ratings.

## Repository structure

```
stage2/    The registered Stage 2 analysis, numbered in run order
R/         Analysis modules (model formulas, priors, ladder, meta-analysis, ...)
scripts/   Runnable scripts of the design analysis (fitting, simulation, aggregation, ...)
hpc/       SLURM batch and array scripts (University of Oxford ARC)
config/    YAML and CSV configuration, including the design grids and the language table
reports/   Report source, rendered PDF and citation style
outputs/   Aggregated design-analysis results (CSV)
tests/     Unit tests
docs/      Pipeline, submission and preregistration documentation
references.bib, renv.lock, targets.R
```

## Priors

Every reported result uses the `primary` prior regime in `R/03_define_priors.R`: zero-centred
normal priors on the focal slopes, so that neither sign is favoured. The `weak`,
`literature_centred` and `heavy_tailed` regimes in the same file were used only to check the
design analysis's sensitivity to the priors.

## Compute

The fits run CPU-only on the University of Oxford's ARC cluster. The `hpc/` scripts use
placeholder identifiers (`PROJECT_ID`, `PROJECT_GROUP`) in place of a specific allocation; set
these for your own project before use.
