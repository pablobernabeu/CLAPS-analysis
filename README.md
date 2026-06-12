# CLAPS — Bayesian Ordinal Design and Power Analysis

Public analysis materials for the CLAPS study: the Bayesian cumulative-ordinal
modelling workflow, the Bayes-factor design (power) analysis and the resulting
reports.

> **Data availability.** Raw pilot and participant data are **not** included in this
> repository, for confidentiality. The design-analysis results here come from
> simulated data sets (a Monte-Carlo design analysis), so the summary tables under
> `outputs/` and the rendered reports need no participant data. Scripts that process
> the pilot data (for example `scripts/00_harmonise_pilot_data.R` and
> `scripts/02_fit_pilot_models.R`) are included for transparency but cannot be run
> without the withheld data.

## What the analysis does

- **Model**: Bayesian cumulative (ordered-logit) mixed-effects model in `brms` /
  `cmdstanr`.
- **Inference**: directional Savage–Dickey Bayes factors, with bridge sampling as a
  calibration check.
- **Design analysis**: a Monte-Carlo Bayes-factor design analysis estimating, by
  simulation, the probability of decisive evidence at each sample size, prior regime
  and model-ladder level.

## Headline result

At the assumed effect the recommended per-language sample size is **N = 80**
(English, Norwegian and Turkish); a pooled cross-language analysis reaches the 80%
target at roughly N = 60 per language. See `reports/claps_report.pdf`.

## Repository structure

```
R/         Analysis modules (model formulas, priors, hypothesis tests, ...)
scripts/   Runnable pipeline scripts (grid generation, fitting, aggregation, ...)
hpc/       SLURM batch and array scripts (University of Oxford ARC)
config/    YAML and CSV configuration, including the design grids
reports/   Quarto report sources, rendered PDFs and generated appendices
outputs/   Aggregated design-analysis result tables (CSV)
docs/      Pipeline, submission and preregistration documentation
references.bib, renv.lock, targets.R
```

## Reproducing the reports

The full report renders from the committed summary tables, with no participant data:

```bash
quarto render reports/claps_report.qmd --to pdf
```

The report's one pilot-data figure (the gender robustness check) is included as a
pre-rendered aggregate image (`outputs/figures/raw_gender_comparison.png`); the
script that builds it from the raw pilot data is included for transparency but, like
the other pilot-data scripts, cannot be run without the withheld data.

## Focal hypotheses

| ID    | Parameter                          | Direction          |
|-------|------------------------------------|--------------------|
| H1a   | affectedness slope (passive)       | > 0                |
| H1b   | active × affectedness interaction  | < 0                |
| H2a/b | pseudo-passive × affectedness      | both (secondary)   |

## Prior regimes

| Regime               | Role                                          |
|----------------------|-----------------------------------------------|
| `proposal`           | Moderately informative; design-analysis default |
| `weak`               | Wider-prior sensitivity check                 |
| `literature_centred` | Direction-encoding sensitivity check          |
| `heavy_tailed`       | Student-t robustness check                    |

## Compute

The fits were run CPU-only on the University of Oxford ARC cluster. The `hpc/`
scripts use placeholder identifiers (`PROJECT_ID`, `PROJECT_GROUP`) and absolute
paths in place of a specific cluster allocation; set these for your own project
before use.
