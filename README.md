# CLAPS: Bayesian Ordinal Design and Power Analysis

Public analysis materials for the CLAPS study. CLAPS asks whether the acceptability
of passive and related sentence types tracks semantic affectedness, and whether that
relationship differs by sentence type, across English, Norwegian and Turkish. This
repository holds the Bayesian cumulative-ordinal modelling workflow, the Bayes-factor
design (power) analysis and the resulting reports.

> **Data availability.** The raw pilot and participant data are **not** included in
> this repository. They are held private to the CLAPS project until it reaches a more
> advanced stage. The design-analysis results here come from simulated data sets (a
> Monte-Carlo design analysis), so the summary tables under `outputs/` and the
> rendered reports need no participant data. Scripts that process the pilot data (for
> example `scripts/00_harmonise_pilot_data.R` and `scripts/02_fit_pilot_models.R`) are
> included for transparency but cannot be run without the withheld data.

## What the analysis does

- **Model**: Bayesian cumulative (ordered-logit) mixed-effects model in `brms` /
  `cmdstanr`.
- **Inference**: directional Savage–Dickey Bayes factors, with bridge sampling as a
  calibration check.
- **Design analysis**: a Monte-Carlo Bayes-factor design analysis estimating, by
  simulation, the probability of decisive evidence at each sample size, prior regime
  and model-ladder level.

## Headline result

At the assumed effect, the recommended per-language sample size is **_N_ = 80** for
English, Norwegian and Turkish. The binding constraint is the active × affectedness
interaction (H1b), because the affectedness slope itself (H1a) is already at ceiling.
A pooled cross-language analysis, which borrows strength across the three languages,
reaches the 80% target at roughly *N* = 60 per language. The full operating
characteristics are in the date-stamped report under `reports/`
(`CLAPS_design_analysis_<date>.pdf`).

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

The committed report PDF is date-stamped, `CLAPS_design_analysis_<date>.pdf`.

The report's one pilot-data figure (the gender robustness check) is included as a
pre-rendered aggregate image (`outputs/figures/raw_gender_comparison.png`). The script
that builds it from the raw pilot data is included for transparency but, like the
other pilot-data scripts, cannot be run without the withheld data.

## Focal hypotheses

| ID    | Parameter                         | Direction          |
|-------|-----------------------------------|--------------------|
| H1a   | Affectedness slope (passive)      | > 0                |
| H1b   | Active × affectedness interaction | < 0                |
| H2a/b | Pseudo-passive × affectedness     | Both (secondary)   |

## Prior regimes

| Regime               | Role                                                                |
|----------------------|---------------------------------------------------------------------|
| `primary`            | Weakly-to-moderately informative, the prespecified analysis prior   |
| `weak`               | Wider scales, a sensitivity check                                   |
| `literature_centred` | Direction-encoding, a sensitivity check (shown as LC in the report) |
| `heavy_tailed`       | Student-t scales, a robustness check                                |

## Compute

The fits were run CPU-only on the University of Oxford ARC cluster. The `hpc/` scripts
use placeholder identifiers (`PROJECT_ID`, `PROJECT_GROUP`) and absolute paths in place
of a specific cluster allocation. Set these for your own project before use.
