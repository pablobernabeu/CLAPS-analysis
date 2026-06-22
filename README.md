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

Power is governed mainly by the number of verbs in the materials, not by the number of
participants. At the design's real verb count, a modest per-language sample of at most
about 50 participants is sufficient, with both focal predictions reaching a detection
probability of about 0.96 to 0.98 at 50 participants. At fewer verbs the design is
borderline (40 verbs) or insufficient (12 to 20 verbs). The primary prediction (the
active × affectedness interaction, H1b) is tested one-tailed, supported by the affectedness
slope itself (H1a). The recommendation rests primarily on a data-grounded analysis that
simulates from each language's fitted pilot, currently running on the cluster, with a
literature-anchored analysis as a cross-check. An earlier figure of *N* = 80 was
conditional on a smaller, fixed verb count and is superseded. Full operating
characteristics are in the current report under `reports/`
(`CLAPS_preliminary_sample_size_analysis_<date>.pdf`).

## Repository structure

```
R/         Analysis modules (model formulas, priors, hypothesis tests, ...)
scripts/   Runnable pipeline scripts (grid generation, fitting, aggregation, ...)
hpc/       SLURM batch and array scripts (University of Oxford ARC)
config/    YAML and CSV configuration, including the design grids
reports/   Rendered report PDF and citation style (apa.csl)
outputs/   Aggregated design-analysis result tables (CSV)
docs/      Pipeline, submission and preregistration documentation
references.bib, renv.lock, targets.R
```

## The report

The current report is provided pre-rendered as
`reports/CLAPS_preliminary_sample_size_analysis_<date>.pdf`. It draws on the private pilot
data for its descriptive figure and sample-composition appendix, so it cannot be
regenerated from this public repository. The aggregated result tables it summarises are
committed under `outputs/`.

The report's pilot-data figure (a gender robustness check) is also provided separately as a
pre-rendered image (`outputs/figures/raw_gender_comparison.png`). The script that builds it
from the raw pilot data is included for transparency but cannot be run without the withheld
data.

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
