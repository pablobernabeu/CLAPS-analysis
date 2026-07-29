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
participants. In the data-grounded analysis at the current verb counts, Norwegian reaches
the conventional 80% joint-power target at about 100 participants, while English and
Turkish plateau below it as participants are added. A pooled cross-linguistic analysis is
therefore the preferred feasible design. Its first cell, at 50 participants per language,
has 96% joint power at the moderate Bayes-factor threshold of 3 but 36% at the
strong-evidence threshold of 10; larger pooled cells are still awaited. The independent
literature-anchored cross-check is more optimistic, reaching target at or below 50
participants with the full verb set. Full operating characteristics and the remaining
design choices are in the current report under `reports/`
(`CLAPS_preliminary_sample_size_analysis_<date>.pdf`).

## Repository structure

```
R/         Analysis modules (model formulas, priors, hypothesis tests, ...)
scripts/   Runnable pipeline scripts (grid generation, fitting, aggregation, ...)
hpc/       SLURM batch and array scripts (University of Oxford ARC)
config/    YAML and CSV configuration, including the design grids
reports/   Reproducible public QMD, rendered report PDF and citation style
outputs/   Aggregated design-analysis result tables (CSV)
docs/      Pipeline, submission and preregistration documentation
references.bib, renv.lock, targets.R
```

## The report

The current report is provided as both
`reports/preliminary_sample_size_analysis.qmd` and the pre-rendered
`reports/CLAPS_preliminary_sample_size_analysis_<date>.pdf`. The public QMD retains the
report's text, figures, tables and appendices, while remaining independent of
participant-level pilot data. It reads the simulation summaries committed under
`outputs/`, embeds the disclosure-reviewed report figure at
`outputs/figures/pilot_acceptability_by_gender_report.png`, and reads aggregate
sample-composition counts from
`outputs/design_summary_pilot/pilot_sample_composition.csv`.

A second and separate image, `outputs/figures/raw_gender_comparison.png`, is the standalone
gender robustness check produced by `scripts/12_plot_raw_gender_comparison.R`. It does not
appear in the report. The scripts that build both images from the raw pilot data are
included for transparency but cannot be run without the withheld data. The private report
source generates the report figure and the sample counts directly from the pilot data,
whereas the public source uses only the committed derivatives.

### Rebuilding the public report

With R, Quarto and a LaTeX distribution installed, restore the locked R environment and
render from the repository root:

```sh
Rscript -e "renv::restore(prompt = FALSE)"
quarto render reports/preliminary_sample_size_analysis.qmd
```

The rebuild requires no participant-level data. It should reproduce the content, tables and
pagination of the committed PDF alongside it, which is the reference to compare against. It
is not expected to be byte-identical. The source carries the render date, TeX records
creation metadata, and the public build embeds a disclosure-reviewed raster of Figure 1
where the private build generates that figure directly from pilot rows.

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
