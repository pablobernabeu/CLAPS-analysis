# Stage 2: the registered analysis

These scripts run the confirmatory analysis registered in the Stage 1 Registered Report,
in this order, from `design_analysis/`. Each implements one part of the report's Methods,
and every setting they use is fixed in the `stage2` block of `config/analysis_config.yaml`,
with the priors in `R/03_define_priors.R` (primary regime). The functions they call are in
`R/14_stage2.R`, and the model ladder is walked by `R/09_model_ladder.R`.

| Order | Script | What it does | Stage 1 report |
|---|---|---|---|
| 1 | `01_check_and_exclude.R` | Checks the teams' data, applies the three exclusion criteria to both tasks, keeps 80 participants and 20 raters per variety | Exclusion criteria; Sampling plan |
| 2 | `02_affectedness.R` | Takes each verb's mean affectedness rating, runs the outcome-neutral check, centres the means and divides them by two SDs | Pre-processing; outcome-neutral check |
| 3 | `03_fit_language_models.R` | Fits each variety's model by the ladder, from the maximal model down, and extracts both interactions per SD of affectedness | Individual-language models |
| 4 | `04_meta_analysis.R` | Averages each language's varieties and tests H1 and H2 by the Bayesian random-effects meta-analysis and its directional Bayes factors | Meta-analytic synthesis |
| 5 | `05_robustness_single_model.R` | Fits the single model across languages by its ladder | Robustness check |
| 6 | `06_robustness_families.R` | Repeats the meta-analysis with languages nested within Glottolog families | Second robustness check |

Steps 3 and 5 are the long ones and run on the cluster: `hpc/submit_stage2_language_fits.sh`
(one array task per analysed variety) and `hpc/submit_stage2_single_model.sh`. Step 5 needs
only steps 1 and 2, so it can run alongside step 3. Steps 1, 2, 4 and 6 take minutes on a
laptop.

```bash
Rscript stage2/01_check_and_exclude.R
Rscript stage2/02_affectedness.R
sbatch --account=<project> --array=1-<number of analysed varieties> hpc/submit_stage2_language_fits.sh
sbatch --account=<project> hpc/submit_stage2_single_model.sh
Rscript stage2/04_meta_analysis.R        # after every step-3 task has finished
Rscript stage2/06_robustness_families.R
```

## How the decisions are made

- **Exclusions.** A participant is excluded for a missing rating or trial, for answering
  fewer than 8 of the 10 attention checks with their target, or for giving the same rating
  on 90% or more of the test trials. Teams replace excluded participants, so a variety is
  analysed once it has 80 usable participants and 20 usable raters, and at least 60 verbs.
  If it has more, the first 80 (or 20) in order of data collection are kept, which needs
  `sessions.csv` (below).
- **Outcome-neutral check.** A variety whose verbs' mean affectedness ratings have an SD
  below 0.5 scale points is left out and reported, since its affectedness slope cannot be
  estimated. The check reads only the ratings.
- **Convergence.** A fit is accepted only with R-hat below 1.01, bulk and tail effective
  sample sizes of at least 400 on every population-level coefficient, threshold,
  random-effect SD and correlation, no divergent transitions and no saturated tree depth
  (`R/07_extract_diagnostics.R`). The ladder descends only when a level errors, runs out
  of memory or wall time, or fails these criteria. A variety with no converging level is
  left out of the meta-analysis and listed in `results/left_out.csv`. A level killed by the
  wall time is resumed at the next level with `--start-level` (or `START_LEVEL=` in the
  `hpc/` scripts).
- **Varieties.** Varieties of one language, as defined by Glottolog in
  `config/stage2_languages.csv`, are analysed separately and averaged into one estimate,
  with the mean of their posterior SDs as its standard error. Each language counts once.
- **The tests.** The per-language interactions, per SD of affectedness (half the model's
  coefficient), enter a random-effects meta-analysis with a normal(0, 0.25) prior on the
  mean for H1 and normal(0, 0.3) for H2, and a half-normal(0.5) prior on the
  between-language SD, computed by the bayesmeta package (Röver, 2020). H1 uses every
  language and H2 the languages whose data include a pseudo-passive. A directional Bayes
  factor of 10 or more is strong evidence for the hypothesis, 1/10 or less strong evidence
  that the interaction runs the other way, and anything in between is inconclusive.
- **Robustness.** The single model and the family model are reported beside the tests and
  decide neither hypothesis.

## Input format

Each team's data arrive as two CSV files in `data/stage2/raw/`, one with `acceptability`
and one with `ratings` in its name, one row per trial, attention checks included:

| Column | Acceptability task | Ratings task |
|---|---|---|
| `Variety` | Team label, as in `config/stage2_languages.csv` | Same |
| `Participant` | Participant code within the variety | Rater code within the variety |
| `Trial` | Trial number | Same |
| `Trial_Type` | `test` or `attention` | Same |
| `Verb` | Verb from the pool (empty on attention trials) | Same |
| `Agent` | `Man` or `Woman` (empty on attention trials) | Same |
| `S_Type` | `Active`, `Passive` or `Pseudo_Passive` (empty on attention trials) | Not used |
| `Attention_Target` | 1 or 7 on attention trials, empty otherwise | Same |
| `Response` | Rating, 1 to 7; empty if not given | Same |

An optional `sessions.csv` in the same directory, with `Variety`, `Participant` and
`Session_Order`, gives the order in which participants completed the task. The formr
exports are converted to this format before step 1; the participant codes must carry no
Prolific or formr identifiers.

## Outputs

`outputs/stage2/processed/` holds the analysed trials and the affectedness scores,
`outputs/stage2/fits/` the selected fits and the diagnostics of every level tried, and
`outputs/stage2/results/` the tables: the exclusions and samples, the outcome-neutral
check, the per-variety estimates, the per-language estimates, the two tests
(`meta_analysis.csv`) and the two robustness checks.

## Dry run

`stage2/dry_run/simulate_data.R` writes small simulated data sets that exercise every rule
above: exclusions of each kind, a team that has not finished, a variety that fails the
outcome-neutral check, a language tested in two varieties and one without a pseudo-passive.
Every script then runs with `--dry-run`, which reads and writes `outputs/stage2_dry_run/`
and uses smaller samples, 12 verbs and a shorter sampler. Nothing else changes.

```bash
Rscript stage2/dry_run/simulate_data.R
Rscript stage2/01_check_and_exclude.R --dry-run
Rscript stage2/02_affectedness.R --dry-run
Rscript stage2/03_fit_language_models.R --all --dry-run
Rscript stage2/04_meta_analysis.R --dry-run
Rscript stage2/05_robustness_single_model.R --dry-run
Rscript stage2/06_robustness_families.R --dry-run
```

The unit tests in `tests/testthat/test-stage2.R` check each rule on constructed data.
With `CLAPS_SLOW_TESTS=true` they also check that the bayesmeta meta-analysis agrees with
the engine of the power simulations to within 0.001 in the posterior probability that the
mean is negative.
