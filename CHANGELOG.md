# Changelog

This file records every change to the analysis that could alter a reported result: the
analysis code, its constants and seeds, its inputs and its dependencies. Design decisions
are argued in `docs/preregistration_decisions.md`, which also carries the dated decisions
from before this file was started on 8 October 2026. Entries here point to it and do not
repeat it. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- The registered two-stage test of H1 and H2 (`docs/preregistration_decisions.md`,
  amendment to section 1). `two_stage_p_negative()` and `detection_rate()` in
  `R/13_simulate_klanguage_pooled.R` apply it, and `scripts/run_klanguage_registered.R`
  writes `outputs/design_summary_pilot/klanguage_registered.csv`: with 20 languages, a
  Bayes factor of 10 in 88% [87, 89] of simulated studies for H1, and in 89% [88, 90] for
  H2 across the 15 languages with a pseudo-passive (2026-10-08).
- `simulate_partly_known_estimates()`, the generation step of `run_partly_known_cell()`,
  so that a second analysis of the reported cells' studies draws exactly the same
  studies. Every K-language table regenerated with it is byte-identical (2026-10-08).
- A 15-of-20 row in `klanguage_pseudo_passive.csv` for the PI's count of 6 October, from
  `scripts/run_klanguage_pseudo_passive.R` (2026-10-07).
- `scripts/build_registered_arc_check.R`, which carries the registered test to ARC's
  R 4.4.2 without writing there. On 2026-10-09 it reproduced `klanguage_registered.csv`
  exactly, and with it the six reported cells it re-creates.
- Run logs for `run_klanguage_registered.R`, `run_klanguage_pseudo_passive.R` and
  `estimate_per_language_approximations.R`, recording the git state, the R version and
  the MD5 of every input, and `config/pilot_models_sha256.csv` for the gitignored pilot
  extracts (2026-10-09).
- Unit tests for the generation step, the two-stage analysis and the shared constants in
  `tests/testthat/test-klanguage-pooled.R`, including the reported H1 rate (2026-10-09).
- `scripts/run_klanguage_registered.R --se-scale=calibrated`, which writes
  `klanguage_registered_se_scaled.csv`: the registered test on the same simulated studies,
  with every language's standard error widened by the calibration factor of the
  per-language approximation (1.165, read from `per_language_approximations.csv`). No
  rate falls by more than 0.4 points (H1 from 87.85% to 87.48%, H2 from 88.85% to
  88.60%), and the report says so. `klanguage_registered.csv` is unchanged (2026-10-09).
- `scripts/run_klanguage_registered.R --n-verbs=72,60,48,36,24`, which writes
  `klanguage_registered_verbs.csv`: the registered test on the same simulated studies,
  with every language using that many verbs, so that the verb part of its standard error
  grows by the square root of 72 over the count, as `se_for()` scales it. With 24 verbs in
  every language, H1 falls from 87.85% to 86.0% and H2 from 88.85% to 88.0%. The rows for
  72 verbs are the registered rates, and `klanguage_registered.csv` is unchanged
  (2026-10-13).
- `scripts/run_klanguage_registered.R --known=none`, which writes
  `klanguage_registered_none_known.csv`: the registered test when none of the measured
  languages takes part, so that every language is drawn from the population estimated
  from them, with the same seeds. H1 falls from 87.85% to 81.2% [80.0, 82.4] and H2
  from 88.85% to 82.3% [81.1, 83.5]. Analysed the earlier way, the H1 studies reproduce
  the "none known" rows of `klanguage_known7.csv` exactly, and `klanguage_registered.csv`
  is unchanged (2026-10-21).
- The same option with `--n-verbs=60,72` and with `--se-scale=calibrated`, which write
  `klanguage_registered_verbs_none_known.csv` and
  `klanguage_registered_se_scaled_none_known.csv`. With 60 verbs per language, H1 falls
  from 81.225% to 81.125% and H2 stays at 82.3%. With widened standard errors, H1 falls
  to 80.8% and H2 to 82.175% (2026-10-22).
- Three options of `scripts/run_klanguage_registered.R`: `--mode=point` holds each
  population's mean and between-language SD at the posterior mean and median, as the
  reported "point" cells held them; `--k=20,25,30` runs those numbers of languages, three
  quarters of them with a pseudo-passive; `--bf-bands` adds the shares of studies whose
  Bayes factor ends at or below 1/10 and between 1/10 and 10. Each adds a word to the
  table's name. Cells that were reported are checked against them, and the default and
  `--known=none` tables are byte-identical (2026-10-24).
- `klanguage_registered_point_bands_none_known.csv` and
  `klanguage_registered_languages_bands_none_known.csv`, from those options with
  `--known=none`. At the estimated mean and spread, H1 reaches a Bayes factor of 10 in
  95.45% of studies and H2 in 97.5%, and no study ends at or below 1/10. Averaged over the
  uncertainty, 0.875% (H1) and 1.225% (H2) end at or below 1/10. With 25 and 30 languages,
  H1 reaches 83.725% at both, and H2 84.9% with 19 and 87.15% with 23 (2026-10-24).
- `scripts/run_klanguage_registered.R --known=recruited`, under which only the measured
  languages with a team in `config/stage2_languages.csv` take part, all but Balinese for
  H1 and English, Turkish, Indonesian and Mandarin for H2, and
  `klanguage_registered_recruited_known.csv` from it: H1 90.65% [89.7, 91.5] and H2
  86.5% [85.4, 87.5]. Both reproduce exactly the check behind the reply of 16 October
  (`reports/correspondence_checks/no_balinese.R`), and the default and `--known=none`
  tables are byte-identical (2026-10-28).
- Three columns in `glossa_tau_refit_summary.csv`, from
  `scripts/summarise_glossa_tau_refit.R`: the share of each verb term's variance common
  to all languages in the refit's verb-split arms, 0.16 (M1) and 0.18 (M2) for the
  active contrast, 0.08 for the intercept and 0.02 for the pseudo-passive contrast. The
  earlier columns are unchanged (2026-10-09).

### Changed

- The report gives the registered rates three ways (`docs/preregistration_decisions.md`,
  section 8, amendment of 2026-10-28): at the estimated effect (H1 95%, H2 97.5%), and
  allowing for uncertainty about it, with no measured language taking part (81%, 82%) or
  with the measured languages that have recruited teams (91%, 86.5%). It adds the share
  of studies ending at or below 1/10 and the H1 rate with 25 or 30 languages, and gives a
  rate that falls exactly on a half to one decimal (2026-10-28).
- The report quotes the registered rates with no measured language taking part, H1 81%
  and H2 82% (`docs/preregistration_decisions.md`, section 8, amendment of 2026-10-22),
  and no longer gives the rates under the other two spreads between languages, at the
  PI's request of 16 and 19 October. Its robustness-check row gives no rate, its
  limitation on widened standard errors is quoted from the matching table, and its
  stopping rule closes data collection six months after the 20th language (2026-10-22).
- Constants that several scripts typed separately are defined once in
  `R/13_simulate_klanguage_pooled.R`: `RHO_EFF_VERB` (0.164), `CLAPS_PILOT_SDX` and
  `REGISTERED_TAU_PRIOR_SCALE` (0.5). The Glossa pseudo-passive terms are read from
  `ALT_PASSIVE_CONTRASTS` and the per-language prior SDs from `config/prior_grid.yaml`.
  No value changed: the regenerated tables, among them `klanguage_known7.csv` and
  `klanguage_pseudo_passive.csv`, are byte-identical (2026-10-09).
- `run_klanguage_pseudo_passive.R` and `estimate_per_language_approximations.R` take
  `--outdir`, and with the other K-language scripts they refuse to write into
  `outputs/design_summary_pilot` without `--overwrite-reported` (2026-10-09).
- The preliminary report presents the pooled tests as confirmatory, analysed by the
  two-stage synthesis, with the single model across languages as a planned robustness
  check (2026-10-08).
- The CI lint gate covers the four scripts behind the registered test as well as `R/`
  (2026-10-09).
- The preliminary report says that 13 of the 120 studies the three-language pooled arm
  planned at 50 to 100 participants ended without a result, 11 of them at 100, and that
  counting them as failures would lower its rate from 57% to 51%. The reported rate is
  unchanged (2026-10-09).
- `docs/preregistration_decisions.md` section 4 describes the Bayes factor as it is
  computed, the posterior odds of the predicted sign over the prior odds, where it named
  a Savage–Dickey density ratio, and withdraws the bridge-sampling calibration, which
  could not calibrate that ratio. The computation is unchanged (2026-10-09).

### Fixed

- The derivation given for the verb correlation now reproduces the value in use and says
  which input was not recorded (`RHO_EFF_VERB`, audit item 8OCT-3). The value is unchanged
  (2026-10-09).
- The verb covariance is indexed by name wherever a verb floor is read, so a reordered
  matrix cannot supply the wrong slope (2026-10-09).
- `compute_bf_exceedance()` in `R/08_summarise_design.R` no longer averages the gender
  arm with the baseline arm at the same design point (finding SEP-3 of the September
  audit). The arm columns join the grouping wherever the cells carry them, and the
  maximal feasible model and the recommended sample size are read from the baseline arm.
  The audit traced the pooling to `outputs/design_summary/`, which fed only the
  archived June report. That table was not regenerated (2026-10-09).
