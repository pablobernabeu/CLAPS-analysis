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
- Three columns in `glossa_tau_refit_summary.csv`, from
  `scripts/summarise_glossa_tau_refit.R`: the share of each verb term's variance common
  to all languages in the refit's verb-split arms, 0.16 (M1) and 0.18 (M2) for the
  active contrast, 0.08 for the intercept and 0.02 for the pseudo-passive contrast. The
  earlier columns are unchanged (2026-10-09).

### Changed

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
