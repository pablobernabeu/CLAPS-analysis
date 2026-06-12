# Preregistration Decisions — CLAPS Bayesian Ordinal Workflow

This document records every major analysis decision made before data collection and before any confirmatory analysis. Each decision is justified by reference to published methodological or substantive literature; full bibliographic details are in `references.bib`. The companion report `reports/claps_report.qmd` develops the rationale in narrative form.

All decisions below are **locked**. Changes require explicit justification, a version-controlled commit, and documentation in this file.

---

## 1. Primary model family

**Decision**: Bayesian cumulative-ordinal mixed-effects model.

**Formula**:
```
Response ~ S_Type * Semantics_scaled
  + (1 + S_Type * Semantics_scaled | Participant)
  + (1 + S_Type | Verb)
```

**Family**: `cumulative(link = "logit", threshold = "flexible")` in `brms` [@burkner2017; @burkner2018].

**Literature basis**: Treating Likert responses as metric is known to bias regression coefficients and inflate both Type I and Type II error rates [@liddellKruschke2018]. The cumulative-logit model is the principled alternative [@burknerVuorre2019] and has been shown to recover psycholinguistic effects more accurately than Gaussian LMMs in independent empirical comparisons [@verissimo2021; @taylorRousseletScheepersSereno2023].

**Locked**: Yes.

---

## 2. Factor coding

**Decision**: Treatment coding for `S_Type`. Reference level = `Passive`.

**Key estimands**:
- `Semantics_scaled`: affectedness slope for canonical passives (H1a).
- `S_TypeActive:Semantics_scaled`: difference in affectedness slope for actives vs passives (H1b).
- `S_TypePseudo_Passive:Semantics_scaled`: difference for pseudo-passives vs passives (H2; languages with pseudo-passives only).

**Literature basis**: Treatment coding preserves direct interpretability of focal coefficients under planned directional contrasts [@brehmAlday2022].

**Locked**: Yes.

---

## 3. Hypotheses and direction

| Hypothesis | Parameter | Direction | BF threshold |
|-----------|-----------|-----------|-------------|
| H1a | `b_Semantics_scaled` | > 0 | BF > 10 (primary), BF > 3 (secondary) |
| H1b | `b_S_TypeActive:Semantics_scaled` | < 0 | BF > 10 (primary) |
| H2a | `b_S_TypePseudo_Passive:Semantics_scaled` | > 0 | BF > 10 (secondary) |
| H2b | `b_S_TypePseudo_Passive:Semantics_scaled` | < 0 | BF > 10 (secondary) |

**Literature basis (substantive)**: The directional H1a–H1b predictions are grounded in independent cross-language passive-affectedness findings [@ambridgeBidgoodPineRowlandFreudenthal2016; @ambridgeArnonBekman2023; @aryawibawaAmbridge2018; @darmasetiyawanAmbridge2022; @liuAmbridge2021], converging eye-tracking and priming evidence [@paolazziGrilloCeraKarageorgouBullmanChowSanti2022; @bidgoodPineRowlandAmbridge2020; @darmasetiyawanMessengerAmbridge2022], and developmental work showing semantic constraints on the passive [@maratsosFoxBeckerChalkley1985; @pinkerLebeauxFrost1987; @nguyenPearl2021; @agostinhoGavarroSantos2025]. H2 direction is unresolved; both directions are reported.

**Literature basis (BF thresholds)**: Thresholds follow the descriptive grades of @wagenmakersLodewyckxKuriyalGrasman2010; full continuous posterior summaries are reported alongside.

**Locked**: Yes (except H2 direction, which is exploratory).

---

## 4. Bayes-factor computation method

**Primary**: Savage–Dickey directional density ratio [@wagenmakersLodewyckxKuriyalGrasman2010; @schadNicenboimBurknerBetancourtVasishth2023], computed in `brms` with `sample_prior = "yes"`.

**Calibration only**: Bridge sampling [@gronauSarafoglouEtAl2017; @gronauSingmannWagenmakers2020] on a prespecified subset of cells (`scripts/05_bf_calibration_cell.R`).

Bridge sampling is **not** the primary route because the Savage–Dickey ratio is more efficient under directional point-null restrictions and avoids the additional variance contributed by marginal-likelihood estimation [@schadNicenboimBurknerBetancourtVasishth2023].

**Locked**: Yes.

---

## 5. Prior regimes

| Regime | Status | Literature basis |
|--------|--------|------------------|
| `primary` | Primary | @gelmanJakulinPittauSu2008; @gelmanSimpsonBetancourt2017; @burknerVuorre2019 |
| `weak` | Sensitivity | @schadNicenboimBurknerBetancourtVasishth2023 |
| `literature_centred` | Sensitivity-only (not used for primary BF) | Pooled prior means from @ambridgeArnonBekman2023; @aryawibawaAmbridge2018; @darmasetiyawanAmbridge2022; @liuAmbridge2021 |
| `heavy_tailed` | Robustness | @gelmanJakulinPittauSu2008 |

**Variance components**: half-Student-t(3, 0, 1) on group-level standard deviations [@gelman2006; @chungGelmanRabeHeskethLiuDorie2015; @simpsonRueRieblerMartinsSorbye2017]. **Group-level correlations**: LKJ(η = 2) [@lewandowskiKurowickaJoe2009], with η = 1 in the `weak` sensitivity regime.

**Threshold modes**: `broad` (primary; Student-t(3, 0, 2.5) on each ordinal intercept) and `ceiling_calibrated` (per-threshold normals centred on logit cumulative proportions from independent pilot data; for languages with documented ceiling effects only) [@burknerVuorre2019; @schadBetancourtVasishth2021].

**Locked**: Yes.

---

## 6. Random-effects structure and model ladder

Descent order: L5 → L4 → L3 → L2 → L1 → L0.

**Literature basis**: The maximal-effects starting point follows @barrLevyScheepersTily2013; the principled-fallback discipline follows @matuschekKlieglVasishthBaayenBates2017. Generalisation beyond the sampled stimuli and participants is assessed under @westfallYarkoni2016.

Fallback triggers (each must be documented in `outputs/design_summary/ladder_selection.csv`):
- Compilation failure
- OOM (out of memory)
- SLURM timeout
- $\hat{R} \geq 1.01$ on any focal parameter [@vehtariGelmanSimpsonCarpenterBurkner2021] (see §9 amendment)
- any post-warmup divergent transition or maximum-treedepth saturation, or bulk/tail ESS $< 400$ on focal coefficients [@gabrySimpsonVehtariBetancourtGelman2019; @nicenboimSchadVasishth2025]

**Forbidden trigger**: choosing a simpler model because it is easier to fit or produces a more convenient result.

**Locked**: Yes.

---

## 7. Pilot–confirmatory data separation

Pilot data are used **only** for:
- Threshold prior calibration (`compute_ceiling_calibrated_thresholds()`)
- `literature_centred` sensitivity-only prior derivation

They are **never included** in the confirmatory analysis dataset. The `split_pilot_confirmatory()` function in `R/01_read_validate_data.R` enforces this.

**Literature basis**: Including pilot data in the confirmatory model and basing prior centres on the same data inflates the apparent precision of confirmatory estimates, biases effect-size-based design analyses upward [@albersLakens2018], and is incompatible with the prior-conditioning logic of Bayes factors [@schadNicenboimBurknerBetancourtVasishth2023].

**Locked**: Yes.

---

## 8. Design analysis success criteria

**Primary**: $P(BF_{10} > 10)$ exceeds a project-agreed threshold (e.g., 80 %) for H1a and H1b simultaneously.

**Secondary**: $P(BF_{10} > 3)$.

**Sensitivity requirement**: The above probabilities must be robust across the `primary`, `weak`, and `heavy_tailed` prior regimes (no regime reversal of the BF category) following the sensitivity workflow of @schadNicenboimBurknerBetancourtVasishth2023.

**Literature basis**: @albersLakens2018 caution against design analyses driven by single point estimates from small pilots; the prespecified four-regime sensitivity grid mitigates this risk.

**Locked**: Yes (specific probability threshold to be agreed with statistical lead before confirmatory run).

**Amendment (2026-06-09):** the design analysis is now run as a Monte-Carlo *power* analysis with `n_simulations_per_cell = 200` independently seeded replicates per single-language and gender design point (and 50 per cross-language point), replacing the earlier single-seed feasibility run in which $P(BF_{10} > 10)$ could only be a 0/1 indicator. With $B = 200$ the Monte-Carlo standard error of an exceedance proportion near 0.8 is $\approx 0.028$. Per-language power curves are estimated over $N \in \{30, 40, 50, 60\}$ participants at the maximal model under the `primary` prior, with prior-sensitivity replicates (`weak`, `literature_centred`, `heavy_tailed`) at $N = 50$, and the gender variation over the full $N$ sweep. The cross-language (pooled) analysis runs at reduced replication ($B = 50$) on the **L4 cross-uncorrelated** model, which converges far more reliably than the L5 cross-maximal model and costs ~3–5 h rather than 8–12 h per fit; the L5 cross-maximal convergence/feasibility is documented by the one-off heavy run (§9). *Implementation:* `scripts/generate_design_grid.R` expands each condition into seeded replicates (`config/design_grid_{single,gender,cross}.csv`); one SLURM array task per replicate.

---

## 9. Convergence and posterior-pathology diagnostics

All confirmatory fits must pass (publication-grade acceptance — a fit failing any
criterion triggers a ladder fallback rather than being reported):

* $\hat{R} < 1.01$ on all parameters using the rank-normalised, folded definition of @vehtariGelmanSimpsonCarpenterBurkner2021;
* bulk and tail ESS $\geq 400$ on focal coefficients [@vehtariGelmanSimpsonCarpenterBurkner2021];
* zero post-warmup divergent transitions and zero maximum-treedepth saturations [@gabrySimpsonVehtariBetancourtGelman2019; @nicenboimSchadVasishth2025];
* a passed prior-predictive check [@gabrySimpsonVehtariBetancourtGelman2019; @schadBetancourtVasishth2021].

**Locked**: Yes.

**Amendment (2026-06-08):** tightened from $\hat{R} < 1.05$ / $<50$ divergences to the publication-grade $\hat{R} < 1.01$ / zero divergences / zero treedepth saturations of @vehtariGelmanSimpsonCarpenterBurkner2021. The model ladder (`R/09_model_ladder.R`) now selects only fits meeting all criteria; previously "marginal" fits (e.g. $\hat{R}$ up to 1.05) were accepted. *Justification:* ensure robust, publication-grade convergence. *Implementation:* `R/07_extract_diagnostics.R` (`classify_convergence`, `convergence_ok`). Sampling settings (16 chains × 5000 iterations, `adapt_delta = 0.99`, `max_treedepth = 12`) are unchanged — they already provide ~48,000 post-warmup draws.

**Amendment (2026-06-09):** the heavy 16-chain × 5000-iteration sampler is retained only for the one-off maximal-model **convergence demonstration** (establishing the maximal feasible model and runtimes; archived under `outputs/design_analysis_feasibility_v1`). The Monte-Carlo **power** replicates (§8 amendment) use the lighter `replication_sampler` — 4 chains × 3000 iterations, 1000 warmup, `adapt_delta = 0.99`, `max_treedepth = 12` — standard BFDA per-fit settings that yield a reliable single Bayes factor at roughly a fifth of the cost and free cores for ~4× more concurrent fits. Every replicate is still held to the full §9 criteria; non-converging replicates are flagged in the per-cell diagnostics and reported as a convergence rate (column `p_convergence_ok`), never silently included. *Justification:* a stable power estimate requires many replicates, which is infeasible at 48,000 draws/fit; per-fit reliability of a single BF does not. *Implementation:* `config/analysis_config.yaml` (`replication_sampler`), consumed by `scripts/generate_design_grid.R`.

---

## 10. Computational constraints and final model selection

The selected confirmatory model is the **highest feasible model level** from the design analysis. "Feasible" means: converges, finishes within ARC wall time, requires no more than 256 GB RAM (HTC partition), and shows no material posterior pathologies. The principle of preferring the most complex model that converges follows @barrLevyScheepersTily2013, tempered by the parsimony argument of @matuschekKlieglVasishthBaayenBates2017.

Selection is documented in `outputs/design_summary/ladder_selection.csv`.

**Locked**: Yes.

---

## 11. Open issues (to be resolved before confirmatory data collection)

| Issue | Status | Resolution needed by |
|-------|--------|---------------------|
| Exact R module name on ARC | Unresolved | Before first production submission |
| H2 direction for Turkish | Unresolved (both ways reported) | After pilot analysis |
| Specific $P(BF > 10)$ success threshold | To be agreed with statistical lead | Before confirmatory pre-reg |
| `literature_centred` prior centres for languages without pooled estimates | Use weakly informative defaults until single-language data available | Before pilot fits |
| Which languages have pseudo-passives | Partially resolved (Turkish: yes; English, German: no) | Before confirmatory data collection |

---

## References

Full bibliographic details for every cited source are in `references.bib`. The reference list is audited automatically by `scripts/00_verify_references.R`, which fails the build if any DOI, title or author field is missing or unverifiable against Crossref.
