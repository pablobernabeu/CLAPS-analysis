# Preregistration Decisions — CLAPS Bayesian Ordinal Workflow

This document records every major analysis decision made before data collection and before any confirmatory analysis. Each decision is justified by reference to published methodological or substantive literature; full bibliographic details are in `references.bib`. The companion report `reports/preliminary_sample_size_analysis.qmd` develops the rationale in narrative form.

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

**Amendment (2026-10-08), proposed to the PI with his draft Stage 1 protocol:** the hypotheses are tested across languages, as the PI's draft states, and the confirmatory test is a two-stage synthesis. Each language is analysed with the model above. The per-language posterior means and SDs of each interaction are then combined in a Bayesian random-effects meta-analysis with Language as the random effect, a normal(0, *P*) prior on the mean across languages (*P* the per-language prior SD per SD of affectedness: 0.25 for H1, 0.30 for H2) and a half-normal(0, 0.5) prior on the between-language SD. The one-stage model across languages in the PI's draft (by-Language slopes on all fixed effects, participants and verbs nested within language) becomes a planned robustness check, reported but deciding neither hypothesis. *Why*: the one-stage model saturated tree depth on three quarters of its iterations when fitted to the pilot data (`outputs/logs/pooledfit_8614954.*`), and in the feasibility runs at 20 languages no rung above `L2_cross_stype_participant_verb` (two of three data sets) converged; the two-stage synthesis is the analysis the power simulations evaluate (`scripts/run_klanguage_registered.R`), and one- and two-stage analyses agree when each study is large and both are specified alike [@burkeEnsorRiley2017].

**Locked**: Yes (the amendment awaits the PI's agreement).

---

## 2. Factor coding

**Decision**: Treatment coding for `S_Type`. Reference level = `Passive`.

**Key estimands**:
- `Semantics_scaled`: affectedness slope for passives (reported; no hypothesis, see §3).
- `S_TypeActive:Semantics_scaled`: difference in affectedness slope for actives vs passives (H1).
- `S_TypePseudo_Passive:Semantics_scaled`: difference for pseudo-passives vs passives (H2; languages with pseudo-passives only).

**Literature basis**: Treatment coding preserves direct interpretability of focal coefficients under planned directional contrasts [@brehmAlday2022].

**Locked**: Yes.

---

## 3. Hypotheses and direction

| Hypothesis | Parameter | Direction | BF threshold |
|-----------|-----------|-----------|-------------|
| H1 (primary) | mean across languages of `b_S_TypeActive:Semantics_scaled` (§1 amendment) | < 0 | BF ≥ 10 |
| H2 (secondary) | mean across the languages with a pseudo-passive of `b_S_TypePseudo_Passive:Semantics_scaled` | < 0 | BF ≥ 10 |

**Amendment (2026-10-06), agreed with the PI:** only these two hypotheses are registered, and they are named H1 and H2. The passive slope (formerly H1a) carries no prediction and need not be positive, since the theory requires only that affectedness matter more for passives than for the comparison sentence type. The reversed pseudo-passive direction (formerly H2a) is not a hypothesis. A Bayes factor of 1/10 or less counts as evidence against a hypothesis. Turkish is tested on both, like every language. The code and output files keep the earlier labels (`h1b`, `h2b`; `H1b_active_interaction_negative`, `H2b_pseudo_negative`), since they are written into the design grids and the replicate files.

**Literature basis (substantive)**: The directional predictions are grounded in independent cross-language passive-affectedness findings [@ambridgeBidgoodPineRowlandFreudenthal2016; @ambridgeArnonBekman2023; @aryawibawaAmbridge2018; @darmasetiyawanAmbridge2022; @liuAmbridge2021], converging eye-tracking and priming evidence [@paolazziGrilloCeraKarageorgouBullmanChowSanti2022; @bidgoodPineRowlandAmbridge2020; @darmasetiyawanMessengerAmbridge2022], and developmental work showing semantic constraints on the passive [@maratsosFoxBeckerChalkley1985; @pinkerLebeauxFrost1987; @nguyenPearl2021; @agostinhoGavarroSantos2025]. H2 is registered in the predicted direction only.

**Literature basis (BF thresholds)**: Thresholds follow the descriptive grades of @wagenmakersLodewyckxKuriyalGrasman2010; full continuous posterior summaries are reported alongside.

**Locked**: Yes.

---

## 4. Bayes-factor computation method

**Primary**: a directional Bayes factor, the posterior odds that the coefficient has the predicted sign divided by the same odds under the prior [@moreyWagenmakers2014]. The hypothesis is the half of the zero-centred prior on the predicted side, and its competitor is the half on the other side. Within a language, the Bayes factor is computed in `brms` from the posterior draws and from prior draws sampled alongside them (`sample_prior = "yes"`). In the two-stage synthesis, it is computed on the mean across languages from the meta-analytic posterior, where the prior odds are exactly one because the prior on the mean is centred at zero. The functions in `R/05_hypothesis_tests.R` keep the name `savage_dickey_directional_bf()` from an earlier description, and they compute this ratio.

**Amendment (2026-10-09):** the description above replaces one that called the Bayes factor a Savage–Dickey density ratio under directional point-null restrictions, with bridge sampling [@gronauSarafoglouEtAl2017; @gronauSingmannWagenmakers2020] as a calibration on a subset of cells (`scripts/05_bf_calibration_cell.R`). The computation is unchanged. The calibration is withdrawn: bridge sampling compares the full model with an intercepts-only model, which is a different quantity from the sign-odds ratio, so it cannot calibrate it.

**Misleading evidence (recorded 2026-10-07, for reviewers):** because the Bayes factor compares the predicted sign with its reverse, a true effect of exactly zero still yields a Bayes factor of 10 in the predicted direction in about 6% to 8.5% of studies per language and hypothesis, and as often in the opposite direction. The figure comes from a normal approximation calibrated against the decision arm (`scripts/estimate_per_language_approximations.R`, `outputs/design_summary_pilot/per_language_approximations.csv`), not from simulation, and the preliminary report does not quote it.

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

**Amendment (2026-10-12):** the list above omits the design analysis, which takes its effects and each language's precision from the pilot models. The published per-language estimates of @ambridgeArnonBekman2023 serve the design analysis as well, and nothing else. Neither the pilot data nor the published data enter the confirmatory analyses, which use only new data collected under the registered protocol.

**Literature basis**: Including pilot data in the confirmatory model and basing prior centres on the same data inflates the apparent precision of confirmatory estimates, biases effect-size-based design analyses upward [@albersLakens2018], and is incompatible with the prior-conditioning logic of Bayes factors [@schadNicenboimBurknerBetancourtVasishth2023].

**Locked**: Yes.

---

## 8. Design analysis success criteria

**Primary**: $P(BF_{10} \geq 10)$ for H1 and for H2, each in its own right, in the two-stage synthesis at the planned sample: 88% for H1 with 20 languages and 89% for H2 with the 15 of them that have a pseudo-passive, on the spread across the measured languages (`outputs/design_summary_pilot/klanguage_registered.csv`, amended 2026-10-08 from per-language rates when the hypotheses became tests across languages).

**Amendment (2026-10-06):** the joint H1a-and-H1b criterion and the secondary $P(BF_{10} > 3)$ are withdrawn with H1a (§3). The proposal to the PI of 1 October 2026 replaces the power threshold with the route *Nature Human Behaviour* offers Bayesian plans with limited resources: a maximum feasible sample of 80 participants per language, fixed in advance, at which a test that has reached neither 10 nor 1/10 is reported as inconclusive. The design analysis then reports detection rates as information rather than as a requirement. The PI asked on 6 October what the route means, and his draft protocol of 8 October takes it; the closing date it needs is in §11.

**Sensitivity requirement**: The above probabilities must be robust across the `primary`, `weak`, and `heavy_tailed` prior regimes (no regime reversal of the BF category) following the sensitivity workflow of @schadNicenboimBurknerBetancourtVasishth2023.

**Amendment (2026-10-07):** prior sensitivity moves to the confirmatory analysis. At Stage 2 each confirmatory Bayes factor is also reported under the `weak` and `heavy_tailed` regimes, and the BF category is checked for reversal there. The pilot-grounded design analysis uses the `primary` regime only. The summaries in `outputs/design_summary/prior_sensitivity.csv` come from the earlier literature-anchored grid and do not bear on it.

**Literature basis**: @albersLakens2018 caution against design analyses driven by single point estimates from small pilots; the prespecified four-regime sensitivity grid mitigates this risk.

**Locked**: Yes.

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
| Exact R module name on ARC | Resolved 2026-07-30: `R/4.4.2-gfbf-2024a`, verified on ARC and pinned in `config/arc_modules.yaml` | Settled |
| H2 direction | Resolved 2026-10-06 with the PI: one-tailed and negative in every language with a pseudo-passive, Turkish included although its pilot points the other way (§3). The simulated rates are in `outputs/design_summary_pilot/pseudo_passive_decision_arm.csv` | Settled |
| Specific $P(BF > 10)$ success threshold | Replaced by the maximum-feasible route (§8), which the PI's draft protocol of 8 October takes; only the closing date remains (below) | Settled |
| `literature_centred` prior centres for languages without pooled estimates | Resolved: the regime centres every language on the same pooled values (`R/03_define_priors.R`) and is sensitivity-only, never used for the primary Bayes factors | Settled |
| Which languages have pseudo-passives | Resolved (English: yes; Turkish: yes; Norwegian: no; Balinese: yes), as set by `has_pseudo_passive` in `config/analysis_config.yaml`. Indonesian and Mandarin also have one in the Ambridge, Arnon & Bekman (2023) data. Each sentence type is treated as the same construction in every language that has it, the Turkish pseudo-passive included (PI, 6 October 2026) | Settled |
| Which construction counts as the passive, language by language | Resolved 2026-10-01 with the PI: the coding of Ambridge, Arnon & Bekman (2023), with the passive as the reference level. The ka-, ma- and a- passives (Balinese), the di- passive (Indonesian) and the bei passive (Mandarin) are the passive, the basic, non-canonical and notional passives are the pseudo-passives, and the Mandarin BA active is set aside. Both interactions are predicted negative in every language. Every reported design-analysis figure uses this coding | Settled |
| Minimum number of language teams, and whether it gates Stage 2 | Resolved 2026-10-06: asked whether he still wanted a minimum number of languages before the Stage 2 paper goes ahead, the PI promised at least 20, 15 of them with a pseudo-passive and all seven measured languages among them. Every language whose team completes data collection enters the synthesis (§12) | Settled |
| Outcome-neutral quality checks for the Stage 1 protocol | Asked of the PI on 8 October 2026; his draft has attention checks and exclusions but no outcome-neutral criteria section. With H1a no longer a requirement (§3) there is no positive control; the candidates put to him are a ceiling check on the actives, attention checks and the reliability of the affectedness norms | Before Stage 1 submission |
| Verb class as an alternative explanation | Settled 2026-10-08: an exploratory Stage 2 model with verb class as a covariate, reported under Exploratory analyses, since the journal's template keeps exploratory analyses out of the Stage 1 protocol. The lists are already matched on verb class (report, Appendix B) | Settled |
| Affectedness norming scale across languages | Settled: affectedness is standardised within each language before analysis (centred and divided by twice its SD), so the estimates do not depend on each language's norming scale. The pilot normed Norwegian on 0 to 100 and English and Turkish on 1 to 7 | Settled |
| Verb term in the pooled model | Resolved 2026-10-06: each language keeps its own by-verb effects (verb IDs prefixed by language), as in the pilot and `R/11_simulate_pooled_v2.R`. The PI asked for whichever option has more power. One shared effect per verb would lower power by 0.3 to 0.75 points at the cross-language correlation of verb effects the design assumes (0.164, the observed 0.22 [0.08, 0.46] diluted for partly different verb lists and norming), and by about 4 points were verbs fully shared, and in the verb-split arms of the Glossa refit, only 0.16 to 0.18 of the active contrast's verb variance is common to all languages (`outputs/design_summary_pilot/glossa_tau_refit_summary.csv`; audit item klanguage-engine#15) | Settled |
| One agent-gender version per participant, counterbalanced across lists | The PI's draft protocol (8 October 2026) uses eight lists, each with the man and the woman as agent on half the trials and the same agent across a verb's sentence types. Balance needs each verb to appear with each agent in half the lists; confirmation of that asked of the PI on 8 October. The power evidence and the two-list design are in `reports/gender_versions_note.md` | Before Stage 1 submission |
| Affectedness rating task: one or both agent-gender versions per rater | The PI's draft says 72 ratings per rater (one version per verb) but its instructions assume 144. Both versions recommended to the PI on 8 October, because one version per rater leaves each version's score on ten raters and attenuates every slope (report, Appendix B) | Before Stage 1 submission |
| Data-collection closing date | Needed because the draft commits to no maximum number of languages, and the journal's maximum-feasible route requires a point at which collection stops regardless of the Bayes factors. Asked of the PI on 8 October | Before Stage 1 submission |
| Which entries count as languages | The draft's list of 43 omits Balinese, which the PI's email of 6 October counted among the 20, and lists varieties separately (four Spanish, two Mandarin, Afrikaans with Kaaps and Namakwalands, three Arabic). The pooled test treats languages as independent draws, so asked of the PI on 8 October whether varieties count separately towards the 20 | Before Stage 1 submission |

---

## 12. The cross-linguistic synthesis and the number of language teams

Added 2026-09-10, in response to the project lead's proposal of 7 September that the Stage 1
report carry a rule of the form "we proceed to the Stage 2 paper only if we have a minimum of
Y languages", with Y at least 6 and possibly 8 or 10.

**The design analysis behind this section.** `scripts/run_klanguage_design_analysis.R`, on the
engine in `R/13_simulate_klanguage_pooled.R`. It is a two-stage Bayesian random-effects
meta-analysis evaluated by grid integration, calibrated at $K = 3$ against the 135 completed
`L5_cross_maximal` replicates in `outputs/design_pooled_v2`, and its inputs are the seven
languages harmonised by `scripts/extract_cross_language_effects.R`.

**What the analysis establishes.** The probability that the pooled test reaches a directional
Bayes factor of 10 is governed by the between-language SD of the interaction, $\tau$, and only
secondarily by $K$. The seven available languages do not form one population: split by the
protocol that produced them, the CLAPS pilot three give $\tau \approx 0.16$ with Cochran's
$Q$ non-significant, and the Ambridge, Arnon & Bekman (2023) four give $\tau \approx 0.53$
with $Q$ significant at $10^{-7}$. The means barely differ. Conditional on each $\tau$, the pooled test
reaches 90% at six languages on the first reading and has not reached it at twenty on the
second, where it stands at 0.89.

**Decisions taken.**

1. The pooled cross-linguistic model remains the preregistered **synthesis**, reported at
   whatever $K$ is achieved. The confirmatory claim continues to rest on the per-language
   tests, which do not depend on $K$ at all (§8).

2. Any minimum-$K$ rule is written as an **interpretation band over an always-run analysis**,
   not as a publication gate. A gate would make the per-language confirmatory tests hostage to
   a recruitment outcome nobody controls, and PCI RR does not permit partial withdrawal of a
   registered plan.

3. The language set is **frozen on a named date by a named role**, before any pooled analysis
   is run and before any data are inspected. Without that, a rule conditioned on the number of
   languages is optional stopping on a between-language quantity.

4. The design analysis reports the three $\tau$ regimes side by side and never a single
   headline figure, together with the expected number of languages whose own test reaches the
   threshold in the direction **opposite** to the prediction. Balinese is already such a
   language, so this is not hypothetical. The rules agreed on 6 October give the prespecified
   reading: such a language is reported as evidence against the hypothesis in that language
   (a Bayes factor of 1/10 or less, §3) and stays in the synthesis.

5. The two-stage synthesis is preregistered as the pooled model's own **fallback rung**. The
   27 August refit of the maximal pooled model to the real pilot data saturated tree depth on
   three quarters of its iterations and so failed the §9 criteria, meaning the branch has
   already occurred once on real data.

**Amendment (2026-10-08):** item 1 is reversed by the PI's draft protocol, which registers H1 and H2 as tests across languages. The pooled analysis is now confirmatory, run as the two-stage synthesis of §1, with the one-stage model as a robustness check; the per-language estimates are reported but not tested. Item 5's fallback becomes the confirmatory analysis itself.

**Locked**: Yes, as of 2026-10-07, with item 1 amended on 8 October. Item 1 was the plan the preliminary report stated until then. Item 2
is overtaken by the PI's minimum of 20 languages, $Y = 20$ (6 October 2026). Item 3 is met by
including every language whose team completes data collection, so the set of languages
follows completion and is never chosen after the data are seen. Items 4 and 5 are implemented.
The corrected estimate of $\tau$ from `scripts/refit_glossa_pooled_tau.R` returned on
25 September at 0.47, against the published 0.45
(`outputs/design_summary_pilot/glossa_tau_refit_summary.csv`), so it does not narrow the fork.

**Update (2026-10-06).** The PI expects at least 20 languages, 15 of them with a pseudo-passive, all seven measured languages among them, and the preliminary report now plans on that count. With the seven measured languages drawn from their estimates and the rest from a population with the spread of the measured languages, the pooled test reaches $BF_{10} \geq 10$ in 89% of simulated studies for H1 (`klanguage_known7.csv`) and in 90% for H2 over the 15 languages with a pseudo-passive (`klanguage_pseudo_passive.csv`, from `scripts/run_klanguage_pseudo_passive.R`). Across the three spreads the rates run from 76% to 92% for H1 and from 58% to 95% for H2. These cells analyse each study as the one-stage model would; analysed as the registered two-stage synthesis (§1), the same simulated studies give 88% for H1 and 89% for H2, and from 74% to 92% and from 56% to 94% across the spreads (`klanguage_registered.csv`, 2026-10-08).


---

## References

Full bibliographic details for every cited source are in `references.bib`. The reference list is audited automatically by `scripts/00_verify_references.R`, which fails the build if any DOI, title or author field is missing or unverifiable against Crossref.
