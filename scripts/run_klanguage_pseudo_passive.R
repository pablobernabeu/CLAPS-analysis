#!/usr/bin/env Rscript
# scripts/run_klanguage_pseudo_passive.R
# ---------------------------------------------------------------------------
# The pooled test of the pseudo-passive prediction (H2b) at 10 and 20 languages,
# built the way the PI asked for the main interaction on 23 September: the
# languages that already have data enter at their estimates, and the remaining
# places are simulated. It answers his note of 1 October, which confirmed the
# coding of the Glossa analysis (Ambridge, Arnon & Bekman, 2023): Passive is the
# reference level, and the non-canonical (Indonesian), basic (Balinese) and
# notional (Mandarin) passives are the pseudo-passives.
#
# Everything statistical is R/13_simulate_klanguage_pooled.R, unchanged. Only the
# inputs differ from scripts/run_klanguage_design_analysis.R, section 5d:
#
#   Known languages. Five have a pseudo-passive and data: English and Turkish
#     from the CLAPS pilots, Balinese, Indonesian and Mandarin from the Glossa
#     single-language models. Norwegian and Hebrew have no pseudo-passive.
#     Per SD of affectedness: English -0.383, Turkish +0.158, Balinese -0.534,
#     Indonesian -0.342, Mandarin -1.100.
#
#   Per-language SE. Built as the main interaction's is, from the verb floor and
#     the pilot posterior SD, but with the pseudo-passive term's own values and
#     rescaled to 80 participants (the main-interaction run uses the main effect's
#     SD at the pilot samples, audit item klanguage-engine#7). The pseudo-passive
#     term is about a quarter more precise than the main interaction, which the
#     pilot posteriors, the verb floors and the decision-arm replicates all agree on.
#
#   How many of the K languages have a pseudo-passive. Five of the seven measured
#     languages do. Two readings are reported: every new language has one, or new
#     languages have one at the measured rate of five in seven.
#
#   Spread. Three populations, as for the main interaction: all five languages,
#     the three Glossa languages, and the two CLAPS pilots.
#
# The analysis prior on tau stays at the scale calibrated for the main
# interaction (1.1), with the bracket 1.0 to 1.2 as a sensitivity; no pooled
# brms run of H2b exists to recalibrate it against, beyond the three-language
# runs in which English and Turkish cancel (0 of 135 reach a Bayes factor of 10).
#
# Inputs: outputs/design_summary_pilot/{cross_language_effects_harmonised,
#   pilot_params_ceilings,klanguage_known7}.csv and the pilot posteriors in
#   outputs/pilot_models (gitignored).
# Output: outputs/design_summary_pilot/klanguage_pseudo_passive.csv
#   and klanguage_pseudo_passive_populations.csv.
# Run from design_analysis/:  Rscript scripts/run_klanguage_pseudo_passive.R [--reps=4000]
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
aval <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
REPS <- as.integer(aval("--reps", 4000L))
IN <- "outputs/design_summary_pilot"
OUT <- aval("--outdir", IN)
PILOT_MODELS <- aval("--pilot-models", "outputs/pilot_models")

source("R/13_simulate_klanguage_pooled.R")

H  <- read.csv(file.path(IN, "cross_language_effects_harmonised.csv"))
PP <- read.csv(file.path(IN, "pilot_params_ceilings.csv"))
RHO_EFF <- 0.164  # as in run_klanguage_design_analysis.R

# ---------------------------------------------------------------------------
# 1. Check: this script calls the engine exactly as the main run does
# ---------------------------------------------------------------------------
# Rebuild the main-interaction inputs as run_klanguage_design_analysis.R does and
# reproduce one reported cell of klanguage_known7.csv to the replicate.
SE_VERB <- mean(PP$se_floor_h1b  * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
SE_PPT  <- sqrt(max(SE_POST^2 - SE_VERB^2, 0))
sp_h1b  <- split_se(SE_VERB, SE_PPT, RHO_EFF)
CLAPS_LANGS <- c("English", "Turkish", "Norwegian")
# Written out as the main run writes them, so that the check reproduces it exactly.
CLAPS_SDX   <- c(0.5035005, 0.5035005, 0.5038147)
stopifnot(isTRUE(all.equal(CLAPS_SDX, PP$affectedness_sd[match(CLAPS_LANGS, PP$language)], tolerance = 1e-6)))
claps_h1b   <- read_claps_interaction_draws(PILOT_MODELS, CLAPS_LANGS)
claps_se    <- setNames(vapply(seq_along(CLAPS_LANGS), function(i) sd(claps_h1b[[i]]) * CLAPS_SDX[i],
                               numeric(1)), CLAPS_LANGS)
P_h1b <- build_populations(prepare_language_table(H, claps_se), "canonical")
ref <- read.csv(file.path(IN, "klanguage_known7.csv"))
ref <- ref[ref$reading_variant == "canonical" & ref$population == "seven" & ref$mode == "assurance" &
             ref$known_set == "seven, drawn, 2023 widened 0.15" & ref$K == 10, ]
if (REPS == ref$n_reps) {
  chk <- run_partly_known_cell(K = 10, known = P_h1b$seven, n_reps = REPS, se_lang = sp_h1b$indep,
                               pop = P_h1b$seven, shared_verb_sd = sp_h1b$shared,
                               transfer_sd = transfer_sd_by_source(P_h1b$seven, 0.15),
                               mode = "assurance", known_truth = "drawn", seed = 7010L)
  cat(sprintf("check: main interaction, seven, K = 10: %.4f here, %.4f reported\n", chk$power, ref$power))
  stopifnot(isTRUE(all.equal(chk$power, ref$power)))
} else {
  cat("check skipped: it needs --reps =", ref$n_reps, "\n")
}

# ---------------------------------------------------------------------------
# 2. Inputs for the pseudo-passive term
# ---------------------------------------------------------------------------
PP_LANGS <- c("English", "Turkish")
term <- "S_TypePseudo_Passive:Semantics_scaled"
pil <- lapply(PP_LANGS, function(L) readRDS(file.path(PILOT_MODELS, paste0("pilot_dgp_v2_pilot_", L, ".rds"))))
names(pil) <- PP_LANGS
sdx <- PP$affectedness_sd[match(PP_LANGS, PP$language)]
sxx <- PP$sxx[match(PP_LANGS, PP$language)]
pp_draws <- lapply(seq_along(PP_LANGS), function(i) as.numeric(as.matrix(pil[[i]]$fixef_draws)[, term]) * sdx[i])
names(pp_draws) <- PP_LANGS
stopifnot(vapply(pil, function(p) identical(p$s_types, c("Passive", "Active", "Pseudo_Passive")), logical(1)))

# Verb floor for the term: the by-verb SD of the pseudo-passive contrast over
# sqrt(sxx), as se_floor_h1b is built from the Active contrast's (Sigma_verb[2, 2]).
floor_pp <- vapply(seq_along(PP_LANGS), function(i) sqrt(pil[[i]]$Sigma_verb[3, 3] / sxx[i]) * sdx[i], numeric(1))
post_pp  <- vapply(pp_draws, sd, numeric(1))
N_PILOT <- c(English = 70, Turkish = 70)
SE_VERB2 <- mean(floor_pp)
SE_POST2 <- mean(post_pp)
# The participant part shrinks with N; the verb part does not.
SE_PPT2  <- sqrt(max(SE_POST2^2 - SE_VERB2^2, 0)) * sqrt(mean(N_PILOT) / 80)
SE_LANG2 <- sqrt(SE_VERB2^2 + SE_PPT2^2)
sp <- split_se(SE_VERB2, SE_PPT2, RHO_EFF)

# Glossa single-language models (OSF dumps), Passive as reference: the
# pseudo-passive-by-Semantics term, rescaled per SD with the harmonised sd_x.
glossa_raw <- data.frame(language = c("Balinese", "Indonesian", "Mandarin"),
                         construction = c("Pseudo_Passive (Passive_basic)", "Non_Canonical_Passive",
                                          "Notional_Passive"),
                         file = c("Balinese_Sem_Only_Bayes_redux.txt", "Indonesian_Sem_Only_Bayes.txt",
                                  "Mandarin_Sem_Only_Bayes.txt"),
                         raw_est = c(-0.53, -0.34, -1.09), raw_se = c(0.13, 0.10, 0.18))
# These are ALT_PASSIVE_CONTRASTS' b_other, typed from the same dumps.
stopifnot(isTRUE(all.equal(glossa_raw$raw_est, ALT_PASSIVE_CONTRASTS$b_other[match(glossa_raw$language,
                                                                               ALT_PASSIVE_CONTRASTS$language)])))
gx <- H[H$source == "Glossa2023", ]
glossa_raw$sd_x <- gx$sd_x[match(glossa_raw$language, gx$language)]

known5 <- rbind(
  data.frame(language = PP_LANGS, source = "CLAPS_pilot",
             y = vapply(pp_draws, mean, numeric(1)), s = post_pp),
  data.frame(language = glossa_raw$language, source = "Glossa2023",
             y = glossa_raw$raw_est * glossa_raw$sd_x, s = glossa_raw$raw_se * glossa_raw$sd_x)
)
rownames(known5) <- NULL
POPS2 <- list(five            = known5,
              glossa_protocol = known5[known5$source == "Glossa2023", ],
              claps_protocol  = known5[known5$source == "CLAPS_pilot", ])

cat("\nPseudo-passive interaction per SD of affectedness, the five languages that have one:\n")
print(known5, row.names = FALSE, digits = 3)
cat(sprintf(paste("\nper-language SE at 80 participants and 72 verbs: total %.4f =",
                  "verb %.4f + participant %.4f; split %.4f independent, %.4f shared\n"),
            SE_LANG2, SE_VERB2, SE_PPT2, sp$indep, sp$shared))
cat(sprintf("(main interaction, as run: total %.4f, split %.4f / %.4f)\n",
            sqrt(SE_VERB^2 + SE_PPT^2), sp_h1b$indep, sp_h1b$shared))

pop_tab <- do.call(rbind, lapply(names(POPS2), function(nm) population_summary(POPS2[[nm]], nm)))
cat("\nPopulations (design prior on tau, half-normal 0.5):\n")
print(pop_tab[, c("population", "K", "mean_effect", "mu_mean", "mu_sd", "tau_median", "tau_lo", "tau_hi",
                  "p_mu_negative")], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 3. Engine check against the three-language brms runs
# ---------------------------------------------------------------------------
# In those runs only English and Turkish carry a pseudo-passive, each drawn from
# its pilot posterior, at the verb-floor SEs, with separate verb lists; 0 of 135
# pooled fits reached a Bayes factor of 10 for H2b.
set.seed(20261001L)
nchk <- 4000L
Ychk <- cbind(sample(pp_draws$English, nchk, TRUE) + rnorm(nchk, 0, floor_pp[1]),
              sample(pp_draws$Turkish, nchk, TRUE) + rnorm(nchk, 0, floor_pp[2]))
p_chk <- ra_meta_batch(Ychk, floor_pp, prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE)$p_mu_negative
cat(sprintf("\nengine check, English and Turkish pooled: P(BF >= 10) = %.4f (brms: 0/135)\n",
            mean(bf_from_p_negative(p_chk) >= 10)))

# ---------------------------------------------------------------------------
# 4. The cells
# ---------------------------------------------------------------------------
# K_pp is the number of languages with a pseudo-passive, the only ones H2b is
# tested in. It maps onto the K languages of the main-interaction figures in two ways.
K_TOTAL <- c(7, 8, 10, 12, 16, 20)
kmap <- rbind(data.frame(K = K_TOTAL, new_rule = "every new language has one", K_pp = K_TOTAL - 2L),
              data.frame(K = K_TOTAL, new_rule = "new languages at five in seven",
                         K_pp = 5L + as.integer(round((K_TOTAL - 7) * 5 / 7))))
KPP <- sort(unique(c(kmap$K_pp, 2L)))

known_sets <- list(
  list(known_set = "five, drawn, Glossa widened 0.15", known = known5, tsd = 0.15),
  list(known_set = "five, drawn, Glossa widened 0.40", known = known5, tsd = 0.40),
  list(known_set = "English and Turkish only",         known = POPS2$claps_protocol, tsd = 0)
)
scales <- c(CALIBRATED_TAU_PRIOR_SCALE, CALIBRATION_BRACKET)

res <- do.call(rbind, lapply(names(POPS2), function(pn) {
  post <- ra_meta_posterior(POPS2[[pn]]$y, POPS2[[pn]]$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
  do.call(rbind, lapply(c("assurance", "point"), function(md) {
    do.call(rbind, lapply(known_sets, function(v) {
      do.call(rbind, lapply(KPP[KPP >= nrow(v$known)], function(k) {
        do.call(rbind, lapply(scales, function(sc) {
          # The bracket is only needed for the headline rows.
          if (sc != CALIBRATED_TAU_PRIOR_SCALE && (md != "assurance" || v$tsd == 0.40)) return(NULL)
          cbind(population = pn, known_set = v$known_set, tsd_glossa = v$tsd, prior_tau_scale = sc,
                K_pp = k,
                run_partly_known_cell(K = k, known = v$known, n_reps = REPS, se_lang = sp$indep,
                                      pop = POPS2[[pn]], shared_verb_sd = sp$shared,
                                      transfer_sd = transfer_sd_by_source(v$known, v$tsd),
                                      mode = md, mu_point = post$mu_mean, tau_point = post$tau_median,
                                      known_truth = "drawn", prior_tau_scale = sc,
                                      seed = 8000L + k))
        }))
      }))
    }))
  }))
}))
# The number of draws must not depend on the population when nothing is drawn.
k5 <- res[res$K_pp == 5 & res$known_set != "English and Turkish only" & res$prior_tau_scale == 1.1, ]
stopifnot(all(tapply(k5$power, k5$known_set, function(x) diff(range(x))) == 0))

res$se_lang_total <- SE_LANG2
write.csv(res, file.path(OUT, "klanguage_pseudo_passive.csv"), row.names = FALSE)
write.csv(cbind(pop_tab, se_lang_total = SE_LANG2, se_indep = sp$indep, se_shared = sp$shared),
          file.path(OUT, "klanguage_pseudo_passive_populations.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Summary on the K of the main-interaction figures
# ---------------------------------------------------------------------------
head_rows <- res[res$mode == "assurance" & res$prior_tau_scale == 1.1, ]
for (ks in unique(head_rows$known_set)) {
  cat("\n", ks, ", assurance, P(BF >= 10) by languages with a pseudo-passive:\n", sep = "")
  d <- head_rows[head_rows$known_set == ks, ]
  print(round(tapply(d$power, list(d$population, d$K_pp), function(x) x[1]), 3))
}
cat("\nOn the K of the main-interaction figures (five known, Glossa widened 0.15, assurance):\n")
h <- head_rows[head_rows$known_set == "five, drawn, Glossa widened 0.15", ]
for (r in seq_len(nrow(kmap))) {
  x <- h[h$K_pp == kmap$K_pp[r], ]
  cat(sprintf("  K = %2d, %-31s K_pp = %2d: five %.3f | Glossa %.3f | CLAPS %.3f\n", kmap$K[r], kmap$new_rule[r],
              kmap$K_pp[r], x$power[x$population == "five"], x$power[x$population == "glossa_protocol"],
              x$power[x$population == "claps_protocol"]))
}
br <- res[res$mode == "assurance" & res$known_set == "five, drawn, Glossa widened 0.15" &
            res$population == "five", ]
cat("\nCalibration bracket, five-language spread:\n")
print(round(tapply(br$power, list(br$prior_tau_scale, br$K_pp), function(x) x[1]), 3))
cat("\nMonte-Carlo error at most", sprintf("%.3f", max(res$mcse)), "\n")
