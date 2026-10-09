#!/usr/bin/env Rscript
# scripts/run_klanguage_registered.R
# ---------------------------------------------------------------------------
# Detection rates for the confirmatory pooled tests as registered: a two-stage
# synthesis, proposed to the PI on 8 October 2026 (docs/preregistration_decisions.md,
# amendment to section 1).
#
#   Stage 1. Each language is analysed on its own with the planned per-language
#            model, whose priors shrink its estimate of each interaction towards
#            zero.
#   Stage 2. The per-language posterior means and SDs are combined in a Bayesian
#            random-effects meta-analysis with Language as the random effect, a
#            normal(0, P) prior on the mean across languages (so each hypothesis is
#            the half-normal on its predicted side) and a half-normal(0, 0.5) prior
#            on the between-language SD, tau. H1 and H2 are directional Bayes factors
#            on that mean.
#
# Both stages are two_stage_p_negative() in R/13_simulate_klanguage_pooled.R. The
# simulated studies are those of the reported 20-language cells, drawn by the
# engine's simulate_partly_known_estimates() with the same seeds, 7020 for H1
# (klanguage_known7.csv, K = 20) and 8015 for H2 (klanguage_pseudo_passive.csv, 15
# languages with a pseudo-passive), so that only the analysis differs. The earlier
# cells analysed each study with a flat prior on the mean and the tau prior scale
# calibrated to mimic the one-stage model (1.1). Here the analysis is the registered
# one, so no calibration is involved. Each set of measured languages keeps its
# transfer doubt as in the report: 0.40 for the CLAPS-pilot spread of H1, 0.15
# otherwise. Before anything is written, the same studies are re-analysed the earlier
# way, and the script stops unless every reported rate comes back exactly, which
# confirms that the inputs are the ones the reported cells used.
#
# Inputs: outputs/design_summary_pilot/{cross_language_effects_harmonised,
#   pilot_params_ceilings,klanguage_known7,klanguage_pseudo_passive}.csv,
#   config/prior_grid.yaml and the pilot posteriors in outputs/pilot_models
#   (gitignored; checksums in config/pilot_models_sha256.csv).
# Output: klanguage_registered.csv in --outdir, and a run log. Writing into
#   outputs/design_summary_pilot replaces the reported table, so it needs
#   --overwrite-reported, and the log then goes to outputs/logs.
# Run from design_analysis/:
#   Rscript scripts/run_klanguage_registered.R --outdir=DIR [--reps=4000] [--log=FILE]
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
aval <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
REPS <- suppressWarnings(as.integer(aval("--reps", "4000")))
if (is.na(REPS) || REPS < 1L) stop("--reps must be a positive whole number")
IN <- "outputs/design_summary_pilot"
OUT <- aval("--outdir", IN)
PILOT_MODELS <- aval("--pilot-models", "outputs/pilot_models")

source("R/13_simulate_klanguage_pooled.R")
writes_reported <- guard_reported_dir(OUT, args)
LOG <- aval("--log", file.path(if (writes_reported) "outputs/logs" else OUT, "klanguage_registered.log"))
INPUTS <- c(file.path(IN, c("cross_language_effects_harmonised.csv", "pilot_params_ceilings.csv",
                            "klanguage_known7.csv", "klanguage_pseudo_passive.csv")),
            "config/prior_grid.yaml",
            file.path(PILOT_MODELS, paste0("pilot_dgp_v2_pilot_", c("English", "Turkish", "Norwegian"), ".rds")))
log_con <- start_run_log(LOG, INPUTS, replicates = REPS, out = OUT)

H  <- read.csv(file.path(IN, "cross_language_effects_harmonised.csv"))
PP <- read.csv(file.path(IN, "pilot_params_ceilings.csv"))

# The planned number of languages, from the PI's count of 6 October 2026: at least 20,
# 15 of them with a pseudo-passive. The report checks that this table carries the same K.
K_PLAN <- 20L
K_PLAN_PP <- 15L

# The per-language model's prior SDs for the two interactions (primary regime). They
# are read from the config that R/00_reference_audit.R keeps in step with
# R/03_define_priors.R, and expressed per SD of affectedness, which the 2-SD scaling
# halves.
normal_sd <- function(spec) {
  x <- suppressWarnings(as.numeric(sub("^normal\\(0, *([0-9.]+)\\)$", "\\1", spec)))
  if (is.na(x)) stop("expected a zero-centred normal prior, got ", spec)
  x
}
PRIMARY <- yaml::read_yaml("config/prior_grid.yaml")$prior_regimes$primary$priors
P_H1 <- normal_sd(PRIMARY$b_active_int) * mean(PP$affectedness_sd)
P_H2 <- normal_sd(PRIMARY$b_pseudo_int) * mean(PP$affectedness_sd)

#' Stop unless the earlier analysis of these studies returns the reported rate.
#'
#' The earlier analysis is the one run_partly_known_cell() applies: a flat prior on the
#' mean, the calibrated tau prior and the independent standard error. Agreement shows
#' that the tables, pilot extracts and constants are the ones the reported cells used.
#'
#' @param Y Simulated per-language estimates from simulate_partly_known_estimates().
#' @param sp Standard-error split from split_se(), a list with indep and shared.
#' @param ref The one reported row these studies should reproduce.
#' @param label What the row is, for the message.
check_reported <- function(Y, sp, ref, label) {
  if (nrow(ref) != 1L) stop("[check] expected one reported row for ", label, ", found ", nrow(ref))
  if (REPS != ref$n_reps) {
    cat("check skipped for", label, ": it needs --reps =", ref$n_reps, "\n")
    return(invisible(NULL))
  }
  got <- detection_rate(ra_meta_batch(Y, sp$indep, prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE)$p_mu_negative)
  cat(sprintf("check: %s, %.4f here, %.4f reported\n", label, got$power, ref$power))
  if (!isTRUE(all.equal(got$power, ref$power))) {
    stop("[check] ", label, " no longer reproduces the reported cell")
  }
}

# --- H1 ---------------------------------------------------------------------------
# Inputs rebuilt as run_klanguage_design_analysis.R builds them, including its use of
# the main effect's posterior SD for the interaction's (audit item klanguage-engine#7).
SE_VERB <- mean(PP$se_floor_h1b * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
sp1 <- split_se(SE_VERB, sqrt(max(SE_POST^2 - SE_VERB^2, 0)), RHO_EFF_VERB)
CL <- c("English", "Turkish", "Norwegian")
SDX <- claps_pilot_sdx(PP, CL)  # the seven-digit scaling the reported cells used
raw <- read_claps_interaction_draws(PILOT_MODELS, CL)
cse <- setNames(vapply(seq_along(CL), function(i) sd(raw[[i]]) * SDX[i], numeric(1)), CL)
P1 <- build_populations(prepare_language_table(H, cse), "canonical")
TSD1 <- c(claps_protocol = 0.40, seven = 0.15, glossa_protocol = 0.15)
ref1 <- read.csv(file.path(IN, "klanguage_known7.csv"))

h1 <- do.call(rbind, lapply(names(TSD1), function(pn) {
  Y <- simulate_partly_known_estimates(K_PLAN, P1$seven, REPS, sp1$indep, pop = P1[[pn]],
                                       shared_verb_sd = sp1$shared,
                                       transfer_sd = transfer_sd_by_source(P1$seven, TSD1[[pn]]),
                                       seed = 7000L + K_PLAN)$Y
  check_reported(Y, sp1,
                 ref1[ref1$reading_variant == "canonical" & ref1$population == pn &
                        ref1$mode == "assurance" & ref1$K == K_PLAN &
                        ref1$known_set == sprintf("seven, drawn, 2023 widened %.2f", TSD1[[pn]]), ],
                 paste("H1,", pn))
  cbind(hypothesis = "H1", population = pn, K = K_PLAN, transfer_sd = TSD1[[pn]],
        prior_mu_sd = P_H1, prior_tau_scale = REGISTERED_TAU_PRIOR_SCALE,
        detection_rate(two_stage_p_negative(Y, sqrt(sp1$indep^2 + sp1$shared^2), P_H1)))
}))

# --- H2 ---------------------------------------------------------------------------
# Inputs rebuilt as run_klanguage_pseudo_passive.R builds them.
PP_LANGS <- c("English", "Turkish")
term <- "S_TypePseudo_Passive:Semantics_scaled"
pl <- lapply(PP_LANGS, function(L) {
  readRDS(file.path(PILOT_MODELS, paste0("pilot_dgp_v2_pilot_", L, ".rds")))
})
stopifnot(vapply(pl, function(p) identical(p$s_types, c("Passive", "Active", "Pseudo_Passive")),
                 logical(1)))
sdx <- PP$affectedness_sd[match(PP_LANGS, PP$language)]
sxx <- PP$sxx[match(PP_LANGS, PP$language)]
stopifnot(!anyNA(sdx), !anyNA(sxx))
d2 <- lapply(seq_along(PP_LANGS), function(i) {
  as.numeric(as.matrix(pl[[i]]$fixef_draws)[, term]) * sdx[i]
})
# Verb floor: the by-verb SD of the pseudo-passive contrast over sqrt(sxx).
fl <- vapply(seq_along(PP_LANGS), function(i) {
  sqrt(pl[[i]]$Sigma_verb["S_TypePseudo_Passive", "S_TypePseudo_Passive"] / sxx[i]) * sdx[i]
}, numeric(1))
# The participant part of the standard error shrinks from the pilot samples to the
# study's 80 participants, and the verb part does not.
N_PILOT <- c(English = 70, Turkish = 70)  # harmonised pilot samples (report, Appendix A)
N_STUDY <- 80
se_ppt2 <- sqrt(max(mean(vapply(d2, sd, numeric(1)))^2 - mean(fl)^2, 0)) *
  sqrt(mean(N_PILOT) / N_STUDY)
sp2 <- split_se(mean(fl), se_ppt2, RHO_EFF_VERB)

# The Glossa single-language pseudo-passive terms, rescaled per SD of affectedness with
# the harmonised sd_x, as run_klanguage_pseudo_passive.R does. The engine's
# ALT_PASSIVE_CONTRASTS holds them, typed from the model dumps.
GL <- c("Balinese", "Indonesian", "Mandarin")
gx <- H[H$source == "Glossa2023", ]
sx <- gx$sd_x[match(GL, gx$language)]
alt <- ALT_PASSIVE_CONTRASTS[match(GL, ALT_PASSIVE_CONTRASTS$language), ]
stopifnot(!anyNA(sx), !anyNA(alt$b_other), !anyNA(alt$se_other))
k5 <- rbind(data.frame(language = PP_LANGS, source = "CLAPS_pilot",
                       y = vapply(d2, mean, numeric(1)), s = vapply(d2, sd, numeric(1))),
            data.frame(language = GL, source = "Glossa2023",
                       y = alt$b_other * sx, s = alt$se_other * sx))
P2 <- list(claps_protocol = k5[k5$source == "CLAPS_pilot", ], five = k5,
           glossa_protocol = k5[k5$source == "Glossa2023", ])
ref2 <- read.csv(file.path(IN, "klanguage_pseudo_passive.csv"))

h2 <- do.call(rbind, lapply(names(P2), function(pn) {
  Y <- simulate_partly_known_estimates(K_PLAN_PP, k5, REPS, sp2$indep, pop = P2[[pn]],
                                       shared_verb_sd = sp2$shared,
                                       transfer_sd = transfer_sd_by_source(k5, 0.15),
                                       seed = 8000L + K_PLAN_PP)$Y
  check_reported(Y, sp2,
                 ref2[ref2$population == pn & ref2$known_set == "five, drawn, Glossa widened 0.15" &
                        ref2$mode == "assurance" & ref2$K_pp == K_PLAN_PP &
                        ref2$prior_tau_scale == CALIBRATED_TAU_PRIOR_SCALE, ],
                 paste("H2,", pn))
  cbind(hypothesis = "H2", population = pn, K = K_PLAN_PP, transfer_sd = 0.15,
        prior_mu_sd = P_H2, prior_tau_scale = REGISTERED_TAU_PRIOR_SCALE,
        detection_rate(two_stage_p_negative(Y, sqrt(sp2$indep^2 + sp2$shared^2), P_H2)))
}))

out <- rbind(h1, h2)
print(out, row.names = FALSE, digits = 3)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(OUT, "klanguage_registered.csv"), row.names = FALSE)
cat("wrote", file.path(OUT, "klanguage_registered.csv"), "\n")
stop_run_log(log_con)
