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
# Each simulated language's estimate has the SE of a normal approximation to its
# analysis. scripts/estimate_per_language_approximations.R calibrates the same
# approximation against the decision arm by widening that SE by a factor k. With
# --se-scale=calibrated the script reads k from per_language_approximations.csv and
# widens every simulated language's SE by it (or by the number given, as in
# --se-scale=1.2). The studies keep their seeds, so the result differs from the
# registered table only through the wider SEs, and it goes to a table of its own. The
# reproduction check runs on the unwidened studies either way.
#
# With --n-verbs=72,48,36 (any counts up to 72) every simulated language uses that
# many verbs. The verb part of its SE scales as sqrt(72 / n_verbs), as se_for() in
# run_klanguage_design_analysis.R scales it, and the participant part is unchanged.
# The rates go to klanguage_registered_verbs.csv, one set per count; at 72 verbs they
# are the registered table's. This option and --se-scale cannot be combined.
#
# The registered table assumes that the measured languages take part, seven for H1
# and five for H2, each with its effect drawn from its estimate. With --known=none
# none of them does, and every one of the K languages is drawn from the population
# estimated from them, as the "none known" rows of klanguage_known7.csv were drawn.
# The seeds are unchanged and the rates go to a table with "_none_known" added to its
# name. The H1 studies are checked against those reported rows; no reported H2 cell
# draws every language, so H2 has no check.
#
# The registered table is an assurance: each study draws the population's mean and
# between-language SD from their posterior. With --mode=point both are held at the
# posterior mean and median, as the reported "point" cells held them, which gives the
# power at the estimated effect. With --k=20,25,30 the studies have those numbers of
# languages, three quarters of them (rounded half up) with a pseudo-passive, each with
# its own seed, 7000 + K or 8000 + the H2 count, as the reported cells seeded theirs.
# --bf-bands adds the shares of studies whose Bayes factor ends at or below 1/10 and
# between 1/10 and 10. Each option adds a word to the table's name ("_point",
# "_languages", "_bands"), so the registered tables are never overwritten. A study
# whose cell was reported is checked against it; the others say that no check exists.
#
# Inputs: outputs/design_summary_pilot/{cross_language_effects_harmonised,
#   pilot_params_ceilings,klanguage_known7,klanguage_pseudo_passive}.csv,
#   config/prior_grid.yaml and the pilot posteriors in outputs/pilot_models
#   (gitignored; checksums in config/pilot_models_sha256.csv), and with
#   --se-scale=calibrated, per_language_approximations.csv.
# Output: klanguage_registered.csv in --outdir, or klanguage_registered_se_scaled.csv
#   when the SEs are widened, or klanguage_registered_verbs.csv for other verb counts,
#   each with the words above and "_none_known" (under --known=none) added, and a run
#   log. Writing into outputs/design_summary_pilot replaces the reported table, so it
#   needs --overwrite-reported, and the log then goes to outputs/logs.
# Run from design_analysis/:
#   Rscript scripts/run_klanguage_registered.R --outdir=DIR [--reps=4000]
#     [--se-scale=calibrated | --n-verbs=72,48,36 | --k=20,25,30] [--mode=point]
#     [--bf-bands] [--known=none] [--log=FILE]
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

# The factor applied to every simulated language's SE: 1 for the registered table.
SE_SCALE_ARG <- aval("--se-scale", "1")
APPROX_TABLE <- file.path(IN, "per_language_approximations.csv")
SE_SCALE <- if (identical(SE_SCALE_ARG, "calibrated")) {
  unique(read.csv(APPROX_TABLE)$k)
} else {
  suppressWarnings(as.numeric(SE_SCALE_ARG))
}
if (length(SE_SCALE) != 1L || !is.finite(SE_SCALE) || SE_SCALE <= 0) {
  stop("--se-scale must be 'calibrated', with one k in ", APPROX_TABLE, ", or a positive number")
}
# The number of verbs every simulated language uses: 72 for the registered table.
N_VERBS <- suppressWarnings(as.integer(strsplit(aval("--n-verbs", "72"), ",", fixed = TRUE)[[1]]))
if (!length(N_VERBS) || anyNA(N_VERBS) || any(N_VERBS < 1L | N_VERBS > 72L) || anyDuplicated(N_VERBS)) {
  stop("--n-verbs must list distinct whole numbers from 1 to 72")
}
VERBS_RUN <- !identical(N_VERBS, 72L)
if (VERBS_RUN && SE_SCALE != 1) stop("--n-verbs and --se-scale cannot be combined")
# The numbers of languages for H1: 20 for the registered table.
K_LIST <- suppressWarnings(as.integer(strsplit(aval("--k", "20"), ",", fixed = TRUE)[[1]]))
if (!length(K_LIST) || anyNA(K_LIST) || any(K_LIST < 2L) || anyDuplicated(K_LIST)) {
  stop("--k must list distinct whole numbers of at least 2")
}
K_RUN <- !identical(K_LIST, 20L)
if (K_RUN && (VERBS_RUN || SE_SCALE != 1)) stop("--k cannot be combined with --n-verbs or --se-scale")
# How each study's population mean and SD are set: drawn ("assurance") for the registered table.
MODE <- aval("--mode", "assurance")
if (!MODE %in% c("assurance", "point")) stop("--mode must be 'assurance' or 'point'")
BANDS <- "--bf-bands" %in% args
OUT_NAME <- if (VERBS_RUN) {
  "klanguage_registered_verbs"
} else if (SE_SCALE == 1) {
  "klanguage_registered"
} else {
  "klanguage_registered_se_scaled"
}
if (K_RUN) OUT_NAME <- paste0(OUT_NAME, "_languages")
if (MODE == "point") OUT_NAME <- paste0(OUT_NAME, "_point")
if (BANDS) OUT_NAME <- paste0(OUT_NAME, "_bands")
# Whether the measured languages take part: "measured" for the registered table.
KNOWN <- aval("--known", "measured")
if (!KNOWN %in% c("measured", "none")) stop("--known must be 'measured' or 'none'")
NONE_KNOWN <- KNOWN == "none"
if (NONE_KNOWN) OUT_NAME <- paste0(OUT_NAME, "_none_known")

source("R/13_simulate_klanguage_pooled.R")
writes_reported <- guard_reported_dir(OUT, args)
LOG <- aval("--log", file.path(if (writes_reported) "outputs/logs" else OUT, paste0(OUT_NAME, ".log")))
INPUTS <- c(file.path(IN, c("cross_language_effects_harmonised.csv", "pilot_params_ceilings.csv",
                            "klanguage_known7.csv", "klanguage_pseudo_passive.csv")),
            "config/prior_grid.yaml",
            file.path(PILOT_MODELS, paste0("pilot_dgp_v2_pilot_", c("English", "Turkish", "Norwegian"), ".rds")),
            if (identical(SE_SCALE_ARG, "calibrated")) APPROX_TABLE)
log_con <- start_run_log(LOG, INPUTS, replicates = REPS, se_scale = SE_SCALE,
                         n_verbs = paste(N_VERBS, collapse = ","), known = KNOWN,
                         k = paste(K_LIST, collapse = ","), mode = MODE, bf_bands = BANDS, out = OUT)

H  <- read.csv(file.path(IN, "cross_language_effects_harmonised.csv"))
PP <- read.csv(file.path(IN, "pilot_params_ceilings.csv"))

# The planned number of languages, from the PI's count of 6 October 2026: at least 20,
# 15 of them with a pseudo-passive. The report checks that this table carries the same K.
K_PLAN <- 20L
K_PLAN_PP <- 15L
# The H2 count for any number of languages: three quarters with a pseudo-passive, the
# protocol's assumption, rounded half up.
k_pp <- function(K) as.integer(floor(0.75 * K + 0.5))
stopifnot(k_pp(K_PLAN) == K_PLAN_PP)

# The per-language model's prior SDs for the two interactions (primary regime). They
# are read from config/prior_grid.yaml, which a unit test in tests/testthat/test-priors.R
# checks against R/03_define_priors.R, and expressed per SD of affectedness, which the
# 2-SD scaling halves.
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

#' check_reported() where a reported cell exists, and otherwise a note that none does.
#'
#' The registered cell (20 languages for H1, 15 for H2, assurance) must always be found,
#' so check_reported() is called for it even when the row is missing, and stops.
#'
#' @param registered Whether this is the registered cell.
#' @inheritParams check_reported
check_if_reported <- function(Y, sp, ref, label, registered) {
  if (!registered && nrow(ref) == 0L) {
    cat("check skipped for", label, ": no reported cell\n")
    return(invisible(NULL))
  }
  check_reported(Y, sp, ref, label)
}

#' The point values of a population's mean and SD, as the reported "point" cells set them.
#'
#' @param pop Anchor set with columns y and s.
#' @return list(mu, tau): the posterior mean of the mean and the posterior median of tau.
point_values <- function(pop) {
  post <- ra_meta_posterior(pop$y, pop$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
  list(mu = post$mu_mean, tau = post$tau_median)
}

#' A standard-error split from split_se() with both parts multiplied by one factor.
scale_se <- function(sp, k) list(shared = sp$shared * k, indep = sp$indep * k)

#' The SE splits to analyse: one per verb count, or the one SE scale.
#'
#' @param se_verb,se_ppt The verb and participant parts of a language's SE at 72 verbs.
#' @return A list of list(n_verbs, sp); n_verbs is NA unless this is a verb-count run.
se_variants <- function(se_verb, se_ppt) {
  if (!VERBS_RUN) {
    return(list(list(n_verbs = NA_integer_, sp = scale_se(split_se(se_verb, se_ppt, RHO_EFF_VERB), SE_SCALE))))
  }
  lapply(N_VERBS, function(n) {
    list(n_verbs = n, sp = split_se(se_verb * sqrt(72 / n), se_ppt, RHO_EFF_VERB))
  })
}

# --- H1 ---------------------------------------------------------------------------
# Inputs rebuilt as run_klanguage_design_analysis.R builds them, including its use of
# the main effect's posterior SD for the interaction's, kept by decision on 9 October
# 2026 (audit item klanguage-engine#7).
SE_VERB <- mean(PP$se_floor_h1b * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
SE_PPT1 <- sqrt(max(SE_POST^2 - SE_VERB^2, 0))
sp1 <- split_se(SE_VERB, SE_PPT1, RHO_EFF_VERB)
CL <- c("English", "Turkish", "Norwegian")
SDX <- claps_pilot_sdx(PP, CL)  # the seven-digit scaling the reported cells used
raw <- read_claps_interaction_draws(PILOT_MODELS, CL)
cse <- setNames(vapply(seq_along(CL), function(i) sd(raw[[i]]) * SDX[i], numeric(1)), CL)
P1 <- build_populations(prepare_language_table(H, cse), "canonical")
TSD1 <- c(claps_protocol = 0.40, seven = 0.15, glossa_protocol = 0.15)
ref1 <- read.csv(file.path(IN, "klanguage_known7.csv"))

#' Rates for one set of simulated studies under each SE variant.
#'
#' The studies drawn at the registered SEs are reused wherever a variant leaves the SEs
#' unchanged, so the registered rows are reproduced exactly.
#'
#' @param Y0 Studies drawn at the registered SEs.
#' @param sp0 The registered SE split.
#' @param variants Output of se_variants().
#' @param draw Function drawing the same studies at another SE split.
#' @param lead One-row data frame of identifying columns for the rows.
#' @param prior_sd Prior SD on the mean across languages.
variant_rows <- function(Y0, sp0, variants, draw, lead, prior_sd) {
  do.call(rbind, lapply(variants, function(v) {
    Y <- if (identical(v$sp, sp0)) Y0 else draw(v$sp)
    p <- two_stage_p_negative(Y, sqrt(v$sp$indep^2 + v$sp$shared^2), prior_sd)
    row <- cbind(lead, detection_rate(p))
    if (VERBS_RUN) row$n_verbs <- v$n_verbs
    if (BANDS) {
      bf <- bf_from_p_negative(p)
      row$p_bf_low <- mean(bf <= 1 / 10)
      row$p_inconclusive <- mean(bf > 1 / 10 & bf < 10)
    }
    row
  }))
}

variants1 <- se_variants(SE_VERB, SE_PPT1)
known1 <- if (NONE_KNOWN) P1$seven[0, ] else P1$seven

#' H1 rows for one number of languages and one population.
#'
#' @param K Number of languages.
#' @param pn Name of the population in P1.
h1_cell <- function(K, pn) {
  tsd <- if (NONE_KNOWN) 0 else TSD1[[pn]]
  pv <- point_values(P1[[pn]])
  draw <- function(sp) {
    simulate_partly_known_estimates(K, known1, REPS, sp$indep, pop = P1[[pn]],
                                    shared_verb_sd = sp$shared,
                                    transfer_sd = transfer_sd_by_source(known1, tsd),
                                    mode = MODE, mu_point = pv$mu, tau_point = pv$tau,
                                    seed = 7000L + K)$Y
  }
  Y0 <- draw(sp1)
  ks <- if (NONE_KNOWN) "none known" else sprintf("seven, drawn, 2023 widened %.2f", tsd)
  check_if_reported(Y0, sp1,
                    ref1[ref1$reading_variant == "canonical" & ref1$population == pn &
                           ref1$mode == MODE & ref1$K == K & ref1$known_set == ks, ],
                    sprintf("H1, %s, K = %d", pn, K), registered = K == K_PLAN && MODE == "assurance")
  variant_rows(Y0, sp1, variants1, draw,
               data.frame(hypothesis = "H1", population = pn, K = K, transfer_sd = tsd,
                          prior_mu_sd = P_H1, prior_tau_scale = REGISTERED_TAU_PRIOR_SCALE),
               P_H1)
}
h1 <- do.call(rbind, lapply(K_LIST, function(K) do.call(rbind, lapply(names(TSD1), h1_cell, K = K))))

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

variants2 <- se_variants(mean(fl), se_ppt2)
known2 <- if (NONE_KNOWN) k5[0, ] else k5
tsd2 <- if (NONE_KNOWN) 0 else 0.15

#' H2 rows for one number of languages with a pseudo-passive and one population.
#'
#' @param kp Number of languages with a pseudo-passive.
#' @param pn Name of the population in P2.
h2_cell <- function(kp, pn) {
  pv <- point_values(P2[[pn]])
  draw <- function(sp) {
    simulate_partly_known_estimates(kp, known2, REPS, sp$indep, pop = P2[[pn]],
                                    shared_verb_sd = sp$shared,
                                    transfer_sd = transfer_sd_by_source(known2, tsd2),
                                    mode = MODE, mu_point = pv$mu, tau_point = pv$tau,
                                    seed = 8000L + kp)$Y
  }
  Y0 <- draw(sp2)
  label <- sprintf("H2, %s, K = %d", pn, kp)
  if (NONE_KNOWN) {
    cat("check skipped for", label, ": no reported cell draws every language\n")
  } else {
    check_if_reported(Y0, sp2,
                      ref2[ref2$population == pn & ref2$known_set == "five, drawn, Glossa widened 0.15" &
                             ref2$mode == MODE & ref2$K_pp == kp &
                             ref2$prior_tau_scale == CALIBRATED_TAU_PRIOR_SCALE, ],
                      label, registered = kp == K_PLAN_PP && MODE == "assurance")
  }
  variant_rows(Y0, sp2, variants2, draw,
               data.frame(hypothesis = "H2", population = pn, K = kp, transfer_sd = tsd2,
                          prior_mu_sd = P_H2, prior_tau_scale = REGISTERED_TAU_PRIOR_SCALE),
               P_H2)
}
h2 <- do.call(rbind, lapply(k_pp(K_LIST), function(kp) do.call(rbind, lapply(names(P2), h2_cell, kp = kp))))

out <- rbind(h1, h2)
# The registered table keeps its columns; a widened run records its factor.
if (SE_SCALE != 1) out$se_scale <- SE_SCALE
if (MODE == "point") out$mode <- MODE
if (NONE_KNOWN) out$known <- KNOWN
print(out, row.names = FALSE, digits = 3)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT, paste0(OUT_NAME, ".csv"))
write.csv(out, out_file, row.names = FALSE)
cat("wrote", out_file, "\n")
stop_run_log(log_con)
