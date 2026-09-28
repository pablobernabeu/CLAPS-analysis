#!/usr/bin/env Rscript
# scripts/run_klanguage_design_analysis.R
# ---------------------------------------------------------------------------
# The K-language pooled design analysis the PI asked for on 7 September 2026:
# how often does the pooled cross-linguistic test clear a directional Bayes
# factor of 10 on the affectedness-by-sentence-type interaction when K language
# teams each contribute 80 participants and 72 verbs?
#
# THE ANSWER IS NOT A CURVE IN K. IT IS A FORK IN tau.
#   The seven available languages do not form one population. Split them by the
#   protocol that produced them and the heterogeneity is almost entirely inside
#   the Ambridge, Arnon & Bekman (2023) set:
#
#     CLAPS pilot (English, Turkish, Norwegian)      tau ~ 0.16   Q p = 0.13
#     Glossa 2023 (Balinese, Hebrew, Indonesian,
#                  Mandarin)                          tau ~ 0.53   Q p = 1.4e-07
#     all seven pooled                                tau ~ 0.39
#
#   The MEANS barely differ (-0.29 against -0.34; posterior means -0.27 and -0.31).
#   The SPREADS differ by a factor of about three and a half. Stage 2 will run the
#   CLAPS protocol in every language, so which of these is the right tau is a
#   question about whether the Glossa spread is a property of the languages or of
#   the measurement. With tau known, 90% is reached at K = 6 under the CLAPS
#   reading and at K = 16 over all seven, and not by K = 20 under the Glossa
#   reading (klanguage_conditional.csv).
#
#   This script therefore reports THREE regimes side by side and never a single
#   headline. A single number here would be the most misleading thing the project
#   could put in a Stage 1 protocol.
#
# WHAT ELSE IT PRODUCES, AND WHY
#   A design analysis that stops at "power at K" leaves the PI without the things
#   a Registered Report actually has to say. So the script also computes:
#
#     - the expected number of languages whose OWN test reaches the threshold in
#       the direction OPPOSITE to the prediction. Balinese already is one. At
#       K = 12 the expectation is close to two, and the protocol needs a
#       preregistered reading of that before it happens rather than after.
#     - the probability that the pooled test fails while most individual languages
#       succeed, which is the awkward cell nobody has costed.
#     - the trade the PI will propose the moment he sees a K curve: participants
#       against languages at a fixed budget, and verbs against languages. The
#       second is the useful one, because it turns out a language that can only
#       field 47 verbs costs the pooled test almost nothing.
#     - the recruitment model. K is a random variable and Y is a rule about it.
#     - the PI's construction of 23 September (section 5d): the seven languages
#       already measured all take part at their estimated effects, and only the
#       other K - 7 places are drawn, for K = 7 to 20.
#     - that construction again under the other reading of which construction is
#       the passive in Balinese, Indonesian and Mandarin, a ruling the project has
#       not made and one that moves the seven-language mean from -0.32 to -0.03.
#
# INPUTS, read from --inputs and never from --out
#   cross_language_effects_harmonised.csv   written by scripts/extract_cross_language_effects.R
#   pilot_params_ceilings.csv               written by scripts/aggregate_pilot.R
#   and, because they are gitignored and so absent from a fresh clone or the public
#   mirror, two directories the run stops without:
#   outputs/pilot_models/pilot_dgp_v2_pilot_{English,Turkish,Norwegian}.rds
#                                           the pilot posteriors (--pilot-models)
#   outputs/design_pooled_v2/*.rds          the three-language brms replicates, or
#                                           the table a previous run wrote from them
#                                           (--pooled-replicates; section calib only)
#
# OUTPUTS, in --out
#   klanguage_populations.csv                what each language set says about (mu, tau)
#   klanguage_populations_by_reading.csv     the same under both readings of the passive
#   klanguage_calibration.csv                the K = 3 check against the brms runs
#   klanguage_calibration_target_replicates.csv  the replicates the target comes from
#   klanguage_power.csv                      the three tau regimes, assurance
#   klanguage_conditional.csv                the same with tau known
#   klanguage_sensitivity.csv                assumption sweeps
#   klanguage_minimum_k.csv                  smallest K per target
#   klanguage_named_subsets.csv              every subset of the seven, no synthetic language
#   klanguage_set_variance.csv               where the outcome variance comes from
#   klanguage_partly_known.csv               sets with some languages already measured
#   klanguage_known7.csv                     the PI's construction, both readings
#   klanguage_outcomes.csv                   wrong-sign languages, contradiction
#   klanguage_tradeoffs.csv                  N-vs-K, verbs-vs-K
#   klanguage_recruitment.csv                P(K < Y) by completion rate
#   klanguage_design_analysis.log            this run's output with its provenance; it
#                                            goes to outputs/logs/ instead when --out is
#                                            the reported directory, which the public
#                                            mirror copies wholesale
#
# COST
#   The full run is dominated by the one-replicate-at-a-time cells of sections 3, 4,
#   6 and 7 and takes one and a half to three and a half hours on one core. The
#   sections that use the batch posterior are cheaper: calib and partly about a
#   minute each, known7 (432 cells) about eleven minutes, and varshare (a million
#   simulated studies per K at the default size) about half an hour, measured on a
#   four-core laptop at 4,000 replicates a cell. The expensive part was already paid
#   for: the 135 three-language brms replicates that calibrate the engine.
#
# USAGE, from design_analysis/
#   Rscript scripts/run_klanguage_design_analysis.R --out=<dir> [options]
#     --out=<dir>                 where to write. Writing into
#                                 outputs/design_summary_pilot replaces the reported
#                                 CSVs, so it also needs --overwrite-reported.
#     --inputs=<dir>              default outputs/design_summary_pilot
#     --pilot-models=<dir>        default outputs/pilot_models
#     --pooled-replicates=<dir|csv>  default outputs/design_pooled_v2
#     --sections=<list>           comma-separated, default all of: calib, power, cond,
#                                 sens, mink, named, varshare, partly, known7,
#                                 outcomes, trades, recruit (mink and recruit need power)
#     --passive=canonical|alternative  reading of the passive for every section but
#                                 known7, which runs both; default canonical
#     --alt-se=independent|canonical   SE rule for the alternative reading
#     --reps=4000   --quick (800 replicates)   --decomp-par=<n>   --log=<file>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(stats) })

# Every path below is relative to design_analysis/, so stop at once from anywhere else
# instead of failing half way with a missing-file error that names the wrong thing.
if (!file.exists("R/13_simulate_klanguage_pooled.R")) {
  stop("run this script from design_analysis/ (R/13_simulate_klanguage_pooled.R not found from ",
       getwd(), ")")
}
source("R/13_simulate_klanguage_pooled.R")

args  <- commandArgs(trailingOnly = TRUE)
aval  <- function(f, d) { h <- grep(paste0("^", f, "="), args, value = TRUE)
                          if (!length(h)) d else sub(paste0("^", f, "="), "", h[1]) }
REPS  <- as.integer(aval("--reps", if ("--quick" %in% args) 800L else 4000L))
REPORTED_DIR  <- "outputs/design_summary_pilot"
IN            <- aval("--inputs", REPORTED_DIR)
OUT           <- aval("--out", REPORTED_DIR)
PILOT_MODELS  <- aval("--pilot-models", "outputs/pilot_models")
POOLED_K3     <- aval("--pooled-replicates", "outputs/design_pooled_v2")
READING       <- aval("--passive", "canonical")
ALT_SE        <- aval("--alt-se", "independent")

# The reported CSVs are quoted in the note and the reply to the PI, so replacing them
# is a decision, not a side effect of running the script with its defaults.
same_dir <- function(a, b) {
  normalizePath(a, winslash = "/", mustWork = FALSE) == normalizePath(b, winslash = "/", mustWork = FALSE)
}
WRITES_REPORTED <- same_dir(OUT, REPORTED_DIR)
if (WRITES_REPORTED && !("--overwrite-reported" %in% args)) {
  stop("--out is ", REPORTED_DIR, ", which holds the reported CSVs. Pass --overwrite-reported ",
       "to replace them deliberately, or --out=<dir> for a run elsewhere.")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ALL_SECTIONS <- c("calib", "power", "cond", "sens", "mink", "named", "varshare",
                  "partly", "known7", "outcomes", "trades", "recruit")
SECTIONS <- strsplit(aval("--sections", paste(ALL_SECTIONS, collapse = ",")), ",", fixed = TRUE)[[1]]
if (length(setdiff(SECTIONS, ALL_SECTIONS))) {
  stop("unknown --sections: ", paste(setdiff(SECTIONS, ALL_SECTIONS), collapse = ", "))
}
want <- function(x) x %in% SECTIONS
# A minimum-K or recruitment table built on a power table from another run would pair
# two sets of random numbers without saying so, so both need power in the same run.
for (dep in c("mink", "recruit")) {
  if (want(dep) && !want("power")) stop("--sections: ", dep, " needs power in the same run")
}

LOG <- aval("--log", if (WRITES_REPORTED) "outputs/logs/klanguage_design_analysis.log"
                     else file.path(OUT, "klanguage_design_analysis.log"))
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
log_con <- file(LOG, open = "wt")
sink(log_con, split = TRUE)

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
KS   <- c(3, 6, 8, 10, 12, 16, 20)

written <- character(0)
write_out <- function(obj, name) {
  write.csv(obj, file.path(OUT, name), row.names = FALSE)
  written <<- c(written, name)
}

# Seeds. Every seed that already existed keeps the value the reported CSVs were
# produced with, so the arms the 28 September audit left alone regenerate
# unchanged. A matching seed does not make a matching table, though. varshare
# (5100 + K) keeps its seed, but its method changed from a two-way to a nested
# three-way decomposition, so its CSV will differ. The audit also changed which
# rows partly (5200) widens, how outcomes scales its own test, and how the
# calibration target is computed, so those tables move too, and the populations
# and partly-known tables gain columns. known7 (7000 + K) is new on 28 September.
#   Common random numbers are intended. Within an arm one seed serves every tau
#   regime, rho, prior scale, threshold and mode, so that differences between rows
#   are not blurred by independent noise. The pairing is partial, because
#   draw_population_params() uses a different amount of the random stream for
#   different populations (sample.int switches to Walker's alias method when many
#   grid points carry weight, and a tau of exactly zero draws nothing), so rows are
#   correlated, not paired.
#   Arms that loop over KS add the POSITION of K in KS, not K, so editing KS moves
#   every seed after the edit. Add new K values in a new arm instead.
SEEDS <- list(calibration = 20260910L, power = 4100L, conditional = 4150L,
              verb_sharing = 4200L, leave_one_out = 4300L, calibration_axis = 4400L,
              threshold = 4500L, mode = 4600L, outcomes = 4700L, budget = 4800L,
              verbs = 4900L, named = 5000L, varshare = 5100L, partly = 5200L,
              known7 = 7000L)

git_state <- function() {
  quiet <- function(expr) tryCatch(expr, error = function(e) character(0),
                                   warning = function(w) character(0))
  sha <- quiet(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE))
  if (!length(sha)) return("not a git checkout, or git not on the PATH")
  dirty <- quiet(system2("git", c("status", "--porcelain", "--", "R", "scripts", "tests"),
                         stdout = TRUE, stderr = FALSE))
  paste(sha, if (length(dirty)) "with uncommitted changes under R/, scripts/ or tests/" else "(clean)")
}

rule("PROVENANCE")
cat("started:  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("git:      ", git_state(), "\n")
cat("R:        ", R.version.string, "| replicates per cell:", REPS, "\n")
cat("sections: ", paste(SECTIONS, collapse = ", "), "\n")
cat("reading:  ", READING, if (READING == "alternative") paste("(SE rule:", ALT_SE, ")"), "\n")
cat("inputs:   ", IN, "\nout:      ", OUT, "\nlog:      ", LOG, "\n")

# ---------------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------------

H  <- read.csv(file.path(IN, "cross_language_effects_harmonised.csv"))
PP <- read.csv(file.path(IN, "pilot_params_ceilings.csv"))

# Per-language SE, split into the part that can be shared across languages and the
# part that cannot. The verb-information floor is a property of the verb count
# rather than of the sample, so it is flat in N above roughly 80 and is what the
# design actually delivers; it is the shareable part, because the teams translate
# one verb list. What is left once the floor is removed from the total per-language
# posterior SD is participant-driven, and participants are never shared.
#
# SE_POST is the posterior SD of the affectedness MAIN effect at the pilot samples
# (N = 70, 70, 57), standing in for the interaction's at N = 80. The two errors
# partly cancel: the interaction's SD rescaled to N = 80 and 72 verbs gives a total
# of 0.0871 against the 0.0885 used here, which moves no power figure by more than
# 0.0025 (audit item klanguage-engine#7). Rebuilding it waits for a sanctioned
# regeneration, because it rewrites every table.
SE_VERB <- mean(PP$se_floor_h1b  * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
SE_PPT  <- sqrt(max(SE_POST^2 - SE_VERB^2, 0))
SE_LANG <- sqrt(SE_VERB^2 + SE_PPT^2)

# How much of the verb-level noise is common across languages. Estimated two ways
# from the pilot and the Glossa data, both giving rho about 0.22, then diluted by
# the fraction of the verb list two languages actually share (0.83) and by the
# correlation of their affectedness ratings for a shared verb (0.89), since each
# team norms affectedness with its own raters. rho_eff = 0.22 * 0.83 * 0.89.
# The published Glossa model's rho = 1 is rejected by its own data at dBIC 2310.
RHO_EFF <- 0.164
split_at <- function(rho = RHO_EFF) split_se(SE_VERB, SE_PPT, rho)
sp_eff <- split_at(RHO_EFF)

# The pilot posteriors give the interaction's own posterior SD, which the harmonised
# table does not carry (it has the main effect's). The per-SD scaling is written out
# because the reported numbers were computed with these values; the check stops a
# changed pilot table from being paired with them.
CLAPS_LANGS <- c("English", "Turkish", "Norwegian")
CLAPS_SDX   <- c(0.5035005, 0.5035005, 0.5038147)
stopifnot(isTRUE(all.equal(CLAPS_SDX, PP$affectedness_sd[match(CLAPS_LANGS, PP$language)],
                           tolerance = 1e-6)))
claps_raw <- read_claps_interaction_draws(PILOT_MODELS, CLAPS_LANGS)
claps_se  <- setNames(vapply(seq_along(CLAPS_LANGS), function(i) sd(claps_raw[[i]]) * CLAPS_SDX[i],
                             numeric(1)), CLAPS_LANGS)
claps_pd  <- lapply(seq_along(CLAPS_LANGS), function(i) claps_raw[[i]] * CLAPS_SDX[i])
Hc <- prepare_language_table(H, claps_se)

input_files <- c(file.path(IN, c("cross_language_effects_harmonised.csv", "pilot_params_ceilings.csv")),
                 file.path(PILOT_MODELS, paste0("pilot_dgp_v2_pilot_", CLAPS_LANGS, ".rds")))
cat("\ninput md5:\n")
for (f in input_files) cat(sprintf("  %s  %s\n", tools::md5sum(f), f))

# The populations under each reading of the passive. The main sections use the one
# named by --passive; section 5d runs all three.
READING_VARIANTS <- list(
  canonical                = list(reading = "canonical",   alt_se = "independent"),
  alternative              = list(reading = "alternative", alt_se = "independent"),
  alternative_canonical_se = list(reading = "alternative", alt_se = "canonical")
)
POPS_BY_READING <- lapply(READING_VARIANTS, function(v) build_populations(Hc, v$reading, v$alt_se))
POPS <- build_populations(Hc, READING, ALT_SE)
# The Glossa regime stays the two-stage meta-analysis of the four 2023 languages. The
# one-stage refit of the published model (scripts/refit_glossa_pooled_tau.R, finished
# 25 September 2026) puts the by-language SD of the interaction at 0.452 as published,
# 0.525 with the verb factor split and 0.470 with per-language thresholds as well, so
# the corrected specification does not shrink the spread (0.470 / 0.452 = 1.04). That
# SD is a different quantity, estimated in one stage over five languages, so it
# confirms the regime's tau of about 0.53 and is not substituted for it.
REGIMES <- c("claps_protocol", "seven", "glossa_protocol")

rule("INPUTS")
cat(sprintf(paste("per-language SE, per SD of affectedness: total %.4f =",
                  "verb-driven %.4f (shareable) + participant-driven %.4f (never shared)\n"),
            SE_LANG, SE_VERB, SE_PPT))
cat(sprintf("cross-language correlation of verb-level effects, after dilution: rho = %.3f\n", RHO_EFF))
cat(sprintf("so each language's SE splits into %.4f independent and %.4f shared\n\n",
            sp_eff$indep, sp_eff$shared))
# The DESIGN prior here, not the calibrated analysis prior: this table describes what
# we believe about the languages, not what the confirmatory model assumes.
pop_tab <- do.call(rbind, lapply(names(POPS), function(nm) population_summary(POPS[[nm]], nm)))
print(pop_tab[, 1:9], row.names = FALSE, digits = 3)
cat("\nThe means barely differ. The spreads differ by a factor of",
    round(pop_tab$tau_median[pop_tab$population == "glossa_protocol"] /
          pop_tab$tau_median[pop_tab$population == "claps_protocol"], 1), "\n")
write_out(cbind(pop_tab, reading = READING, alt_se = if (READING == "alternative") ALT_SE else NA),
          "klanguage_populations.csv")

pop_by_reading <- do.call(rbind, lapply(names(READING_VARIANTS), function(rv) {
  v <- READING_VARIANTS[[rv]]
  P <- POPS_BY_READING[[rv]]
  cbind(reading_variant = rv, reading = v$reading,
        alt_se = if (v$reading == "alternative") v$alt_se else NA,
        do.call(rbind, lapply(names(P), function(nm) population_summary(P[[nm]], nm))))
}))
cat("\nThe same populations under the other reading of the passive (Balinese, Indonesian\n",
    "and Mandarin measured against their second passive-like construction):\n\n", sep = "")
print(pop_by_reading[pop_by_reading$population %in% REGIMES,
                     c("reading_variant", "population", "mean_effect", "mu_mean", "mu_sd",
                       "tau_median", "tau_lo", "tau_hi", "p_mu_negative")],
      row.names = FALSE, digits = 3)
write_out(pop_by_reading, "klanguage_populations_by_reading.csv")

# ---------------------------------------------------------------------------
# 2. Calibration against the three-language brms replicates
# ---------------------------------------------------------------------------
if (want("calib")) {
  rule("CALIBRATION AGAINST THE THREE-LANGUAGE brms REPLICATES")

  k3 <- read_pooled_k3_replicates(POOLED_K3)
  tg <- pooled_k3_target(k3)
  cat(sprintf(paste("%d of %d replicates finished and are scored (%d failed, all at N >= %d);",
                    "%d of the %d pass convergence_ok. Mean prior probability %.4f.\n"),
              tg$n_success, nrow(k3), tg$n_error,
              min(k3$n_participants[k3$status != "success"], Inf),
              tg$n_convergence_ok, tg$n_success, tg$mean_prior_prob))
  cat("The target pools every sample size of that run:\n")
  print(tg$by_n, row.names = FALSE, digits = 3)

  # As calibrated: each language at its verb-information floor, no shared term,
  # because each brms replicate gave its three languages separate verb lists.
  ix <- match(CLAPS_LANGS, PP$language)
  s_floor <- PP$se_floor_h1b[ix] * PP$affectedness_sd[ix]
  floor_scales <- c(CALIBRATION_BRACKET[1], CALIBRATED_TAU_PRIOR_SCALE, CALIBRATION_BRACKET[2])
  floor_rows <- do.call(rbind, lapply(floor_scales, function(tsc)
    calibrate_against_pooled_k3(claps_pd, s_floor, n_reps = REPS, prior_tau_scale = tsc,
                                seed = SEEDS$calibration)))
  # At the SEs every production cell uses (audit item klanguage-engine#4): the
  # independent part generated and analysed, plus the shared part generated only.
  prod_scales <- c(0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2)
  prod_rows <- do.call(rbind, lapply(prod_scales, function(tsc)
    calibrate_against_pooled_k3(claps_pd, rep(sp_eff$indep, 3), n_reps = REPS,
                                shared_sd = sp_eff$shared, se_analysis = rep(sp_eff$indep, 3),
                                prior_tau_scale = tsc, seed = SEEDS$calibration)))
  target_row <- data.frame(tau_prior_scale = NA, t(tg$target), n = tg$n_success)
  calib <- rbind(target_row, floor_rows, prod_rows)
  calib$row <- c(sprintf("brms target (%d replicates, N = %s, exact scoring)", tg$n_success,
                         paste(range(tg$by_n$N), collapse = "-")),
                 paste0("engine, verb-floor SEs as calibrated, tau prior scale ", floor_scales),
                 paste0("engine, production SEs, tau prior scale ", prod_scales))
  calib$se_config <- c("brms", rep("verb_floor", length(floor_scales)),
                       rep("production", length(prod_scales)))
  calib <- calib[, c("tau_prior_scale", "mean_z", "sd_z", "p_bf10", "p_bf6", "p_bf3", "row",
                     "se_config", "n")]
  print(calib[, c("row", "mean_z", "sd_z", "p_bf10", "p_bf6", "p_bf3")], row.names = FALSE, digits = 3)
  pr <- calib[calib$se_config == "production", ]
  best <- pr$tau_prior_scale[which.min(abs(pr$p_bf10 - tg$target[["p_bf10"]]))]
  cat(sprintf(paste("\nAt the production SEs the scale closest to the target P(BF >= 10) of %.3f is %.2f,",
                    "and sd z stays between %.3f and %.3f against %.3f.\n"),
              tg$target[["p_bf10"]], best, min(pr$sd_z), max(pr$sd_z), tg$target[["sd_z"]]))
  cat("The analysis keeps CALIBRATED_TAU_PRIOR_SCALE =", CALIBRATED_TAU_PRIOR_SCALE,
      "until the maintainer chooses between the two SE configurations.\n")
  write_out(calib, "klanguage_calibration.csv")
  write_out(k3, "klanguage_calibration_target_replicates.csv")
}

# ---------------------------------------------------------------------------
# 3. The three regimes
# ---------------------------------------------------------------------------
# The analysis is told each language's independent SE only and never sees the common
# verb offset the data carry. A pooled model with verb IDs prefixed by language, which
# is what R/11 and the pilot use, treats each language's whole error as independent;
# passing that total instead of the independent part moves power by at most 0.002. If
# Stage 2 shares translated verb IDs across languages, the exact analysis adds
# shared^2 to the variance of mu, which lowers power by 0.003 to 0.0075 at rho_eff and
# by about 0.04 at rho = 1 (audit item klanguage-engine#15). Which applies depends on
# how Stage 2 codes its verbs.
cell <- function(K, pop_name, mode = "assurance", rho = RHO_EFF,
                 tsc = CALIBRATED_TAU_PRIOR_SCALE, thr = 10, seed) {
  sp <- split_at(rho)
  cbind(population = pop_name, rho = rho,
        run_klanguage_cell(K = K, n_reps = REPS, se_lang = sp$indep, mode = mode,
                           pop = POPS[[pop_name]], shared_verb_sd = sp$shared,
                           prior_tau_scale = tsc, bf_threshold = thr, seed = seed))
}

if (want("power")) {
  rule(sprintf("P(BF >= 10) FOR THE POOLED INTERACTION, BY tau REGIME (%d replicates a cell)", REPS))

  power <- do.call(rbind, lapply(REGIMES, function(pn)
    do.call(rbind, lapply(seq_along(KS), function(i) cell(KS[i], pn, seed = SEEDS$power + i)))))

  tb <- tapply(power$power, list(power$population, power$K), function(x) x[1])
  cat("\nASSURANCE: (mu, tau) drawn from their joint posterior given that regime's\n",
      "languages, so the figure carries our uncertainty about tau as well as the\n",
      "sampling variation. This is the honest headline.\n\n", sep = "")
  print(round(tb[REGIMES, ], 3))
  cat("\nMonte-Carlo error on each entry is at most",
      sprintf("%.3f", max(power$mcse)), "\n")
  write_out(power, "klanguage_power.csv")
}

# Conditional power: tau treated as KNOWN at each regime's posterior median. The
# gap between this table and the one above is the price of not knowing tau, and it
# is the larger of the two effects at every K. Separating them matters because they
# have different remedies: a small conditional power is fixed by recruiting more
# languages, and a large gap is not.
if (want("cond")) {
  cond <- do.call(rbind, lapply(REGIMES, function(pn) {
    post <- ra_meta_posterior(POPS[[pn]]$y, POPS[[pn]]$s,
                              prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
    do.call(rbind, lapply(seq_along(KS), function(i) {
      r <- run_klanguage_cell(K = KS[i], n_reps = REPS, se_lang = sp_eff$indep, mode = "point",
                              mu_point = post$mu_mean, tau_point = post$tau_median,
                              shared_verb_sd = sp_eff$shared, bf_threshold = 10,
                              seed = SEEDS$conditional + i)
      cbind(population = pn, mu = post$mu_mean, tau = post$tau_median, r)
    }))
  }))
  cat("\nCONDITIONAL on tau being known at that regime's posterior median:\n\n")
  print(round(tapply(cond$power, list(cond$population, cond$K), function(x) x[1])[REGIMES, ], 3))
  cat("\ntau assumed:", paste(sprintf("%s %.2f", REGIMES,
      vapply(REGIMES, function(p) cond$tau[cond$population == p][1], numeric(1))), collapse = "  "), "\n")
  cat("The gap between the two tables is what not knowing tau costs. It does not\n",
      "close as K rises, because six to twelve languages cannot pin tau down either.\n", sep = "")
  write_out(cond, "klanguage_conditional.csv")
}

# ---------------------------------------------------------------------------
# 4. Sensitivities
# ---------------------------------------------------------------------------
if (want("sens")) {
  rule("SENSITIVITIES (seven-language regime unless stated)")

  sens <- list()
  add <- function(tag, df) sens[[tag]] <<- cbind(axis = tag, df)

  add("verb_sharing", do.call(rbind, lapply(c(0, RHO_EFF, 0.22, 1), function(rho)
    do.call(rbind, lapply(seq_along(KS), function(i)
      cell(KS[i], "seven", rho = rho, seed = SEEDS$verb_sharing + i))))))

  add("leave_one_out", do.call(rbind, lapply(c("seven", "drop_balinese", "drop_hebrew"), function(pn)
    do.call(rbind, lapply(seq_along(KS), function(i) cell(KS[i], pn, seed = SEEDS$leave_one_out + i))))))

  add("calibration", do.call(rbind, lapply(c(CALIBRATION_BRACKET[1], CALIBRATED_TAU_PRIOR_SCALE,
                                             CALIBRATION_BRACKET[2]), function(tsc)
    do.call(rbind, lapply(seq_along(KS), function(i)
      cell(KS[i], "seven", tsc = tsc, seed = SEEDS$calibration_axis + i))))))

  add("threshold", do.call(rbind, lapply(c(3, 6, 10), function(thr)
    do.call(rbind, lapply(seq_along(KS), function(i)
      cell(KS[i], "seven", thr = thr, seed = SEEDS$threshold + i))))))

  add("mode", do.call(rbind, lapply(c("assurance", "safeguard"), function(m)
    do.call(rbind, lapply(seq_along(KS), function(i)
      cell(KS[i], "seven", mode = m, seed = SEEDS$mode + i))))))

  sensitivity <- do.call(rbind, sens)
  KEY <- c(verb_sharing = "rho", leave_one_out = "population", calibration = "prior_tau_scale",
           threshold = "bf_threshold", mode = "mode")
  for (a in unique(sensitivity$axis)) {
    d <- sensitivity[sensitivity$axis == a, ]
    cat("\n", a, " (rows = ", KEY[[a]], ", columns = K)\n", sep = "")
    print(round(tapply(d$power, list(as.character(d[[KEY[[a]]]]), d$K), function(x) x[1]), 3))
  }
  write_out(sensitivity, "klanguage_sensitivity.csv")
}

# ---------------------------------------------------------------------------
# 5. Smallest K reaching each target
# ---------------------------------------------------------------------------
if (want("mink")) {
  rule("SMALLEST K REACHING EACH TARGET, BY REGIME")

  min_k <- do.call(rbind, lapply(c(0.80, 0.90, 0.95), function(tgt)
    do.call(rbind, lapply(REGIMES, function(pn) {
      d <- power[power$population == pn, ]
      d <- d[order(d$K), ]
      hit <- d$K[d$power >= tgt]
      data.frame(target = tgt, regime = pn,
                 min_K = if (length(hit)) min(hit) else NA_integer_,
                 at_6 = d$power[d$K == 6], at_8 = d$power[d$K == 8],
                 at_10 = d$power[d$K == 10], at_12 = d$power[d$K == 12])
    }))))
  print(min_k, row.names = FALSE, digits = 3)
  cat("\nNA means no K up to", max(KS), "reaches the target.\n")
  write_out(min_k, "klanguage_minimum_k.csv")
}

# ---------------------------------------------------------------------------
# 5b. The named-language arm: which languages join, not how many
# ---------------------------------------------------------------------------
# Everything above draws languages exchangeably from a population, which is the
# only way to say anything about teams that have not been recruited, and is also
# the step a reviewer will press hardest: at K = 12 it means nine of the twelve
# languages are inventions. For K <= 7 that step is avoidable, and the PI named
# seven specific languages, so this arm answers what he literally asked.
if (want("named")) {
  rule("NAMED LANGUAGES: EVERY SUBSET OF THE SEVEN, NO SYNTHETIC LANGUAGE ANYWHERE")

  subsets <- do.call(rbind, lapply(c(5, 6, 7), function(K)
    enumerate_language_subsets(POPS$seven, K, sp_eff$indep, sp_eff$shared,
                               n_reps = REPS, seed = SEEDS$named + K)))
  cat("\nAt K = 6, one language is left out. Power by which one:\n")
  s6 <- subsets[subsets$K == 6, c("omitted", "realised_tau", "mean_effect", "power", "mcse")]
  print(s6[order(-s6$power), ], row.names = FALSE, digits = 3)
  cat("\nAll seven together (K = 7):",
      sprintf("%.3f\n", subsets$power[subsets$K == 7]))
  cat("Range across the seven K = 6 subsets:",
      sprintf("%.3f to %.3f", min(s6$power), max(s6$power)),
      "- a", sprintf("%.1f", max(s6$power) / max(min(s6$power), 0.001)),
      "fold spread at FIXED K.\n")
  if (want("power")) {
    cat("Exchangeable arm at K = 6, for comparison:",
        sprintf("%.3f\n", power$power[power$population == "seven" & power$K == 6]))
  }
  cat("\nThe gap between the named and exchangeable arms at K = 7 is the\n",
      "fixed-versus-random-language estimand difference. It is not a discrepancy\n",
      "to be reconciled: the two arms answer different questions.\n", sep = "")
  write_out(subsets, "klanguage_named_subsets.csv")
}

# Where the variation in the outcome comes from, at a fixed number of languages:
# uncertainty about (mu, tau), which languages join, and sampling noise. K is fixed
# within a row, so how much the NUMBER of languages matters is read across rows.
if (want("varshare")) {
  rule("WHERE THE OUTCOME VARIANCE COMES FROM, AT FIXED K")

  DECOMP <- list(n_par = as.integer(aval("--decomp-par", max(200L, REPS %/% 4L))),
                 n_sets = 10L, n_reps = 100L, n_boot = 200L)
  cat(sprintf("%d draws of (mu, tau) under the design prior, %d language sets each, %d studies a set\n\n",
              DECOMP$n_par, DECOMP$n_sets, DECOMP$n_reps))
  vshare <- do.call(rbind, lapply(c(6, 8, 10, 12, 20), function(K)
    decompose_outcome_variance(K, POPS$seven, sp_eff$indep, sp_eff$shared,
                               n_par = DECOMP$n_par, n_sets = DECOMP$n_sets,
                               n_reps = DECOMP$n_reps, n_boot = DECOMP$n_boot,
                               seed = SEEDS$varshare + K)))
  print(vshare[, c("K", "mean_power", "mcse_power", "share_population_params", "share_composition",
                   "share_noise", "share_composition_given_population", "share_true_effects")],
        row.names = FALSE, digits = 3)
  cat("\nshare_population_params is our uncertainty about how large and how variable the\n",
      "effect is across languages; share_composition is which languages join, given the\n",
      "population; share_noise is sampling error. share_true_effects adds the first two,\n",
      "the quantity an earlier version of this table reported as composition alone.\n", sep = "")
  write_out(vshare, "klanguage_set_variance.csv")
}

# ---------------------------------------------------------------------------
# 5c. Sets in which some of the languages have already been measured
# ---------------------------------------------------------------------------
# The exchangeable arm treats every recruited language as unknown. That is right
# for a language nobody has measured and too pessimistic for one we have piloted,
# and the English, Turkish and Norwegian teams are the likeliest members of any
# Stage 2 set. What is known about them is their true effect, not their data: no
# pilot observation enters the Stage 2 analysis, and none enters here. The pilot
# supplies a design prior on where that language sits, which is what the rest of
# this package already uses it for.
#
# The Ambridge, Arnon & Bekman (2023) languages are a weaker case and are treated
# as one. If such a team joins, the 2023 estimate bears on what CLAPS would find
# there only as far as the protocol transfers, so the prior on that language, and
# on no CLAPS pilot language, is widened by the transfer SD. At a transfer SD equal
# to tau the earlier measurement is worth nothing and the language is an unknown
# draw again.
claps3  <- POPS$claps_protocol
glossa4 <- POPS$glossa_protocol
g <- function(...) glossa4[match(c(...), glossa4$language), ]

if (want("partly")) {
  rule("SETS WITH SOME LANGUAGES ALREADY MEASURED")

  # The scenarios name their languages. An earlier version took "two of the 2023
  # four" by position, which silently meant Balinese and Hebrew, and so reported the
  # cost of recruiting the one wrong-signed language as though it were the general
  # effect of knowing two more. Which languages are known matters here for exactly
  # the reason it matters in section 5b, so both ends are shown.
  known_scenarios <- list(
    list(tag = "a) nothing known: every language drawn",     known = claps3[0, ], tsd = 0),
    list(tag = "b) the 3 CLAPS pilot teams",                 known = claps3,      tsd = 0),
    list(tag = "c) b + Hebrew and Mandarin, transfer good",  known = rbind(claps3, g("Hebrew", "Mandarin")),   tsd = 0.15),
    list(tag = "d) b + Balinese and Indonesian, transfer good", known = rbind(claps3, g("Balinese", "Indonesian")), tsd = 0.15),
    list(tag = "e) as d, transfer poor",                     known = rbind(claps3, g("Balinese", "Indonesian")), tsd = 0.40)
  )

  partly <- do.call(rbind, lapply(seq_along(known_scenarios), function(j) {
    sc <- known_scenarios[[j]]
    do.call(rbind, lapply(c(6, 8, 10, 12), function(K)
      cbind(scenario = sc$tag,
            run_partly_known_cell(K = K, known = sc$known, n_reps = REPS,
                                  se_lang = sp_eff$indep, pop = POPS$seven,
                                  shared_verb_sd = sp_eff$shared,
                                  transfer_sd = transfer_sd_by_source(sc$known, sc$tsd),
                                  seed = SEEDS$partly + 10L * j + K))))
  }))
  stopifnot(!anyDuplicated(partly[c("scenario", "K")]))
  print(round(tapply(partly$power, list(partly$scenario, partly$K), function(x) x[1]), 3))
  cat("\nRow a is the exchangeable arm restated. The others take the named teams' effects\n",
      "from their earlier estimates and draw only the remaining places. Rows c and d differ\n",
      "only in WHICH two of the 2023 languages take part. Row e widens the prior on those\n",
      "two 2023 languages alone; the CLAPS teams are never widened.\n", sep = "")
  write_out(partly, "klanguage_partly_known.csv")
}

# ---------------------------------------------------------------------------
# 5d. All seven measured languages known, the rest drawn: the PI's construction
# ---------------------------------------------------------------------------
# On 23 September the PI asked for the pooled test at 10 and 20 languages when the
# seven languages already measured all take part at their estimated effects and only
# the other three or thirteen are new. This is section 5c with every measured
# language known. Four versions of "known" are run, because the phrase admits more
# than one reading and the answer depends on which:
#   drawn, 0.15   each known language's effect is drawn from its earlier posterior,
#                 the three CLAPS ones as they are and the four 2023 ones widened by
#                 0.15 for doubt about protocol transfer.
#   drawn, 0.40   the same with transfer doubted.
#   fixed         each known language held exactly at its estimate, with no
#                 uncertainty about it at all (the named arm's convention at K = 7).
#                 This is the literal reading of "kept at their estimated effects",
#                 so it is reported beside "drawn", not relegated to a footnote. The
#                 two differ by about 19 points at K = 7 and 4 to 5 points at K = 10.
#   none known    every place drawn, for reference.
# Not every pairing of a transfer doubt with a tau regime is coherent. The CLAPS
# regime assumes that languages measured under the CLAPS protocol spread by only
# about 0.155, so most of the 2023 languages' spread of about 0.53 would come from
# their protocol. That same premise says a 2023 estimate transfers poorly: the
# protocol-induced deviation it implies is roughly sqrt(0.53^2 - 0.155^2), about
# 0.5, far more than 0.15. The coherent row for the CLAPS regime is therefore the
# 0.40 one. The 0.15 doubt fits the seven and Glossa regimes, which do not
# attribute the 2023 spread to protocol. Every cell is still computed, so that the
# pairing is a choice made in the reporting, in the open.
# Each runs in assurance and conditional mode under each tau regime, and under both
# readings of the passive; the alternative reading changes the known Balinese,
# Indonesian and Mandarin effects as well as the populations. At K = 7 nothing is
# drawn, so the regime and the mode cannot matter, and the check below makes sure
# they do not. Every cell uses seed 7000 + K, so readings, regimes and variants are
# compared on common random numbers.
if (want("known7")) {
  rule("ALL SEVEN MEASURED LANGUAGES KNOWN, THE REMAINING PLACES DRAWN")

  KS_KNOWN7 <- c(7, 8, 10, 12, 16, 20)
  known7 <- do.call(rbind, lapply(names(READING_VARIANTS), function(rv) {
    P <- POPS_BY_READING[[rv]]
    k7 <- P$seven
    variants <- list(
      list(known_set = "seven, drawn, 2023 widened 0.15", known = k7, truth = "drawn", tsd = 0.15),
      list(known_set = "seven, drawn, 2023 widened 0.40", known = k7, truth = "drawn", tsd = 0.40),
      list(known_set = "seven, fixed at estimates",       known = k7, truth = "fixed", tsd = 0),
      list(known_set = "none known",                      known = k7[0, ], truth = "drawn", tsd = 0)
    )
    do.call(rbind, lapply(REGIMES, function(pn) {
      post <- ra_meta_posterior(P[[pn]]$y, P[[pn]]$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
      do.call(rbind, lapply(c("assurance", "point"), function(md) {
        do.call(rbind, lapply(variants, function(v) {
          do.call(rbind, lapply(KS_KNOWN7, function(K) {
            cbind(reading_variant = rv, reading = READING_VARIANTS[[rv]]$reading,
                  alt_se = if (READING_VARIANTS[[rv]]$reading == "alternative") READING_VARIANTS[[rv]]$alt_se else NA,
                  population = pn, known_set = v$known_set, tsd_glossa = v$tsd,
                  run_partly_known_cell(K = K, known = v$known, n_reps = REPS,
                                        se_lang = sp_eff$indep, pop = P[[pn]],
                                        shared_verb_sd = sp_eff$shared,
                                        transfer_sd = transfer_sd_by_source(v$known, v$tsd),
                                        mode = md, mu_point = post$mu_mean, tau_point = post$tau_median,
                                        known_truth = v$truth, seed = SEEDS$known7 + K))
          }))
        }))
      }))
    }))
  }))

  # Nothing is drawn at K = 7 when all seven are known, so every regime and mode must
  # give the same answer; a difference would mean the population leaked into it.
  k7_all <- known7[known7$K == 7 & known7$n_known == 7, ]
  spread <- tapply(k7_all$power, list(k7_all$reading_variant, k7_all$known_set), function(x) diff(range(x)))
  stopifnot(all(spread == 0, na.rm = TRUE))

  for (rv in names(READING_VARIANTS)) {
    cat("\nReading:", rv, "\n")
    for (md in c("assurance", "point")) {
      d <- known7[known7$reading_variant == rv & known7$mode == md, ]
      tab <- tapply(d$power, list(paste(d$population, "|", d$known_set), d$K), function(x) x[1])
      cat(if (md == "assurance") "  assurance, (mu, tau) drawn for the new places:\n"
          else "  conditional, (mu, tau) at the regime's posterior mean and median:\n")
      print(round(tab, 3))
    }
  }
  cat("\nMonte-Carlo error on each entry is at most", sprintf("%.3f", max(known7$mcse)), "\n")
  write_out(known7, "klanguage_known7.csv")
}

# ---------------------------------------------------------------------------
# 6. The two outcomes the protocol has no sentence for
# ---------------------------------------------------------------------------
# A cross-linguistic claim can fail in ways a power figure does not describe. Both
# of these are computed on the SAME simulated studies as the power above, so the
# pooled and per-language results are the joint outcome of one dataset rather than
# two independent calculations.
outcome_cell <- function(K, pop_name, seed) {
  set.seed(seed)
  sp   <- split_at(RHO_EFF)
  pars <- draw_population_params(REPS, POPS[[pop_name]]$y, POPS[[pop_name]]$s)
  z_thr <- qnorm(10 / 11)
  # A single-language analysis reports the whole of its verb-sampling error, shared
  # or not, so its own test is standardised by the total SE. Dividing by the
  # independent part alone inflated every |z| by 6% (audit item klanguage-engine#6).
  se_own <- sqrt(sp$indep^2 + sp$shared^2)
  res <- vapply(seq_len(REPS), function(i) {
    theta  <- rnorm(K, pars$mu[i], pars$tau[i])
    common <- rnorm(1, 0, sp$shared)
    y      <- theta + common + rnorm(K, 0, sp$indep)
    # Each language's own test, at its own precision. The prediction is that the
    # interaction is NEGATIVE, so a language that supports it has a negative
    # estimate and a language that contradicts it at the same threshold has a
    # positive one. Getting these two the wrong way round reports the count of
    # confirming languages as the count of contradicting ones.
    z_lang <- y / se_own
    right  <- sum(z_lang <= -z_thr)          # reaches threshold, predicted direction
    wrong  <- sum(z_lang >=  z_thr)          # reaches threshold, opposite direction
    # The pooled half is analysed as in section 3 (see the note on cell()).
    p      <- ra_meta_posterior(y, rep(sp$indep, K))
    pooled <- directional_bf_meta(p) >= 10
    c(wrong = wrong, right = right, pooled = pooled)
  }, numeric(3))
  data.frame(
    population = pop_name, K = K,
    e_wrong_sign      = mean(res["wrong", ]),
    p_any_wrong_sign  = mean(res["wrong", ] >= 1),
    # The pooled test fails while at least half the languages individually succeed.
    p_contradiction   = mean(res["pooled", ] == 0 & res["right", ] >= K / 2),
    p_pooled          = mean(res["pooled", ] == 1))
}

if (want("outcomes")) {
  rule("OUTCOMES THE PROTOCOL NEEDS A PRESPECIFIED READING OF")

  outcomes <- do.call(rbind, lapply(REGIMES, function(pn)
    do.call(rbind, lapply(c(6, 8, 10, 12), function(K)
      outcome_cell(K, pn, seed = SEEDS$outcomes + K)))))
  print(outcomes, row.names = FALSE, digits = 3)
  cat("\ne_wrong_sign is the expected NUMBER of languages whose own test reaches a\n",
      "directional Bayes factor of 10 in the direction OPPOSITE to the prediction.\n",
      "Balinese, at +0.33 per SD, is already such a language.\n", sep = "")
  write_out(outcomes, "klanguage_outcomes.csv")
}

# ---------------------------------------------------------------------------
# 7. The trades the PI will propose
# ---------------------------------------------------------------------------
# Per-language SE as a function of participants and verbs. The verb term is the
# floor and scales as 1/sqrt(n_verbs); the participant term scales as 1/sqrt(N).
# The anchor is nominal: SE_VERB and SE_PPT are treated as the values at N = 80 and
# 72 verbs, although they were measured on pilots of N = 70, 70 and 57 with 72, 72
# and 66 verbs (see the note on SE_POST above). The mis-anchoring is about 2% of the
# total SE at N = 40, too small to change a trade.
se_for <- function(N = 80, n_verbs = 72) {
  v <- SE_VERB * sqrt(72 / n_verbs)
  p <- SE_PPT  * sqrt(80 / N)
  list(verb = v, total = sqrt(v^2 + p^2),
       indep = sqrt(v^2 * (1 - RHO_EFF) + p^2), shared = v * sqrt(RHO_EFF))
}
trade_cell <- function(K, N, n_verbs, pop_name, seed) {
  s <- se_for(N, n_verbs)
  r <- run_klanguage_cell(K = K, n_reps = REPS, se_lang = s$indep, mode = "assurance",
                          pop = POPS[[pop_name]], shared_verb_sd = s$shared,
                          bf_threshold = 10, seed = seed)
  data.frame(population = pop_name, K = K, N = N, n_verbs = n_verbs,
             total_participants = K * N, se_lang = s$total, power = r$power, mcse = r$mcse)
}

if (want("trades")) {
  rule("TRADES: PARTICIPANTS AGAINST LANGUAGES, VERBS AGAINST LANGUAGES")

  # a. Fixed participant budget of 640: buy languages or buy participants?
  budget <- do.call(rbind, lapply(list(c(8, 80), c(10, 64), c(12, 53), c(16, 40)), function(kn)
    trade_cell(kn[1], kn[2], 72, "seven", seed = SEEDS$budget + kn[1])))
  cat("\na. Fixed budget of 640 participants (seven-language regime)\n")
  print(budget, row.names = FALSE, digits = 3)

  # b. Verb count. Three of the five Glossa languages could not field 72 verbs
  #    (Balinese 47, Hebrew 56, Mandarin 57), yet the design assumes 72 everywhere.
  verbs <- do.call(rbind, lapply(c(40, 47, 56, 72), function(v)
    do.call(rbind, lapply(c(6, 8, 12), function(K)
      trade_cell(K, 80, v, "seven", seed = SEEDS$verbs + v + K)))))
  cat("\nb. Verb count, at N = 80 (seven-language regime)\n")
  print(round(tapply(verbs$power, list(verbs$n_verbs, verbs$K), function(x) x[1]), 3))
  cat("\nA team that can only field 47 verbs costs the pooled test very little, while\n",
      "forfeiting its own confirmatory test. That is a recruitment-policy fact.\n", sep = "")

  tradeoffs <- rbind(cbind(trade = "participant_budget", budget),
                     cbind(trade = "verb_count", verbs))
  write_out(tradeoffs, "klanguage_tradeoffs.csv")
}

# ---------------------------------------------------------------------------
# 8. Recruitment: K is a random variable and Y is a rule about it
# ---------------------------------------------------------------------------
if (want("recruit")) {
  rule("RECRUITMENT: WHAT A MINIMUM OF Y LANGUAGES COSTS")

  pw <- setNames(power$power[power$population == "seven"], power$K[power$population == "seven"])
  pw_at <- function(k) {
    ks <- as.numeric(names(pw))
    if (k < min(ks)) return(pw[[1]])
    if (k > max(ks)) return(pw[[length(pw)]])
    approx(ks, as.numeric(pw), xout = k)$y
  }
  recruit <- do.call(rbind, lapply(c(10, 12, 14, 16), function(invited)
    do.call(rbind, lapply(c(0.5, 0.6, 0.7, 0.8), function(p)
      do.call(rbind, lapply(c(6, 8, 10), function(Y) {
        k    <- 0:invited
        pk   <- dbinom(k, invited, p)
        data.frame(invited = invited, completion = p, Y = Y,
                   E_K = invited * p,
                   p_below_Y = sum(pk[k < Y]),
                   # Power averaged over the K values that clear the rule, i.e. the
                   # operating characteristic of the study that actually goes ahead.
                   power_given_proceed = if (sum(pk[k >= Y]) > 0)
                     sum(pk[k >= Y] * vapply(k[k >= Y], pw_at, numeric(1))) / sum(pk[k >= Y])
                   else NA_real_)
      }))))))
  print(recruit[recruit$completion %in% c(0.6, 0.8), ], row.names = FALSE, digits = 3)
  cat("\np_below_Y is the chance the rule bites and the project does not proceed to Stage 2.\n")
  write_out(recruit, "klanguage_recruitment.csv")
}

# ---------------------------------------------------------------------------
rule("FILES WRITTEN")
for (f in written) cat(" ", file.path(OUT, f), "\n")
cat("\nfinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    "| R:", R.version.string, "| replicates per cell:", REPS, "\n")
sink()
close(log_con)
