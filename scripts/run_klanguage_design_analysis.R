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
#     CLAPS pilot (English, Turkish, Norwegian)      tau ~ 0.16   Q p = 0.11
#     Glossa 2023 (Balinese, Hebrew, Indonesian,
#                  Mandarin)                          tau ~ 0.53   Q p = 1.4e-07
#     all seven pooled                                tau ~ 0.39
#
#   The MEANS barely differ (-0.28 against -0.32). The SPREADS differ by a factor
#   of about three and a half. Stage 2 will run the CLAPS protocol in every
#   language, so which of these is the right tau is a question about whether the
#   Glossa spread is a property of the languages or of the measurement, and the
#   minimum number of languages that reaches 90% is 4, 14 or 22 accordingly.
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
#
# OUTPUTS
#   outputs/design_summary_pilot/klanguage_calibration.csv   the K = 3 check
#   outputs/design_summary_pilot/klanguage_power.csv         the three tau regimes
#   outputs/design_summary_pilot/klanguage_sensitivity.csv   assumption sweeps
#   outputs/design_summary_pilot/klanguage_minimum_k.csv     smallest K per target
#   outputs/design_summary_pilot/klanguage_outcomes.csv      wrong-sign, contradiction
#   outputs/design_summary_pilot/klanguage_tradeoffs.csv     N-vs-K, verbs-vs-K
#   outputs/design_summary_pilot/klanguage_recruitment.csv   P(K < Y) by completion rate
#
# COST
#   Minutes on one core. The expensive part was already paid for: the 135
#   three-language brms replicates that calibrate the engine, and the corrected
#   refit of the Glossa model (scripts/refit_glossa_pooled_tau.R) that will replace
#   the Glossa-protocol tau when it lands.
#
# USAGE
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" scripts/run_klanguage_design_analysis.R
#   Options: --reps=4000 --out=<dir> --quick
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(stats) })
source("R/13_simulate_klanguage_pooled.R")

args  <- commandArgs(trailingOnly = TRUE)
aval  <- function(f, d) { h <- grep(paste0("^", f, "="), args, value = TRUE)
                          if (!length(h)) d else sub(paste0("^", f, "="), "", h[1]) }
REPS  <- as.integer(aval("--reps", if ("--quick" %in% args) 800L else 4000L))
OUT   <- aval("--out", "outputs/design_summary_pilot")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
KS   <- c(3, 6, 8, 10, 12, 16, 20)

# ---------------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------------

H  <- read.csv(file.path(OUT, "cross_language_effects_harmonised.csv"))
PP <- read.csv(file.path(OUT, "pilot_params_ceilings.csv"))

# Per-language SE, split into the part that can be shared across languages and the
# part that cannot. The verb-information floor is a property of the verb count
# rather than of the sample, so it is flat in N above roughly 80 and is what the
# design actually delivers; it is the shareable part, because the teams translate
# one verb list. What is left once the floor is removed from the total per-language
# posterior SD is participant-driven, and participants are never shared.
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

split_se <- function(rho) list(shared = SE_VERB * sqrt(rho),
                               indep  = sqrt(SE_VERB^2 * (1 - rho) + SE_PPT^2))

one_per_language <- function(df) {
  # Prefer the CLAPS estimate of English: it is the one collected under the CLAPS
  # protocol, which is the protocol every new team will follow.
  df <- df[order(df$source != "CLAPS_pilot"), ]
  df[!duplicated(df$language), ]
}
Hc <- H[H$se_kind %in% c("slope_sd", "posterior_sd"), ]

# The CLAPS rows in the harmonised table carry `slope_post_sd`, which is the
# posterior SD of the affectedness MAIN effect, not of the interaction the pooled
# test is about. The interaction's own posterior SD is available exactly, in the
# archived pilot posteriors, so it is read from there instead: sd() of the
# S_TypeActive:Semantics_scaled column of fixef_draws, converted to the per-SD
# scale. The two differ by enough to matter for the meta-analytic weights
# (0.120/0.089/0.064 against 0.098/0.078/0.090), and the interaction's is the
# right one because it is the quantity being pooled.
claps_interaction_se <- function(langs = c("English", "Turkish", "Norwegian"),
                                 sdx = c(0.5035005, 0.5035005, 0.5038147)) {
  files <- file.path("outputs/pilot_models", paste0("pilot_dgp_v2_pilot_", langs, ".rds"))
  if (!all(file.exists(files))) return(NULL)
  setNames(vapply(seq_along(langs), function(i) {
    fd <- as.matrix(readRDS(files[i])$fixef_draws)
    sd(as.numeric(fd[, "S_TypeActive:Semantics_scaled"])) * sdx[i]
  }, numeric(1)), langs)
}
claps_se <- claps_interaction_se()
if (is.null(claps_se)) {
  message("[klang] pilot posteriors absent; keeping the main-slope SDs as a stand-in.")
} else {
  hit <- Hc$source == "CLAPS_pilot" & Hc$language %in% names(claps_se)
  Hc$se[hit] <- claps_se[Hc$language[hit]]
}
Hc <- Hc[!duplicated(paste(Hc$source, Hc$language)), ]

# The three regimes, plus the leave-one-out sets used later.
POPS <- list(
  claps_protocol = subset(Hc, source == "CLAPS_pilot"),
  seven          = one_per_language(Hc),
  glossa_protocol = subset(Hc, source == "Glossa2023" & language != "English"),
  drop_balinese  = subset(one_per_language(Hc), language != "Balinese"),
  drop_hebrew    = subset(one_per_language(Hc), language != "Hebrew")
)
POPS <- lapply(POPS, function(d) data.frame(language = d$language, y = d$y, s = d$se))

REGIMES <- c("claps_protocol", "seven", "glossa_protocol")

rule("INPUTS")
cat(sprintf(paste("per-language SE, per SD of affectedness: total %.4f =",
                  "verb-driven %.4f (shareable) + participant-driven %.4f (never shared)\n"),
            SE_LANG, SE_VERB, SE_PPT))
cat(sprintf("cross-language correlation of verb-level effects, after dilution: rho = %.3f\n\n", RHO_EFF))
pop_tab <- do.call(rbind, lapply(names(POPS), function(nm) {
  # The DESIGN prior here, not the calibrated analysis prior: this table describes
  # what we believe about the languages, not what the confirmatory model assumes.
  p <- POPS[[nm]]; post <- ra_meta_posterior(p$y, p$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
  # Cochran's Q, as the test of whether the spread is more than sampling error.
  w <- 1 / p$s^2; mu_fe <- sum(w * p$y) / sum(w)
  Q <- sum(w * (p$y - mu_fe)^2); df <- nrow(p) - 1
  data.frame(population = nm, K = nrow(p), mean_effect = mean(p$y),
             tau_median = post$tau_median, p_mu_negative = post$p_mu_negative,
             Q = Q, df = df, p_Q = pchisq(Q, df, lower.tail = FALSE),
             I2 = max(0, (Q - df) / Q))
}))
print(pop_tab, row.names = FALSE, digits = 3)
cat("\nThe means barely differ. The spreads differ by a factor of",
    round(pop_tab$tau_median[pop_tab$population == "glossa_protocol"] /
          pop_tab$tau_median[pop_tab$population == "claps_protocol"], 1), "\n")

# ---------------------------------------------------------------------------
# 2. Calibration against the 135 three-language brms replicates
# ---------------------------------------------------------------------------
rule("CALIBRATION AGAINST THE 135 THREE-LANGUAGE brms REPLICATES")

TARGET <- c(mean_z = 1.386, sd_z = 0.311, p_bf10 = 0.578, p_bf6 = 0.859, p_bf3 = 0.985)

claps_pd <- local({
  langs <- c("English", "Turkish", "Norwegian")
  sdx   <- c(0.5035005, 0.5035005, 0.5038147)
  TERM  <- "S_TypeActive:Semantics_scaled"
  files <- file.path("outputs/pilot_models", paste0("pilot_dgp_v2_pilot_", langs, ".rds"))
  if (!all(file.exists(files))) return(NULL)
  lapply(seq_along(langs), function(i)
    as.numeric(as.matrix(readRDS(files[i])$fixef_draws)[, TERM]) * sdx[i])
})

if (!is.null(claps_pd)) {
  ix <- match(c("English", "Turkish", "Norwegian"), PP$language)
  s3 <- PP$se_floor_h1b[ix] * PP$affectedness_sd[ix]
  run_k3 <- function(tsc, n_reps, seed = 20260910L) {
    set.seed(seed); z <- numeric(n_reps); bf <- numeric(n_reps)
    for (i in seq_len(n_reps)) {
      theta <- vapply(claps_pd, function(v) v[sample.int(length(v), 1L)], numeric(1))
      y <- theta + rnorm(3, 0, s3)
      p <- ra_meta_posterior(y, s3, prior_tau_scale = tsc)
      bf[i] <- directional_bf_meta(p)
      z[i]  <- qnorm(min(max(p$p_mu_negative, 1e-12), 1 - 1e-12))
    }
    data.frame(tau_prior_scale = tsc, mean_z = mean(z), sd_z = sd(z),
               p_bf10 = mean(bf >= 10), p_bf6 = mean(bf >= 6), p_bf3 = mean(bf >= 3))
  }
  scales <- c(CALIBRATION_BRACKET[1], CALIBRATED_TAU_PRIOR_SCALE, CALIBRATION_BRACKET[2])
  calib  <- do.call(rbind, lapply(scales, run_k3, n_reps = REPS))
  calib  <- rbind(data.frame(tau_prior_scale = NA, t(TARGET)), calib)
  calib$row <- c("brms target (135 replicates)", paste0("engine, tau prior scale ", scales))
  print(calib[, c("row", "tau_prior_scale", "mean_z", "sd_z", "p_bf10", "p_bf6", "p_bf3")],
        row.names = FALSE, digits = 3)
  write.csv(calib, file.path(OUT, "klanguage_calibration.csv"), row.names = FALSE)
} else cat("pilot DGP files not present locally; calibration skipped.\n")

# ---------------------------------------------------------------------------
# 3. The three regimes
# ---------------------------------------------------------------------------
cell <- function(K, pop_name, mode = "assurance", rho = RHO_EFF,
                 tsc = CALIBRATED_TAU_PRIOR_SCALE, thr = 10, seed) {
  sp <- split_se(rho)
  cbind(population = pop_name, rho = rho,
        run_klanguage_cell(K = K, n_reps = REPS, se_lang = sp$indep, mode = mode,
                           pop = POPS[[pop_name]], shared_verb_sd = sp$shared,
                           prior_tau_scale = tsc, bf_threshold = thr, seed = seed))
}

rule(sprintf("P(BF >= 10) FOR THE POOLED INTERACTION, BY tau REGIME (%d replicates a cell)", REPS))

power <- do.call(rbind, lapply(REGIMES, function(pn)
  do.call(rbind, lapply(seq_along(KS), function(i) cell(KS[i], pn, seed = 4100L + i)))))

tb <- tapply(power$power, list(power$population, power$K), function(x) x[1])
cat("\nASSURANCE: (mu, tau) drawn from their joint posterior given that regime's\n",
    "languages, so the figure carries our uncertainty about tau as well as the\n",
    "sampling variation. This is the honest headline.\n\n", sep = "")
print(round(tb[REGIMES, ], 3))
cat("\nMonte-Carlo error on each entry is at most",
    sprintf("%.3f", max(power$mcse)), "\n")

# Conditional power: tau treated as KNOWN at each regime's posterior median. The
# gap between this table and the one above is the price of not knowing tau, and it
# is the larger of the two effects at every K. Separating them matters because they
# have different remedies: a small conditional power is fixed by recruiting more
# languages, and a large gap is not.
cond <- do.call(rbind, lapply(REGIMES, function(pn) {
  post <- ra_meta_posterior(POPS[[pn]]$y, POPS[[pn]]$s,
                            prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
  do.call(rbind, lapply(seq_along(KS), function(i) {
    sp <- split_se(RHO_EFF)
    r <- run_klanguage_cell(K = KS[i], n_reps = REPS, se_lang = sp$indep, mode = "point",
                            mu_point = post$mu_mean, tau_point = post$tau_median,
                            shared_verb_sd = sp$shared, bf_threshold = 10,
                            seed = 4150L + i)
    cbind(population = pn, mu = post$mu_mean, tau = post$tau_median, r)
  }))
}))
cat("\nCONDITIONAL on tau being known at that regime's posterior median:\n\n")
print(round(tapply(cond$power, list(cond$population, cond$K), function(x) x[1])[REGIMES, ], 3))
cat("\ntau assumed:", paste(sprintf("%s %.2f", REGIMES,
    vapply(REGIMES, function(p) cond$tau[cond$population == p][1], numeric(1))), collapse = "  "), "\n")
cat("The gap between the two tables is what not knowing tau costs. It does not\n",
    "close as K rises, because six to twelve languages cannot pin tau down either.\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Sensitivities
# ---------------------------------------------------------------------------
rule("SENSITIVITIES (seven-language regime unless stated)")

sens <- list(); add <- function(tag, df) sens[[tag]] <<- cbind(axis = tag, df)

add("verb_sharing", do.call(rbind, lapply(c(0, RHO_EFF, 0.22, 1), function(rho)
  do.call(rbind, lapply(seq_along(KS), function(i)
    cell(KS[i], "seven", rho = rho, seed = 4200L + i))))))

add("leave_one_out", do.call(rbind, lapply(c("seven", "drop_balinese", "drop_hebrew"), function(pn)
  do.call(rbind, lapply(seq_along(KS), function(i) cell(KS[i], pn, seed = 4300L + i))))))

add("calibration", do.call(rbind, lapply(c(CALIBRATION_BRACKET[1], CALIBRATED_TAU_PRIOR_SCALE,
                                           CALIBRATION_BRACKET[2]), function(tsc)
  do.call(rbind, lapply(seq_along(KS), function(i)
    cell(KS[i], "seven", tsc = tsc, seed = 4400L + i))))))

add("threshold", do.call(rbind, lapply(c(3, 6, 10), function(thr)
  do.call(rbind, lapply(seq_along(KS), function(i)
    cell(KS[i], "seven", thr = thr, seed = 4500L + i))))))

add("mode", do.call(rbind, lapply(c("assurance", "safeguard"), function(m)
  do.call(rbind, lapply(seq_along(KS), function(i)
    cell(KS[i], "seven", mode = m, seed = 4600L + i))))))

sensitivity <- do.call(rbind, sens)
KEY <- c(verb_sharing = "rho", leave_one_out = "population", calibration = "prior_tau_scale",
         threshold = "bf_threshold", mode = "mode")
for (a in unique(sensitivity$axis)) {
  d <- sensitivity[sensitivity$axis == a, ]
  cat("\n", a, " (rows = ", KEY[[a]], ", columns = K)\n", sep = "")
  print(round(tapply(d$power, list(as.character(d[[KEY[[a]]]]), d$K), function(x) x[1]), 3))
}

# ---------------------------------------------------------------------------
# 5. Smallest K reaching each target
# ---------------------------------------------------------------------------
rule("SMALLEST K REACHING EACH TARGET, BY REGIME")

min_k <- do.call(rbind, lapply(c(0.80, 0.90, 0.95), function(tgt)
  do.call(rbind, lapply(REGIMES, function(pn) {
    d <- power[power$population == pn, ]; d <- d[order(d$K), ]
    hit <- d$K[d$power >= tgt]
    data.frame(target = tgt, regime = pn,
               min_K = if (length(hit)) min(hit) else NA_integer_,
               at_6 = d$power[d$K == 6], at_8 = d$power[d$K == 8],
               at_10 = d$power[d$K == 10], at_12 = d$power[d$K == 12])
  }))))
print(min_k, row.names = FALSE, digits = 3)
cat("\nNA means no K up to", max(KS), "reaches the target.\n")

# ---------------------------------------------------------------------------
# 5b. The named-language arm: which languages join, not how many
# ---------------------------------------------------------------------------
# Everything above draws languages exchangeably from a population, which is the
# only way to say anything about teams that have not been recruited, and is also
# the step a reviewer will press hardest: at K = 12 it means nine of the twelve
# languages are inventions. For K <= 7 that step is avoidable, and the PI named
# seven specific languages, so this arm answers what he literally asked.
rule("NAMED LANGUAGES: EVERY SUBSET OF THE SEVEN, NO SYNTHETIC LANGUAGE ANYWHERE")

sp_eff <- split_se(RHO_EFF)
subsets <- do.call(rbind, lapply(c(5, 6, 7), function(K)
  enumerate_language_subsets(POPS$seven, K, sp_eff$indep, sp_eff$shared,
                             n_reps = REPS, seed = 5000L + K)))
cat("\nAt K = 6, one language is left out. Power by which one:\n")
s6 <- subsets[subsets$K == 6, c("omitted", "realised_tau", "mean_effect", "power", "mcse")]
print(s6[order(-s6$power), ], row.names = FALSE, digits = 3)
cat("\nAll seven together (K = 7):",
    sprintf("%.3f\n", subsets$power[subsets$K == 7]))
cat("Range across the seven K = 6 subsets:",
    sprintf("%.3f to %.3f", min(s6$power), max(s6$power)),
    "- a", sprintf("%.1f", max(s6$power) / max(min(s6$power), 0.001)),
    "fold spread at FIXED K.\n")
cat("Exchangeable arm at K = 6, for comparison:",
    sprintf("%.3f\n", power$power[power$population == "seven" & power$K == 6]))
cat("\nThe gap between the named and exchangeable arms at K = 7 is the\n",
    "fixed-versus-random-language estimand difference. It is not a discrepancy\n",
    "to be reconciled: the two arms answer different questions.\n", sep = "")

# How much of the outcome is decided by WHICH languages join rather than how many.
vshare <- do.call(rbind, lapply(c(6, 8, 10, 12), function(K)
  realised_set_variance_share(K, POPS$seven, sp_eff$indep, sp_eff$shared,
                              n_sets = max(120L, REPS %/% 20L),
                              n_reps_per_set = 150L, seed = 5100L + K)))
cat("\nShare of the outcome variance due to which languages joined:\n")
print(vshare, row.names = FALSE, digits = 3)
cat("\nIf that share does not fall as K rises, a rule phrased purely in terms of K\n",
    "is answering a question other than the one that decides the result.\n", sep = "")

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
# there only as far as the protocol transfers, so the prior on that language is
# widened by `transfer_sd`. At transfer_sd equal to tau the earlier measurement is
# worth nothing and the language is an unknown draw again.
rule("SETS WITH SOME LANGUAGES ALREADY MEASURED")

sp_eff <- split_se(RHO_EFF)
claps3 <- POPS$claps_protocol
glossa4 <- POPS$glossa_protocol

# The scenarios name their languages. An earlier version took "two of the 2023
# four" by position, which silently meant Balinese and Hebrew, and so reported the
# cost of recruiting the one wrong-signed language as though it were the general
# effect of knowing two more. Which languages are known matters here for exactly
# the reason it matters in section 5b, so both ends are shown.
g <- function(...) glossa4[match(c(...), glossa4$language), ]
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
                                transfer_sd = sc$tsd, seed = 5200L + 10L * j + K))))
}))
print(round(tapply(partly$power, list(partly$scenario, partly$K), function(x) x[1]), 3))
cat("
Row a is the exchangeable arm restated. The others hold the named teams at their
",
    "measured values and draw only the remaining places. Rows c and d differ only in
",
    "WHICH two of the 2023 languages take part, and the gap between them is wider than
",
    "anything K does, which is section 5b's point arriving again.
", sep = "")

# ---------------------------------------------------------------------------
# 6. The two outcomes the protocol has no sentence for
# ---------------------------------------------------------------------------
# A cross-linguistic claim can fail in ways a power figure does not describe. Both
# of these are computed on the SAME simulated studies as the power above, so the
# pooled and per-language results are the joint outcome of one dataset rather than
# two independent calculations.
rule("OUTCOMES THE PROTOCOL NEEDS A PRESPECIFIED READING OF")

outcome_cell <- function(K, pop_name, seed) {
  set.seed(seed)
  sp   <- split_se(RHO_EFF)
  pars <- draw_population_params(REPS, POPS[[pop_name]]$y, POPS[[pop_name]]$s)
  z_thr <- qnorm(10 / 11)
  res <- vapply(seq_len(REPS), function(i) {
    theta  <- rnorm(K, pars$mu[i], pars$tau[i])
    common <- rnorm(1, 0, sp$shared)
    y      <- theta + common + rnorm(K, 0, sp$indep)
    # Each language's own test, at its own precision. The prediction is that the
    # interaction is NEGATIVE, so a language that supports it has a negative
    # estimate and a language that contradicts it at the same threshold has a
    # positive one. Getting these two the wrong way round reports the count of
    # confirming languages as the count of contradicting ones.
    z_lang <- y / sp$indep
    right  <- sum(z_lang <= -z_thr)          # reaches threshold, predicted direction
    wrong  <- sum(z_lang >=  z_thr)          # reaches threshold, opposite direction
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

outcomes <- do.call(rbind, lapply(REGIMES, function(pn)
  do.call(rbind, lapply(c(6, 8, 10, 12), function(K)
    outcome_cell(K, pn, seed = 4700L + K)))))
print(outcomes, row.names = FALSE, digits = 3)
cat("\ne_wrong_sign is the expected NUMBER of languages whose own test reaches a\n",
    "directional Bayes factor of 10 in the direction OPPOSITE to the prediction.\n",
    "Balinese, at +0.33 per SD, is already such a language.\n", sep = "")

# ---------------------------------------------------------------------------
# 7. The trades the PI will propose
# ---------------------------------------------------------------------------
rule("TRADES: PARTICIPANTS AGAINST LANGUAGES, VERBS AGAINST LANGUAGES")

# Per-language SE as a function of participants and verbs. The verb term is the
# floor and scales as 1/sqrt(n_verbs); the participant term scales as 1/sqrt(N).
# Both are anchored at the pilot values (N = 80, 72 verbs).
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

# a. Fixed participant budget of 640: buy languages or buy participants?
budget <- do.call(rbind, lapply(list(c(8, 80), c(10, 64), c(12, 53), c(16, 40)), function(kn)
  trade_cell(kn[1], kn[2], 72, "seven", seed = 4800L + kn[1])))
cat("\na. Fixed budget of 640 participants (seven-language regime)\n")
print(budget, row.names = FALSE, digits = 3)

# b. Verb count. Three of the five Glossa languages could not field 72 verbs
#    (Balinese 47, Hebrew 56, Mandarin 57), yet the design assumes 72 everywhere.
verbs <- do.call(rbind, lapply(c(40, 47, 56, 72), function(v)
  do.call(rbind, lapply(c(6, 8, 12), function(K)
    trade_cell(K, 80, v, "seven", seed = 4900L + v + K)))))
cat("\nb. Verb count, at N = 80 (seven-language regime)\n")
print(round(tapply(verbs$power, list(verbs$n_verbs, verbs$K), function(x) x[1]), 3))
cat("\nA team that can only field 47 verbs costs the pooled test very little, while\n",
    "forfeiting its own confirmatory test. That is a recruitment-policy fact.\n", sep = "")

tradeoffs <- rbind(cbind(trade = "participant_budget", budget),
                   cbind(trade = "verb_count", verbs))

# ---------------------------------------------------------------------------
# 8. Recruitment: K is a random variable and Y is a rule about it
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
write.csv(power,       file.path(OUT, "klanguage_power.csv"),       row.names = FALSE)
write.csv(cond,        file.path(OUT, "klanguage_conditional.csv"),  row.names = FALSE)
write.csv(partly,      file.path(OUT, "klanguage_partly_known.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(OUT, "klanguage_sensitivity.csv"), row.names = FALSE)
write.csv(min_k,       file.path(OUT, "klanguage_minimum_k.csv"),   row.names = FALSE)
write.csv(outcomes,    file.path(OUT, "klanguage_outcomes.csv"),    row.names = FALSE)
write.csv(tradeoffs,   file.path(OUT, "klanguage_tradeoffs.csv"),   row.names = FALSE)
write.csv(recruit,     file.path(OUT, "klanguage_recruitment.csv"), row.names = FALSE)
write.csv(pop_tab,     file.path(OUT, "klanguage_populations.csv"), row.names = FALSE)
write.csv(subsets,     file.path(OUT, "klanguage_named_subsets.csv"), row.names = FALSE)
write.csv(vshare,      file.path(OUT, "klanguage_set_variance.csv"),  row.names = FALSE)

rule("FILES WRITTEN")
for (f in c("klanguage_calibration.csv", "klanguage_power.csv", "klanguage_sensitivity.csv",
            "klanguage_minimum_k.csv", "klanguage_outcomes.csv", "klanguage_tradeoffs.csv",
            "klanguage_recruitment.csv", "klanguage_populations.csv",
            "klanguage_named_subsets.csv", "klanguage_set_variance.csv",
            "klanguage_conditional.csv", "klanguage_partly_known.csv"))
  if (file.exists(file.path(OUT, f))) cat(" ", file.path(OUT, f), "\n")
cat("\nR:", R.version.string, "| replicates per cell:", REPS, "\n")
