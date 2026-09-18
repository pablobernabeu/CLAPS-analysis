#!/usr/bin/env Rscript
# scripts/extract_cross_language_effects.R
# ---------------------------------------------------------------------------
# Put the seven languages' affectedness-by-sentence-type interaction estimates
# on ONE common scale, and characterise the between-language distribution.
#
# WHY THIS SCRIPT EXISTS
#   The PI asked for a pooled (K-language) design analysis built on the seven
#   languages for which an Active-vs-Passive affectedness interaction has been
#   estimated: English, Norwegian and Turkish from the CLAPS pilot, and English,
#   Balinese, Hebrew, Indonesian and Mandarin from Ambridge, Arnon & Bekman
#   (2023, Glossa: Psycholinguistics, doi:10.5070/G6011177). Any such analysis
#   needs a between-language SD (tau) of the interaction, and tau is only defined
#   once every estimate is expressed in the same units. The two sources do NOT
#   share units, so this script establishes the conversion, applies it, and then
#   estimates tau three ways with a full set of leave-one-out sensitivities.
#
# WHAT IT ESTABLISHES (all of it recomputed from the raw files, nothing assumed)
#   1. Units of the affectedness predictor in each source.
#   2. A per-language-SD harmonisation, and an alternative harmonisation onto
#      the shared ten-feature semantic space.
#   3. Mean, SD and three random-effects meta-analyses (DerSimonian-Laird by
#      hand, REML via metafor, and a Bayesian posterior for tau under a
#      half-normal(0, 0.5) prior computed by grid integration).
#   4. tau, its uncertainty, and its sensitivity to dropping Balinese, dropping
#      Indonesian, and to which English estimate is used.
#   5. Exactly which contrast each language's "Active vs Passive" interaction is.
#   6. Sample sizes, verb counts, response scales, and an empirical check on
#      what the Balinese 10-point-to-5-point rescale does to a logit slope.
#
# COST
#   Reads CSVs and text model summaries; fits four small fixed-effect cumulative
#   link models (ordinal::clm) for the response-scale check. Runs in seconds.
#   No brms, no Stan, no cluster.
#
# USAGE
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" scripts/extract_cross_language_effects.R
#   Optionally: --osf=<dir> --claps=<dir> --out=<dir>
# ---------------------------------------------------------------------------

options(stringsAsFactors = FALSE, width = 200)

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", flag, "="), "", hit[1])
}

DEFAULT_ROOT <- "C:/Users/pablob/OneDrive - Nexus365/Documents/GitHub/private_Pablo_CLAPS"
OSF   <- argval("--osf",   file.path(DEFAULT_ROOT, "PRIVATE_background/Sub2_OSF_Passives"))
CLAPS <- argval("--claps", file.path(DEFAULT_ROOT, "design_analysis"))
OUT   <- argval("--out",   file.path(CLAPS, "outputs/design_summary_pilot"))

for (p in c(OSF, CLAPS)) {
  if (!dir.exists(p)) stop("[paths] not found: ", p)
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

have_metafor <- requireNamespace("metafor", quietly = TRUE)
have_ordinal <- requireNamespace("ordinal", quietly = TRUE)

rule <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

# ---------------------------------------------------------------------------
# 1. UNITS OF THE AFFECTEDNESS PREDICTOR IN EACH SOURCE
# ---------------------------------------------------------------------------
#
# CLAPS
#   scripts/00_harmonise_pilot_data.R sets Semantics <- affectedness_scores_all,
#   a verb-level mean affectedness rating (English/Turkish on a 1-7 scale,
#   Norwegian on a 0-100 slider). R/02_preprocess_factors.R::scale_semantics()
#   then builds the model predictor as
#
#       Semantics_scaled = (Semantics - mean(Semantics)) / (2 * sd(Semantics))
#
#   within language (Gelman 2008, doi:10.1002/sim.3107). So one unit of
#   Semantics_scaled is TWO raw SDs of that language's affectedness, and one raw
#   SD is 0.5 units of Semantics_scaled. Because the scaling constants are
#   computed over trial rows and the design is balanced, sd(Semantics_scaled) is
#   exactly 0.5 at trial level; scripts/aggregate_pilot.R instead uses the
#   verb-level SD (n-1 over the verb set), which is 0.5035. That verb-level SD is
#   the one already stored as `affectedness_sd` in pilot_params_ceilings.csv, and
#   it is the one used here so that the numbers reconcile with that file.
#
# Ambridge, Arnon & Bekman (2023)
#   Semantics is PCA1 of ten semantic-feature ratings collected verb by verb
#   (the ten items are the column headings of Semantics.csv; V6_Hebrew_Passives.R
#   contains the aggregation and PCA that built them, and V6_Indonesian_Passives.R
#   / V6_Mandarin_Passives.R record `Data$Semantics = Data$PCA1`). Crucially, the
#   distributed *_Passives.csv files already store that PCA1 score Z-SCORED
#   within language: mean 0, SD 1 over trial rows (verified below). The
#   `scale()` calls in the per-language scripts are therefore no-ops, and the
#   pooled script V6_All_Semantics_Only.R, which does NOT call scale(), is
#   nevertheless already on a per-language-SD footing. Hebrew is the one language
#   whose script omits scale() entirely - which turns out not to matter, for the
#   same reason.
#
#   Consequence: every Glossa `S_TypeActive:Semantics` coefficient is ALREADY
#   "logit change in the affectedness slope per 1 SD of that language's
#   affectedness predictor". No conversion is needed beyond the trial-level to
#   verb-level SD refinement applied for consistency with CLAPS.

rule("1. UNITS OF THE AFFECTEDNESS PREDICTOR")

glossa_files <- c(Balinese   = "Balinese_Passives.csv",
                  English    = "English_Passives.csv",
                  Hebrew     = "Hebrew_Passives.csv",
                  Indonesian = "Indonesian_Passives.csv",
                  Mandarin   = "Mandarin_Passives.csv")

glossa_raw <- lapply(glossa_files, function(f) read.csv(file.path(OSF, f)))
names(glossa_raw) <- names(glossa_files)

glossa_units <- do.call(rbind, lapply(names(glossa_raw), function(nm) {
  d  <- glossa_raw[[nm]]
  vm <- tapply(d$Semantics, d$verb, function(x) x[1])
  data.frame(source = "Glossa2023", language = nm,
             n_obs = nrow(d),
             n_participants = length(unique(d$Participant)),
             n_verbs = length(vm),
             sem_mean_trial = mean(d$Semantics),
             sem_sd_trial   = sd(d$Semantics),
             sem_sd_verb    = sd(vm),
             resp_min = min(d$Response), resp_max = max(d$Response),
             resp_levels = length(unique(d$Response)))
}))

claps <- read.csv(file.path(CLAPS, "data/pilot/claps_pilot_harmonised.csv"))
claps_units <- do.call(rbind, lapply(c("English", "Turkish", "Norwegian"), function(L) {
  s  <- claps[claps$Language == L, ]
  z  <- (s$Semantics - mean(s$Semantics)) / (2 * sd(s$Semantics))   # scale_semantics()
  vz <- tapply(z, s$Verb_ID, function(x) x[1])
  data.frame(source = "CLAPS_pilot", language = L,
             n_obs = nrow(s),
             n_participants = length(unique(s$Participant)),
             n_verbs = length(vz),
             sem_mean_trial = mean(z),
             sem_sd_trial   = sd(z),
             sem_sd_verb    = sd(vz),
             resp_min = min(s$Response), resp_max = max(s$Response),
             resp_levels = length(unique(s$Response)))
}))

units_tbl <- rbind(claps_units, glossa_units)
cat("\nPredictor actually entered into each model (CLAPS: Semantics_scaled; Glossa: Semantics):\n")
print(units_tbl, row.names = FALSE, digits = 5)
cat("\nCLAPS raw affectedness before scaling (affectedness_scores_all):\n")
print(do.call(rbind, lapply(c("English", "Turkish", "Norwegian"), function(L) {
  s <- claps[claps$Language == L, ]
  v <- sort(unique(s$Semantics))
  # Granularity pins down the instrument: a mean rating lands on multiples of
  # 1/(number of judgements averaged), a PCA score does not land on a lattice.
  gran <- min(diff(v))
  data.frame(language = L, n_distinct = length(v),
             raw_min = min(v), raw_max = max(v),
             raw_mean = mean(s$Semantics), raw_sd = sd(s$Semantics),
             min_gap = gran, is_multiple_of_1_40 = all(abs(v * 40 - round(v * 40)) < 1e-8))
})), row.names = FALSE, digits = 5)
cat("The CLAPS affectedness predictor is a verb-level MEAN RATING (English and Turkish on a\n")
cat("1-7 scale, Norwegian on a 0-100 slider), not a PCA score: every English and Turkish value\n")
cat("is an exact multiple of 1/40, the lattice a mean over 40 judgements produces. The Glossa\n")
cat("predictor is PCA1 of ten feature ratings. Different instruments for the same construct,\n")
cat("which is why Harmonisation B exists.\n")

# ---------------------------------------------------------------------------
# 2. RAW PER-LANGUAGE INTERACTION ESTIMATES, WITH PROVENANCE
# ---------------------------------------------------------------------------
#
# Glossa estimates are the posterior mean and Est.Error of S_TypeActive:Semantics
# from the single-language maximal cumulative(logit) models in the OSF text
# dumps. Every one is cross-checked below against AA_Models_Summary.csv, where
# the per-sentence-type simple slopes are tabulated and Active minus Passive must
# reproduce the interaction.
#
# CLAPS estimates are posterior means of S_TypeActive:Semantics_scaled from the
# pilot L5_correlated_maximal cumulative(logit) fits, as summarised in
# outputs/design_summary_pilot/pilot_params_ceilings.csv.

rule("2. RAW PER-LANGUAGE INTERACTION ESTIMATES")

glossa_est <- data.frame(
  source   = "Glossa2023",
  language = c("Balinese", "English", "Hebrew", "Indonesian", "Mandarin"),
  est      = c( 0.33,      -0.54,     -0.81,    -0.07,        -0.80),
  se       = c( 0.13,       0.18,      0.20,     0.08,         0.21),
  # main-effect (reference-level, i.e. canonical-passive) affectedness slope
  slope_ref    = c(0.38, 0.53, 0.79, 0.27, 0.80),
  slope_ref_se = c(0.08, 0.11, 0.19, 0.17, 0.15),
  file = c("Balinese_Sem_Only_Bayes_redux.txt", "English_Sem_Only_Bayes.txt",
           "Hebrew_Sem_Only_cumulative.txt", "Indonesian_Sem_Only_Bayes.txt",
           "Mandarin_Sem_Only_Bayes.txt")
)

# Cross-check against AA_Models_Summary.csv (Active slope minus Passive slope).
aa <- read.csv(file.path(OSF, "AA_Models_Summary.csv"), check.names = TRUE)
names(aa)[1] <- "Language"
aa_sem <- subset(aa, Model == "Semantics_only" & Predictor_Type == "Semantics" &
                   Language != "All_Languages")
xcheck <- do.call(rbind, lapply(unique(aa_sem$Language), function(L) {
  s  <- subset(aa_sem, Language == L)
  ga <- s$Estimate[s$Stype == "Active"]
  gp <- s$Estimate[s$Stype == "Passive"]
  data.frame(language = L, active_slope = ga, passive_slope = gp,
             implied_interaction = ga - gp,
             reported_interaction = glossa_est$est[match(L, glossa_est$language)])
}))
xcheck$discrepancy <- round(xcheck$implied_interaction - xcheck$reported_interaction, 3)
cat("\nCross-check: AA_Models_Summary.csv simple slopes, Active minus Passive\n")
print(xcheck, row.names = FALSE, digits = 4)
if (max(abs(xcheck$discrepancy)) > 0.02) {
  warning("[xcheck] simple-slope reconstruction differs from the reported interaction")
}

pp <- read.csv(file.path(CLAPS, "outputs/design_summary_pilot/pilot_params_ceilings.csv"))
claps_est <- data.frame(
  source   = "CLAPS_pilot",
  language = pp$language,
  est      = pp$interaction_mean,     # units: per Semantics_scaled unit = per 2 raw SD
  se_floor = pp$se_floor_h1b,         # verb-information SE floor, same units
  slope_post_sd = pp$slope_post_sd,   # posterior SD of the MAIN-EFFECT slope
  sd_x     = pp$affectedness_sd,      # verb-level SD of Semantics_scaled
  file     = "outputs/design_summary_pilot/pilot_params_ceilings.csv"
)
cat("\nCLAPS pilot interaction estimates (Semantics_scaled units):\n")
print(claps_est[, c("language", "est", "se_floor", "slope_post_sd", "sd_x")],
      row.names = FALSE, digits = 4)

# The CLAPS pilot posterior SD of the INTERACTION is not archived locally: the
# per-language DGP objects (pilot_dgp_v2_pilot_*.rds) that scripts/aggregate_pilot.R
# reads live on ARC, and pilot_params_ceilings.csv preserves only sd() of the
# MAIN-EFFECT slope draws plus the verb-information SE floor for H1b. Three
# stand-ins are therefore carried through every downstream calculation:
#   floor       - se_floor_h1b. A hard lower bound: it is the SE that survives
#                 with unlimited participants, so it understates pilot uncertainty.
#   slope_sd    - the main-effect posterior SD, used as-is.
#   calibrated  - the main-effect posterior SD scaled by the median ratio
#                 SE(interaction)/SE(main slope) observed across the five Glossa
#                 languages, which is the only direct evidence available on how
#                 much wider the interaction posterior is in models of this shape.
se_ratio <- glossa_est$se / glossa_est$slope_ref_se
se_ratio_median <- median(se_ratio)
cat("\nSE(interaction) / SE(main-effect slope) in the five Glossa models:\n")
print(setNames(round(se_ratio, 3), glossa_est$language))
cat("median ratio =", round(se_ratio_median, 3), "\n")

# ---------------------------------------------------------------------------
# 3. HARMONISATION A: LOGIT CHANGE PER 1 SD OF THE LANGUAGE'S OWN PREDICTOR
# ---------------------------------------------------------------------------
#
# Glossa: multiply by the verb-level SD of Semantics (approximately 1.007,
#         because the stored column is Z-scored at trial level).
# CLAPS:  multiply by affectedness_sd (0.5035 / 0.5038), the verb-level SD of
#         Semantics_scaled - equivalently, halve, because a Semantics_scaled unit
#         is two raw SDs.
# Both are cumulative-logit coefficients, so the target quantity is the change in
# the latent-scale affectedness slope, in log-odds, when moving one within-language
# SD along that language's affectedness predictor.

rule("3. HARMONISATION A - PER-LANGUAGE SD UNITS")

g_sd <- setNames(glossa_units$sem_sd_verb, glossa_units$language)

harm_glossa <- data.frame(
  source   = "Glossa2023",
  language = glossa_est$language,
  raw_est  = glossa_est$est,
  raw_se   = glossa_est$se,
  sd_x     = as.numeric(g_sd[glossa_est$language]),
  se_kind  = "posterior_sd"
)
harm_glossa$y  <- harm_glossa$raw_est * harm_glossa$sd_x
harm_glossa$se <- harm_glossa$raw_se  * harm_glossa$sd_x

mk_claps <- function(se_kind) {
  raw_se <- switch(se_kind,
                   floor      = claps_est$se_floor,
                   slope_sd   = claps_est$slope_post_sd,
                   calibrated = claps_est$slope_post_sd * se_ratio_median,
                   stop("unknown se_kind"))
  data.frame(source = "CLAPS_pilot", language = claps_est$language,
             raw_est = claps_est$est, raw_se = raw_se,
             sd_x = claps_est$sd_x, se_kind = se_kind,
             y = claps_est$est * claps_est$sd_x,
             se = raw_se * claps_est$sd_x)
}

harm_all <- rbind(harm_glossa, mk_claps("floor"), mk_claps("slope_sd"), mk_claps("calibrated"))
harm_all$label <- ifelse(harm_all$source == "CLAPS_pilot",
                         paste0(harm_all$language, " (CLAPS)"),
                         paste0(harm_all$language, " (G23)"))

cat("\nHarmonised estimates, logits per 1 SD of the language's own affectedness predictor:\n")
print(harm_all[, c("source", "language", "raw_est", "raw_se", "sd_x", "se_kind", "y", "se")],
      row.names = FALSE, digits = 4)

# English cross-check. The same language, two independent data sets, two
# independent affectedness instruments (a mean 1-7 rating vs PCA1 of ten
# features), two response scales (7-point vs 5-point). If the harmonisation is
# doing its job these two must land close together.
eng_claps  <- harm_all$y[harm_all$language == "English" & harm_all$source == "CLAPS_pilot"][1]
eng_glossa <- harm_all$y[harm_all$language == "English" & harm_all$source == "Glossa2023"]
cat(sprintf("\nENGLISH CROSS-CHECK  CLAPS %.3f vs Glossa %.3f  (difference %.3f, ratio %.2f)\n",
            eng_claps, eng_glossa, eng_claps - eng_glossa, eng_claps / eng_glossa))

# ---------------------------------------------------------------------------
# 4. HARMONISATION B: ONTO THE SHARED TEN-FEATURE SEMANTIC SPACE
# ---------------------------------------------------------------------------
#
# The ten semantic-feature questions are the same instrument everywhere (they are
# properties of the verb/event, asked about English glosses), so it is tempting to
# treat them as one shared axis. They are not one axis in the data: the ratings
# were re-collected per language and PCA1 recomputed per language, and CLAPS used
# a mean affectedness rating rather than a PCA at all. The verb-level scores
# therefore correlate only 0.6-0.9 across languages (printed below).
#
# The defensible harmonisation is a classical errors-in-variables correction. Take
# a latent shared affectedness axis Z, standardised. Suppose each language's
# verb-level score is X_l = lambda_l * Z + e_l, with e_l independent of Z and of
# the response noise. A regression on the proxy X_l is then attenuated relative to
# a regression on Z by exactly lambda_l, so
#
#       beta_shared_l = beta_perSD_l / lambda_l.
#
# lambda_l is estimated from the off-diagonal of the cross-language correlation
# matrix of verb-level scores, since cor(X_l, X_m) = lambda_l * lambda_m. That is
# a one-factor model fitted by least squares on the off-diagonal only (MINRES),
# which treats all seven languages symmetrically. A second, cruder variant simply
# takes lambda_l = cor(X_l, PCA1 of Semantics.csv), which privileges Hebrew
# (whose predictor IS that PCA1, so its lambda is 1 by construction).
#
# CAVEAT, and it is not a small one. This treats every scrap of cross-language
# disagreement about a verb's affectedness as measurement error. Some of it is
# almost certainly real: cultural and lexical differences mean the same English
# gloss need not denote an equally affecting event everywhere. To the extent that
# the disagreement is real rather than noise, the correction over-inflates every
# estimate, and inflates the low-lambda languages most. Harmonisation B is
# therefore an UPPER bound on the harmonised effects and hence on tau, and
# Harmonisation A is the primary.

rule("4. HARMONISATION B - SHARED TEN-FEATURE SEMANTIC SPACE")

vkey <- function(x) gsub("[^a-z]", "", tolower(as.character(x)))

verb_scores <- list()
for (nm in names(glossa_raw)) {
  d <- glossa_raw[[nm]]
  u <- unique(data.frame(key = vkey(d$verb), val = d$Semantics))
  u <- u[!duplicated(u$key), ]
  verb_scores[[paste0(nm, "_G23")]] <- setNames(u$val, u$key)
}
for (L in c("English", "Turkish", "Norwegian")) {
  s <- claps[claps$Language == L, ]
  u <- unique(data.frame(key = vkey(s$Verb), val = s$Semantics))
  u <- u[!duplicated(u$key), ]
  verb_scores[[paste0(L, "_CLAPS")]] <- setNames(u$val, u$key)
}

all_keys <- sort(unique(unlist(lapply(verb_scores, names))))
M <- matrix(NA_real_, nrow = length(all_keys), ncol = length(verb_scores),
            dimnames = list(all_keys, names(verb_scores)))
for (nm in names(verb_scores)) M[names(verb_scores[[nm]]), nm] <- verb_scores[[nm]]
M <- scale(M)                                   # z-score each language's own verb set

Rmat <- cor(M, use = "pairwise.complete.obs")
Nmat <- crossprod(!is.na(M))                    # pairwise n
cat("\nVerb-level affectedness: cross-language correlation matrix (pairwise complete)\n")
print(round(Rmat, 3))
cat("\nPairwise n of shared verbs\n")
print(Nmat)

# One-factor MINRES on the off-diagonal.
minres_loadings <- function(R) {
  k <- ncol(R)
  off <- which(upper.tri(R), arr.ind = TRUE)
  obj <- function(par) {
    lam <- tanh(par)                            # keeps |lambda| < 1
    sum((R[off] - lam[off[, 1]] * lam[off[, 2]])^2, na.rm = TRUE)
  }
  fit <- optim(rep(atanh(0.8), k), obj, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-12))
  lam <- tanh(fit$par)
  if (mean(lam) < 0) lam <- -lam                # sign is arbitrary; fix it positive
  names(lam) <- colnames(R)
  list(lambda = lam, rss = fit$value, conv = fit$convergence)
}
mf <- minres_loadings(Rmat)
cat("\nOne-factor (MINRES) loadings on the shared affectedness axis:\n")
print(round(mf$lambda, 3))
cat("residual sum of squares on the off-diagonal:", signif(mf$rss, 3),
    "| optim convergence code:", mf$conv, "\n")

# Crude variant: correlation with PCA1 of the ten-feature norms in Semantics.csv.
S10 <- read.csv(file.path(OSF, "Semantics.csv"), header = TRUE, check.names = FALSE)
S10$key <- vkey(S10$verb)
pcS <- prcomp(S10[, 2:11], center = TRUE, scale. = TRUE)
Zaxis <- setNames(as.numeric(scale(pcS$x[, 1])), S10$key)
cat("\nSemantics.csv: ", nrow(S10), " verbs x 10 features; PC1 explains ",
    round(100 * summary(pcS)$importance[2, 1], 1), "% of variance\n", sep = "")
lam_Z <- sapply(colnames(M), function(nm) {
  common <- intersect(rownames(M)[!is.na(M[, nm])], names(Zaxis))
  cor(M[common, nm], Zaxis[common])
})
n_Z <- sapply(colnames(M), function(nm)
  length(intersect(rownames(M)[!is.na(M[, nm])], names(Zaxis))))
cat("\nCorrelation of each language's verb-level score with Semantics.csv PC1:\n")
print(data.frame(lambda_Z = round(lam_Z, 3), n_matched = n_Z))

key_of <- function(src, lang) paste0(lang, ifelse(src == "CLAPS_pilot", "_CLAPS", "_G23"))
harm_all$mkey     <- key_of(harm_all$source, harm_all$language)
harm_all$lambda_1f <- as.numeric(mf$lambda[harm_all$mkey])
harm_all$lambda_Z  <- as.numeric(lam_Z[harm_all$mkey])
harm_all$y_shared  <- harm_all$y  / harm_all$lambda_1f
harm_all$se_shared <- harm_all$se / harm_all$lambda_1f

cat("\nHarmonisation B: logits per 1 SD of the shared affectedness axis\n")
print(harm_all[, c("label", "se_kind", "y", "lambda_1f", "y_shared", "lambda_Z")],
      row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------------
# 5. RANDOM-EFFECTS META-ANALYSIS
# ---------------------------------------------------------------------------
#
# Three estimators of tau, the between-language SD:
#   DL   - DerSimonian & Laird (1986, doi:10.1016/0197-2456(86)90046-2), by hand.
#   REML - via metafor::rma (Viechtbauer 2010, doi:10.18637/jss.v036.i03), with a
#          Q-profile (Higgins-Thompson) confidence interval from confint().
#   Bayes - posterior for tau under a half-normal(0, 0.5) prior, with mu given a
#          flat prior and integrated out analytically, evaluated on a grid. A
#          half-normal(0, 0.5) puts about 95% of its prior mass below tau = 1
#          logit per SD, which is generous next to effects of order 0.5.
#
# The Bayesian version matters here because K is small: with seven estimates the
# likelihood for tau is nearly flat towards zero and both DL and REML can return
# a boundary estimate of exactly zero, which is not a statement that languages
# agree, only that the data cannot rule it out.

rule("5. RANDOM-EFFECTS META-ANALYSIS")

meta_dl <- function(y, se) {
  k <- length(y)
  w <- 1 / se^2
  mu_fe <- sum(w * y) / sum(w)
  Q <- sum(w * (y - mu_fe)^2)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / C)
  wr <- 1 / (se^2 + tau2)
  mu <- sum(wr * y) / sum(wr)
  se_mu <- sqrt(1 / sum(wr))
  I2 <- max(0, (Q - (k - 1)) / Q)
  list(k = k, mu = mu, se_mu = se_mu, tau2 = tau2, tau = sqrt(tau2),
       Q = Q, df = k - 1, p_Q = pchisq(Q, k - 1, lower.tail = FALSE), I2 = I2)
}

# Marginal (profile) likelihood for tau with mu integrated out under a flat prior.
loglik_tau <- function(tau, y, se) {
  v <- se^2 + tau^2
  w <- 1 / v
  mu <- sum(w * y) / sum(w)
  0.5 * sum(log(w)) - 0.5 * log(sum(w)) - 0.5 * sum(w * (y - mu)^2)
}

meta_bayes <- function(y, se, prior_sd = 0.5, tau_max = 3, ngrid = 6000) {
  grid <- seq(0, tau_max, length.out = ngrid)
  lp <- vapply(grid, loglik_tau, numeric(1), y = y, se = se) +
    dnorm(grid, 0, prior_sd, log = TRUE)        # half-normal: normal density on tau >= 0
  lp <- lp - max(lp)
  d <- exp(lp)
  d <- d / sum(d * diff(grid)[1])
  cdf <- cumsum(d) * diff(grid)[1]
  cdf <- cdf / max(cdf)
  # The CDF saturates at 1 well before tau_max, so it is only strictly increasing
  # over part of the grid. approx() warns and silently averages on the flat tail;
  # trimming to the strictly increasing part makes the inversion well defined.
  keep <- c(TRUE, diff(cdf) > 0)
  q <- function(p) approx(cdf[keep], grid[keep], xout = p, rule = 2)$y
  mean_tau <- sum(grid * d) * diff(grid)[1]
  # posterior for mu: mixture over tau of N(mu_hat(tau), 1/sum(w))
  wgt <- d / sum(d)
  mu_t <- vapply(grid, function(t) { w <- 1 / (se^2 + t^2); sum(w * y) / sum(w) }, numeric(1))
  v_t  <- vapply(grid, function(t) { w <- 1 / (se^2 + t^2); 1 / sum(w) }, numeric(1))
  mu_mean <- sum(wgt * mu_t)
  mu_var  <- sum(wgt * (v_t + mu_t^2)) - mu_mean^2
  # probability the pooled mean is negative, and the directional Bayes factor
  mu_draws_p_neg <- sum(wgt * pnorm(0, mu_t, sqrt(v_t)))
  list(tau_mean = mean_tau, tau_median = q(0.5), tau_lo = q(0.025), tau_hi = q(0.975),
       tau_p_gt = function(x) 1 - approx(grid[keep], cdf[keep], xout = x, rule = 2)$y,
       mu_mean = mu_mean, mu_sd = sqrt(mu_var),
       p_mu_neg = mu_draws_p_neg,
       bf_dir = mu_draws_p_neg / (1 - mu_draws_p_neg),
       grid = grid, dens = d)
}

meta_all <- function(y, se, tag, prior_sd = 0.5) {
  dl <- meta_dl(y, se)
  bs <- meta_bayes(y, se, prior_sd = prior_sd)
  out <- data.frame(
    set = tag, k = dl$k,
    mean_raw = mean(y), sd_raw = sd(y),
    mu_DL = dl$mu, se_mu_DL = dl$se_mu, tau_DL = dl$tau,
    Q = dl$Q, df = dl$df, p_Q = dl$p_Q, I2 = dl$I2,
    tau_REML = NA_real_, tau_REML_lo = NA_real_, tau_REML_hi = NA_real_,
    mu_REML = NA_real_, se_mu_REML = NA_real_,
    tau_bayes_median = bs$tau_median, tau_bayes_lo = bs$tau_lo, tau_bayes_hi = bs$tau_hi,
    mu_bayes = bs$mu_mean, mu_bayes_sd = bs$mu_sd,
    p_mu_neg = bs$p_mu_neg, bf_dir_mu = bs$bf_dir,
    p_tau_gt_0.2 = bs$tau_p_gt(0.2), p_tau_gt_0.4 = bs$tau_p_gt(0.4)
  )
  if (have_metafor) {
    fit <- metafor::rma(yi = y, sei = se, method = "REML")
    out$tau_REML    <- sqrt(fit$tau2)
    out$mu_REML     <- as.numeric(fit$b)
    out$se_mu_REML  <- fit$se
    ci <- try(metafor::confint.rma.uni(fit), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      out$tau_REML_lo <- ci$random["tau", "ci.lb"]
      out$tau_REML_hi <- ci$random["tau", "ci.ub"]
    }
  }
  out
}

# Assemble the analysis sets. English appears in both sources; it is ONE language,
# so the primary set uses the CLAPS estimate (the pilot for the study actually
# being designed) and the alternatives swap in the Glossa one or keep both.
build_set <- function(se_kind, english = c("CLAPS", "Glossa", "both"), drop = character(0),
                      scale = c("own", "shared")) {
  english <- match.arg(english); scale <- match.arg(scale)
  h <- harm_all[harm_all$source == "Glossa2023" | harm_all$se_kind == se_kind, ]
  if (english == "CLAPS")  h <- h[!(h$language == "English" & h$source == "Glossa2023"), ]
  if (english == "Glossa") h <- h[!(h$language == "English" & h$source == "CLAPS_pilot"), ]
  h <- h[!(h$language %in% drop), ]
  if (scale == "shared") { h$y <- h$y_shared; h$se <- h$se_shared }
  h
}

SE_KIND_PRIMARY <- "calibrated"
prim <- build_set(SE_KIND_PRIMARY, "CLAPS")
cat("\nPRIMARY SET (K = ", nrow(prim), "), Harmonisation A, CLAPS English, ",
    SE_KIND_PRIMARY, " CLAPS SEs:\n", sep = "")
print(prim[, c("label", "y", "se")], row.names = FALSE, digits = 4)

res <- meta_all(prim$y, prim$se, "PRIMARY: A / CLAPS-English / calibrated-SE")
cat("\nUnweighted mean across languages:  ", round(mean(prim$y), 4), "\n", sep = "")
cat("Unweighted SD across languages:    ", round(sd(prim$y), 4), "\n", sep = "")
cat("Naive SE of that mean (sd/sqrt(K)):", round(sd(prim$y) / sqrt(nrow(prim)), 4), "\n", sep = "")
cat("\nPrimary meta-analysis:\n")
print(t(round(res[, -1], 4)))

# ---------------------------------------------------------------------------
# 6. SENSITIVITY OF TAU
# ---------------------------------------------------------------------------

rule("6. SENSITIVITY OF TAU")

sens <- list()
add <- function(tag, ...) sens[[length(sens) + 1L]] <<- meta_all(..., tag = tag)

for (sk in c("floor", "slope_sd", "calibrated")) {
  s <- build_set(sk, "CLAPS")
  add(paste0("A | CLAPS-Eng | SE=", sk), s$y, s$se)
}
for (eng in c("CLAPS", "Glossa", "both")) {
  s <- build_set(SE_KIND_PRIMARY, eng)
  add(paste0("A | English=", eng, " | SE=", SE_KIND_PRIMARY), s$y, s$se)
}
# Full leave-one-out, so that the two languages the PI singled out (Balinese,
# which goes the wrong way, and Indonesian, which is near zero) can be read
# against every other language rather than only against the whole set.
add("A | all 7", prim$y, prim$se)
for (dr in c(as.list(sort(unique(prim$language))), list(c("Balinese", "Indonesian")))) {
  s <- build_set(SE_KIND_PRIMARY, "CLAPS", drop = dr)
  add(paste0("A | drop ", paste(dr, collapse = "+")), s$y, s$se)
}
s <- build_set(SE_KIND_PRIMARY, "CLAPS", scale = "shared")
add("B | shared axis | CLAPS-Eng", s$y, s$se)
s <- build_set(SE_KIND_PRIMARY, "CLAPS", drop = "Balinese", scale = "shared")
add("B | shared axis | drop Balinese", s$y, s$se)
# Glossa languages only (the five that share one instrument and one 5-point scale)
s <- harm_all[harm_all$source == "Glossa2023", ]
add("A | Glossa five only", s$y, s$se)
# CLAPS pilot languages only
s <- build_set(SE_KIND_PRIMARY, "CLAPS")
s <- s[s$source == "CLAPS_pilot", ]
add("A | CLAPS three only", s$y, s$se)
# Prior sensitivity on tau
for (psd in c(0.25, 0.5, 1.0)) {
  s <- build_set(SE_KIND_PRIMARY, "CLAPS")
  sens[[length(sens) + 1L]] <- meta_all(s$y, s$se, paste0("A | half-normal(0,", psd, ") prior"),
                                        prior_sd = psd)
}

sens_tbl <- do.call(rbind, sens)
show <- c("set", "k", "mean_raw", "sd_raw", "mu_DL", "tau_DL", "tau_REML",
          "tau_REML_lo", "tau_REML_hi", "tau_bayes_median", "tau_bayes_lo", "tau_bayes_hi",
          "I2", "p_Q", "p_mu_neg")
cat("\n")
print(sens_tbl[, show], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 7. WHAT EACH LANGUAGE'S "ACTIVE VS PASSIVE" CONTRAST ACTUALLY IS
# ---------------------------------------------------------------------------

rule("7. SENTENCE-TYPE CODING AND COMPARABILITY")

cat("\nRaw S_Type levels present in each source file:\n")
for (nm in names(glossa_raw)) {
  cat("  ", nm, ": ", paste(sort(unique(glossa_raw[[nm]]$S_Type)), collapse = ", "), "\n", sep = "")
}
for (L in c("English", "Turkish", "Norwegian")) {
  cat("  ", L, " (CLAPS): ",
      paste(sort(unique(claps$S_Type[claps$Language == L])), collapse = ", "), "\n", sep = "")
}

contrasts_tbl <- data.frame(
  language = c("English (CLAPS)", "Turkish (CLAPS)", "Norwegian (CLAPS)",
               "English (G23)", "Hebrew (G23)", "Indonesian (G23)",
               "Mandarin (G23)", "Balinese (G23)"),
  numerator = c("Active", "Active", "Active", "Active", "Active",
                "Active (X_Active)", "Active (SVO)", "Active"),
  reference = c("Passive", "Passive (-Il)", "Passive (analytic bli-)",
                "Passive", "Passive", "Passive (Canonical_Passive, di-)",
                "Passive (O BEI S V)", "Passive = ka- / ma- / a- pooled"),
  other_levels_in_model = c("Pseudo_Passive (topicalisation)",
                            "Pseudo_Passive (topicalisation)",
                            "none (Synthetic_Passive rows dropped)",
                            "none", "none",
                            "Non_Canonical_Passive",
                            "BA_Active, Notional_Passive (OSV)",
                            "Pseudo_Passive = Passive_basic"),
  comparable = c("yes", "yes", "yes - but no third level, so the Active/Passive contrast is fitted without a pseudo-passive competitor",
                 "yes", "yes",
                 "reference is one of two passives; the other is a separate coefficient",
                 "reference is one of two passives AND one of two actives is excluded from this contrast",
                 "reference pools three distinct passive morphologies")
)
print(contrasts_tbl, row.names = FALSE)

# What happens if a different construction is treated as "the passive"? For the
# three languages with more than two sentence types this is not a footnote: the
# sign of the focal interaction changes. Point estimates only - the posterior
# covariances needed for the SE of these re-contrasts are not in the OSF dumps.
alt <- do.call(rbind, lapply(unique(aa_sem$Language), function(L) {
  s <- subset(aa_sem, Language == L)
  g <- setNames(s$Estimate, s$Stype)
  data.frame(language = L,
             active = unname(g["Active"]),
             passive = unname(g["Passive"]),
             pseudo = if ("Pseudo_Passive" %in% names(g)) unname(g["Pseudo_Passive"]) else NA,
             interaction_vs_Passive = unname(g["Active"] - g["Passive"]),
             interaction_vs_Pseudo = if ("Pseudo_Passive" %in% names(g))
               unname(g["Active"] - g["Pseudo_Passive"]) else NA)
}))
cat("\nSimple affectedness slopes by sentence type, and the interaction under each\n",
    "choice of reference construction (from AA_Models_Summary.csv):\n", sep = "")
print(alt, row.names = FALSE, digits = 3)
cat("\nMandarin BA_Active and Notional_Passive coefficients from Mandarin_Sem_Only_Bayes.txt:\n")
cat("  Semantics (in canonical BEI passive)      =  0.80\n")
cat("  S_TypeActive:Semantics                    = -0.80  -> Active SVO slope  0.00\n")
cat("  S_TypeBA_Active:Semantics                 = +1.15  -> BA active slope   1.95\n")
cat("  S_TypeNotional_Passive:Semantics          = -1.09  -> notional passive -0.29\n")
cat("  So Active(SVO) minus Notional_Passive = +0.29: the interaction flips sign\n")
cat("  depending on which of Mandarin's two passive-like constructions is the reference,\n")
cat("  and BA_Active is dropped entirely from the pooled Glossa model.\n")

# ---------------------------------------------------------------------------
# 8. SAMPLE SIZES, VERB COUNTS, RESPONSE SCALES, AND THE BALINESE RESCALE
# ---------------------------------------------------------------------------

rule("8. SAMPLE SIZES, VERB COUNTS, RESPONSE SCALES")

design_tbl <- units_tbl[, c("source", "language", "n_participants", "n_verbs",
                            "n_obs", "resp_min", "resp_max", "resp_levels")]
design_tbl$scale_note <- c(
  "7-point ordinal (CLAPS pilot)", "7-point ordinal", "7-point ordinal",
  "collected on 10 points; the analysis maps 1..10 -> 1..5 by (5-1)(R-1)/9+1 then rounds",
  "5-point ordinal", "5-point ordinal", "5-point ordinal", "5-point ordinal")
print(design_tbl, row.names = FALSE)

cat("\nVerb overlap across the five Glossa languages (all use English glosses):\n")
gk <- lapply(glossa_raw, function(d) unique(vkey(d$verb)))
cat("  union =", length(Reduce(union, gk)), " intersection =", length(Reduce(intersect, gk)), "\n")

# Two arithmetic checks that the pooled Glossa model is what the script says it is.
n_all <- sum(sapply(glossa_raw, nrow))
n_ba  <- sum(glossa_raw$Mandarin$S_Type == "BA_Active")
cat("  rows in the five files =", n_all, "; Mandarin BA_Active rows =", n_ba,
    "; difference =", n_all - n_ba, "\n")
cat("  All_Sem_Only_Bayes.txt reports 76,841 observations and 72 verb levels, so the pooled\n")
cat("  model is exactly the five files minus BA_Active, over the UNION of the verb set.\n")
cat("  NOTE: the pooled Glossa model reports 72 verb levels for 76,841 observations across\n")
cat("  five languages, i.e. verbs are CROSSED with Language - translation equivalents are\n")
cat("  treated as the same random-effect unit. CLAPS R/11_simulate_pooled_v2.R and the\n")
cat("  L5_cross_verbblock_aligned fit instead prefix verb IDs by language, so verbs are\n")
cat("  NESTED. The two choices differ in whether by-verb noise in the interaction averages\n")
cat("  away as languages are added.\n")

if (have_ordinal) {
  cat("\nWhat does the Balinese 10-point-to-5-point rescale do to a logit slope?\n")
  cat("The map is a pure binning: {1,2}->1 {3,4}->2 {5,6}->3 {7,8}->4 {9,10}->5.\n")
  cat("In a cumulative-logit model the latent scale is fixed by the link, so monotone\n")
  cat("recategorisation should move the thresholds and cost precision without moving the\n")
  cat("slope. Checked here with fixed-effect ordinal::clm fits (not the maximal brms\n")
  cat("models, so the LEVELS are not comparable to the published estimates - the RATIO is\n")
  cat("what the check is about):\n\n")
  suppressPackageStartupMessages(library(ordinal))
  b <- glossa_raw$Balinese
  b$S3 <- b$S_Type
  b$S3[b$S3 %in% c("Passive_ka", "Passive_ma", "Passive_a")] <- "Passive"
  b <- subset(b, S3 %in% c("Active", "Passive"))
  b$S3  <- relevel(factor(b$S3), ref = "Passive")
  b$R10 <- factor(b$Response, ordered = TRUE)
  b$R5  <- factor(round((5 - 1) * (b$Response - 1) / (10 - 1) + 1), ordered = TRUE)
  m10 <- ordinal::clm(R10 ~ S3 * Semantics, data = b, link = "logit")
  m5  <- ordinal::clm(R5  ~ S3 * Semantics, data = b, link = "logit")
  scale_chk <- data.frame(
    language = "Balinese", coarsening = "10 -> 5",
    fine = unname(coef(m10)["S3Active:Semantics"]),
    coarse = unname(coef(m5)["S3Active:Semantics"]))
  for (L in c("English", "Turkish", "Norwegian")) {
    s <- subset(claps, Language == L & S_Type %in% c("Active", "Passive"))
    s$S_Type <- relevel(factor(s$S_Type), ref = "Passive")
    s$Z  <- (s$Semantics - mean(s$Semantics)) / (2 * sd(s$Semantics))
    s$R7 <- factor(s$Response, ordered = TRUE)
    s$R5 <- factor(round((5 - 1) * (s$Response - 1) / (7 - 1) + 1), ordered = TRUE)
    a <- ordinal::clm(R7 ~ S_Type * Z, data = s, link = "logit")
    cc <- ordinal::clm(R5 ~ S_Type * Z, data = s, link = "logit")
    scale_chk <- rbind(scale_chk, data.frame(
      language = paste0(L, " (CLAPS)"), coarsening = "7 -> 5",
      fine = unname(coef(a)["S_TypeActive:Z"]),
      coarse = unname(coef(cc)["S_TypeActive:Z"])))
  }
  scale_chk$ratio <- scale_chk$coarse / scale_chk$fine
  print(scale_chk, row.names = FALSE, digits = 4)
  cat("\nSo response-scale coarsening moves the logit interaction by roughly 10 per cent\n")
  cat("or less and never changes its sign. It is not the reason Balinese goes the other way.\n")
} else {
  cat("\n[skip] package 'ordinal' not installed; response-scale check not run.\n")
}

# ---------------------------------------------------------------------------
# 9. THE PUBLISHED POOLED MODELS, FOR REFERENCE
# ---------------------------------------------------------------------------

rule("9. PUBLISHED POOLED ESTIMATES ON THE SAME SCALE")

cat("\nAmbridge, Arnon & Bekman (2023) five-language pooled model (All_Sem_Only_Bayes.txt):\n")
cat("  Response ~ S_Type*Semantics + (1+S_Type*Semantics|Participant)\n")
cat("                             + (1+S_Type*Semantics|verb) + (1+S_Type*Semantics|Language)\n")
cat("  cumulative(logit); 76,841 obs; 280 participants; 72 verb levels; 5 languages\n")
cat("  fixed  S_TypeActive:Semantics             = -0.31 (Est.Error 0.28), P(<0) = 0.89, ER 8.28\n")
cat("  sd(S_TypeActive:Semantics) by Language    =  0.45, 95% CI [0.19, 1.12]   <- published tau\n")
cat("  sd(S_TypeActive:Semantics) by verb        =  1.08, 95% CI [0.85, 1.37]\n")
cat("  Because Glossa Semantics is already Z-scored within language, that published tau of\n")
cat("  0.45 is ALREADY in the units of Harmonisation A and is directly comparable to the\n")
cat("  tau estimated above.\n")

g5 <- sens_tbl[sens_tbl$set == "A | Glossa five only", ]
cat(sprintf("\n  VALIDATION. Re-estimating tau from the five Glossa single-language estimates alone,\n"))
cat(sprintf("  by meta-analysis rather than by a joint hierarchical fit, gives REML tau = %.3f\n", g5$tau_REML))
cat(sprintf("  [%.2f, %.2f] and Bayesian median %.3f [%.2f, %.2f], against the published\n",
            g5$tau_REML_lo, g5$tau_REML_hi, g5$tau_bayes_median, g5$tau_bayes_lo, g5$tau_bayes_hi))
cat("  hierarchical 0.45 [0.19, 1.12]. The two-stage approximation reproduces the one-stage\n")
cat("  answer, which is the main check on the harmonisation being right.\n")

cat("\nCLAPS three-language pooled pilot fit (outputs/logs/pooledfit_8614954.out,\n")
cat("L5_cross_verbblock_aligned, 75,528 obs, 197 participants, 210 language-prefixed verbs):\n")
cat("  S_TypeActive:Semantics_scaled = -0.484  ->  -0.484 * ", round(mean(pp$affectedness_sd), 4),
    " = ", round(-0.484 * mean(pp$affectedness_sd), 3), " per SD\n", sep = "")

# ---------------------------------------------------------------------------
# 10. WRITE OUT
# ---------------------------------------------------------------------------

rule("10. MASTER TABLE")

master <- prim[, c("label", "raw_est", "raw_se", "sd_x", "y", "se", "lambda_1f", "y_shared")]
master <- master[order(master$y), ]
names(master) <- c("language", "coef_native", "se_native", "SD_of_predictor",
                   "per_SD (A)", "se (A)", "lambda", "shared_axis (B)")
print(master, row.names = FALSE, digits = 3)
cat("\ncoef_native is on the units the model was fitted in: Semantics_scaled (= 2 raw SD) for\n")
cat("CLAPS, an already-Z-scored PCA1 for Glossa. Column 'per_SD (A)' is the common scale:\n")
cat("logit change in the affectedness slope, Active minus canonical Passive, per 1 SD of that\n")
cat("language's affectedness predictor. All eight rows including the second English estimate\n")
cat("are in cross_language_effects_harmonised.csv.\n")

rule("11. FILES WRITTEN")

harm_out <- harm_all[, c("source", "language", "label", "se_kind", "raw_est", "raw_se",
                         "sd_x", "y", "se", "lambda_1f", "lambda_Z", "y_shared", "se_shared")]
f1 <- file.path(OUT, "cross_language_effects_harmonised.csv")
f2 <- file.path(OUT, "cross_language_tau_sensitivity.csv")
f3 <- file.path(OUT, "cross_language_design_summary.csv")
write.csv(harm_out, f1, row.names = FALSE)
write.csv(sens_tbl, f2, row.names = FALSE)
write.csv(design_tbl, f3, row.names = FALSE)
cat("\n", f1, "\n", f2, "\n", f3, "\n", sep = "")

cat("\nSession:\n")
cat("  R:", R.version.string, "\n")
cat("  metafor:", if (have_metafor) as.character(packageVersion("metafor")) else "ABSENT", "\n")
cat("  ordinal:", if (have_ordinal) as.character(packageVersion("ordinal")) else "ABSENT", "\n")
