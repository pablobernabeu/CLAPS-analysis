# R/13_simulate_klanguage_pooled.R
# ---------------------------------------------------------------------------
# Operating characteristics of the pooled cross-linguistic test at K languages.
#
# WHAT THIS ANSWERS
#   The PI asked (7 September 2026) for the probability that the pooled test
#   clears a directional Bayes factor of 10 on the affectedness-by-sentence-type
#   interaction when K = 6, 8, 10 or 12 language teams each contribute 80
#   participants and 72 verbs, so that the Stage 1 Registered Report can carry a
#   rule of the form "we proceed to Stage 2 only if at least Y languages join".
#   He later asked (23 September) for the same test at 10 and 20 languages when the
#   seven languages already measured all take part at their estimated effects and
#   only the remaining places are unknown; run_partly_known_cell() covers that.
#
# WHY IT IS NOT THE POOLED SIMULATOR IN R/11
#   R/11_simulate_pooled_v2.R simulates participant-level data for three FIXED
#   languages and fits L5_cross_maximal to it. Two things rule that route out here.
#
#   First cost. Measured over the 135 usable replicates in
#   outputs/design_pooled_v2, one such fit takes a median of 34 to 170 hours
#   depending on sample size. Twelve languages is four times the data and a larger
#   by-Language covariance block, so a grid over K = 6, 8, 10, 12 at any useful
#   replication is out of reach.
#
#   Second, and more important, the estimand changes. The three-language run treats
#   English, Turkish and Norwegian as fixed and asks how well the model recovers
#   their average. "At least Y languages, we do not yet know which" makes the K
#   languages draws from a population of languages, which adds a source of variation
#   the fixed-language run does not carry. A K-language answer is therefore not an
#   extrapolation of the existing curve, and would not be one however cheap the fits
#   became.
#
# THE REDUCTION THAT MAKES A CHEAP ENGINE LEGITIMATE
#   A directional Savage-Dickey Bayes factor is a monotone function of a single
#   number. Writing p for the posterior probability that the coefficient has the
#   predicted sign and p0 for the same probability under the prior,
#
#       BF = [p / (1 - p)] / [p0 / (1 - p0)]
#
#   so BF >= B is exactly z := qnorm(p) >= qnorm(B * o0 / (1 + B * o0)) with
#   o0 = p0 / (1 - p0). Over the 135 completed three-language replicates the
#   measured prior probability is 0.4997, so the threshold for B = 10 is z = 1.335,
#   and the realised z is very close to normal: N(1.386, 0.311), which returns
#   P(BF >= 10) = 0.566 against the 0.578 actually observed. The whole operating
#   characteristic of a sixty-five-hour fit is therefore carried by two numbers, and
#   an engine that reproduces those two numbers at K = 3 has earned the right to be
#   run at K = 12. calibrate_against_pooled_k3() in this file is that check, and
#   pooled_k3_target() computes its target from the replicates themselves. The
#   target pools the six sample sizes of that run (N = 50 to 150; z shows no trend in
#   N) and counts every fit that finished, of which only 13 pass the stricter
#   convergence flag; pooled_k3_target() reports both, so neither is hidden.
#
# THE ENGINE
#   Two-stage: simulate each language's interaction estimate, then combine them in a
#   Bayesian random-effects meta-analysis and read the directional Bayes factor off
#   the posterior for the population mean. The two-stage step is not assumed to be
#   adequate; it is checked twice. Against the published one-stage hierarchical fit,
#   where meta-analysing the five Ambridge, Arnon & Bekman (2023) single-language
#   estimates returns tau = 0.48 against the published 0.45 (scripts/extract_cross_language_effects.R,
#   section 9). And against the three-language brms runs, by the calibration above.
#
#   The posterior is computed by grid integration over tau rather than by MCMC.
#   With one parameter to integrate out and a conjugate normal conditional for mu
#   this is exact to the width of the grid, costs about a millisecond, and removes
#   sampler noise from a quantity that is itself an operating characteristic.
#
# THE PART THAT MATTERS MOST
#   Power here is governed by tau, the between-language SD, far more than by K.
#   tau is estimated from seven languages, one of which (Balinese) shows the effect
#   in the opposite direction, so it is itself uncertain: the posterior median is
#   0.39 with a 95% interval of roughly [0.20, 0.76]. Every entry point below
#   therefore runs in one of three modes, and the report is expected to carry all
#   three rather than the first alone:
#
#     "point"     — fix (mu, tau) at a stated value. Conditional power.
#     "assurance" — draw (mu, tau) from their joint posterior given the seven
#                   observed languages, once per replicate. This is the probability
#                   of success in the sense of O'Hagan, Stevens & Campbell (2005,
#                   doi:10.1002/pst.175), and it is the honest headline.
#     "safeguard" — fix tau at a stated upper quantile of its posterior and mu at a
#                   lower quantile of its own, mirroring the safeguard convention
#                   already used elsewhere in this package (Perugini, Gallucci &
#                   Costantini, 2014, doi:10.1027/1614-2241/a000081).
#
# WHAT IT CANNOT ESTABLISH
#   Nothing here can tell you which languages will join. The seven in hand are a
#   convenience sample, they are not phylogenetically or areally independent, and
#   the design analysis inherits that. What the engine gives is the operating
#   characteristic of the pooled test under a stated population of languages; the
#   population is an assumption, and it is the assumption the report has to defend.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stats)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# 0. The calibrated tau prior
# ---------------------------------------------------------------------------
#
# A standalone random-effects meta-analysis of seven languages would use a weakly
# informative half-normal(0, 0.5) on tau. That is NOT what belongs here, because
# the engine has to imitate a particular analysis model, not a generic one.
#
# L5_cross_maximal carries a six-term correlated by-Language block, estimated from
# as few as three language groups. Almost none of that block is identified at that
# K, and the resulting uncertainty propagates into the fixed interaction: over the
# 135 completed three-language replicates the posterior SD of the pooled
# interaction is 0.206 per SD of affectedness, where a plug-in random-effects
# calculation on the same three languages gives 0.103. The confirmatory model is
# twice as uncertain as a one-coefficient meta-analysis of the same data, and an
# engine that ignores that would overstate power at every K.
#
# So the tau prior is treated as the engine's single calibration parameter and
# fitted to the three-language runs. The fit holds each language's SE at its
# verb-information floor for the interaction (0.104, 0.076 and 0.042 per SD for
# English, Turkish and Norwegian) with no verb term shared across languages,
# because each brms run gave its three languages separate verb lists. Scored
# exactly, with every Bayes factor compared with the threshold unrounded, over 4,000
# engine replicates a row:
#
#              mean z   sd z   P(BF>=10)  P(BF>=6)  P(BF>=3)
#   brms        1.386  0.311     0.578     0.852     0.978
#   scale 1.0   1.430  0.310     0.607     0.878     0.992
#   scale 1.1   1.408  0.308     0.582     0.864     0.991
#   scale 1.2   1.389  0.307     0.557     0.848     0.990
#
# 1.1 comes closest on P(BF >= 10), the quantity the design analysis reports; no
# scale matches every moment at once. (The brms row once read 0.859 and 0.985 for
# the last two columns, which is what rounding each Bayes factor to one decimal
# before thresholding gives.)
#
# The production cells do not use the floor SEs. They give every language the SE a
# team at 80 participants and 72 verbs would achieve, split into an independent part
# (0.083) and a part shared through the common verb list (0.030). Rerun at those SEs
# the K = 3 check gives P(BF >= 10) = 0.536 at scale 1.1, the scale that would
# restore 0.578 there is about 0.95, and no scale brings sd z down to 0.311. A refit
# to the production SEs would raise power at six or more languages by under a point
# (seven-language assurance at K = 6: 0.562 to 0.569) and at K = 3 by two to four
# points. Which SE configuration the calibration should use is a decision for the
# maintainer (audit item klanguage-engine#4). Until it is made the scale stays at
# 1.1, and the calibration section of scripts/run_klanguage_design_analysis.R
# reports both configurations side by side.
#
# Holding the PRIOR fixed while K grows is the conservative choice that also
# happens to be the principled one: a prior is overwhelmed by data as groups
# accumulate, so the inflation it represents shrinks on its own as K rises, exactly
# as the real model's by-Language block becomes identified. Nothing has to be
# assumed about how fast.
CALIBRATED_TAU_PRIOR_SCALE <- 1.1

# The bracketing values the calibration cannot separate. Every headline figure is
# reported across this range, because a single scale would read as more settled
# than five moments from 135 replicates can make it.
CALIBRATION_BRACKET <- c(1.0, 1.2)

# The DESIGN prior on tau is a different object from the analysis prior above, and
# conflating them is easy to do and wrong in a way that quietly moves every number.
#
#   The analysis prior asks: when the real K-language model is fitted, how much
#   uncertainty about the by-Language SD does it carry? That is a property of
#   L5_cross_maximal and its LKJ-plus-half-t block, and it is what the 1.1 above was
#   calibrated to reproduce.
#
#   The design prior asks: given the seven languages we can see, what do we believe
#   tau is? That is a statement about the world, and it should use the weakly
#   informative half-normal(0, 0.5) that the meta-analysis literature recommends for
#   few groups, not a scale chosen to imitate a particular sampler's behaviour.
#
# Using the calibrated 1.1 for both inflates the design prior's tau (0.42 rather than
# 0.39 on the seven languages) and so understates power at every K.
DESIGN_TAU_PRIOR_SCALE <- 0.5

# ---------------------------------------------------------------------------
# 1. The Bayesian random-effects meta-analysis, by grid integration
# ---------------------------------------------------------------------------

#' Posterior for the population mean of a random-effects meta-analysis.
#'
#' Model: y_k ~ N(theta_k, s_k^2), theta_k ~ N(mu, tau^2), with a normal prior on
#' mu and a half-normal prior on tau. Conditional on tau the posterior for mu is
#' normal in closed form, so only tau is integrated numerically.
#'
#' The half-normal prior on tau follows the recommendation for random-effects
#' meta-analysis with few studies, where a flat or half-Cauchy prior lets tau run
#' away on the strength of two or three groups and an inverse-gamma prior is
#' sharply informative near zero without anyone intending it (Gelman, 2006,
#' doi:10.1214/06-BA117A; Roever, Bender, Dias et al., 2021,
#' doi:10.1002/jrsm.1475). The default scale of 0.5 is on the same scale as the
#' effect itself, which is weakly informative here: the seven-language posterior
#' median for tau moves from 0.33 to 0.42 as the prior scale moves from 0.25 to 1
#' (outputs/design_summary_pilot/cross_language_tau_sensitivity.csv), so the
#' choice is reported as a sensitivity rather than defended as correct.
#'
#' @param y,s Numeric vectors of per-language estimates and their standard errors.
#' @param prior_mu_sd Normal prior SD on mu, centred at zero. The CLAPS prior on
#'   the interaction is centred at zero and this mirrors it; set to Inf for a flat
#'   prior, which is what the reference two-stage meta-analyses use.
#' @param prior_tau_scale Half-normal scale for tau. The default is NOT the
#'   weakly informative 0.5 that a standalone meta-analysis would use. It is
#'   CALIBRATED_TAU_PRIOR_SCALE, set so that the engine reproduces what the
#'   confirmatory model actually does at K = 3 (see calibrate_against_pooled_k3).
#' @param tau_grid Grid of tau values to integrate over.
#' @return List with the posterior mean, SD and sign probability for mu, the
#'   posterior median of tau, and the grid weights.
ra_meta_posterior <- function(y, s,
                              prior_mu_sd = Inf,
                              prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                              tau_grid = seq(0, 3, length.out = 601)) {
  stopifnot(length(y) == length(s), all(is.finite(y)), all(s > 0))

  # log marginal likelihood of the data given tau, integrating mu out analytically.
  # For each tau: v_k = s_k^2 + tau^2, w_k = 1/v_k.
  # With a flat prior on mu the marginal is prod N(y_k | mu, v_k) integrated over mu.
  # With a proper normal prior N(0, p^2) the same integral is available in closed form.
  p2 <- if (is.finite(prior_mu_sd)) prior_mu_sd^2 else Inf

  # Vectorised over the grid: V is G x K. The loop version of this ran at about
  # forty replicates a second, which is too slow for a twenty-thousand-replicate
  # calibration; the matrix form runs the same arithmetic in one pass.
  V  <- outer(tau_grid^2, s^2, `+`)          # G x K
  W  <- 1 / V
  Sw   <- rowSums(W)
  Swy  <- as.vector(W %*% y)
  Swy2 <- as.vector(W %*% (y^2))
  base <- -0.5 * (rowSums(log(2 * pi * V)) + Swy2)
  prec <- if (is.infinite(p2)) Sw else Sw + 1 / p2
  log_marg <- base + 0.5 * (Swy^2 / prec) + 0.5 * log(2 * pi / prec) -
    if (is.infinite(p2)) 0 else 0.5 * log(2 * pi * p2)

  log_prior_tau <- dnorm(tau_grid, 0, prior_tau_scale, log = TRUE)  # half-normal, up to 2x
  lp <- log_marg + log_prior_tau
  lp <- lp - max(lp)
  wt <- exp(lp)
  # Trapezoidal weights, so the grid spacing is honoured and an unevenly spaced
  # grid integrates correctly. A one-point grid pins tau at that value, which is
  # how a conditional-on-tau calculation is expressed here; it has no spacing, so
  # it takes unit weight rather than a difference of an empty vector.
  trap <- trapezoid_weights(tau_grid)
  stopifnot(length(trap) == length(tau_grid))
  wt <- wt * trap
  wt <- wt / sum(wt)

  # Conditional posterior for mu at each tau, from the same matrices.
  m  <- Swy / prec
  sd <- sqrt(1 / prec)

  post_mean <- sum(wt * m)
  post_var  <- sum(wt * (sd^2 + m^2)) - post_mean^2
  # P(mu < 0) as a mixture over the tau grid
  p_neg <- sum(wt * pnorm(0, mean = m, sd = sd))

  cum <- cumsum(wt)
  tau_median <- tau_grid[which.min(abs(cum - 0.5))]

  list(mu_mean = post_mean, mu_sd = sqrt(max(post_var, 0)),
       p_mu_negative = p_neg, tau_median = tau_median,
       tau_grid = tau_grid, weights = wt)
}

#' Quadrature weights for a tau grid, shared by the one-at-a-time and batch posteriors.
trapezoid_weights <- function(tau_grid) {
  if (length(tau_grid) == 1L) {
    return(1)
  }
  dx <- diff(tau_grid)
  c(dx[1], (utils::head(dx, -1) + utils::tail(dx, -1)), utils::tail(dx, 1)) / 2
}

#' The same posterior for many simulated studies at once.
#'
#' A design cell analyses thousands of simulated studies that share one set of
#' standard errors, so everything on the grid that depends only on s is computed
#' once and the rest runs as matrix arithmetic over replicates. The arithmetic is
#' that of ra_meta_posterior(), step for step and in the same order, so each
#' replicate's p_mu_negative is the number the one-at-a-time function returns
#' (tests/testthat/test-klanguage-pooled.R checks this), at about a sixth of the cost.
#' That is what makes cells of 4,000 replicates affordable at K = 20 and the
#' million-replicate variance decomposition affordable at all.
#'
#' @param Y Matrix of per-language estimates, one row per replicate and one column
#'   per language.
#' @param s Standard errors, one per column of Y, or a single value for all.
#' @param want_tau Also return each replicate's posterior median of tau.
#' @param chunk Replicates per pass; memory is about 60 bytes per grid point per
#'   replicate in a pass.
#' @return List with p_mu_negative, and tau_median when asked for, one per row of Y.
ra_meta_batch <- function(Y, s,
                          prior_mu_sd = Inf,
                          prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                          tau_grid = seq(0, 3, length.out = 601),
                          want_tau = FALSE, chunk = 500L) {
  Y <- as.matrix(Y)
  K <- ncol(Y)
  s <- if (length(s) == 1L) rep(s, K) else s
  stopifnot(length(s) == K, all(is.finite(Y)), all(s > 0))
  p2 <- if (is.finite(prior_mu_sd)) prior_mu_sd^2 else Inf

  V <- outer(tau_grid^2, s^2, `+`)
  W <- 1 / V
  Sw <- rowSums(W)
  log_det <- rowSums(log(2 * pi * V))
  prec <- if (is.infinite(p2)) Sw else Sw + 1 / p2
  norm_const <- 0.5 * log(2 * pi / prec)
  prior_const <- if (is.infinite(p2)) 0 else 0.5 * log(2 * pi * p2)
  log_prior_tau <- dnorm(tau_grid, 0, prior_tau_scale, log = TRUE)
  trap <- trapezoid_weights(tau_grid)
  sd_mu <- sqrt(1 / prec)

  n <- nrow(Y)
  p_neg <- numeric(n)
  tau_med <- if (want_tau) numeric(n) else NULL
  for (start in seq(1L, max(n, 1L), by = chunk)) {
    if (n == 0L) break
    idx <- start:min(n, start + chunk - 1L)
    Yt <- t(Y[idx, , drop = FALSE])
    Swy <- W %*% Yt                                   # G x chunk
    Swy2 <- W %*% (Yt^2)
    base <- -0.5 * (log_det + Swy2)
    log_marg <- base + 0.5 * (Swy^2 / prec) + norm_const - prior_const
    lp <- log_marg + log_prior_tau
    lp <- sweep(lp, 2L, apply(lp, 2L, max))
    wt <- exp(lp) * trap
    wt <- sweep(wt, 2L, colSums(wt), "/")
    m <- Swy / prec
    p_neg[idx] <- colSums(wt * pnorm(0, mean = m, sd = sd_mu))
    if (want_tau) {
      # apply() drops to a vector when the grid has one point, so reshape.
      cum <- matrix(apply(wt, 2L, cumsum), nrow = length(tau_grid))
      tau_med[idx] <- tau_grid[apply(abs(cum - 0.5), 2L, which.min)]
    }
  }
  list(p_mu_negative = p_neg, tau_median = tau_med)
}

#' Directional Bayes factor for the population mean having the predicted sign.
#'
#' Uses the same odds-ratio definition as R/05_hypothesis_tests.R, so a Bayes
#' factor from this engine and one from a fitted brms model mean the same thing.
#' The prior sign probability is 0.5 whenever the prior on mu is symmetric about
#' zero, which it is here and is in the CLAPS priors; it is kept as an argument so
#' that an asymmetric design prior can be substituted without editing the arithmetic.
directional_bf_meta <- function(post, direction = "negative", prior_sign_prob = 0.5) {
  p <- if (direction == "negative") post$p_mu_negative else 1 - post$p_mu_negative
  p <- min(max(p, 1e-12), 1 - 1e-12)
  (p / (1 - p)) / (prior_sign_prob / (1 - prior_sign_prob))
}

#' directional_bf_meta() for a vector of P(mu < 0), element by element, clamp included.
bf_from_p_negative <- function(p, prior_sign_prob = 0.5) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  (p / (1 - p)) / (prior_sign_prob / (1 - prior_sign_prob))
}

#' The z on which the calibration is stated, qnorm of P(mu < 0) with the same clamp.
z_from_p_negative <- function(p) qnorm(pmin(pmax(p, 1e-12), 1 - 1e-12))

# ---------------------------------------------------------------------------
# 2. The population of languages
# ---------------------------------------------------------------------------

#' Joint posterior draws of (mu, tau) given the observed languages.
#'
#' Draws are taken from the same grid posterior used above: tau is drawn from its
#' marginal grid weights and mu from its conditional normal at that tau. This keeps
#' the design analysis and the estimation of the language population on one set of
#' assumptions instead of two.
#'
#' @param y,s Observed per-language estimates and SEs, on the common per-SD scale
#'   (outputs/design_summary_pilot/cross_language_effects_harmonised.csv).
draw_population_params <- function(n, y, s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE,
                                   prior_mu_sd = Inf, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  post <- ra_meta_posterior(y, s, prior_mu_sd = prior_mu_sd,
                            prior_tau_scale = prior_tau_scale)
  idx <- sample.int(length(post$tau_grid), n, replace = TRUE, prob = post$weights)
  tau <- post$tau_grid[idx]
  # The conditional mean and SD of mu at every grid point are already available
  # from the same posterior, so the draw is an index lookup rather than a refit.
  p2   <- if (is.finite(prior_mu_sd)) prior_mu_sd^2 else Inf
  V    <- outer(post$tau_grid^2, s^2, `+`)
  W    <- 1 / V
  prec <- rowSums(W) + if (is.infinite(p2)) 0 else 1 / p2
  m    <- as.vector(W %*% y) / prec
  sdm  <- sqrt(1 / prec)
  mu   <- rnorm(n, m[idx], sdm[idx])
  data.frame(mu = mu, tau = tau)
}

# ---------------------------------------------------------------------------
# 2b. The observed languages, and which construction counts as the passive
# ---------------------------------------------------------------------------
#
# These turn the harmonised table (outputs/design_summary_pilot/
# cross_language_effects_harmonised.csv) into the language sets the design draws
# from. They live here rather than in the run script so that the tests exercise the
# code that produced the reported numbers, not a copy of it.

#' Split a language's SE into the part shared across languages and the part that is not.
#'
#' The verb-information floor is shareable, because the teams translate one verb
#' list, and a fraction rho of it is common to any two languages. The participant
#' term is never shared. The total is unchanged by the split.
split_se <- function(se_verb, se_ppt, rho) {
  list(shared = se_verb * sqrt(rho),
       indep  = sqrt(se_verb^2 * (1 - rho) + se_ppt^2))
}

#' One row per language, preferring the CLAPS estimate of a language both sources measured.
#'
#' The CLAPS estimate of English is the one collected under the protocol every new
#' team will follow. The sort is stable, so languages keep their order within source.
one_per_language <- function(df) {
  df <- df[order(df$source != "CLAPS_pilot"), ]
  df[!duplicated(df$language), ]
}

#' Posterior draws of the pilot interaction, in Semantics_scaled units.
#'
#' The archived pilot posteriors are the only exact source of the interaction's own
#' posterior SD; the harmonised table carries the main-effect SD in its place. They
#' are gitignored, so a fresh clone lacks them, and the run must stop instead of
#' quietly substituting the main-effect SDs, which moves the CLAPS regime's power by
#' about two points (audit item klanguage-engine#8).
read_claps_interaction_draws <- function(dir = "outputs/pilot_models",
                                         langs = c("English", "Turkish", "Norwegian"),
                                         term = "S_TypeActive:Semantics_scaled") {
  files <- file.path(dir, paste0("pilot_dgp_v2_pilot_", langs, ".rds"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("[klang] pilot posteriors not found: ", paste(missing, collapse = ", "),
         ". They are gitignored; copy them into ", dir, " or pass --pilot-models=<dir>.")
  }
  stats::setNames(lapply(files, function(f) {
    as.numeric(as.matrix(readRDS(f)$fixef_draws)[, term])
  }), langs)
}

#' The harmonised estimates the engine uses, one row per source and language.
#'
#' Keeps the rows whose SE is a posterior SD and swaps each CLAPS row's SE for the
#' interaction's own posterior SD, `claps_se`, named by language and already on the
#' per-SD scale. The table's CLAPS SEs are those of the affectedness MAIN effect,
#' which is the wrong quantity to weight a meta-analysis of the interaction by.
prepare_language_table <- function(h, claps_se) {
  hc <- h[h$se_kind %in% c("slope_sd", "posterior_sd"), ]
  hit <- hc$source == "CLAPS_pilot" & hc$language %in% names(claps_se)
  hc$se[hit] <- claps_se[hc$language[hit]]
  hc[!duplicated(paste(hc$source, hc$language)), ]
}

# For three of the seven languages the published model has a pseudo-passive as well as
# the passive. The project follows the published coding, with the passive as the
# reference (settled with the PI on 1 October 2026), and every reported figure uses
# it. The "alternative" reading measures the interaction against the pseudo-passive
# instead, a sensitivity the outputs still carry.
# Both coefficients come from the same published single-language model, fitted with
# the canonical passive as the reference level, so the interaction measured against
# the other construction is the difference of the two S_Type-by-Semantics terms:
#
#   Balinese    Balinese_Sem_Only_Bayes_redux.txt  Active  0.33 (0.13), Passive_basic         -0.53 (0.13)
#   Indonesian  Indonesian_Sem_Only_Bayes.txt       Active -0.07 (0.08), Non_Canonical_Passive -0.34 (0.10)
#   Mandarin    Mandarin_Sem_Only_Bayes.txt         Active -0.80 (0.21), Notional_Passive      -1.09 (0.18)
#
# The differences, +0.86, +0.27 and +0.29, reproduce Active minus Pseudo_Passive in
# the simple slopes of AA_Models_Summary.csv (Ambridge, Arnon & Bekman, 2023, OSF
# materials). Hebrew has only two sentence types and the Glossa English row is not
# in any population, so every other language is the same under both readings.
ALT_PASSIVE_CONTRASTS <- data.frame(
  language     = c("Balinese", "Indonesian", "Mandarin"),
  construction = c("Pseudo_Passive (Passive_basic)", "Non_Canonical_Passive", "Notional_Passive (OSV)"),
  b_active     = c(0.33, -0.07, -0.80),
  se_active    = c(0.13, 0.08, 0.21),
  b_other      = c(-0.53, -0.34, -1.09),
  se_other     = c(0.13, 0.10, 0.18),
  stringsAsFactors = FALSE
)

#' Re-express the three forked languages against the alternative passive construction.
#'
#' The OSF dumps carry no posterior covariances, so the SE of a difference of two
#' coefficients cannot be recovered exactly. "independent" treats the two as
#' uncorrelated. Two treatment contrasts against the same reference level are more
#' often positively correlated than not, which would make the true SE smaller, so
#' this is the cautious choice unless a pair happens to correlate negatively.
#' "canonical" keeps the canonical contrast's SE, as a sensitivity. An exact SE needs
#' the three models refitted with the other construction as the reference.
#'
#' @param hc Table from prepare_language_table().
#' @param reading "canonical" returns hc unchanged.
#' @param alt_se "independent" or "canonical", as above.
apply_passive_reading <- function(hc, reading = c("canonical", "alternative"),
                                  alt_se = c("independent", "canonical"),
                                  alt = ALT_PASSIVE_CONTRASTS) {
  reading <- match.arg(reading)
  alt_se <- match.arg(alt_se)
  if (reading == "canonical") {
    return(hc)
  }
  i <- which(hc$source == "Glossa2023" & hc$language %in% alt$language)
  a <- alt[match(hc$language[i], alt$language), ]
  # The contrasts are typed from the model dumps, while the canonical estimates reach
  # the harmonised table by another route, so check they agree before subtracting.
  if (!isTRUE(all.equal(hc$raw_est[i], a$b_active)) ||
        !isTRUE(all.equal(hc$raw_se[i], a$se_active))) {
    stop("[passive] ALT_PASSIVE_CONTRASTS no longer matches the harmonised Active coefficients")
  }
  raw_se <- if (alt_se == "independent") sqrt(a$se_active^2 + a$se_other^2) else a$se_active
  hc$raw_est[i] <- a$b_active - a$b_other
  hc$raw_se[i] <- raw_se
  hc$y[i] <- hc$raw_est[i] * hc$sd_x[i]
  hc$se[i] <- hc$raw_se[i] * hc$sd_x[i]
  # The shared-axis harmonisation was never computed against the other construction.
  hc$y_shared[i] <- NA_real_
  hc$se_shared[i] <- NA_real_
  hc
}

#' The language sets the design draws from, under one reading of the passive.
#'
#' claps_protocol, seven and glossa_protocol are the three tau regimes; the two
#' leave-one-out sets are sensitivities. Each keeps `source`, which is how a
#' partly-known run tells a CLAPS pilot language from a 2023 one.
build_populations <- function(hc, reading = "canonical", alt_se = "independent") {
  hc <- apply_passive_reading(hc, reading, alt_se)
  per_language <- one_per_language(hc)
  pops <- list(
    claps_protocol  = hc[hc$source == "CLAPS_pilot", ],
    seven           = per_language,
    glossa_protocol = hc[hc$source == "Glossa2023" & hc$language != "English", ],
    drop_balinese   = per_language[per_language$language != "Balinese", ],
    drop_hebrew     = per_language[per_language$language != "Hebrew", ]
  )
  lapply(pops, function(d) data.frame(language = d$language, source = d$source, y = d$y, s = d$se))
}

#' What a language set says about the population, under the DESIGN prior on tau.
#'
#' The first nine columns are those of klanguage_populations.csv as first reported.
population_summary <- function(pop, name = NA_character_, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE) {
  post <- ra_meta_posterior(pop$y, pop$s, prior_tau_scale = prior_tau_scale)
  # Cochran's Q, as the test of whether the spread is more than sampling error.
  w <- 1 / pop$s^2
  mu_fe <- sum(w * pop$y) / sum(w)
  q_stat <- sum(w * (pop$y - mu_fe)^2)
  df <- nrow(pop) - 1
  cum <- cumsum(post$weights)
  tau_at <- function(a) post$tau_grid[which.min(abs(cum - a))]
  data.frame(population = name, K = nrow(pop), mean_effect = mean(pop$y),
             tau_median = post$tau_median, p_mu_negative = post$p_mu_negative,
             Q = q_stat, df = df, p_Q = pchisq(q_stat, df, lower.tail = FALSE),
             I2 = max(0, (q_stat - df) / q_stat),
             sd_effect = stats::sd(pop$y), mu_mean = post$mu_mean, mu_sd = post$mu_sd,
             tau_lo = tau_at(0.025), tau_hi = tau_at(0.975))
}

#' Per-row protocol-transfer doubt: none for a CLAPS pilot language, `transfer_sd` otherwise.
#'
#' A CLAPS pilot language was measured under the protocol Stage 2 will use, so its
#' earlier estimate transfers by construction. Only a language measured under
#' another protocol carries the extra doubt.
transfer_sd_by_source <- function(known, transfer_sd, pilot_source = "CLAPS_pilot") {
  if (is.null(known$source)) {
    stop("[transfer] `known` needs a source column to tell CLAPS pilot rows from the others")
  }
  transfer_sd * (known$source != pilot_source)
}

# ---------------------------------------------------------------------------
# 3. One replicate of a K-language study
# ---------------------------------------------------------------------------

#' Simulate one K-language pooled study and score it.
#'
#' @param K Number of languages contributing.
#' @param mu,tau Population mean and between-language SD of the interaction, on
#'   the per-SD-of-affectedness scale.
#' @param se_lang Per-language standard error of the interaction that a CLAPS team
#'   will achieve at the planned sample. Either one number applied to every
#'   language, or a vector of length K.
#' @param shared_verb_sd Standard deviation of a verb-sample offset COMMON to every
#'   language. Nonzero when the K teams translate one verb list, so that the
#'   idiosyncrasy of the 72 chosen verbs biases every language the same way and does
#'   not average out as languages are added. Zero when each team draws its own verbs.
#'   The truth lies between; the grid runs both ends and the middle.
#' @param theta_fixed Optional vector of K true per-language effects. When supplied,
#'   languages are treated as FIXED rather than sampled, which is what the existing
#'   three-language runs assume.
#' @return One-row data frame with the estimates, the Bayes factor and z.
simulate_klanguage_replicate <- function(K, mu, tau, se_lang,
                                         shared_verb_sd = 0,
                                         prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                         prior_mu_sd = Inf,
                                         theta_fixed = NULL) {
  s <- if (length(se_lang) == 1L) rep(se_lang, K) else se_lang
  stopifnot(length(s) == K)

  theta <- if (!is.null(theta_fixed)) {
    stopifnot(length(theta_fixed) == K)
    theta_fixed
  } else {
    rnorm(K, mu, tau)
  }

  # One offset shared by every language, then independent estimation error.
  common <- if (shared_verb_sd > 0) rnorm(1, 0, shared_verb_sd) else 0
  y <- theta + common + rnorm(K, 0, s)

  post <- ra_meta_posterior(y, s, prior_mu_sd = prior_mu_sd,
                            prior_tau_scale = prior_tau_scale)
  bf <- directional_bf_meta(post, direction = "negative")

  data.frame(K = K, mu_true = mu, tau_true = tau,
             mu_hat = post$mu_mean, mu_sd = post$mu_sd,
             tau_hat = post$tau_median,
             p_negative = post$p_mu_negative,
             z = qnorm(min(max(post$p_mu_negative, 1e-12), 1 - 1e-12)),
             bf = bf)
}

# ---------------------------------------------------------------------------
# 4. A cell: many replicates at one K
# ---------------------------------------------------------------------------

#' Run one cell of the K-language design analysis.
#'
#' @param mode "point", "assurance" or "safeguard" — see the file header.
#' @param pop Data frame with columns y and s: the observed languages that define
#'   the population. Required in assurance and safeguard modes.
#' @param mu_point,tau_point Used in point mode.
#' @param safeguard_q Quantile pair c(mu, tau): mu is taken at the conservative end
#'   of its posterior and tau at the pessimistic end of its own.
#' @param bf_threshold Evidence threshold. 10 is the preregistered criterion.
run_klanguage_cell <- function(K, n_reps, se_lang,
                               mode = c("assurance", "point", "safeguard"),
                               pop = NULL,
                               mu_point = NULL, tau_point = NULL,
                               safeguard_q = c(0.25, 0.75),
                               shared_verb_sd = 0,
                               prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                               design_tau_prior_scale = DESIGN_TAU_PRIOR_SCALE,
                               prior_mu_sd = Inf,
                               bf_threshold = 10,
                               seed = 1L) {
  mode <- match.arg(mode)
  set.seed(seed)

  if (mode == "point") {
    stopifnot(!is.null(mu_point), !is.null(tau_point))
    params <- data.frame(mu = rep(mu_point, n_reps), tau = rep(tau_point, n_reps))
  } else if (mode == "assurance") {
    stopifnot(!is.null(pop))
    params <- draw_population_params(n_reps, pop$y, pop$s,
                                     prior_tau_scale = design_tau_prior_scale,
                                     prior_mu_sd = prior_mu_sd)
  } else {
    stopifnot(!is.null(pop))
    d <- draw_population_params(20000L, pop$y, pop$s,
                                prior_tau_scale = design_tau_prior_scale,
                                prior_mu_sd = prior_mu_sd)
    # mu towards zero, tau towards its upper tail: the unfavourable corner.
    mu_sg  <- quantile(d$mu,  1 - safeguard_q[1], names = FALSE)  # least negative
    tau_sg <- quantile(d$tau, safeguard_q[2],     names = FALSE)
    params <- data.frame(mu = rep(mu_sg, n_reps), tau = rep(tau_sg, n_reps))
  }

  reps <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
    simulate_klanguage_replicate(K, params$mu[i], params$tau[i], se_lang,
                                 shared_verb_sd = shared_verb_sd,
                                 prior_tau_scale = prior_tau_scale,
                                 prior_mu_sd = prior_mu_sd)
  }))

  p <- mean(reps$bf >= bf_threshold)
  data.frame(
    K = K, mode = mode, n_reps = n_reps, se_lang = se_lang[1],
    shared_verb_sd = shared_verb_sd, prior_tau_scale = prior_tau_scale,
    design_tau_prior_scale = design_tau_prior_scale,
    bf_threshold = bf_threshold,
    power = p,
    mcse = sqrt(p * (1 - p) / n_reps),
    power_lo = binom.test(round(p * n_reps), n_reps)$conf.int[1],
    power_hi = binom.test(round(p * n_reps), n_reps)$conf.int[2],
    mean_z = mean(reps$z), sd_z = sd(reps$z),
    median_bf = median(reps$bf),
    mean_tau_hat = mean(reps$tau_hat),
    # The share of replicates whose pooled posterior favours the UNPREDICTED sign,
    # P(BF < 1), whether or not either threshold was reached. It is not the
    # conditional Type S rate of Gelman & Carlin (2014, doi:10.1177/1745691614551642),
    # which conditions on the threshold being cleared and is several times smaller;
    # renaming the column and adding that rate waits for the next regeneration
    # (audit item klanguage-engine#10).
    p_wrong_sign = mean(reps$p_negative < 0.5),
    seed = seed
  )
}

# ---------------------------------------------------------------------------
# 5. The calibration check against the completed three-language brms runs
# ---------------------------------------------------------------------------

#' Does the cheap engine reproduce the expensive one at K = 3?
#'
#' The three-language brms runs draw each language's true effect from its own pilot
#' posterior, once per replicate, so this does the same: each replicate takes one
#' draw per language from `pilot_draws`, adds sampling error at the stated SEs, and
#' scores the pooled Bayes factor. The random-number stream is used in the order the
#' reported calibration was produced with (three sample.int calls, an optional shared
#' offset, then three normals), which is why the reported rows reproduce exactly.
#'
#' A pass is not a proof that the engine is right at K = 12; it is evidence that the
#' two-stage reduction and the tau prior together reproduce what the full model does
#' at the only K where both can be run. That is the strongest check available.
#'
#' @param pilot_draws List of three numeric vectors, per-SD posterior draws of the
#'   interaction for English, Turkish and Norwegian.
#' @param se Per-language SEs used to generate the estimates.
#' @param shared_sd SD of an offset common to the three languages; zero in the
#'   calibration as run, because the brms replicates gave each language its own verbs.
#' @param se_analysis SEs the meta-analysis is told, by default the generating ones.
#' @return One-row data frame of the moments the target is stated in.
calibrate_against_pooled_k3 <- function(pilot_draws, se, n_reps = 4000L, shared_sd = 0,
                                        se_analysis = se,
                                        prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                        seed = 20260910L) {
  stopifnot(length(pilot_draws) == 3L, length(se) == 3L, length(se_analysis) == 3L)
  set.seed(seed)
  Y <- matrix(NA_real_, n_reps, 3L)
  for (i in seq_len(n_reps)) {
    theta <- vapply(pilot_draws, function(v) v[sample.int(length(v), 1L)], numeric(1))
    common <- if (shared_sd > 0) rnorm(1, 0, shared_sd) else 0
    Y[i, ] <- theta + common + rnorm(3, 0, se)
  }
  p <- ra_meta_batch(Y, se_analysis, prior_tau_scale = prior_tau_scale)$p_mu_negative
  bf <- bf_from_p_negative(p)
  z <- z_from_p_negative(p)
  data.frame(tau_prior_scale = prior_tau_scale, mean_z = mean(z), sd_z = sd(z),
             p_bf10 = mean(bf >= 10), p_bf6 = mean(bf >= 6), p_bf3 = mean(bf >= 3),
             n = n_reps)
}

#' One row per three-language brms replicate, for the calibration target.
#'
#' @param path Either the directory of per-replicate .rds files pulled from ARC
#'   (outputs/design_pooled_v2, gitignored) or a CSV this function wrote before,
#'   which lets a machine without the .rds files reproduce the target.
read_pooled_k3_replicates <- function(path = "outputs/design_pooled_v2",
                                      hypothesis = "H1b_active_interaction_negative") {
  if (grepl("\\.csv$", path)) {
    if (!file.exists(path)) stop("[calibration] replicate table not found: ", path)
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }
  files <- sort(list.files(path, pattern = "\\.rds$", full.names = TRUE))
  if (!length(files)) {
    stop("[calibration] no per-replicate .rds files in ", path,
         ". They are gitignored; pull them from ARC or pass --pooled-replicates=<csv>.")
  }
  do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    sm <- x$summary
    row <- if (is.null(x$bf_results)) NULL else x$bf_results[x$bf_results$hypothesis == hypothesis, ]
    has_bf <- !is.null(row) && nrow(row) == 1L
    # A failed fit's summary has no model_level and its record has no diagnostics.
    data.frame(cell_id = sm$cell_id, n_participants = sm$n_participants, status = sm$status,
               model_level = if ("model_level" %in% names(sm)) sm$model_level else NA_character_,
               convergence_ok = if (is.null(x$diagnostics)) NA else x$diagnostics$convergence_ok,
               prior_prob = if (has_bf) row$prior_prob else NA_real_,
               posterior_prob = if (has_bf) row$posterior_prob else NA_real_,
               BF_10 = if (has_bf) row$BF_10 else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' The calibration target, computed from the replicates instead of typed in.
#'
#' Every replicate whose fit finished counts, scored exactly: z = qnorm of the
#' posterior probability of the predicted sign, and each Bayes factor compared with
#' its threshold unrounded, as the engine's own rows are.
#' @return List: target (the five moments), counts, and a by-sample-size table.
pooled_k3_target <- function(reps) {
  ok <- reps$status == "success" & is.finite(reps$posterior_prob)
  z <- z_from_p_negative(reps$posterior_prob[ok])
  bf <- reps$BF_10[ok]
  moments <- function(z, bf) {
    c(mean_z = mean(z), sd_z = stats::sd(z), p_bf10 = mean(bf >= 10),
      p_bf6 = mean(bf >= 6), p_bf3 = mean(bf >= 3))
  }
  n_all <- reps$n_participants
  by_n <- do.call(rbind, lapply(sort(unique(n_all)), function(N) {
    here <- ok & n_all == N
    data.frame(N = N, n_fits = sum(n_all == N), n_success = sum(here),
               n_error = sum(n_all == N & reps$status != "success"),
               p_bf10 = if (any(here)) mean(reps$BF_10[here] >= 10) else NA_real_)
  }))
  list(target = moments(z, bf), n_success = sum(ok), n_error = sum(reps$status != "success"),
       n_convergence_ok = sum(reps$convergence_ok[ok] %in% TRUE),
       mean_prior_prob = mean(reps$prior_prob[ok]), by_n = by_n)
}

# ---------------------------------------------------------------------------
# 6. The named-language arm: which languages join, not how many
# ---------------------------------------------------------------------------
#
# Everything above treats the K languages as exchangeable draws from a population,
# which is the only way to say anything about languages that have not been
# recruited. It is also the step a reviewer will press hardest on, because at
# K = 12 it means nine of the twelve languages are inventions.
#
# For K <= 7 that step is avoidable. The seven languages in hand can be enumerated
# directly: hold each one's interaction at its harmonised observed value, add only
# the sampling error a CLAPS team at 80 participants and 72 verbs would incur, and
# ask how often the pooled test succeeds. No synthetic language enters anywhere.
#
# The two arms answer different questions and must be printed side by side rather
# than blended, because the gap between them at K = 7 IS the fixed-versus-random
# language estimand difference, and it is large. The named arm is also the one that
# answers what the PI literally asked, since he named seven specific languages.

#' Power over every subset of the observed languages at a given size.
#'
#' @param pop Data frame with language, y (harmonised effect) and s columns.
#' @param K Subset size; must not exceed nrow(pop).
#' @param se_indep,se_shared Per-language sampling error, split as elsewhere.
#' @return Data frame, one row per subset, plus the languages it omits.
enumerate_language_subsets <- function(pop, K, se_indep, se_shared,
                                       n_reps = 4000L, bf_threshold = 10,
                                       prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                       seed = 1L) {
  stopifnot(K <= nrow(pop))
  set.seed(seed)
  combs <- utils::combn(nrow(pop), K, simplify = FALSE)
  do.call(rbind, lapply(combs, function(ix) {
    theta <- pop$y[ix]
    s     <- rep(se_indep, K)
    bf <- vapply(seq_len(n_reps), function(i) {
      y <- theta + rnorm(1, 0, se_shared) + rnorm(K, 0, se_indep)
      directional_bf_meta(ra_meta_posterior(y, s, prior_tau_scale = prior_tau_scale))
    }, numeric(1))
    p <- mean(bf >= bf_threshold)
    data.frame(K = K,
               languages = paste(pop$language[ix], collapse = "+"),
               omitted   = paste(setdiff(pop$language, pop$language[ix]), collapse = "+"),
               realised_tau = if (K > 1) sd(theta) else NA_real_,
               mean_effect  = mean(theta),
               power = p, mcse = sqrt(p * (1 - p) / n_reps))
  }))
}

#' Where the variation in the pooled outcome comes from, at a fixed number of languages.
#'
#' A study's success depends on three nested things, and the law of total variance
#' splits the Bernoulli variance pbar (1 - pbar) of its outcome among them:
#'
#'   population   Var over (mu, tau) of the expected power given (mu, tau): our
#'                uncertainty about how large and how variable the effect is across
#'                languages, which no recruitment rule touches.
#'   composition  the expected variance, given (mu, tau), of a set's power across the
#'                sets of K languages that could join: which languages they are.
#'   noise        the expected p (1 - p) within a set: sampling error in the study.
#'
#' It is estimated from a nested simulation, n_par draws of (mu, tau), n_sets
#' language sets under each and n_reps studies of each set, with the finite-sample
#' corrections that make each component unbiased. An earlier version drew one set per
#' (mu, tau), so its "which languages" share also held all of the population term, and
#' it drew (mu, tau) under the analysis prior instead of the design prior (audit items
#' klanguage-engine#2 and #3). K is fixed within a call, so none of this measures how
#' much the NUMBER of languages matters; comparing calls across K does that.
#'
#' @param pop The language set defining the population, drawn under the DESIGN prior.
#' @param n_boot Bootstrap resamples of the (mu, tau) draws, for the 90% intervals.
#' @return One row: mean power and its Monte Carlo SE, the three shares with 90%
#'   bootstrap intervals, and share_true_effects (population plus composition),
#'   the quantity the earlier version reported under the composition label.
decompose_outcome_variance <- function(K, pop, se_indep, se_shared,
                                       n_par = 1000L, n_sets = 10L, n_reps = 100L,
                                       bf_threshold = 10,
                                       prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                       design_tau_prior_scale = DESIGN_TAU_PRIOR_SCALE,
                                       n_boot = 200L, seed = 1L) {
  stopifnot(n_par >= 2L, n_sets >= 2L, n_reps >= 2L)
  set.seed(seed)
  pars <- draw_population_params(n_par, pop$y, pop$s, prior_tau_scale = design_tau_prior_scale)
  s <- rep(se_indep, K)
  P <- matrix(NA_real_, n_par, n_sets)
  for (a in seq_len(n_par)) {
    Y <- do.call(rbind, lapply(seq_len(n_sets), function(b) {
      theta <- rnorm(K, pars$mu[a], pars$tau[a])          # the set, drawn once
      common <- rnorm(n_reps, 0, se_shared)
      matrix(theta, n_reps, K, byrow = TRUE) + common +
        matrix(rnorm(n_reps * K, 0, se_indep), n_reps, K)
    }))
    p_neg <- ra_meta_batch(Y, s, prior_tau_scale = prior_tau_scale, chunk = nrow(Y))$p_mu_negative
    hit <- bf_from_p_negative(p_neg) >= bf_threshold
    P[a, ] <- colMeans(matrix(hit, n_reps, n_sets))
  }

  shares <- function(P) {
    pbar <- mean(P)
    total <- pbar * (1 - pbar)
    noise <- mean(P * (1 - P)) * n_reps / (n_reps - 1)
    comp <- max(mean(apply(P, 1L, stats::var)) - noise / n_reps, 0)
    popn <- max(stats::var(rowMeans(P)) - comp / n_sets - noise / (n_reps * n_sets), 0)
    if (total <= 0) {
      return(c(population = NA_real_, composition = NA_real_, noise = NA_real_, given = NA_real_))
    }
    c(population = popn / total, composition = comp / total, noise = noise / total,
      given = if (comp + noise > 0) comp / (comp + noise) else NA_real_)
  }
  est <- shares(P)
  boot <- vapply(seq_len(n_boot), function(r) {
    shares(P[sample.int(n_par, n_par, replace = TRUE), , drop = FALSE])
  }, numeric(4))
  lo <- apply(boot, 1L, stats::quantile, probs = 0.05, na.rm = TRUE, names = FALSE)
  hi <- apply(boot, 1L, stats::quantile, probs = 0.95, na.rm = TRUE, names = FALSE)

  data.frame(K = K, mean_power = mean(P), mcse_power = stats::sd(rowMeans(P)) / sqrt(n_par),
             share_population_params = est[["population"]],
             share_composition = est[["composition"]],
             share_noise = est[["noise"]],
             share_true_effects = est[["population"]] + est[["composition"]],
             share_composition_given_population = est[["given"]],
             population_lo = lo[1], population_hi = hi[1],
             composition_lo = lo[2], composition_hi = hi[2],
             noise_lo = lo[3], noise_hi = hi[3],
             n_par = n_par, n_sets = n_sets, n_reps = n_reps, n_boot = n_boot,
             prior_tau_scale = prior_tau_scale, design_tau_prior_scale = design_tau_prior_scale,
             bf_threshold = bf_threshold, seed = seed)
}

# ---------------------------------------------------------------------------
# 7. Partly known language sets
# ---------------------------------------------------------------------------
#
# The exchangeable arm treats every one of the K languages as an unknown draw from
# the population. That is right for a language nobody has measured. It is too
# pessimistic for a language we have already piloted, and the CLAPS teams for
# English, Turkish and Norwegian are the likeliest members of any Stage 2 set.
#
# What is known about those languages is their TRUE EFFECT, not their data. No
# pilot observation enters the Stage 2 dataset, and none enters here either: what
# the pilot supplies is a design prior on where that language's effect sits, which
# is what the rest of this package already uses it for. A recruited language whose
# effect is pinned down in advance contributes a known quantity to the pooled mean
# instead of a draw from a distribution with SD tau, and with tau three to five
# times the per-language standard error that is a real gain.
#
# The distinction is visible in the code below and is worth stating, because it is
# the thing that makes this arm legitimate. An already-measured language enters the
# GENERATIVE step with its effect drawn from the pilot posterior instead of from
# N(mu, tau). It enters the ANALYSIS step exactly like every other language: a fresh
# estimate with the same standard error, meta-analysed with a flat prior on mu. The
# earlier measurement never appears in the likelihood. So the gain here is a
# reduction in the variance of the truth being estimated, not prior information
# smuggled into the test.
#
# The languages measured only under the Ambridge, Arnon & Bekman (2023) protocol
# are a weaker case, and deliberately so. If a Hebrew or a Mandarin team joins
# CLAPS, the 2023 estimate for that language is evidence about what CLAPS would
# find there only to the extent that the protocol transfers, which is the open
# question of section 2 of the note. `transfer_sd` widens the prior on such a
# language by an amount representing that doubt; at transfer_sd equal to tau it
# is worth nothing at all, and the language is back to being an unknown draw.
# A CLAPS pilot language takes no such widening, so the widening is set row by row.

#' Expand transfer_sd to one value per known language, refusing to recycle a doubt.
#'
#' A single non-zero value would widen every known language, the CLAPS pilot ones
#' included, which is the defect that put partly-known rows c to e out by up to
#' seven points (audit item klanguage-engine#1). So a non-zero doubt must be given
#' per row, usually by transfer_sd_by_source(); a single zero means none anywhere.
#' With one known language a single value is already one per row, so this function
#' cannot tell it from a recycled scalar. run_partly_known_cell() therefore also
#' checks the rows themselves, and refuses any doubt on a CLAPS pilot row.
per_row_transfer_sd <- function(transfer_sd, m) {
  if (m == 0L) {
    stopifnot(length(transfer_sd) <= 1L)
    return(numeric(0))
  }
  if (length(transfer_sd) == 1L && transfer_sd == 0) {
    return(rep(0, m))
  }
  if (length(transfer_sd) != m) {
    stop("[transfer] transfer_sd needs one value per known language (", m, "), got ",
         length(transfer_sd), ". A single non-zero value would also widen the CLAPS pilot ",
         "languages; build the vector with transfer_sd_by_source().")
  }
  stopifnot(all(is.finite(transfer_sd)), all(transfer_sd >= 0))
  transfer_sd
}

#' Power when some, or all, of the K languages have already been measured.
#'
#' @param K Total number of languages contributing.
#' @param known Data frame of the already-measured languages that will take part,
#'   with columns y (their harmonised effect) and s (its standard error). At most K
#'   rows; the remaining K - nrow(known) places are drawn from the population.
#' @param pop The anchor set defining the population the unknown places are drawn
#'   from, as elsewhere. Not consulted when every place is known.
#' @param mode "assurance" draws (mu, tau) for the unknown places from pop's
#'   posterior under the design prior, once per replicate; "point" holds them at
#'   mu_point and tau_point.
#' @param transfer_sd Extra SD added to a known language's prior, representing
#'   doubt that its earlier measurement transfers to the CLAPS protocol: one value
#'   per row of `known` (see per_row_transfer_sd), zero for a CLAPS pilot language.
#' @param known_truth "drawn" takes each known language's true effect from
#'   N(y, s^2 + transfer_sd^2), once per replicate, which is how the pilot supplies a
#'   design prior. "fixed" holds it at y exactly, the literal reading of "kept at
#'   their estimated effects"; there is then no draw for a transfer doubt to widen.
#' @return One row. Its first eleven columns are those of the first reported table;
#'   transfer_sd there is the largest per-row value and n_widened counts the rows.
run_partly_known_cell <- function(K, known, n_reps, se_lang,
                                  pop = NULL, shared_verb_sd = 0,
                                  transfer_sd = 0,
                                  mode = c("assurance", "point"),
                                  mu_point = NULL, tau_point = NULL,
                                  known_truth = c("drawn", "fixed"),
                                  prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                  design_tau_prior_scale = DESIGN_TAU_PRIOR_SCALE,
                                  bf_threshold = 10, seed = 1L) {
  mode <- match.arg(mode)
  known_truth <- match.arg(known_truth)
  m <- nrow(known)
  n_new <- K - m
  stopifnot(m <= K, length(se_lang) == 1L)
  tsd <- per_row_transfer_sd(transfer_sd, m)
  # The length check above lets a lone value through when exactly one language is
  # known, so a CLAPS pilot language on its own could still be widened. Wherever the
  # rows say where they came from, the rule is checked on the rows.
  if (!is.null(known$source) && any(tsd[known$source == "CLAPS_pilot"] > 0)) {
    stop("[transfer] a CLAPS pilot language was measured under the Stage 2 protocol and ",
         "takes no transfer doubt; build transfer_sd with transfer_sd_by_source()")
  }
  if (known_truth == "fixed" && any(tsd > 0)) {
    stop("[transfer] a transfer doubt widens a drawn effect; with known_truth = 'fixed' there is none")
  }
  set.seed(seed)
  # With every place known nothing is drawn from the population, so it is not
  # consulted at all. Drawing it anyway would leave the answer unchanged in law but
  # not in practice: the draw uses up a different amount of the random stream for
  # different populations, and cells that must agree would differ by Monte Carlo
  # noise.
  pars <- NULL
  if (n_new > 0L && mode == "assurance") {
    stopifnot(!is.null(pop))
    pars <- draw_population_params(n_reps, pop$y, pop$s,
                                   prior_tau_scale = design_tau_prior_scale)
  } else if (n_new > 0L) {
    stopifnot(!is.null(mu_point), !is.null(tau_point))
    pars <- data.frame(mu = rep(mu_point, n_reps), tau = rep(tau_point, n_reps))
  }
  sd_known <- sqrt(known$s^2 + tsd^2)

  # Generation stays one replicate at a time so the random stream is consumed in the
  # order the reported rows were produced with; only the analysis is batched.
  Y <- matrix(NA_real_, n_reps, K)
  for (i in seq_len(n_reps)) {
    theta_known <- if (m == 0L) {
      numeric(0)
    } else if (known_truth == "fixed") {
      known$y
    } else {
      rnorm(m, known$y, sd_known)
    }
    theta_new <- if (n_new > 0L) rnorm(n_new, pars$mu[i], pars$tau[i]) else numeric(0)
    Y[i, ] <- c(theta_known, theta_new) + rnorm(1, 0, shared_verb_sd) + rnorm(K, 0, se_lang)
  }
  p_neg <- ra_meta_batch(Y, se_lang, prior_tau_scale = prior_tau_scale)$p_mu_negative
  bf <- bf_from_p_negative(p_neg)
  z <- z_from_p_negative(p_neg)

  p <- mean(bf >= bf_threshold)
  ci <- binom.test(round(p * n_reps), n_reps)$conf.int
  data.frame(K = K, n_known = m, transfer_sd = if (m > 0L) max(tsd) else 0, n_reps = n_reps,
             power = p, mcse = sqrt(p * (1 - p) / n_reps),
             power_lo = ci[1], power_hi = ci[2],
             mean_z = mean(z), median_bf = median(bf), seed = seed,
             n_widened = sum(tsd > 0), mode = mode, known_truth = known_truth,
             mu_point = if (mode == "point") mu_point %||% NA_real_ else NA_real_,
             tau_point = if (mode == "point") tau_point %||% NA_real_ else NA_real_)
}
