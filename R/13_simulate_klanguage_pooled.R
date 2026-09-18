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
#   P(BF >= 10) = 0.565 against the 0.578 actually observed. The whole operating
#   characteristic of a sixty-five-hour fit is therefore carried by two numbers, and
#   an engine that reproduces those two numbers at K = 3 has earned the right to be
#   run at K = 12. calibrate_against_pooled_k3() in this file is that check.
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
#   0.39 with a 95% interval of roughly [0.20, 0.77]. Every entry point below
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
# fitted to the three-language runs. With the per-language SEs held at the
# verb-information floors the per-language design analysis already establishes
# (0.104, 0.076, 0.042 per SD for English, Turkish and Norwegian), a half-normal
# scale of 1.1 reproduces every moment of the brms runs:
#
#            mean z   sd z   P(BF>=10)  P(BF>=6)  P(BF>=3)
#   brms      1.386  0.311     0.578     0.859     0.985
#   scale 1.0 1.429  0.311     0.604     0.874     0.992
#   scale 1.2 1.389  0.308     0.556     0.845     0.990
#
# Two things are worth saying plainly about this. The per-language SE was fixed in
# advance from the per-language design analysis rather than tuned, and the value the
# calibration would have chosen for it independently is 0.06 to 0.08, against the
# 0.074 the design analysis gives. So the two routes agree on the input that was not
# fitted, which is the reason to believe the one that was.
#
# And holding the PRIOR fixed while K grows is the conservative choice that also
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
  if (length(tau_grid) == 1L) {
    trap <- 1
  } else {
    dx   <- diff(tau_grid)
    trap <- c(dx[1], (utils::head(dx, -1) + utils::tail(dx, -1)), utils::tail(dx, 1)) / 2
  }
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
#'   three-language runs assume and is how the calibration check is set up.
#' @return One-row data frame with the estimates, the Bayes factor and z.
simulate_klanguage_replicate <- function(K, mu, tau, se_lang,
                                         shared_verb_sd = 0,
                                         prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                         prior_mu_sd = Inf,
                                         theta_fixed = NULL) {
  s <- if (length(se_lang) == 1L) rep(se_lang, K) else se_lang
  stopifnot(length(s) == K)

  theta <- if (!is.null(theta_fixed)) {
    stopifnot(length(theta_fixed) == K); theta_fixed
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

  params <- switch(mode,
    point = {
      stopifnot(!is.null(mu_point), !is.null(tau_point))
      data.frame(mu = rep(mu_point, n_reps), tau = rep(tau_point, n_reps))
    },
    assurance = {
      stopifnot(!is.null(pop))
      draw_population_params(n_reps, pop$y, pop$s,
                             prior_tau_scale = design_tau_prior_scale,
                             prior_mu_sd = prior_mu_sd)
    },
    safeguard = {
      stopifnot(!is.null(pop))
      d <- draw_population_params(20000L, pop$y, pop$s,
                                  prior_tau_scale = design_tau_prior_scale,
                                  prior_mu_sd = prior_mu_sd)
      # mu towards zero, tau towards its upper tail: the unfavourable corner.
      mu_sg  <- quantile(d$mu,  1 - safeguard_q[1], names = FALSE)  # least negative
      tau_sg <- quantile(d$tau, safeguard_q[2],     names = FALSE)
      data.frame(mu = rep(mu_sg, n_reps), tau = rep(tau_sg, n_reps))
    })

  reps <- do.call(rbind, lapply(seq_len(n_reps), function(i)
    simulate_klanguage_replicate(K, params$mu[i], params$tau[i], se_lang,
                                 shared_verb_sd = shared_verb_sd,
                                 prior_tau_scale = prior_tau_scale,
                                 prior_mu_sd = prior_mu_sd)))

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
    # Type S: the sign is wrong given that the threshold was cleared in EITHER
    # direction. Reported because a directional Bayes factor can be met by an
    # effect that runs the wrong way, and Balinese shows that is not hypothetical
    # (Gelman & Carlin, 2014, doi:10.1177/1745691614551642).
    p_wrong_sign = mean(reps$p_negative < 0.5),
    seed = seed
  )
}

# ---------------------------------------------------------------------------
# 5. The calibration check against the completed three-language brms runs
# ---------------------------------------------------------------------------

#' Does the cheap engine reproduce the expensive one at K = 3?
#'
#' The three-language runs hold the languages fixed and draw each language's true
#' effect from its own pilot posterior, so the check mirrors that: theta is fixed
#' at the three pilot interaction values on the common scale, and the engine is run
#' with the per-language SEs those fits achieve. The target, measured over the 135
#' usable replicates in outputs/design_pooled_v2, is mean z 1.386, sd z 0.311 and
#' P(BF >= 10) 0.578.
#'
#' A pass is not a proof that the engine is right at K = 12; it is evidence that the
#' two-stage reduction and the tau prior together reproduce what the full model does
#' at the only K where both can be run. That is the strongest check available.
calibrate_against_pooled_k3 <- function(theta_claps, se_claps, n_reps = 20000L,
                                        prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                        seed = 20260910L) {
  set.seed(seed)
  reps <- do.call(rbind, lapply(seq_len(n_reps), function(i)
    simulate_klanguage_replicate(3L, mu = mean(theta_claps), tau = sd(theta_claps),
                                 se_lang = se_claps,
                                 prior_tau_scale = prior_tau_scale,
                                 theta_fixed = theta_claps)))
  list(
    engine = c(mean_z = mean(reps$z), sd_z = sd(reps$z),
               p_bf10 = mean(reps$bf >= 10), p_bf6 = mean(reps$bf >= 6),
               p_bf3 = mean(reps$bf >= 3)),
    target = c(mean_z = 1.386, sd_z = 0.311, p_bf10 = 0.578, p_bf6 = NA, p_bf3 = NA),
    reps = reps
  )
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

#' How much of the outcome depends on WHICH languages join rather than how many.
#'
#' Draws language sets, computes each set's realised success probability, and
#' returns the share of the total outcome variance attributable to the set:
#'
#'     Var_set(p) / [ Var_set(p) + E_set(p (1 - p)) ]
#'
#' The first term is variation between language sets; the second is the residual
#' Bernoulli variation within a set. A value near zero would mean the identity of
#' the recruited languages hardly matters and a count-based rule is adequate. A
#' value that stays high as K grows means the opposite, and means a rule phrased
#' purely in terms of K is answering the wrong question.
realised_set_variance_share <- function(K, pop, se_indep, se_shared,
                                        n_sets = 400L, n_reps_per_set = 200L,
                                        bf_threshold = 10,
                                        prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                        seed = 1L) {
  set.seed(seed)
  pars <- draw_population_params(n_sets, pop$y, pop$s, prior_tau_scale = prior_tau_scale)
  ps <- vapply(seq_len(n_sets), function(j) {
    theta <- rnorm(K, pars$mu[j], pars$tau[j])   # the set, drawn once
    s <- rep(se_indep, K)
    bf <- vapply(seq_len(n_reps_per_set), function(i) {
      y <- theta + rnorm(1, 0, se_shared) + rnorm(K, 0, se_indep)
      directional_bf_meta(ra_meta_posterior(y, s, prior_tau_scale = prior_tau_scale))
    }, numeric(1))
    mean(bf >= bf_threshold)
  }, numeric(1))
  v_between <- var(ps)
  v_within  <- mean(ps * (1 - ps))
  data.frame(K = K, mean_power = mean(ps),
             var_between_sets = v_between, mean_within_var = v_within,
             share_due_to_which_languages = v_between / (v_between + v_within),
             p10 = unname(quantile(ps, 0.10)), p90 = unname(quantile(ps, 0.90)))
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

#' Power when some of the K languages have already been measured.
#'
#' @param K Total number of languages contributing.
#' @param known Data frame of the already-measured languages that will take part,
#'   with columns y (their harmonised effect) and s (its standard error). Fewer
#'   than K rows; the remaining K - nrow(known) places are drawn from the population.
#' @param pop The anchor set defining the population the unknown places are drawn
#'   from, as elsewhere.
#' @param transfer_sd Extra SD added to a known language's prior, representing
#'   doubt that its earlier measurement transfers to the CLAPS protocol. Zero for a
#'   CLAPS pilot language, positive for one measured under another protocol.
#' @return One row, in the same shape as run_klanguage_cell().
run_partly_known_cell <- function(K, known, n_reps, se_lang,
                                  pop = NULL, shared_verb_sd = 0,
                                  transfer_sd = 0,
                                  prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE,
                                  design_tau_prior_scale = DESIGN_TAU_PRIOR_SCALE,
                                  bf_threshold = 10, seed = 1L) {
  m <- nrow(known)
  stopifnot(m <= K, !is.null(pop))
  set.seed(seed)
  pars <- draw_population_params(n_reps, pop$y, pop$s,
                                 prior_tau_scale = design_tau_prior_scale)
  s <- rep(se_lang, K)

  reps <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
    theta_known <- if (m > 0) {
      rnorm(m, known$y, sqrt(known$s^2 + transfer_sd^2))
    } else numeric(0)
    theta_new <- if (K > m) rnorm(K - m, pars$mu[i], pars$tau[i]) else numeric(0)
    theta <- c(theta_known, theta_new)
    y <- theta + rnorm(1, 0, shared_verb_sd) + rnorm(K, 0, se_lang)
    post <- ra_meta_posterior(y, s, prior_tau_scale = prior_tau_scale)
    data.frame(bf = directional_bf_meta(post),
               z  = qnorm(min(max(post$p_mu_negative, 1e-12), 1 - 1e-12)))
  }))

  p <- mean(reps$bf >= bf_threshold)
  data.frame(K = K, n_known = m, transfer_sd = transfer_sd, n_reps = n_reps,
             power = p, mcse = sqrt(p * (1 - p) / n_reps),
             power_lo = binom.test(round(p * n_reps), n_reps)$conf.int[1],
             power_hi = binom.test(round(p * n_reps), n_reps)$conf.int[2],
             mean_z = mean(reps$z), median_bf = median(reps$bf), seed = seed)
}
