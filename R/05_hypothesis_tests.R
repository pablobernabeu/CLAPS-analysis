# R/05_hypothesis_tests.R
#
# Purpose
#   Convert a fitted brms model into Bayes factors for the CLAPS focal
#   hypotheses. This file is the only place where evidence is quantified; the
#   simulation, ladder and summary modules all route through compute_all_bf(),
#   so a change to the testing rule takes effect everywhere at once.
#
# Entry points
#   compute_all_bf(fit)               All focal hypotheses for one fitted model.
#   savage_dickey_directional_bf(...) One hypothesis, one coefficient.
#   directional_bf_from_draws(...)    The arithmetic alone, for testing.
#   classify_bf(bf)                   Evidence category from a BF value.
#
# Why Savage-Dickey rather than bridge sampling
#   Every hypothesis here is nested: H0 fixes one coefficient to zero while the
#   rest of the model is unchanged, which is exactly the case the Savage-Dickey
#   density ratio is derived for (Verdinelli & Wasserman, 1995,
#   doi:10.1080/01621459.1995.10476554; the psychology-facing tutorial is
#   Wagenmakers et al., 2010, doi:10.1016/j.cogpsych.2009.12.001). It is computed
#   from draws the model has already produced, so it costs nothing beyond the fit,
#   which matters when the design analysis fits thousands of models. Bridge
#   sampling estimates the marginal likelihood directly and so does not need
#   nesting, but it requires refitting under each hypothesis and is markedly more
#   expensive. It is therefore kept only as an occasional calibration check
#   (scripts/05_bf_calibration_cell.R) on a subset of cells, to confirm the
#   cheaper route agrees with it, and never as the primary route.
#
# Requirement
#   Every fit passed in must have been sampled with sample_prior = "yes". The
#   directional Savage-Dickey ratio compares a posterior probability with the
#   corresponding prior probability, so without stored prior draws it cannot be
#   formed. The functions below raise an error rather than substituting an
#   analytic approximation.
#
# PREREGISTERED TEST DIRECTIONS, lettered as the coefficients enter the model:
#   H1a — ONE-TAILED, positive direction:
#         Semantics > 0 for passives, the reference level.
#         Supporting prediction, held as a prespecified positive control.
#
#   H1b — ONE-TAILED, negative direction:
#         S_TypeActive:Semantics < 0  (smaller affectedness slope for actives than passives).
#         Primary prediction; the confirmatory decision rests on this one alone.
#
#   H2  — TWO-TAILED (both directions reported as H2a + H2b):
#         S_TypePseudo_Passive:Semantics.
#         Secondary prediction; Turkish pilot showed larger effect for pseudo-passives
#         (opposite direction).

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# ---------------------------------------------------------------------------
# Savage-Dickey Bayes factors (primary route)
# ---------------------------------------------------------------------------
# The Savage-Dickey density ratio is valid for nested hypotheses where the
# restricted model is obtained by setting a parameter to a point value
# (H0: theta = 0) and the prior is continuous at that point. Both conditions hold
# here: the priors in R/03_define_priors.R are normal or Student-t, hence
# continuous everywhere, and no hypothesis fixes a variance component to a
# boundary.
#
# For directional hypotheses the one-sided variant is used. It compares the
# probability that the coefficient has the predicted sign before and after seeing
# the data, as an odds ratio:
#
#   BF_10 = [ p(theta > 0 | data)  / p(theta < 0 | data)  ]
#           -------------------------------------------------
#           [ p(theta > 0 | prior) / p(theta < 0 | prior) ]
#
# Reading it as an odds ratio explains the reference-prior division: a prior that
# already places most of its mass on the predicted side must not be credited as
# evidence. A symmetric prior gives prior odds of 1 and the expression reduces to
# the posterior odds, but the sensitivity regimes include deliberately asymmetric
# priors, so the denominator is not redundant.
#
# The equivalent test on the *unsigned* coefficient is not used, because a
# two-sided BF near 1 conflates "no effect" with "an effect in the unpredicted
# direction", and the CLAPS predictions are directional.
#
# Requires sample_prior = "yes" in brms, which stores the prior draws the
# denominator needs.

#' Directional Savage-Dickey Bayes factor from draws.
#'
#' The arithmetic core of savage_dickey_directional_bf(), split out so it can be
#' exercised on plain numeric vectors without fitting a model. Given posterior
#' and prior draws of one coefficient, it compares the probability that the
#' coefficient has the predicted sign under the posterior with the same
#' probability under the prior, as an odds ratio.
#'
#' @param post_vals,prior_vals Numeric draws of the coefficient.
#' @param direction "positive" or "negative"; the sign the hypothesis predicts.
#' @return List with posterior_prob, prior_prob, bf_10 and log_bf.
directional_bf_from_draws <- function(post_vals, prior_vals,
                                      direction = "positive") {
  stopifnot(direction %in% c("positive", "negative"))

  if (direction == "positive") {
    post_prob  <- mean(post_vals  > 0)
    prior_prob <- mean(prior_vals > 0)
  } else {
    post_prob  <- mean(post_vals  < 0)
    prior_prob <- mean(prior_vals < 0)
  }

  # Clamp away from 0 and 1, where the odds ratio would be infinite. With the
  # draw counts used here the clamp binds only when every draw falls on one
  # side, which is itself a signal that the sample is too small to resolve.
  post_prob  <- pmin(pmax(post_prob,  1e-6), 1 - 1e-6)
  prior_prob <- pmin(pmax(prior_prob, 1e-6), 1 - 1e-6)

  bf_10 <- (post_prob / (1 - post_prob)) / (prior_prob / (1 - prior_prob))
  list(posterior_prob = post_prob, prior_prob = prior_prob,
       bf_10 = bf_10, log_bf = log(bf_10))
}

#' Compute the Savage-Dickey BF for a directional hypothesis from a fitted model.
#'
#' @param fit A brmsfit with sample_prior = "yes". Without that argument brms
#'   stores no prior draws and the ratio cannot be formed, which is why the
#'   missing-prior case below is an error rather than a fallback.
#' @param param Character; Stan parameter name (e.g. "b_Semantics_scaled").
#' @param direction "positive" or "negative"; the direction H1 predicts.
#' @return One-row tibble with BF_10, BF_01, posterior_prob, prior_prob,
#'   log_BF10, the evidence category and the method label.
#' @details A thin wrapper: it locates the posterior and prior draws for one
#'   coefficient, hands the arithmetic to directional_bf_from_draws(), and labels
#'   the result. Both error messages name what was looked for and what was
#'   available, because a wrong parameter name is by far the commonest way to
#'   misuse this function and the brms naming scheme is not obvious.
savage_dickey_directional_bf <- function(fit, param, direction = "positive") {
  stopifnot(direction %in% c("positive", "negative"))

  # brms stores prior draws under the parameter name prefixed with "prior_".
  post_draws  <- brms::as_draws_df(fit)
  prior_draws <- brms::as_draws_df(fit, variable = paste0("prior_", param))

  if (!param %in% names(post_draws)) {
    stop("[BF] Parameter '", param, "' not found in posterior draws. ",
         "Available: ", paste(head(names(post_draws), 20), collapse = ", "))
  }
  prior_param <- paste0("prior_", param)
  if (!prior_param %in% names(prior_draws)) {
    stop("[BF] '", prior_param, "' not found. Was sample_prior = 'yes' set?")
  }

  post_vals  <- post_draws[[param]]
  prior_vals <- prior_draws[[prior_param]]

  bf <- directional_bf_from_draws(post_vals, prior_vals, direction)

  tibble::tibble(
    param          = param,
    direction      = direction,
    posterior_prob = bf$posterior_prob,
    prior_prob     = bf$prior_prob,
    BF_10          = bf$bf_10,
    BF_01          = 1 / bf$bf_10,
    log_BF10       = bf$log_bf,
    bf_category    = classify_bf(bf$bf_10),
    method         = "savage_dickey_directional"
  )
}

#' Classify a Bayes factor using the prespecified CLAPS thresholds.
#' Primary threshold: BF > 10 (strong evidence for H1).
#' Secondary threshold: BF > 3 (moderate evidence).
classify_bf <- function(bf_10) {
  # Conventional evidence bands (Lee & Wagenmakers, 2013). A Bayes factor of
  # exactly 1 is no evidence either way; the >= ordering below places that
  # boundary case on the H1 side, so the bands are half-open upwards.
  dplyr::case_when(
    bf_10 >= 100    ~ "extreme_H1",
    bf_10 >= 30     ~ "very_strong_H1",
    bf_10 >= 10     ~ "strong_H1",
    bf_10 >= 3      ~ "moderate_H1",
    bf_10 >= 1      ~ "anecdotal_H1",
    bf_10 >= 1 / 3  ~ "anecdotal_H0",
    bf_10 >= 1 / 10 ~ "moderate_H0",
    TRUE            ~ "strong_or_more_H0"
  )
}

#' Compute all focal hypothesis tests for a fitted CLAPS model.
#'
#' @param fit A brmsfit object with sample_prior = "yes".
#' @param has_pseudo_passive Logical; FALSE for languages lacking the
#'   Pseudo_Passive level (Norwegian, Balinese per config/analysis_config.yaml),
#'   where the H2 rows are simply absent from the result rather than NA.
#' @param semantics_var Name of the affectedness predictor as it appears in the
#'   model formula, used to build the Stan coefficient names.
#' @return Tibble with one row per hypothesis, identified by the `hypothesis`
#'   column: H1a_semantics_positive, H1b_active_interaction_negative and, when
#'   pseudo-passives are present, H2a_pseudo_positive and H2b_pseudo_negative.
#' @details H2 contributes two rows because it is preregistered as two-tailed, and
#'   a two-tailed directional Savage-Dickey test is reported as the pair of
#'   one-tailed tests rather than collapsed into one number. The two are
#'   reciprocal by construction, so reporting both makes the direction of the
#'   evidence explicit instead of leaving it implied by a threshold.
#'
#'   H1a is a context check rather than a claim under test on its own: it
#'   confirms that affectedness predicts acceptability in the passive, the
#'   reference level against which the H1b interaction is interpreted.
compute_all_bf <- function(fit, has_pseudo_passive = TRUE,
                           semantics_var = "Semantics_scaled") {
  # Population-level coefficients carry a "b_" prefix in brms.
  s_param <- paste0("b_", semantics_var)
  # The interaction separator is matched rather than assumed: brms writes ":" but
  # some downstream draw containers report "." instead, so the patterns below
  # accept either via the [:.] class rather than hard-coding one.
  post_names <- names(brms::as_draws_df(fit))
  # Match the gender-AVERAGED focal interactions ONLY. The end-anchor ($) excludes
  # the higher-order gender terms (e.g. b_S_TypeActive:Semantics_scaled:Gender1) that
  # gender_spec = "interaction" adds; the old unanchored ".*" pattern also matched
  # those and [1] returned an order-dependent wrong term. The length()==1 guard turns
  # any future coefficient-naming drift into a loud failure instead of a silent one.
  ia_active_hits   <- grep(paste0("^b_S_TypeActive[:.]", semantics_var, "$"),
                           post_names, value = TRUE)
  if (length(ia_active_hits) != 1L) {
    stop("[BF] expected exactly one Active:Semantics coefficient, found ",
         length(ia_active_hits), ": ", paste(ia_active_hits, collapse = ", "))
  }
  ia_active_param  <- ia_active_hits[1]
  ia_pseudo_param  <- if (has_pseudo_passive) {
    ia_pseudo_hits <- grep(paste0("^b_S_TypePseudo_Passive[:.]", semantics_var, "$"),
                           post_names, value = TRUE)
    if (length(ia_pseudo_hits) != 1L) {
      stop("[BF] expected exactly one Pseudo_Passive:Semantics coefficient, found ",
           length(ia_pseudo_hits), ": ", paste(ia_pseudo_hits, collapse = ", "))
    }
    ia_pseudo_hits[1]
  } else {
    NA_character_
  }

  results <- list()

  # H1a: Semantics > 0 for passives (reference level)
  if (!is.na(s_param) && s_param %in% post_names) {
    results[["H1a_semantics_positive"]] <-
      savage_dickey_directional_bf(fit, s_param, "positive")
  }

  # H1b: Active:Semantics < 0  — ONE-TAILED (negative).
  # Preregistered direction: smaller affectedness slope for actives than passives.
  if (!is.na(ia_active_param)) {
    results[["H1b_active_interaction_negative"]] <-
      savage_dickey_directional_bf(fit, ia_active_param, "negative")
  }

  # H2: Pseudo_Passive:Semantics — TWO-TAILED: report both directions.
  # Secondary prediction; Turkish pilot showed opposite direction to prediction.
  if (has_pseudo_passive && !is.na(ia_pseudo_param)) {
    results[["H2a_pseudo_positive"]] <-
      savage_dickey_directional_bf(fit, ia_pseudo_param, "positive")
    results[["H2b_pseudo_negative"]] <-
      savage_dickey_directional_bf(fit, ia_pseudo_param, "negative")
  }

  dplyr::bind_rows(results, .id = "hypothesis")
}

#' Select and round the BF columns for presentation in the report.
#'
#' @param bf_df Output of compute_all_bf().
#' @return The same rows with the reporting columns only, rounded.
#' @details Rounding happens here, at the point of display, and not in
#'   compute_all_bf(), so that stored results keep full precision for any later
#'   re-aggregation. Posterior probabilities are kept to four decimal places
#'   rather than three because they sit close to 1 for well-supported hypotheses,
#'   where the third decimal alone would hide the difference between 0.9991 and
#'   0.9999.
summarise_bf_table <- function(bf_df) {
  bf_df |>
    dplyr::select(hypothesis, param, direction, posterior_prob,
                  BF_10, BF_01, bf_category, method) |>
    dplyr::mutate(
      BF_10 = round(BF_10, 3),
      BF_01 = round(BF_01, 3),
      posterior_prob = round(posterior_prob, 4)
    )
}
