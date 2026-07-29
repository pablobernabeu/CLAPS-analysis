# R/05_hypothesis_tests.R
# Bayes-factor hypothesis tests for the CLAPS focal estimands.
# Uses Savage-Dickey density ratio for nested directional hypotheses where valid.
# Bridge sampling retained only as calibration check when explicitly requested.
#
# PREREGISTERED TEST DIRECTIONS:
#   H1 — ONE-TAILED, negative direction:
#         S_TypeActive:Semantics < 0  (smaller affectedness slope for actives than passives).
#
#   H2 — TWO-TAILED (both directions reported as H2a + H2b):
#         S_TypePseudo_Passive:Semantics.
#         Secondary prediction; Turkish pilot showed larger effect for pseudo-passives
#         (opposite direction).
#
#   H1a — Semantics > 0 for passives (reference level) — supporting context check.

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
# restricted model is obtained by setting a parameter to a point value (H0: theta = 0)
# and the prior is continuous at that point.
# Reference: Wagenmakers et al. (2010); Verdinelli & Wasserman (1995).
# For directional hypotheses, we use the one-sided variant:
#   BF_10 = p(theta > 0 | data) / p(theta > 0 | prior)
# This requires sample_prior = "yes" in brms.

#' Compute Savage-Dickey BF for a directional hypothesis from a brmsfit object.
#' @param fit A brmsfit with sample_prior = "yes".
#' @param param Character; Stan parameter name (e.g. "b_Semantics_scaled").
#' @param direction "positive" or "negative"; the direction of H1.
#' @return Tibble with BF_10, BF_01, posterior_prob, prior_prob, log_BF10.
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

savage_dickey_directional_bf <- function(fit, param, direction = "positive") {
  stopifnot(direction %in% c("positive", "negative"))

  # Extract posterior samples
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
    bf_10 >= 100  ~ "extreme_H1",
    bf_10 >= 30   ~ "very_strong_H1",
    bf_10 >= 10   ~ "strong_H1",
    bf_10 >= 3    ~ "moderate_H1",
    bf_10 >= 1    ~ "anecdotal_H1",
    bf_10 >= 1/3  ~ "anecdotal_H0",
    bf_10 >= 1/10 ~ "moderate_H0",
    TRUE          ~ "strong_or_more_H0"
  )
}

#' Compute all focal hypothesis tests for a fitted CLAPS model.
#' @param fit A brmsfit object with sample_prior = "yes".
#' @param has_pseudo_passive Logical.
#' @param semantics_var Name suffix used in Stan parameter, e.g. "Semantics_scaled".
#' @return Tibble with one row per hypothesis.
compute_all_bf <- function(fit, has_pseudo_passive = TRUE,
                            semantics_var = "Semantics_scaled") {
  s <- gsub("_", "", semantics_var)  # brms replaces _ with nothing in param names sometimes
  # brms uses underscores; actual param name constructed here
  s_param <- paste0("b_", semantics_var)
  ia_active_param  <- paste0("b_S_TypeActive:", semantics_var)
  # Fix: brms uses '.' or ':' — check actual posterior names
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
  } else NA_character_

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

#' Summarise BF results as a wide table for reporting.
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
