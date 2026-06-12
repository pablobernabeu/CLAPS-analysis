# R/04_model_formulas.R
# Prespecified model ladder for the CLAPS Bayesian ordinal mixed-effects model.
# Single-language ladder: L5 (maximal correlated) down to L0 (intercepts only).
# For languages without Pseudo_Passive, terms are dropped before returning formulas.
# S_Type must be treatment-coded with Passive as reference.

suppressPackageStartupMessages({
  library(brms)
})

#' Return the ordered model ladder as a named list of brmsformula objects.
#' @param has_pseudo_passive Logical; if FALSE, pseudo-passive terms are absent.
#' @param response_var Character; name of the response variable.
#' @param semantics_var Character; name of the Semantics predictor (scaled or raw).
#' @param include_gender Logical; if TRUE, add a Gender fixed-effect covariate
#'   (referent gender, Man/Woman; see R/02_preprocess_factors.R::derive_gender()).
#'   Used by the gender model variation; non-focal, so it takes the default
#'   weakly-informative "b" prior and is not part of any focal hypothesis.
#' @return Named list of brmsformula objects, in order L5 … L0.
build_model_ladder <- function(has_pseudo_passive = TRUE,
                                response_var  = "Response",
                                semantics_var = "Semantics_scaled",
                                include_gender = FALSE) {

  r <- response_var
  s <- semantics_var

  # Fixed-effects structure is identical across ladder levels (all include the
  # full S_Type * Semantics interaction). When include_gender is TRUE the
  # referent-gender covariate is added as a FIXED EFFECT ONLY.
  #
  # Random-effects rationale (verified against the pilot data): Gender is fully
  # crossed with the grouping factors — every Participant, every Verb, and even
  # every Participant x Verb cell contains both Man and Woman items; only Item
  # itself determines gender. The fixed Gender effect is therefore cleanly
  # identified by the existing crossed by-Participant / by-Verb random structure.
  # Gender is deliberately NOT added to any random slope, and no by-Item or
  # Participant:Verb random effect is introduced: a random Gender slope would add
  # many parameters (overfitting), and an Item-level random effect would be
  # collinear with Gender. This is the "keep the random effects crossed" case.
  fe <- paste0(r, " ~ S_Type * ", s)
  if (isTRUE(include_gender)) fe <- paste0(fe, " + Gender")

  ladder <- list(

    L5_correlated_maximal = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type * ", s, " | Participant)",
        " + (1 + S_Type | Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    L4_uncorrelated_maximal = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type * ", s, " || Participant)",
        " + (1 + S_Type || Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    L3_no_participant_interaction_slope = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type + ", s, " || Participant)",
        " + (1 + S_Type || Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    L2_sentence_type_slopes_only = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type || Participant)",
        " + (1 + S_Type || Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    L1_random_intercepts_plus_participant_semantics = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + ", s, " || Participant)",
        " + (1 | Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    L0_random_intercepts_only = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 | Participant)",
        " + (1 | Verb)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    )
  )

  # For languages without pseudo-passives, S_Type has only two levels (Active, Passive).
  # The formula terms remain the same; the pseudo-passive contrast simply does not
  # appear in the posterior because the level is absent from the data.
  # We annotate the ladder to signal this.
  attr(ladder, "has_pseudo_passive") <- has_pseudo_passive
  attr(ladder, "semantics_var")      <- semantics_var
  attr(ladder, "response_var")       <- response_var
  attr(ladder, "include_gender")     <- include_gender

  ladder
}

#' Return the names of the model ladder in descending order of complexity.
ladder_names <- function() {
  c(
    "L5_correlated_maximal",
    "L4_uncorrelated_maximal",
    "L3_no_participant_interaction_slope",
    "L2_sentence_type_slopes_only",
    "L1_random_intercepts_plus_participant_semantics",
    "L0_random_intercepts_only"
  )
}

#' Given a ladder name, return the next fallback level.
#' Returns NA if already at L0.
next_fallback <- function(level_name) {
  nms <- ladder_names()
  idx <- match(level_name, nms)
  if (is.na(idx) || idx >= length(nms)) return(NA_character_)
  nms[idx + 1L]
}

#' Standard brms control arguments for production fits.
production_control <- function(adapt_delta = 0.99, max_treedepth = 12) {
  list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
}

# ---------------------------------------------------------------------------
# Cross-language (multi-language) model ladder
# ---------------------------------------------------------------------------
# Reference model from the OSF study
# (V6_All_Semantics_Only.R in Sub2_OSF_Passives.zip):
#
#   SemOnly = brm(
#     Response ~ S_Type*Semantics
#       + (1 + S_Type*Semantics | Participant)
#       + (1 + S_Type*Semantics | Verb)
#       + (1 + S_Type*Semantics | Language),
#     family = cumulative(), ...
#   )
#
# Note: S_Type is treatment-coded with Passive as reference level.
# Norwegian rows with S_Type == "Synthetic_Passive" are excluded upstream.
# Verb labels are language-specific (Language_VERB format).

#' Return the cross-language model ladder as a named list of brmsformula objects.
#' The L5 level matches the OSF reference model exactly; lower levels are
#' computational fallbacks.
#' @param response_var Character; name of the response variable.
#' @param semantics_var Character; name of the scaled Semantics predictor.
#' @return Named list of brmsformula objects, in order L5_cross … L0_cross.
build_multilanguage_ladder <- function(response_var  = "Response",
                                        semantics_var = "Semantics_scaled",
                                        include_gender = FALSE) {
  r <- response_var
  s <- semantics_var

  # Fixed effects are the same across all ladder levels. Gender (when included)
  # is a fixed-effect-only covariate — see build_model_ladder() for the crossed
  # random-effects rationale (no Gender random slope; no Item-level term).
  fe <- paste0(r, " ~ S_Type * ", s)
  if (isTRUE(include_gender)) fe <- paste0(fe, " + Gender")

  # The Language grouping factor appears in all levels; complexity of
  # by-Language slopes decreases as we descend the ladder.
  ladder <- list(

    # L5: OSF reference model — maximal correlated random effects for all
    # three grouping factors (Participant, Verb, Language).
    L5_cross_maximal = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type * ", s, " | Participant)",
        " + (1 + S_Type * ", s, " | Verb)",
        " + (1 + S_Type * ", s, " | Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L4: uncorrelated random effects (faster; use if L5 OOMs or diverges).
    L4_cross_uncorrelated = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type * ", s, " || Participant)",
        " + (1 + S_Type * ", s, " || Verb)",
        " + (1 + S_Type * ", s, " || Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L3: drop interaction slope for Participant; keep for Verb and Language.
    L3_cross_no_participant_interaction = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type + ", s, " || Participant)",
        " + (1 + S_Type * ", s, " || Verb)",
        " + (1 + S_Type * ", s, " || Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L2: S_Type slopes only for Participant/Verb; full by-Language slopes.
    L2_cross_stype_participant_verb = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + S_Type || Participant)",
        " + (1 + S_Type || Verb)",
        " + (1 + S_Type * ", s, " || Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L1: random intercepts for Participant/Verb; S_Type slope for Language.
    L1_cross_intercepts_only_ppt_verb = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 | Participant)",
        " + (1 | Verb)",
        " + (1 + S_Type || Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L0: random intercepts only — minimal cross-language model.
    L0_cross_intercepts_only = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 | Participant)",
        " + (1 | Verb)",
        " + (1 | Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    )
  )

  attr(ladder, "model_type")   <- "cross_language"
  attr(ladder, "semantics_var") <- semantics_var
  attr(ladder, "response_var")  <- response_var
  attr(ladder, "include_gender") <- include_gender

  ladder
}

#' Return the names of the cross-language ladder in descending complexity order.
multilanguage_ladder_names <- function() {
  c(
    "L5_cross_maximal",
    "L4_cross_uncorrelated",
    "L3_cross_no_participant_interaction",
    "L2_cross_stype_participant_verb",
    "L1_cross_intercepts_only_ppt_verb",
    "L0_cross_intercepts_only"
  )
}

#' Given a cross-language ladder level name, return the next fallback.
#' Returns NA_character_ if already at L0.
next_multilanguage_fallback <- function(level_name) {
  nms <- multilanguage_ladder_names()
  idx <- match(level_name, nms)
  if (is.na(idx) || idx >= length(nms)) return(NA_character_)
  nms[idx + 1L]
}

#' Standard sampling arguments for production fits.
production_sampling <- function(iter = 4000, warmup = 2000, chains = 4,
                                 cores = as.integer(Sys.getenv("STAN_NUM_THREADS", "4")),
                                 seed = 12345) {
  list(iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed)
}
