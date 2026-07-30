# R/04_model_formulas.R
#
# Purpose
#   Define, once and in advance, every model the analysis is permitted to fit.
#   Two ladders live here: a single-language one (L5 down to L0) and a
#   cross-language one (L7 down to L0) that adds Language as a third grouping
#   factor. R/09_model_ladder.R walks a ladder; this file only builds it.
#
# The likelihood
#   All levels use a cumulative ordinal model with a logit link and flexible
#   thresholds: brms::cumulative(link = "logit", threshold = "flexible").
#   "Cumulative" treats a 1-7 acceptability rating as ordered categories rather
#   than as a number, which matters because the intervals between Likert points
#   are not known to be equal. Modelling such ratings as metric can both invent
#   effects and hide real ones (Liddell & Kruschke, 2018,
#   doi:10.1016/j.jesp.2018.08.009). "Flexible" lets the six thresholds take any
#   ordered values instead of being forced onto a parametric spacing, which is
#   what allows the model to represent the ceiling effects these ratings show.
#
# What varies down a ladder
#   The fixed effects are identical at every level. Only the random-effects
#   structure is reduced, and always in a prespecified order. The first step
#   (L5 to L4) replaces "|" with "||", dropping the correlations among random
#   slopes while keeping the slopes themselves; for a term with k random effects
#   that removes k(k-1)/2 correlation parameters, which is usually the largest
#   single saving available and the least costly in interpretive terms. Later
#   steps remove slopes, starting with the by-participant interaction slope,
#   which the pilot showed to be the least well identified.
#
# Assumption this file depends on
#   S_Type must already be treatment-coded with Passive as the reference level
#   (R/02_preprocess_factors.R). The formulas name no coefficients directly, but
#   the priors and hypothesis tests do, and they assume that coding.
#
# Sampler settings
#   production_control() and production_sampling() at the foot of this file supply
#   defaults for interactive use. The authoritative values live in
#   config/analysis_config.yaml, and the pipeline scripts pass them explicitly;
#   the defaults here match neither the heavy convergence-demonstration sampler
#   nor the lighter per-replicate one, so do not read them as the settings used
#   for any reported result.

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
                               gender_spec   = c("none", "main", "interaction"),
                               include_gender = NULL) {

  # Back-compat shim: legacy callers (R/09_model_ladder.R, pilot scripts, tests)
  # pass include_gender = TRUE/FALSE; translate to gender_spec. New callers pass
  # gender_spec directly.
  if (!is.null(include_gender)) {
    gender_spec <- if (isTRUE(include_gender)) "main" else "none"
  }
  gender_spec <- match.arg(gender_spec)

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
  #
  # Gender enters the FIXED effects only (never a random slope, at any level):
  #   none        -> Response ~ S_Type * Semantics
  #   main        -> ... + Gender                          (legacy main-effect covariate)
  #   interaction -> Response ~ S_Type * Semantics * Gender (full 3-way; Gender must be
  #                  sum/deviation-coded upstream so the focal S_Type/Semantics terms
  #                  stay gender-averaged). The 3-way adds only fixed coefficients, so
  #                  the cost-driving by-group random covariance is unchanged.
  fe <- if (gender_spec == "interaction") {
    paste0(r, " ~ S_Type * ", s, " * Gender")
  } else if (gender_spec == "main") {
    paste0(r, " ~ S_Type * ", s, " + Gender")
  } else {
    paste0(r, " ~ S_Type * ", s)
  }

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

  # For languages without pseudo-passives, S_Type has only two levels (Active,
  # Passive). The formula terms remain the same; the pseudo-passive contrast simply
  # does not appear in the posterior because the level is absent from the data.
  # has_pseudo_passive is therefore recorded as an annotation and does not alter
  # any formula built above. The argument matters to the *callers*: build_brms_prior()
  # uses it to omit a prior for a coefficient that will not exist, and
  # compute_all_bf() uses it to decide whether to report H2 at all.
  attr(ladder, "has_pseudo_passive") <- has_pseudo_passive
  attr(ladder, "semantics_var")      <- semantics_var
  attr(ladder, "response_var")       <- response_var
  attr(ladder, "gender_spec")        <- gender_spec
  attr(ladder, "include_gender")     <- (gender_spec != "none")

  ladder
}

#' Return the names of the model ladder in descending order of complexity.
#'
#' @return Character vector, most complex first. This ordering *is* the fallback
#'   order: R/09_model_ladder.R iterates over it, and next_fallback() below reads
#'   the successor from it. It is duplicated in select_highest_feasible_model()
#'   in R/07_extract_diagnostics.R, which must be updated in step with any change
#'   here.
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
#'
#' @param level_name A level name from ladder_names().
#' @return The next simpler level, or NA_character_ at the bottom of the ladder
#'   or when the name is unrecognised. Both cases return NA so that a caller's
#'   loop terminates rather than erroring; an unrecognised name reaching here
#'   means the caller built it wrongly, and R/09_model_ladder.R checks for that
#'   explicitly at its start.
next_fallback <- function(level_name) {
  nms <- ladder_names()
  idx <- match(level_name, nms)
  if (is.na(idx) || idx >= length(nms)) return(NA_character_)
  nms[idx + 1L]
}

#' Standard brms `control` arguments for production fits.
#'
#' @param adapt_delta Target acceptance rate for the NUTS step-size adaptation.
#'   Raised from Stan's default of 0.8 to 0.99, which forces a smaller step size.
#'   Hierarchical ordinal models of this kind have sharply curved posterior
#'   geometry near the variance-component boundaries, where the default step size
#'   produces divergent transitions; a higher adapt_delta trades sampling speed
#'   for the absence of those divergences.
#' @param max_treedepth Ceiling on the NUTS trajectory length, as a power of 2.
#'   Raised from the default of 10 to 12. Note that saturating this limit is a
#'   warning about efficiency rather than about validity, and it is counted as a
#'   convergence failure in R/07_extract_diagnostics.R only because it degrades the
#'   tail resolution the Bayes factors depend on. It is deliberately not raised
#'   further ad hoc: config/analysis_config.yaml marks it "diagnostic escalation
#'   only", because repeatedly raising it masks a model that is misspecified for
#'   the data rather than merely slow.
#' @return A list for brms::brm(control = ...).
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
#'
#' @param response_var Character; name of the response variable.
#' @param semantics_var Character; name of the scaled Semantics predictor.
#' @param gender_spec "none", "main" or "interaction"; how referent gender enters
#'   the fixed effects.
#' @param include_gender Deprecated logical, translated to gender_spec for older
#'   callers.
#' @return Named list of brmsformula objects, in order L7_cross … L0_cross.
#' @details L5_cross_maximal reproduces the published OSF reference model exactly,
#'   which is why it, rather than the topmost rung, is the point of comparison with
#'   previous work. L6 and L7 sit *above* it and exist to test a specific
#'   expectation rather than to be selected: they place the full fixed structure,
#'   including gender, into the by-Language random effects, and with only three to
#'   twenty language groups that covariance block is not expected to be
#'   identifiable. Running them documents where on the complexity gradient
#'   convergence actually breaks, instead of asserting it. Levels below L5 are
#'   computational fallbacks in the usual sense.
#'
#'   Descending this ladder reduces the by-Participant and by-Verb structure first
#'   and keeps the by-Language slopes longest (compare L2, which has intercepts
#'   only for participants and verbs but full slopes for Language). Cross-language
#'   variation in the affectedness effect is the object of interest here, so it is
#'   the last thing given up.
build_multilanguage_ladder <- function(response_var  = "Response",
                                       semantics_var = "Semantics_scaled",
                                       gender_spec   = c("none", "main", "interaction"),
                                       include_gender = NULL) {
  if (!is.null(include_gender)) {
    gender_spec <- if (isTRUE(include_gender)) "main" else "none"
  }
  gender_spec <- match.arg(gender_spec)
  r <- response_var
  s <- semantics_var

  # Fixed effects are the same across all ladder levels. Gender (when included)
  # is a fixed-effect-only covariate — see build_model_ladder() for the crossed
  # random-effects rationale (no Gender random slope; no Item-level term).
  #   none/main/interaction exactly as in build_model_ladder; "interaction" builds
  #   the full 3-way S_Type * Semantics * Gender in the fixed effects only.
  fe <- if (gender_spec == "interaction") {
    paste0(r, " ~ S_Type * ", s, " * Gender")
  } else if (gender_spec == "main") {
    paste0(r, " ~ S_Type * ", s, " + Gender")
  } else {
    paste0(r, " ~ S_Type * ", s)
  }

  # Random-slope term for the maximal rungs (L6/L7). The truly-maximal model puts
  # the FULL fixed structure into every random-slope term (Barr et al. 2013); with
  # the 3-way gender model that includes gender. These rungs are hypothesised NOT to
  # converge (an over-parameterised by-Language covariance from few language groups);
  # the feasibility run verifies the failure. Without gender they reduce to the
  # existing S_Type*Semantics maximal.
  rs_full  <- if (gender_spec == "interaction") paste0("S_Type * ", s, " * Gender") else paste0("S_Type * ", s)
  rs_gmain <- if (gender_spec == "interaction") paste0("S_Type * ", s, " + Gender") else paste0("S_Type * ", s)

  # The Language grouping factor appears in all levels; complexity of
  # by-Language slopes decreases as we descend the ladder.
  ladder <- list(

    # L7: TRULY MAXIMAL — the full fixed structure as correlated random slopes on
    # all three grouping factors (Barr et al. 2013, "keep it maximal"). With the
    # 3-way gender model this is a 12-term by-group covariance; the by-Language
    # block, estimated from only 3-20 language groups, is hypothesised to be
    # unidentifiable, so the model should NOT converge. Run to verify the failure.
    L7_cross_truly_maximal = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + ", rs_full, " | Participant)",
        " + (1 + ", rs_full, " | Verb)",
        " + (1 + ", rs_full, " | Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

    # L6: maximal + gender random MAIN effect only — the cleanly-nested step
    # between L5 (no gender random term) and L7 (full gender interaction). Locates
    # where on the complexity gradient convergence breaks.
    L6_cross_gender_maximal = brms::bf(
      as.formula(paste0(
        fe,
        " + (1 + ", rs_gmain, " | Participant)",
        " + (1 + ", rs_gmain, " | Verb)",
        " + (1 + ", rs_gmain, " | Language)"
      )),
      family = brms::cumulative(link = "logit", threshold = "flexible")
    ),

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
  attr(ladder, "gender_spec")    <- gender_spec
  attr(ladder, "include_gender") <- (gender_spec != "none")

  ladder
}

#' Return the names of the cross-language ladder in descending complexity order.
multilanguage_ladder_names <- function() {
  c(
    "L7_cross_truly_maximal",
    "L6_cross_gender_maximal",
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
#'
#' @param iter,warmup Total and warm-up iterations per chain, so the retained
#'   draws per chain are iter - warmup.
#' @param chains Number of chains. Four is the minimum at which R-hat is
#'   informative about between-chain disagreement.
#' @param cores Chains to run in parallel. Read from the STAN_NUM_THREADS
#'   environment variable so that a SLURM script can set it from
#'   --cpus-per-task and the same code adapts to the allocation it was given,
#'   rather than hard-coding a core count that would either oversubscribe or
#'   waste the node.
#' @param seed Sampler seed. Fixed so a fit is reproducible; the design-analysis
#'   grids override it per cell, since there the seed is what distinguishes one
#'   simulated replicate from another.
#' @return A list of arguments for brms::brm().
#' @details These defaults are for interactive and smoke-test use. The reported
#'   runs pass values from config/analysis_config.yaml explicitly; see the note in
#'   the file header.
production_sampling <- function(iter = 4000, warmup = 2000, chains = 4,
                                cores = as.integer(Sys.getenv("STAN_NUM_THREADS", "4")),
                                seed = 12345) {
  list(iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed)
}
