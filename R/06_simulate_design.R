# R/06_simulate_design.R
# Bayes-factor design analysis: simulate data from the generative model,
# fit ordinal mixed models, and compute Bayes-factor operating characteristics.
# One cell = one combination of (n_participants, n_verbs, prior_regime,
#             threshold_mode, model_level, language, seed).
# This script is called by scripts/04_design_analysis_cell.R with a row index.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ordinal)
})

source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/05_hypothesis_tests.R")

# ---------------------------------------------------------------------------
# Data simulation
# ---------------------------------------------------------------------------

#' Simulate ordinal CLAPS responses from a cumulative-logit model.
#' @param n_participants Integer.
#' @param n_verbs Integer.
#' @param n_items_per_cell Integer; items per participant × condition cell.
#' @param beta_semantics True semantics slope (logit scale).
#' @param beta_active_interaction True Active:Semantics interaction slope.
#' @param beta_pseudo_interaction True Pseudo_Passive:Semantics interaction slope.
#' @param sd_participant_intercept SD for participant random intercepts.
#' @param sd_participant_semantics SD for participant Semantics random slopes.
#' @param sd_verb_intercept SD for verb random intercepts.
#' @param thresholds Length-6 vector of cumulative logit thresholds.
#' @param has_pseudo_passive Logical.
#' @param include_gender Logical; if TRUE, generate a referent-gender covariate
#'   (Man/Woman) and add beta_gender to the linear predictor. The Gender draws
#'   are made ONLY when include_gender is TRUE, so the baseline (FALSE) RNG
#'   stream and simulated responses are unchanged.
#' @param beta_gender True Gender (Woman vs Man) effect on the logit scale;
#'   used only when include_gender is TRUE.
#' @param seed Integer.
#' @return Tibble with simulated CLAPS data (a Gender column is added when
#'   include_gender is TRUE).
simulate_claps_data <- function(
  n_participants           = 30,
  n_verbs                  = 20,
  n_items_per_cell         = 1,
  beta_semantics           = 0.8,
  beta_active_interaction  = -0.5,
  beta_pseudo_interaction  = 0.0,
  sd_participant_intercept = 0.6,
  sd_participant_semantics = 0.3,
  sd_verb_intercept        = 0.5,
  thresholds               = c(-3, -1.5, -0.2, 1.0, 2.2, 3.5),
  has_pseudo_passive       = TRUE,
  include_gender           = FALSE,
  beta_gender              = 0.3,
  seed                     = 1L
) {
  set.seed(seed)

  s_types <- if (has_pseudo_passive) {
    c("Passive", "Active", "Pseudo_Passive")
  } else {
    c("Passive", "Active")
  }

  # Participant-level random effects
  part_ids      <- paste0("P", seq_len(n_participants))
  part_intercepts <- rnorm(n_participants, 0, sd_participant_intercept)
  part_slopes     <- rnorm(n_participants, 0, sd_participant_semantics)
  names(part_intercepts) <- part_ids
  names(part_slopes)     <- part_ids

  # Verb-level random effects
  verb_ids       <- paste0("V", seq_len(n_verbs))
  verb_intercepts <- rnorm(n_verbs, 0, sd_verb_intercept)
  names(verb_intercepts) <- verb_ids

  # Semantics values (continuous, ~ uniform)
  semantics_vals <- seq(-0.5, 0.5, length.out = 10)

  rows <- purrr::map_dfr(part_ids, function(pid) {
    purrr::map_dfr(s_types, function(st) {
      verb_sample <- sample(verb_ids, size = n_verbs, replace = FALSE)
      sem_sample  <- sample(semantics_vals, size = n_verbs, replace = TRUE)
      purrr::map2_dfr(verb_sample, sem_sample, function(vid, sem) {
        # Optional referent-gender covariate (Man/Woman). The draw is made
        # only when include_gender is TRUE, so the baseline RNG stream — and
        # therefore every previously simulated cell — is left unchanged.
        g        <- if (include_gender) sample(c("Man", "Woman"), size = 1) else NA_character_
        b_gender <- if (include_gender) beta_gender * as.numeric(g == "Woman") else 0

        # Linear predictor
        b_st   <- if (st == "Active") 0 else if (st == "Pseudo_Passive") 0 else 0  # S_Type main effect = 0 for now
        b_int  <- if (st == "Active") beta_active_interaction * sem
                  else if (st == "Pseudo_Passive") beta_pseudo_interaction * sem
                  else 0
        eta <- part_intercepts[pid] + part_slopes[pid] * sem +
               verb_intercepts[vid] + beta_semantics * sem + b_int + b_gender

        # Ordinal probabilities from cumulative logit
        cum_probs <- plogis(thresholds - eta)
        probs     <- diff(c(0, cum_probs, 1))
        probs     <- pmax(probs, 0)
        probs     <- probs / sum(probs)

        resp <- sample(seq_len(7), size = 1, prob = probs)

        row <- tibble::tibble(
          Participant = pid,
          Verb        = vid,
          S_Type      = st,
          Semantics   = sem,
          Semantics_scaled = sem,  # already scaled
          Response    = as.integer(resp)
        )
        if (include_gender) {
          row$Gender <- factor(g, levels = c("Man", "Woman"))
        }
        row
      })
    })
  })

  # Treatment-code S_Type with Passive as the reference level so the fitted
  # coefficient names (S_TypeActive[:Semantics_scaled],
  # S_TypePseudo_Passive[:Semantics_scaled]) match the prior coefficient names
  # in build_brms_prior(). Without this, brms factorises the character column
  # alphabetically (reference = "Active") and the interaction priors fail to
  # correspond to any model parameter.
  rows$S_Type <- factor(rows$S_Type, levels = s_types)
  rows
}

#' Simulate ordinal CLAPS data across multiple languages for the cross-language
#' (multilanguage) model ladder. Adds a Language grouping factor with by-language
#' random intercepts and Semantics slopes, on top of by-participant and by-verb
#' effects. Participants and verbs are language-specific (IDs are prefixed with
#' the language) and so are nested within Language but crossed with each other.
#'
#' Simplifying assumption: each language contributes n_participants participants
#' and n_verbs verbs, and (when has_pseudo_passive) all languages include the
#' Pseudo_Passive level. The real Norwegian data lack pseudo-passives; that
#' imbalance is not reproduced here, which is conservative for a power analysis.
#'
#' @param n_participants Integer; participants PER LANGUAGE.
#' @param n_verbs Integer; verbs PER LANGUAGE.
#' @param languages Character vector of language labels.
#' @param sd_language_intercept,sd_language_semantics By-language random-effect SDs.
#' @inheritParams simulate_claps_data
#' @return Tibble with a Language factor (and Gender when include_gender).
simulate_claps_data_multilanguage <- function(
  n_participants           = 30,
  n_verbs                  = 20,
  beta_semantics           = 0.8,
  beta_active_interaction  = -0.5,
  beta_pseudo_interaction  = 0.0,
  sd_participant_intercept = 0.6,
  sd_participant_semantics = 0.3,
  sd_verb_intercept        = 0.5,
  sd_language_intercept    = 0.4,
  sd_language_semantics    = 0.2,
  thresholds               = c(-3, -1.5, -0.2, 1.0, 2.2, 3.5),
  languages                = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive       = TRUE,
  include_gender           = FALSE,
  beta_gender              = 0.3,
  seed                     = 1L
) {
  set.seed(seed)

  s_types <- if (has_pseudo_passive) {
    c("Passive", "Active", "Pseudo_Passive")
  } else {
    c("Passive", "Active")
  }

  # Language-level random effects (intercept + Semantics slope)
  lang_int <- rnorm(length(languages), 0, sd_language_intercept)
  lang_slp <- rnorm(length(languages), 0, sd_language_semantics)
  names(lang_int) <- languages
  names(lang_slp) <- languages

  semantics_vals <- seq(-0.5, 0.5, length.out = 10)

  rows <- purrr::map_dfr(languages, function(lang) {
    part_ids <- paste0(lang, "_P", seq_len(n_participants))
    part_int <- rnorm(n_participants, 0, sd_participant_intercept)
    part_slp <- rnorm(n_participants, 0, sd_participant_semantics)
    names(part_int) <- part_ids
    names(part_slp) <- part_ids

    verb_ids <- paste0(lang, "_V", seq_len(n_verbs))
    verb_int <- rnorm(n_verbs, 0, sd_verb_intercept)
    names(verb_int) <- verb_ids

    purrr::map_dfr(part_ids, function(pid) {
      purrr::map_dfr(s_types, function(st) {
        verb_sample <- sample(verb_ids, size = n_verbs, replace = FALSE)
        sem_sample  <- sample(semantics_vals, size = n_verbs, replace = TRUE)
        purrr::map2_dfr(verb_sample, sem_sample, function(vid, sem) {
          g        <- if (include_gender) sample(c("Man", "Woman"), size = 1) else NA_character_
          b_gender <- if (include_gender) beta_gender * as.numeric(g == "Woman") else 0
          b_int    <- if (st == "Active") beta_active_interaction * sem
                      else if (st == "Pseudo_Passive") beta_pseudo_interaction * sem
                      else 0
          eta <- lang_int[lang] + part_int[pid] + verb_int[vid] +
                 (beta_semantics + lang_slp[lang] + part_slp[pid]) * sem + b_int + b_gender

          cum_probs <- plogis(thresholds - eta)
          probs     <- diff(c(0, cum_probs, 1))
          probs     <- pmax(probs, 0)
          probs     <- probs / sum(probs)
          resp      <- sample(seq_len(7), size = 1, prob = probs)

          row <- tibble::tibble(
            Participant = pid, Verb = vid, Language = lang, S_Type = st,
            Semantics = sem, Semantics_scaled = sem, Response = as.integer(resp)
          )
          if (include_gender) row$Gender <- factor(g, levels = c("Man", "Woman"))
          row
        })
      })
    })
  })

  rows$S_Type   <- factor(rows$S_Type,   levels = s_types)
  rows$Language <- factor(rows$Language, levels = languages)
  rows
}

# ---------------------------------------------------------------------------
# Single design-analysis cell
# ---------------------------------------------------------------------------

#' Run one design-analysis cell: simulate, fit, compute BF.
#' Returns a tibble with BF results and diagnostics.
#' @param cell A single-row tibble or named list from the design grid.
#' @param out_dir Output directory for saving results.
#' @param overwrite Logical.
run_design_cell <- function(cell, out_dir = "outputs/design_analysis", overwrite = FALSE) {
  cell_id <- paste(
    cell$language, cell$model_level, cell$prior_regime, cell$threshold_mode,
    cell$n_participants, cell$n_verbs, cell$seed,
    sep = "_"
  )
  out_file <- file.path(out_dir, paste0(cell_id, ".rds"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(out_file) && !overwrite) {
    message("[design_cell] Skipping existing: ", cell_id)
    return(readRDS(out_file))
  }

  message("[design_cell] Running: ", cell_id)
  t_start <- proc.time()[["elapsed"]]

  has_pp <- isTRUE(cell$has_pseudo_passive)

  # Optional gender model variation (referent-gender covariate). Defaults to
  # FALSE so existing grid cells are unaffected.
  include_gender <- isTRUE(cell$include_gender)
  beta_gender    <- cell$beta_gender %||% 0.3

  # Cross-language cells (language "AllLanguages" / "*_cross_*" levels) use the
  # multilanguage simulator + ladder; everything else uses the single-language path.
  is_cross <- identical(cell$language, "AllLanguages") || grepl("_cross", cell$model_level)

  if (is_cross) {
    sim_data <- simulate_claps_data_multilanguage(
      n_participants          = cell$n_participants,
      n_verbs                 = cell$n_verbs,
      beta_semantics          = cell$beta_semantics,
      beta_active_interaction = cell$beta_active_interaction,
      beta_pseudo_interaction = if (has_pp) cell$beta_pseudo_interaction else 0,
      has_pseudo_passive      = has_pp,
      include_gender          = include_gender,
      beta_gender             = beta_gender,
      seed                    = cell$seed
    )
    ladder <- build_multilanguage_ladder(include_gender = include_gender)
  } else {
    sim_data <- simulate_claps_data(
      n_participants          = cell$n_participants,
      n_verbs                 = cell$n_verbs,
      beta_semantics          = cell$beta_semantics,
      beta_active_interaction = cell$beta_active_interaction,
      beta_pseudo_interaction = if (has_pp) cell$beta_pseudo_interaction else 0,
      has_pseudo_passive      = has_pp,
      include_gender          = include_gender,
      beta_gender             = beta_gender,
      seed                    = cell$seed
    )
    ladder <- build_model_ladder(has_pseudo_passive = has_pp,
                                 include_gender = include_gender)
  }
  formula <- ladder[[cell$model_level]]
  if (is.null(formula)) {
    stop("[design_cell] model_level '", cell$model_level, "' not found in ladder.")
  }

  # Build prior. Pass has_pseudo_passive so the pseudo-passive interaction prior
  # is omitted for languages without that level (e.g. Norwegian); otherwise the
  # prior references a coefficient the model does not have.
  prior_obj <- build_brms_prior(
    regime_name        = cell$prior_regime,
    threshold_mode     = cell$threshold_mode,
    has_pseudo_passive = has_pp
  )
  # Drop priors that don't correspond to a parameter of this model level (e.g.
  # the lkj correlation prior for uncorrelated / intercept-only ladder levels).
  prior_obj <- align_prior_to_model(prior_obj, formula, sim_data)

  samp <- production_sampling(
    iter   = as.integer(cell$iter   %||% 4000),
    warmup = as.integer(cell$warmup %||% 2000),
    chains = as.integer(cell$chains %||% 4),
    seed   = cell$seed
  )
  ctrl <- production_control()

  # Fit model
  fit_result <- tryCatch({
    fit <- brms::brm(
      formula,
      data         = sim_data,
      prior        = prior_obj,
      backend      = "cmdstanr",
      sample_prior = "yes",
      iter         = samp$iter,
      warmup       = samp$warmup,
      chains       = samp$chains,
      cores        = samp$cores,
      seed         = samp$seed,
      control      = ctrl,
      silent       = 2
    )
    list(fit = fit, error = NULL)
  }, error = function(e) {
    list(fit = NULL, error = conditionMessage(e))
  })

  t_end <- proc.time()[["elapsed"]]
  runtime_sec <- t_end - t_start

  if (!is.null(fit_result$error)) {
    result <- tibble::tibble(
      cell_id       = cell_id,
      status        = "error",
      error_message = fit_result$error,
      runtime_sec   = runtime_sec
    )
    tmp <- paste0(out_file, ".tmp")
    saveRDS(result, tmp)
    file.rename(tmp, out_file)
    return(result)
  }

  fit  <- fit_result$fit
  diag <- extract_convergence_diagnostics(fit)
  bf   <- tryCatch(
    compute_all_bf(fit, has_pseudo_passive = has_pp),
    error = function(e) tibble::tibble(error = conditionMessage(e))
  )

  result <- dplyr::bind_cols(
    tibble::tibble(
      cell_id     = cell_id,
      status      = "success",
      runtime_sec = runtime_sec
    ),
    dplyr::select(dplyr::as_tibble(cell), -dplyr::any_of("cell_id"))
  )

  result <- list(
    summary     = result,
    bf_results  = bf,
    diagnostics = diag
  )

  tmp <- paste0(out_file, ".tmp")
  saveRDS(result, tmp)
  file.rename(tmp, out_file)
  message("[design_cell] Done: ", cell_id, " (", round(runtime_sec, 1), "s)")
  result
}

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
