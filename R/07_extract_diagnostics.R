# R/07_extract_diagnostics.R
#
# Purpose
#   Decide whether a fitted model may be believed. Every diagnostic reported here
#   feeds one binary, `convergence_ok`, which the ladder engine uses to accept or
#   reject a model level and which the aggregation scripts use to qualify the
#   reported power.
#
# What is checked, and why these four
#   R-hat        Between- and within-chain agreement. Chains that have not mixed
#                are exploring different parts of the posterior, so any summary
#                of them is meaningless.
#   Bulk ESS     Effective sample size in the body of the distribution, which
#                governs the accuracy of posterior means and probabilities.
#   Tail ESS     Effective sample size in the tails. This is the one that matters
#                most here: the directional Savage-Dickey ratio in
#                R/05_hypothesis_tests.R is a ratio of tail probabilities, so a
#                fit with adequate bulk but poor tail resolution yields an
#                unreliable Bayes factor while looking healthy on other measures.
#   Divergences and treedepth
#                Sampler-level failures. A divergent transition indicates
#                posterior geometry the sampler could not follow, which biases
#                the draws rather than merely thinning them.
#
#   The thresholds are those of Vehtari et al. (2021, doi:10.1214/20-BA1221):
#   R-hat < 1.01 and ESS >= 400. Both are stricter than the older R-hat < 1.1
#   convention, which that paper shows to be too lenient to detect the failure
#   modes that occur in hierarchical models of this kind.
#
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(posterior)
})

# Null-coalescing operator. Base R gained %||% only in 4.4, while
# config/arc_modules.yaml sets min_version 4.3.0 and CI pins 4.3.3, so it cannot be
# assumed present.
#
# It is defined HERE, in the file that uses it, rather than left to the caller.
# Fourteen files source this module, and several of them (for example
# scripts/03_prior_sensitivity.R) define no such operator themselves. Because %||%
# sits inside extract_convergence_diagnostics(), it is resolved when that function
# is CALLED, not when this file is sourced, so the omission surfaced only at the
# moment a fit was diagnosed, on an R older than 4.4. Defining it alongside its use
# removes that dependence on what the caller happens to have in scope.
#
# Note that this plain variant shares its name with the stricter one in
# R/00_reference_audit.R, and whichever is evaluated last wins; see the hazard note
# in that file.
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Extract convergence diagnostics from a brmsfit object.
#'
#' Works under both the rstan and cmdstanr backends: brms exposes an
#' rstan-compatible stanfit at `fit$fit` in either case, so the slot accesses
#' below are valid for cmdstanr-backed fits too.
#'
#' @param fit A brmsfit object.
#' @param max_treedepth Tree-depth ceiling the sampler was run with, used to
#'   count saturating transitions. Only consulted when the value cannot be read
#'   back off the fit itself. The default tracks `production_control()` in
#'   R/04_model_formulas.R; passing the wrong value here does not corrupt any
#'   estimate, but it does mis-count `n_max_treedepth`, so keep the two in step.
#' @return A one-row tibble of diagnostic summaries.
extract_convergence_diagnostics <- function(fit, max_treedepth = 12L) {
  stopifnot(inherits(fit, "brmsfit"))

  draws <- posterior::as_draws_array(fit)
  summ  <- posterior::summarise_draws(
    draws,
    posterior::default_convergence_measures()
  )

  # Population-level coefficients, group-level SDs and correlations. Deliberately
  # broader than the two focal terms: a fit whose variance components have not
  # mixed cannot be trusted for the focal Bayes factors either.
  focal <- summ |>
    dplyr::filter(grepl("^b_|^sd_|^cor_|^Intercept", variable))

  # Sampler diagnostics
  np <- brms::nuts_params(fit)

  # Prefer the ceiling actually used by this fit; fall back to the argument when
  # the slot is absent, which it can be depending on backend and brms version.
  depth_limit <- fit$fit@sim$max_depth %||% max_treedepth

  # A divergence anywhere signals a posterior geometry the sampler could not
  # follow, so these are counted as transitions, not as chains.
  n_divergent   <- sum(np$Value[np$Parameter == "divergent__"] > 0, na.rm = TRUE)
  # Saturating the tree depth is an efficiency failure rather than a validity
  # one, but it inflates Monte-Carlo error in the Savage-Dickey density ratio
  # and so is tracked with the same weight here.
  n_max_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= depth_limit,
                         na.rm = TRUE)

  # Rhat threshold checks
  max_rhat   <- max(focal$rhat,     na.rm = TRUE)
  min_ess_b  <- min(focal$ess_bulk, na.rm = TRUE)
  min_ess_t  <- min(focal$ess_tail, na.rm = TRUE)

  # Publication-grade criteria of Vehtari et al. (2021, doi:10.1214/20-BA1221):
  # R-hat below 1.01, and both bulk and tail ESS at or above 400, which is the
  # point at which their Monte-Carlo error estimates become dependable. This is
  # a strict flag by design, and is recorded rather than used to discard fits;
  # see the aggregation scripts, which report power with and without it.
  convergence_ok <- max_rhat < 1.01 & min_ess_b >= 400 & min_ess_t >= 400 &
    n_divergent == 0 & n_max_treedepth == 0

  tibble::tibble(
    n_params          = nrow(focal),
    max_rhat          = round(max_rhat,  4),
    min_ess_bulk      = round(min_ess_b, 1),
    min_ess_tail      = round(min_ess_t, 1),
    n_divergent       = n_divergent,
    n_max_treedepth   = n_max_treedepth,
    convergence_ok    = convergence_ok,
    n_chains          = fit$fit@sim$chains,
    n_iter            = fit$fit@sim$iter,
    n_warmup          = fit$fit@sim$warmup
  )
}

#' Classify convergence by prespecified, publication-grade criteria.
#' "converged" (the only acceptable outcome) requires R-hat < 1.01 on all focal
#' parameters, bulk and tail ESS >= 400, zero divergent transitions and zero
#' max-treedepth saturations (Vehtari et al., 2021). Both "marginal_*" and
#' "failed_*" indicate a non-publication-grade fit that must trigger a fallback.
#'
#' @param diag_row A one-row tibble from extract_convergence_diagnostics().
#' @return A single status string.
#' @details Two tiers of failure are distinguished even though both are rejected,
#'   because the distinction tells a reader what to do next. "failed_*" (R-hat
#'   >= 1.05, more than 10 divergences, ESS < 100) means the fit is far from
#'   usable and the model is probably too complex for the data, so descending the
#'   ladder is the right response. "marginal_*" means the fit is close, and more
#'   iterations or a higher adapt_delta at the same ladder level would likely
#'   suffice. The ladder engine treats both as a fallback trigger, so this
#'   information is advisory and lives in the logs.
#'
#'   The order of the cases is load-bearing. All "failed_" tests precede all
#'   "marginal_" tests, so a fit that is both severely and mildly deficient is
#'   reported at its worst; case_when() returns the first match. The NA check
#'   comes first because max_rhat is NA when no draws were produced at all, and
#'   every comparison below would then return NA rather than a status.
classify_convergence <- function(diag_row) {
  dplyr::case_when(
    is.na(diag_row$max_rhat)                                       ~ "unknown",
    diag_row$max_rhat        >= 1.05                               ~ "failed_rhat",
    diag_row$n_divergent     >  10                                 ~ "failed_divergences",
    diag_row$min_ess_bulk    <  100 | diag_row$min_ess_tail < 100  ~ "failed_ess",
    diag_row$max_rhat        >= 1.01                               ~ "marginal_rhat",
    diag_row$n_divergent     >  0                                  ~ "marginal_divergences",
    diag_row$min_ess_bulk    <  400 | diag_row$min_ess_tail < 400  ~ "marginal_ess",
    diag_row$n_max_treedepth >  0                                  ~ "marginal_treedepth",
    TRUE                                                           ~ "converged"
  )
}

#' Summarise model-ladder diagnostics across all levels.
#' @param diag_list Named list of diagnostic tibbles (one per model level). A
#'   NULL element denotes a level that was never reached, which happens whenever
#'   a higher level converged and the descent stopped.
#' @return Tibble with one row per model level. Levels that were not run appear
#'   with convergence_status "not_run" rather than being omitted, so the table
#'   shows the whole ladder and makes clear where the descent halted.
summarise_ladder_diagnostics <- function(diag_list) {
  purrr::imap_dfr(diag_list, function(diag, level_name) {
    if (is.null(diag)) {
      tibble::tibble(model_level = level_name, convergence_status = "not_run")
    } else {
      dplyr::mutate(diag,
        model_level        = level_name,
        convergence_status = classify_convergence(diag)
      )
    }
  })
}

#' Determine the highest feasible model level from ladder diagnostics.
#'
#' @param ladder_summary Tibble from summarise_ladder_diagnostics().
#' @return The name of the most complex converged level, or NA_character_ with a
#'   warning when none converged. NA rather than an error, because "no level fits"
#'   is a reportable finding about the design, not a programming fault.
#' @details Only the status "converged" qualifies. The "marginal_*" statuses are
#'   excluded here for the same reason the ladder engine falls back on them.
select_highest_feasible_model <- function(ladder_summary) {
  feasible <- ladder_summary |>
    dplyr::filter(convergence_status == "converged") |>
    dplyr::pull(model_level)

  # Complexity order, most complex first, so the first surviving element below is
  # the answer. NOTE: this repeats the sequence returned by ladder_names() in
  # R/04_model_formulas.R. The duplication is a maintenance hazard: a level
  # renamed or inserted there must be mirrored here, or it will silently never be
  # selected, since a name absent from this vector is dropped by the %in% filter.
  level_order <- c(
    "L5_correlated_maximal",
    "L4_uncorrelated_maximal",
    "L3_no_participant_interaction_slope",
    "L2_sentence_type_slopes_only",
    "L1_random_intercepts_plus_participant_semantics",
    "L0_random_intercepts_only"
  )

  feasible_ordered <- level_order[level_order %in% feasible]

  if (length(feasible_ordered) == 0) {
    warning("[select_model] No feasible model found in ladder.")
    return(NA_character_)
  }

  feasible_ordered[1]
}

#' Save diagnostics as a CSV alongside the model outputs.
#'
#' @param diag_tibble Diagnostics to write.
#' @param out_path Destination CSV; its directory is created if absent.
#' @return Invisibly, `out_path`.
#' @details Written via a temporary file and a rename, for the same reason the fits
#'   are (see R/09_model_ladder.R): the rename is atomic, so a job killed mid-write
#'   cannot leave a half-written CSV that a later aggregation step would read as
#'   real data.
save_diagnostics_csv <- function(diag_tibble, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(out_path, ".tmp")
  readr::write_csv(diag_tibble, tmp)
  file.rename(tmp, out_path)
  invisible(out_path)
}
