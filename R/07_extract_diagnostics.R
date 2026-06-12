# R/07_extract_diagnostics.R
# Extract and summarise convergence diagnostics from brmsfit objects.
# Reports Rhat, Bulk ESS, Tail ESS, divergences, treedepth exceedances, runtime.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(posterior)
})

#' Extract convergence diagnostics from a brmsfit object.
#' @param fit A brmsfit object.
#' @return A tibble with diagnostic summaries.
extract_convergence_diagnostics <- function(fit) {
  stopifnot(inherits(fit, "brmsfit"))

  draws <- posterior::as_draws_array(fit)
  summ  <- posterior::summarise_draws(
    draws,
    posterior::default_convergence_measures()
  )

  # Focal fixed-effect parameters only
  focal <- summ |>
    dplyr::filter(grepl("^b_|^sd_|^cor_|^Intercept", variable))

  # Sampler diagnostics
  np <- brms::nuts_params(fit)

  n_divergent   <- sum(np$Value[np$Parameter == "divergent__"] > 0, na.rm = TRUE)
  n_max_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= fit$fit@sim$max_depth %||% 12,
                         na.rm = TRUE)

  # Rhat threshold checks
  max_rhat   <- max(focal$rhat,     na.rm = TRUE)
  min_ess_b  <- min(focal$ess_bulk, na.rm = TRUE)
  min_ess_t  <- min(focal$ess_tail, na.rm = TRUE)

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
#' @param diag_list Named list of diagnostic tibbles (one per model level).
#' @return Tibble with one row per model level.
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
#' Feasibility = converged + no timeout/OOM flag.
#' @param ladder_summary Tibble from summarise_ladder_diagnostics().
#' @return Character; name of the highest feasible model level.
select_highest_feasible_model <- function(ladder_summary) {
  feasible <- ladder_summary |>
    dplyr::filter(convergence_status == "converged") |>
    dplyr::pull(model_level)

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
save_diagnostics_csv <- function(diag_tibble, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(out_path, ".tmp")
  readr::write_csv(diag_tibble, tmp)
  file.rename(tmp, out_path)
  invisible(out_path)
}
