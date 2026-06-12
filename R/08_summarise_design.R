# R/08_summarise_design.R
# Aggregate design-analysis cells into operating-characteristic summaries.
# Computes BF exceedance probabilities, sensitivity to prior regime,
# model-ladder selection frequency, and failure rates.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
})

source("R/05_hypothesis_tests.R")

#' Load all completed design-analysis cell results from an output directory.
#' Skips errored cells with a warning.
load_design_cells <- function(out_dir = "outputs/design_analysis") {
  files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("[summarise] No .rds files found in ", out_dir)
  }
  message("[summarise] Loading ", length(files), " cell files.")
  purrr::map_dfr(files, function(f) {
    tryCatch({
      res <- readRDS(f)
      if (is.list(res) && "summary" %in% names(res) && "bf_results" %in% names(res)) {
        dplyr::bind_cols(
          res$summary,
          dplyr::select(res$bf_results, hypothesis, BF_10, BF_01, bf_category,
                        posterior_prob, method),
          dplyr::select(res$diagnostics, convergence_ok, max_rhat, n_divergent,
                        min_ess_bulk)
        )
      } else if (is.data.frame(res) && "status" %in% names(res)) {
        # Error/skip tibble from a previous run — preserve its status
        res
      } else {
        # Truly unrecognised format
        tibble::tibble(cell_id = basename(f), status = "malformed")
      }
    }, error = function(e) {
      warning("[summarise] Failed to load ", f, ": ", conditionMessage(e))
      tibble::tibble(cell_id = basename(f), status = "load_error")
    })
  })
}

#' Compute BF exceedance probabilities by design condition.
#' @param df Combined design cell data frame.
#' @param group_vars Character vector of grouping variables.
#' @param bf_threshold Numeric; primary BF threshold (default 10).
#' @param bf_threshold_secondary Numeric; secondary BF threshold (default 3).
compute_bf_exceedance <- function(df, group_vars = c("language", "model_level",
                                                      "n_participants", "n_verbs",
                                                      "prior_regime", "hypothesis"),
                                   bf_threshold = 10,
                                   bf_threshold_secondary = 3) {
  df |>
    dplyr::filter(status == "success") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_sims              = dplyr::n(),
      p_bf_primary        = mean(BF_10 >= bf_threshold,           na.rm = TRUE),
      p_bf_secondary      = mean(BF_10 >= bf_threshold_secondary, na.rm = TRUE),
      median_bf           = median(BF_10, na.rm = TRUE),
      mean_bf             = mean(BF_10,   na.rm = TRUE),
      p_convergence_ok    = mean(convergence_ok, na.rm = TRUE),
      .groups             = "drop"
    )
}

#' Compute sensitivity of BF category to prior regime.
prior_sensitivity_summary <- function(df) {
  df |>
    dplyr::filter(status == "success") |>
    dplyr::group_by(n_participants, n_verbs, hypothesis, seed) |>
    dplyr::summarise(
      n_regimes      = dplyr::n_distinct(prior_regime),
      n_bf_categories = dplyr::n_distinct(bf_category),
      max_bf         = max(BF_10, na.rm = TRUE),
      min_bf         = min(BF_10, na.rm = TRUE),
      bf_ratio       = max_bf / pmax(min_bf, 1e-6),
      category_stable = n_bf_categories == 1,
      .groups         = "drop"
    )
}

#' Compute model-ladder selection frequency.
ladder_selection_summary <- function(df) {
  df |>
    dplyr::filter(status == "success") |>
    dplyr::count(model_level, name = "n_selected") |>
    dplyr::mutate(prop_selected = n_selected / sum(n_selected))
}

#' Compute failure rate summary.
failure_summary <- function(df) {
  tibble::tibble(
    total_cells      = nrow(df),
    n_success        = sum(df$status == "success", na.rm = TRUE),
    n_error          = sum(df$status == "error",   na.rm = TRUE),
    n_timeout        = sum(df$status == "timeout", na.rm = TRUE),
    n_oom            = sum(df$status == "oom",     na.rm = TRUE),
    n_malformed      = sum(df$status == "malformed", na.rm = TRUE),
    failure_rate     = 1 - sum(df$status == "success", na.rm = TRUE) / nrow(df)
  )
}

#' Highest-complexity model level that converged, per language ("maximal feasible
#' model"). Feasible = fit succeeded AND met the publication-grade convergence
#' criteria (convergence_ok). Level rank is the L-number in the level name.
maximal_feasible_model <- function(df) {
  ok <- dplyr::filter(df, status == "success", convergence_ok %in% TRUE)
  if (nrow(ok) == 0) {
    return(tibble::tibble(language = character(), maximal_feasible_model = character(),
                          level_rank = integer(), n_converged_cells = integer()))
  }
  ok |>
    dplyr::mutate(level_rank = suppressWarnings(as.integer(sub("^L([0-9]).*", "\\1", model_level)))) |>
    dplyr::group_by(language) |>
    dplyr::summarise(
      maximal_feasible_model = model_level[which.max(level_rank)],
      level_rank             = max(level_rank, na.rm = TRUE),
      n_converged_cells      = dplyr::n(),
      .groups = "drop"
    )
}

#' Recommended sample size: per language, the smallest n_participants at which
#' BOTH focal hypotheses exceed the primary BF threshold with probability >=
#' target, evaluated at that language's maximal feasible model and the primary
#' (proposal) prior regime. NA if the target is not reached anywhere in the grid.
#' NOTE: a trustworthy exceedance *probability* needs many simulations per design
#' point (config `design_analysis$n_simulations_per_cell`); with one simulation
#' per cell this is an indicative 0/1 estimate, not a stable probability.
recommended_sample_size <- function(exc, mfm, target = 0.80,
                                    focal = c("H1a_semantics_positive",
                                              "H1b_active_interaction_negative"),
                                    regime = "proposal") {
  empty <- tibble::tibble(language = character(), recommended_n_participants = integer(),
                          n_verbs = integer(), meets_target = logical(),
                          target = numeric(), regime = character())
  if (is.null(exc) || nrow(exc) == 0 || nrow(mfm) == 0) return(empty)
  if (!all(c("language", "model_level") %in% names(exc))) return(empty)

  d <- exc |>
    dplyr::inner_join(dplyr::select(mfm, language, model_level = maximal_feasible_model),
                      by = c("language", "model_level")) |>
    dplyr::filter(prior_regime == regime, hypothesis %in% focal)
  if (nrow(d) == 0) return(empty)

  d |>
    dplyr::group_by(language, n_participants, n_verbs) |>
    dplyr::summarise(all_focal_ok = (dplyr::n() == length(focal)) && all(p_bf_primary >= target),
                     .groups = "drop") |>
    dplyr::group_by(language) |>
    dplyr::summarise(
      meets_target = any(all_focal_ok),
      recommended_n_participants = if (any(all_focal_ok)) min(n_participants[all_focal_ok]) else NA_integer_,
      n_verbs = if (any(all_focal_ok)) n_verbs[all_focal_ok][which.min(n_participants[all_focal_ok])] else NA_integer_,
      .groups = "drop"
    ) |>
    dplyr::mutate(target = target, regime = regime)
}

#' Per-cell fit runtimes, summarised by language and model level (in minutes).
runtime_summary <- function(df) {
  ok <- dplyr::filter(df, status == "success", !is.na(runtime_sec))
  if (nrow(ok) == 0) {
    return(tibble::tibble(language = character(), model_level = character(),
                          n_cells = integer(), median_runtime_min = numeric(),
                          max_runtime_min = numeric(), p_converged = numeric()))
  }
  ok |>
    dplyr::group_by(language, model_level) |>
    dplyr::summarise(
      n_cells            = dplyr::n(),
      median_runtime_min = round(stats::median(runtime_sec) / 60, 2),
      max_runtime_min    = round(max(runtime_sec) / 60, 2),
      p_converged        = round(mean(convergence_ok %in% TRUE), 2),
      .groups = "drop"
    ) |>
    dplyr::arrange(language, dplyr::desc(model_level))
}

#' Write aggregated design summary to CSV.
write_design_summary <- function(df, out_dir = "outputs/design_summary") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fail  <- failure_summary(df)
  readr::write_csv(fail, file.path(out_dir, "failure_summary.csv"))
  message("[summarise] Failure summary: ", fail$n_success, "/", fail$total_cells, " succeeded.")

  df_ok <- dplyr::filter(df, status == "success")

  if (nrow(df_ok) == 0) {
    message("[summarise] No successful cells — skipping BF/sensitivity/ladder summaries.")
    message("[summarise] All cells have status: ",
            paste(sort(unique(df$status)), collapse = ", "))
    message("[summarise] Partial summary written to ", out_dir)
    return(invisible(list(fail = fail)))
  }

  exc    <- compute_bf_exceedance(df)
  sens   <- prior_sensitivity_summary(df)
  ladder <- ladder_selection_summary(df)
  mfm    <- maximal_feasible_model(df)
  rec    <- recommended_sample_size(exc, mfm)
  rt     <- runtime_summary(df)

  readr::write_csv(exc,    file.path(out_dir, "bf_exceedance.csv"))
  readr::write_csv(sens,   file.path(out_dir, "prior_sensitivity.csv"))
  readr::write_csv(ladder, file.path(out_dir, "ladder_selection.csv"))
  readr::write_csv(mfm,    file.path(out_dir, "maximal_feasible_model.csv"))
  readr::write_csv(rec,    file.path(out_dir, "recommended_sample_size.csv"))
  readr::write_csv(rt,     file.path(out_dir, "runtime_summary.csv"))

  message("[summarise] Design summary written to ", out_dir,
          " (incl. maximal feasible model, recommended N, runtimes).")
  invisible(list(exc = exc, sens = sens, ladder = ladder,
                 mfm = mfm, rec = rec, rt = rt, fail = fail))
}
