# R/08_summarise_design.R
#
# Purpose
#   Turn the many per-cell .rds files written by the design analysis into the
#   handful of CSV summaries the report reads. A "cell" is one simulated design
#   point: one combination of language, model level, prior regime, threshold
#   mode, sample size, verb count and seed, fitted once. Aggregating over the
#   seeds within a design point converts individual Bayes factors into the
#   operating characteristic of interest, namely the probability that a study of
#   that size returns BF > 10.
#
# Inputs
#   out_dir  A directory of per-cell .rds files (default
#            "outputs/design_analysis"). Each file holds either a list with
#            $summary, $bf_results and $diagnostics for a cell that fitted, or a
#            one-row tibble carrying a status for a cell that did not.
#
# Outputs (written by write_design_summary into out_dir)
#   failure_summary.csv        Counts by status, plus the overall failure rate.
#   bf_exceedance.csv          P(BF > threshold) per design point. The main result.
#   prior_sensitivity.csv      How far a single simulated data set's BF moves
#                              when only the prior regime changes.
#   ladder_selection.csv       How often each model-ladder level was run.
#   maximal_feasible_model.csv Highest convergent ladder level per language.
#   recommended_sample_size.csv Smallest N meeting the target, per language.
#   runtime_summary.csv        Fit times, for planning HPC walltime.
#
# Thresholds
#   The BF thresholds default to 10 (primary) and 3 (secondary), matching
#   config/analysis_config.yaml design_analysis$bf_threshold_primary and
#   $bf_threshold_secondary. They are function arguments rather than constants so
#   a sensitivity check can vary them without editing this file, but the
#   preregistered value is 10.
#
# A note on interpreting the output
#   An exceedance proportion is only as stable as the number of seeds behind it.
#   The number of rows contributing to each estimate is reported as `n_sims` in
#   bf_exceedance.csv, and should be read before the proportion beside it.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
})

source("R/05_hypothesis_tests.R")

#' Read every per-cell .rds file in a directory into one long data frame.
#'
#' @param out_dir Directory of per-cell .rds files.
#' @return A data frame with one row per hypothesis per cell, so a cell testing
#'   four hypotheses contributes four rows. Cells that failed contribute a single
#'   row carrying their status and no BF columns; dplyr::bind_rows() fills the
#'   missing columns with NA. Every downstream summary therefore filters on
#'   status == "success" before computing anything.
#' @details Failures are recorded rather than dropped. A design analysis in which
#'   a third of the cells hit the walltime is a different result from one in
#'   which they all fitted, and that distinction would be invisible if unreadable
#'   files were skipped silently. Three shapes are recognised, in order: a
#'   complete result list, a status-only tibble from a cell that errored, and
#'   anything else, which is labelled "malformed". A file that cannot be
#'   deserialised at all becomes "load_error" with a warning; this happens when a
#'   job was killed midway through writing, although the atomic
#'   write-to-.tmp-then-rename in the cell runners makes it rare.
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

#' Bayes-factor exceedance probability per design point. This is the headline
#' result of the design analysis: for a study of a given size, the proportion of
#' simulated data sets in which the evidence reached the threshold.
#'
#' @param df Combined cell data frame from load_design_cells().
#' @param group_vars Columns defining a design point. Everything not named here
#'   is averaged over, so the default deliberately retains `hypothesis`: H1a and
#'   H1b have different power and must not be pooled.
#' @param bf_threshold Primary BF threshold; 10 as preregistered.
#' @param bf_threshold_secondary Secondary, more lenient threshold; 3.
#' @return One row per design point, with `n_sims` (seeds behind the estimate),
#'   `p_bf_primary` and `p_bf_secondary` (the exceedance proportions),
#'   the median and mean BF, and `p_convergence_ok`.
#' @details Both the median and the mean BF are kept because the BF distribution
#'   over seeds is heavily right-skewed: a handful of simulated data sets produce
#'   enormous Bayes factors that dominate the mean. The median describes the
#'   typical study, and a mean far above it is a signal of that skew rather than
#'   of a stronger effect. `p_convergence_ok` is reported alongside so that an
#'   apparently high exceedance rate resting on poorly converged fits can be
#'   spotted; na.rm = TRUE means these proportions are computed over cells with a
#'   usable value, so a low `n_sims` deserves attention.
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

#' How much the conclusion moves when only the prior changes.
#'
#' @param df Combined cell data frame.
#' @return One row per (n_participants, n_verbs, hypothesis, seed) combination,
#'   with the number of prior regimes compared, the number of distinct evidence
#'   categories they produced, the largest and smallest BF, their ratio, and
#'   `category_stable`.
#' @details Grouping includes `seed`, which is what makes this a sensitivity
#'   analysis rather than a power analysis: holding the simulated data set fixed
#'   and varying only the prior isolates the prior's contribution. `category_stable`
#'   is the interpretable quantity, since a BF moving from 40 to 120 changes the
#'   number but not the conclusion, whereas one moving from 8 to 12 crosses the
#'   preregistered threshold. `pmax(min_bf, 1e-6)` guards the ratio against a
#'   division by zero when a regime yields a BF that underflows.
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

#' How often each model-ladder level appears among the successful cells.
#'
#' @param df Combined cell data frame.
#' @return One row per model level, with counts and proportions.
#' @details This counts the ladder levels that were *run and succeeded*, which is
#'   not the same as an automatic fallback frequency: the design grid specifies
#'   the level for each cell, so the mix reflects both the grid's composition and
#'   which levels survived fitting. Read it as a feasibility profile, not as a
#'   model-selection result. For the selection result see maximal_feasible_model().
ladder_selection_summary <- function(df) {
  df |>
    dplyr::filter(status == "success") |>
    dplyr::count(model_level, name = "n_selected") |>
    dplyr::mutate(prop_selected = n_selected / sum(n_selected))
}

#' Tally cell outcomes by status.
#'
#' @param df Combined cell data frame.
#' @return A one-row tibble of counts and the overall failure rate.
#' @details Reported first by write_design_summary(), because every other summary
#'   below is conditional on the cells that succeeded and is therefore only as
#'   representative as this tally allows. The "oom" and "timeout" statuses come
#'   from run_model_ladder() in R/09_model_ladder.R, which classifies a failure by
#'   matching its error text. They are therefore only populated for cells that
#'   went through the ladder engine and that died with a message R could catch;
#'   a job killed outright by SLURM leaves no .rds at all and is invisible here.
#'   Counting the expected cells against the files present is the way to detect
#'   that case, and hpc/ arrays are sized so this comparison is possible.
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
#'
#' @param df Combined cell data frame.
#' @return One row per language, naming the level and its rank, or a zero-row
#'   tibble of the same shape when nothing converged. Returning the empty tibble
#'   rather than NULL keeps the downstream join in recommended_sample_size()
#'   working without a special case.
#' @details This implements the repository's standing rule that the reported model
#'   is the most complex one that fits, never the one that is most convenient. The
#'   rank is parsed from the level name, so the naming convention "L<digit>_..."
#'   in R/04_model_formulas.R is load-bearing: a level named otherwise yields
#'   NA_integer_ and drops out of the max(). `convergence_ok %in% TRUE` rather
#'   than `== TRUE` because the column is NA for cells that never fitted, and
#'   `%in%` returns FALSE there where `==` would return NA and keep the row.
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
#' target, evaluated at that language's maximal feasible model under the primary
#' prior regime. NA if the target is not reached anywhere in the grid.
#'
#' @param exc Output of compute_bf_exceedance().
#' @param mfm Output of maximal_feasible_model().
#' @param target Required exceedance probability; 0.80 by convention.
#' @param focal The hypotheses that must *both* be powered. H2 is excluded
#'   because it is a secondary, two-tailed prediction and powering the design for
#'   it would inflate the recommendation beyond what is preregistered.
#' @param regime Prior regime the recommendation is read from. The sensitivity
#'   regimes exist to be compared against this one, not to set the sample size.
#' @return One row per language, or a zero-row tibble of the same shape when the
#'   inputs are empty or lack the join columns.
#' @details Both focal hypotheses must clear the target *at the same design
#'   point*, which is stricter than taking the larger of two per-hypothesis
#'   recommendations, and is the reason for the `dplyr::n() == length(focal)`
#'   check: it rejects a design point where one hypothesis is simply missing from
#'   the grid rather than present and underpowered.
#'
#'   Restricting to the maximal feasible model matters because power is reported
#'   for the model that will actually be used. A simpler ladder level would give a
#'   flatteringly small recommendation by ignoring variance components that the
#'   real analysis estimates.
#'
#' @section Interpretation:
#'   A trustworthy exceedance probability needs many simulations per design point
#'   (config `design_analysis$n_simulations_per_cell`, currently 200). Run with one
#'   seed per cell, `p_bf_primary` can only be 0 or 1, and this function then
#'   reports the smallest N at which a single simulated study happened to succeed.
#'   That is an indicative figure for smoke-testing the pipeline, not a sample-size
#'   recommendation. Check `n_sims` in bf_exceedance.csv before quoting the result.
recommended_sample_size <- function(exc, mfm, target = 0.80,
                                    focal = c("H1a_semantics_positive",
                                              "H1b_active_interaction_negative"),
                                    regime = "primary") {
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
    # A design point qualifies only if every focal hypothesis is present in the
    # grid there AND all of them clear the target.
    dplyr::group_by(language, n_participants, n_verbs) |>
    dplyr::summarise(all_focal_ok = (dplyr::n() == length(focal)) && all(p_bf_primary >= target),
                     .groups = "drop") |>
    dplyr::group_by(language) |>
    dplyr::summarise(
      meets_target = any(all_focal_ok),
      recommended_n_participants = if (any(all_focal_ok)) min(n_participants[all_focal_ok]) else NA_integer_,
      # Report the verb count belonging to the recommended N, since participants
      # and verbs trade off against each other; quoting a marginal minimum of each
      # separately would describe a design point that was never simulated.
      n_verbs = if (any(all_focal_ok)) n_verbs[all_focal_ok][which.min(n_participants[all_focal_ok])] else NA_integer_,
      .groups = "drop"
    ) |>
    dplyr::mutate(target = target, regime = regime)
}

#' Per-cell fit runtimes by language and model level, in minutes.
#'
#' @param df Combined cell data frame.
#' @return One row per (language, model_level), or a zero-row tibble of that
#'   shape when no cell recorded a runtime.
#' @details Exists to size the SLURM walltime requests in hpc/. The maximum
#'   matters more than the median for that purpose, because an array job is
#'   killed on the walltime of its slowest task, so both are reported.
#'   `p_converged` sits alongside them because a level that fits quickly but
#'   rarely converges is not the cheap option it appears to be.
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

#' Compute every summary above and write them as CSVs.
#'
#' @param df Combined cell data frame from load_design_cells().
#' @param out_dir Destination directory; created if absent.
#' @return Invisibly, a named list of the summary tables. When no cell succeeded
#'   the list holds only `fail`, so a caller must not assume the other names are
#'   present.
#' @details The failure summary is written before anything else, so that a run in
#'   which every cell failed still leaves a diagnosable artefact on disk instead
#'   of an empty directory. That case then returns early: the summaries below
#'   would otherwise all be empty tables, which read as "no effect" rather than
#'   "no data". The statuses actually observed are printed to make the reason
#'   visible in the job log.
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
