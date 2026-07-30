# R/09_model_ladder.R
#
# Purpose
#   Fit the most complex random-effects structure the data will support, by
#   working down a prespecified ladder of models and stopping at the first level
#   that both fits and converges to publication standard.
#
# Why a ladder rather than one chosen model
#   A maximal random-effects structure protects the Type-I error rate of the
#   fixed effects, but frequently fails to converge on realistic sample sizes
#   (Barr et al., 2013, doi:10.1016/j.jml.2012.11.001; Matuschek et al., 2017,
#   doi:10.1016/j.jml.2017.01.001, on the resulting power cost). Deciding in
#   advance which simplifications are permitted, and in what order, keeps that
#   trade-off from becoming a data-dependent choice made after seeing the
#   results. The levels are defined in R/04_model_formulas.R and listed in the
#   README's model-ladder table.
#
# The rule this file enforces
#   Descend only on failure, never for convenience. A level is abandoned only if
#   it errors, exhausts memory or walltime, or fails the convergence criteria in
#   FALLBACK_TRIGGERS below. A level that converges is selected even if a simpler
#   one would have been faster or would have given a tidier result.
#
# Inputs and outputs
#   Takes one language's preprocessed data frame and a prior object; returns the
#   selected level, its fit and a per-level log. Fits are saved as
#   <out_dir>/<language>_<prior_label>_<level>.rds and diagnostics as
#   <out_dir>/diagnostics/<...>_diag.csv, the latter for every level attempted,
#   including those that failed, so the descent can be audited afterwards.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(purrr)
})

source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/07_extract_diagnostics.R")

# The exhaustive list of reasons a level may be abandoned. Anything not on this
# list is not grounds for descending the ladder. Declared as a named constant so
# the permitted reasons are stated in one place and can be quoted in the
# preregistration rather than being inferred from the control flow below.
FALLBACK_TRIGGERS <- c(
  "timeout",           # SLURM wall-time exceeded
  "oom",               # Out-of-memory error
  "compilation_error", # Stan compilation failed
  # not publication-grade: R-hat >= 1.01, any divergence /
  # treedepth saturation, or bulk/tail ESS < 400
  "convergence_failed"
)

#' Run the model ladder for a single language's data, applying fallback logic.
#'
#' @param df Preprocessed data frame for one language, as returned by
#'   R/02_preprocess_factors.R. Factor contrasts must already be set: this
#'   function does not touch them, and the priors are matched to coefficient names
#'   that depend on them.
#' @param prior_obj brms prior object from build_brms_prior().
#' @param sampling_args List passed on to brms::brm (iter, warmup, chains, cores,
#'   seed). Supplied by the caller rather than defaulted here so that the heavy
#'   convergence-demonstration sampler and the lighter per-replicate sampler in
#'   config/analysis_config.yaml can both drive the same engine.
#' @param control_args List passed to brm's `control` (adapt_delta, max_treedepth).
#' @param has_pseudo_passive Logical; FALSE drops the Pseudo_Passive terms for
#'   languages without that construction.
#' @param start_level Ladder level to begin at, defaulting to the maximal L5.
#'   Lower values exist for resuming a run, not for choosing a simpler model.
#' @param out_dir Directory for saved fits and diagnostics.
#' @param language,prior_label Labels used to build the output filenames.
#' @param include_gender Logical; if TRUE, fit the gender model variation. The
#'   data frame must then contain a Gender column.
#' @param overwrite Logical; if FALSE, an existing converged fit at a level is
#'   reused instead of refitted.
#' @return A list with `selected_model` (the level name, or NA if every level
#'   failed), `fit` (the brmsfit, or NULL) and `ladder_log`, a per-level record of
#'   status, diagnostics and runtime. Note that the returned list has no
#'   `diagnostics` element; per-level diagnostics live inside `ladder_log` and are
#'   also written to CSV.
#' @details Failing at every level produces a warning and a NULL fit rather than
#'   an error, so that a batch covering several languages is not aborted by one
#'   language that will not fit. Callers must check for a NULL fit.
#'
#'   Reusing an existing fit is what makes a resubmitted job cheap after a
#'   walltime kill. A reused fit is re-diagnosed rather than trusted, so a fit
#'   left behind by an earlier run under different convergence criteria cannot be
#'   silently accepted.
run_model_ladder <- function(df,
                             prior_obj,
                             sampling_args,
                             control_args,
                             has_pseudo_passive = TRUE,
                             start_level        = "L5_correlated_maximal",
                             out_dir            = "outputs/models",
                             language           = "unknown",
                             prior_label        = "primary",
                             include_gender     = FALSE,
                             overwrite          = FALSE) {
  ladder      <- build_model_ladder(has_pseudo_passive = has_pseudo_passive,
                                    include_gender = include_gender)
  level_order <- ladder_names()
  start_idx   <- match(start_level, level_order)

  if (is.na(start_idx)) stop("[ladder] Unknown start_level: ", start_level)

  ladder_log <- list()
  selected_fit  <- NULL
  selected_level <- NA_character_

  for (idx in seq(start_idx, length(level_order))) {
    level_name <- level_order[idx]
    formula    <- ladder[[level_name]]

    out_rds  <- file.path(out_dir, paste0(language, "_", prior_label, "_", level_name, ".rds"))
    diag_csv <- file.path(out_dir, "diagnostics",
                          paste0(language, "_", prior_label, "_", level_name, "_diag.csv"))
    dir.create(file.path(out_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)

    # Resume path: a fit already on disk is re-diagnosed against the current
    # criteria, not assumed good. If it does not meet them, the ladder descends
    # exactly as it would have done had the fit just been computed.
    if (file.exists(out_rds) && !overwrite) {
      message("[ladder] Loading existing fit: ", level_name)
      fit <- tryCatch(readRDS(out_rds), error = function(e) NULL)
      if (!is.null(fit)) {
        diag <- extract_convergence_diagnostics(fit)
        status <- classify_convergence(diag)
        ladder_log[[level_name]] <- list(status = status, diag = diag)
        if (status == "converged") {
          selected_fit   <- fit
          selected_level <- level_name
          break
        } else {
          message("[ladder] Existing fit at ", level_name, " has status: ", status,
                  " — falling back.")
          next
        }
      }
    }

    message("[ladder] Fitting: ", level_name)
    t_start <- proc.time()[["elapsed"]]

    fit_attempt <- tryCatch({
      fit <- do.call(brms::brm, c(
        list(
          formula      = formula,
          data         = df,
          # Each ladder level has fewer parameters than the one above, so the
          # full prior object names coefficients that a lower level does not
          # have. brms treats an unmatched prior as an error, so the prior is
          # filtered to this level's parameters at each step.
          prior        = align_prior_to_model(prior_obj, formula, df),
          backend      = "cmdstanr",
          sample_prior = "yes",
          control      = control_args
        ),
        sampling_args
      ))
      list(fit = fit, status = "success", error = NULL)
    }, error = function(e) {
      # Classify the failure from its message, so that the ladder log distinguishes
      # a resource limit from a genuine modelling problem. case_when() takes the
      # first match, and the order is deliberate: memory first, because an
      # allocation failure during compilation is a memory problem rather than a
      # compilation one; a compile failure next, since it means the model was never
      # sampled at all and no amount of extra walltime would help; then timeout.
      # Anything unrecognised stays "error" rather than being forced into a
      # category, which keeps a novel failure visible in the log.
      msg <- conditionMessage(e)
      status <- dplyr::case_when(
        grepl("out of memory|OOM|cannot allocate", msg, ignore.case = TRUE) ~ "oom",
        grepl("compilation|compile", msg, ignore.case = TRUE) ~ "compilation_error",
        grepl("timeout|time.out|timed out", msg, ignore.case = TRUE)       ~ "timeout",
        TRUE                                                                ~ "error"
      )
      list(fit = NULL, status = status, error = msg)
    })

    t_end <- proc.time()[["elapsed"]]
    runtime <- t_end - t_start

    if (fit_attempt$status != "success") {
      ladder_log[[level_name]] <- list(
        status = fit_attempt$status, error = fit_attempt$error,
        runtime_sec = runtime
      )
      message("[ladder] ", level_name, " failed (", fit_attempt$status, "): ",
              substr(fit_attempt$error, 1, 200))
      next
    }

    fit  <- fit_attempt$fit
    diag <- extract_convergence_diagnostics(fit)
    conv_status <- classify_convergence(diag)

    if (conv_status != "converged") {
      # Publication-grade acceptance: only a strictly converged fit (R-hat < 1.01,
      # ESS >= 400, zero divergences/treedepth) is selected. "marginal_*" fits
      # (e.g. R-hat in [1.01, 1.05)) now trigger a fallback rather than being
      # accepted, so the reported model meets the standard in §9 of the prereg.
      ladder_log[[level_name]] <- list(
        status = "convergence_failed", diag = diag, runtime_sec = runtime,
        convergence_status = conv_status
      )
      message("[ladder] ", level_name, " not publication-grade converged (",
              conv_status, ") — falling back.")

      # Keep the diagnostics of the rejected level. They are the evidence that
      # the descent was warranted, and are needed to report why the maximal model
      # was not used.
      save_diagnostics_csv(dplyr::mutate(diag, model_level = level_name,
                                         convergence_status = conv_status), diag_csv)
      next
    }

    # Converged. Write to a temporary name and rename into place, because rename
    # within a filesystem is atomic while saveRDS is not: a job killed partway
    # through a direct save would leave a truncated .rds that the resume path
    # above would later try to read. The rename means the final path either does
    # not exist or holds a complete fit.
    tmp_rds <- paste0(out_rds, ".tmp")
    saveRDS(fit, tmp_rds)
    file.rename(tmp_rds, out_rds)

    save_diagnostics_csv(dplyr::mutate(diag, model_level = level_name,
                                       convergence_status = conv_status), diag_csv)

    ladder_log[[level_name]] <- list(
      status = "converged", diag = diag, runtime_sec = runtime,
      convergence_status = conv_status
    )
    selected_fit   <- fit
    selected_level <- level_name
    message("[ladder] Selected: ", level_name, " (", round(runtime, 1), "s)")
    break
  }

  if (is.null(selected_fit)) {
    warning("[ladder] All model levels failed for ", language, " / ", prior_label)
  }

  list(
    selected_model = selected_level,
    fit            = selected_fit,
    ladder_log     = ladder_log
  )
}
