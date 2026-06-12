# R/09_model_ladder.R
# Model-ladder execution engine.
# Tries each level from L5 to L0, applying prespecified fallback criteria.
# Never skips to a simpler model because it is easier to fit unless the
# higher model has failed by prespecified criteria.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(purrr)
})

source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/07_extract_diagnostics.R")

# Fallback trigger criteria
FALLBACK_TRIGGERS <- c(
  "timeout",          # SLURM wall-time exceeded
  "oom",              # Out-of-memory error
  "compilation_error",# Stan compilation failed
  "convergence_failed"# not publication-grade: R-hat >= 1.01, any divergence /
                      # treedepth saturation, or bulk/tail ESS < 400
)

#' Run the model ladder for a single language-dataset, applying fallback logic.
#' @param df Preprocessed data frame for one language.
#' @param prior_obj brms prior object.
#' @param sampling_args List; passed to brms::brm.
#' @param control_args List; passed to brms::brm control argument.
#' @param has_pseudo_passive Logical.
#' @param start_level Character; starting ladder level (default L5).
#' @param out_dir Character; directory to save fit objects atomically.
#' @param language Character; language label for file naming.
#' @param prior_label Character; prior regime label.
#' @param include_gender Logical; if TRUE, fit the gender model variation
#'   (adds the Gender covariate). The data frame must contain a Gender column.
#' @param overwrite Logical.
#' @return List with slots: selected_model, fit, diagnostics, ladder_log.
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

    # Check for existing fit
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
          prior        = align_prior_to_model(prior_obj, formula, df),
          backend      = "cmdstanr",
          sample_prior = "yes",
          control      = control_args
        ),
        sampling_args
      ))
      list(fit = fit, status = "success", error = NULL)
    }, error = function(e) {
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

      # Save diagnostics even for non-converged fits
      save_diagnostics_csv(dplyr::mutate(diag, model_level = level_name,
                                         convergence_status = conv_status), diag_csv)
      next
    }

    # Converged — save atomically
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
