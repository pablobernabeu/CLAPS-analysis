# targets.R
# {targets} pipeline for the CLAPS Bayesian ordinal workflow.
# Runs reference audit, prior predictive checks, pilot model fitting,
# prior sensitivity analysis, design analysis aggregation, and report rendering.
# NOT used for the SLURM array jobs themselves (those use hpc/ scripts).
# This pipeline is for local/interactive orchestration and CI checks.

library(targets)
library(tarchetypes)

# Source R/ modules
tar_source("R/")

# Global options
tar_option_set(
  packages   = c("brms", "dplyr", "readr", "ggplot2", "yaml", "posterior"),
  format     = "qs",
  error      = "stop",
  workspace_on_error = TRUE
)

cfg <- yaml::read_yaml("config/analysis_config.yaml")

list(

  # ---------------------------------------------------------------------------
  # 0. Reference audit
  # ---------------------------------------------------------------------------
  tar_target(
    reference_audit,
    run_reference_audit(bib_path = "references.bib",
                        out_dir  = "outputs/reference_audit"),
    format = "rds"
  ),

  # ---------------------------------------------------------------------------
  # 1. Load and validate pilot data
  # ---------------------------------------------------------------------------
  tar_target(
    pilot_data_path,
    cfg$pilot_data_path,
    format = "file"
  ),

  tar_target(
    pilot_raw,
    read_raw_data(pilot_data_path)
  ),

  tar_target(
    pilot_validated,
    validate_raw_data(pilot_raw, pilot_data_path)
  ),

  tar_target(
    pilot_split,
    split_pilot_confirmatory(pilot_validated)
  ),

  # ---------------------------------------------------------------------------
  # 2. Prior predictive checks (one per language)
  # ---------------------------------------------------------------------------
  # One branch per language in the config. Note that this iterates over the CONFIG,
  # not over the languages present in the data, so a language listed there before its
  # data arrive still gets a branch. Such a branch yields NULL rather than failing,
  # for the reason given below.
  #
  # This file is the one public copy that scrub_to_public.sh does not regenerate (it
  # is listed in PUBLIC_SPECIFIC, because the public pipeline is deliberately smaller
  # than the private one). That exemption is about scope, not correctness, so the two
  # fixes below were applied here by hand to match the private source.
  tar_map(
    values = tibble::tibble(
      language = names(cfg$languages),
      has_pp   = purrr::map_lgl(cfg$languages, "has_pseudo_passive")
    ),
    tar_target(
      pilot_lang_df,
      {
        d <- dplyr::filter(pilot_split$pilot, Language == language)
        # A config language with no data yet filters to zero rows. Return NULL
        # instead of calling preprocess_data(), which stops in code_s_type()
        # because no Passive level is present. Combined with error = "stop", that
        # aborted the ENTIRE pipeline rather than the one branch with nothing to do.
        if (nrow(d) == 0L) NULL else preprocess_data(d, has_pseudo_passive = has_pp)
      }
    ),
    tar_target(
      threshold_params,
      # Ceiling calibration is a property of how participants used the response
      # scale, so it applies to any language with data. The gate here was previously
      # has_pseudo_passive, which is unrelated: whether a language has
      # pseudo-passives says nothing about whether its ratings pile up at the top of
      # the scale. That wrongly excluded Norwegian. The only real precondition is
      # that the branch has data.
      if (is.null(pilot_lang_df)) NULL
      else compute_ceiling_calibrated_thresholds(pilot_lang_df, language)
    )
  ),

  # ---------------------------------------------------------------------------
  # 3. Design grid
  # ---------------------------------------------------------------------------
  tar_target(
    design_grid,
    readr::read_csv("config/design_grid.csv", show_col_types = FALSE)
  ),

  # ---------------------------------------------------------------------------
  # 4. Aggregate design results (run after SLURM array completes)
  # ---------------------------------------------------------------------------
  tar_target(
    design_cells_loaded,
    load_design_cells("outputs/design_analysis"),
    format = "rds"
  ),

  tar_target(
    design_summary,
    write_design_summary(design_cells_loaded, out_dir = "outputs/design_summary")
  ),

  # ---------------------------------------------------------------------------
  # 5. Manifest
  # ---------------------------------------------------------------------------
  tar_target(
    manifest,
    write_manifest("outputs/manifest.csv")
  )

  # The report is produced outside this pipeline. The current report draws on the
  # private pilot data, so it is rendered locally (reports/render_report.sh) and
  # committed as a PDF rather than built here.

)
