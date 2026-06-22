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
  tar_map(
    values = tibble::tibble(
      language = names(cfg$languages),
      has_pp   = purrr::map_lgl(cfg$languages, "has_pseudo_passive")
    ),
    tar_target(
      pilot_lang_df,
      preprocess_data(
        dplyr::filter(pilot_split$pilot, Language == language),
        has_pseudo_passive = has_pp
      )
    ),
    tar_target(
      threshold_params,
      if (isTRUE(cfg$languages[[language]]$has_pseudo_passive))
        compute_ceiling_calibrated_thresholds(pilot_lang_df, language)
      else NULL
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
