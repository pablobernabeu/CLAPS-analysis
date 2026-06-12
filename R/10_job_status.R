# R/10_job_status.R
# Monitor SLURM array job progress, collect output manifests, and
# report which design cells are complete, pending, failed or missing.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
})

#' List completed design-analysis output files and extract cell metadata from filenames.
list_completed_cells <- function(out_dir = "outputs/design_analysis") {
  files <- list.files(out_dir, pattern = "\\.qs$", full.names = FALSE)
  if (length(files) == 0) {
    message("[status] No completed cells in ", out_dir)
    return(tibble::tibble())
  }
  tibble::tibble(
    filename    = files,
    cell_id     = tools::file_path_sans_ext(files),
    completed   = TRUE,
    mtime       = file.mtime(file.path(out_dir, files))
  )
}

#' Compare completed cells against the design grid to find gaps.
#' @param design_grid Tibble; the full design grid.
#' @param completed Tibble from list_completed_cells().
#' @return Tibble with completion status per grid row.
check_grid_completion <- function(design_grid, completed) {
  design_grid <- dplyr::mutate(design_grid,
    expected_cell_id = paste(language, model_level, prior_regime,
                             threshold_mode, n_participants, n_verbs, seed,
                             sep = "_")
  )
  design_grid |>
    dplyr::left_join(
      dplyr::select(completed, cell_id, completed, mtime),
      by = c("expected_cell_id" = "cell_id")
    ) |>
    dplyr::mutate(
      completed = dplyr::coalesce(completed, FALSE),
      status    = dplyr::if_else(completed, "done", "pending")
    )
}

#' Query SLURM job status for a given job ID via sacct.
#' Returns NULL silently if sacct is not available (e.g., on Windows dev machine).
query_slurm_status <- function(job_id) {
  if (!nchar(Sys.which("sacct"))) {
    message("[status] sacct not available; skipping SLURM query.")
    return(NULL)
  }
  result <- tryCatch(
    system2("sacct",
            args    = c("-j", job_id, "--format=JobID,State,ExitCode,Elapsed,MaxRSS",
                        "--noheader", "--parsable2"),
            stdout  = TRUE,
            stderr  = FALSE),
    error = function(e) NULL
  )
  if (is.null(result) || length(result) == 0) return(NULL)
  readr::read_delim(
    paste(result, collapse = "\n"),
    delim = "|",
    col_names = c("JobID", "State", "ExitCode", "Elapsed", "MaxRSS"),
    show_col_types = FALSE
  )
}

#' Write a status report CSV.
write_status_report <- function(grid_status, out_path = "outputs/job_status_report.csv") {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(grid_status, out_path)
  n_done    <- sum(grid_status$status == "done",    na.rm = TRUE)
  n_pending <- sum(grid_status$status == "pending", na.rm = TRUE)
  message("[status] Done: ", n_done, " | Pending: ", n_pending,
          " | Total: ", nrow(grid_status))
  invisible(out_path)
}

#' Write outputs/manifest.csv with Git SHA, datetime, software versions.
write_manifest <- function(out_path = "outputs/manifest.csv",
                            additional_cols = list()) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  git_sha <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "unknown"
  )

  r_version <- paste0(R.version$major, ".", R.version$minor)

  brms_version <- tryCatch(
    as.character(utils::packageVersion("brms")), error = function(e) "unknown"
  )
  cmdstanr_version <- tryCatch(
    as.character(utils::packageVersion("cmdstanr")), error = function(e) "unknown"
  )

  manifest <- tibble::tibble(
    git_sha          = git_sha,
    datetime         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version        = r_version,
    brms_version     = brms_version,
    cmdstanr_version = cmdstanr_version
  )

  for (nm in names(additional_cols)) {
    manifest[[nm]] <- additional_cols[[nm]]
  }

  readr::write_csv(manifest, out_path)
  message("[manifest] Written to ", out_path)
  invisible(manifest)
}
