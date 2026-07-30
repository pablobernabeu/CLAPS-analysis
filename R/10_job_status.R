# R/10_job_status.R
#
# Purpose
#   Answer "how far has this run got?" for a SLURM array job, and record the
#   provenance of whatever it produced. Progress is inferred from the output
#   files on disk rather than from SLURM's own accounting, because a task can be
#   reported as COMPLETED by SLURM while its cell errored inside R, and because
#   sacct records expire from the accounting database after a site-configured
#   retention window whereas the .rds files persist.
#
# Entry points
#   list_completed_cells()  What is on disk.
#   check_grid_completion() On disk, compared with what the grid asked for.
#   query_slurm_status()    Optional sacct cross-check for one job ID.
#   write_status_report()   Progress CSV.
#   write_manifest()        Provenance CSV (git SHA, R and package versions).
#
# Used by
#   scripts/08_submit_status_report.R, and the polling helpers
#   scripts/poll_arc_status.sh and scripts/parse_arc_status.sh.
#
# Caveat on judging success
#   A cell's .rds file existing means the cell *finished*, not that it fitted.
#   A cell that errored still writes an .rds carrying status = "error". Read the
#   status column via R/08_summarise_design.R for that distinction; this module
#   deliberately reports only completion.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
})

#' List the design-analysis cells that have finished, from the files on disk.
#'
#' @param out_dir Directory the cell runners write into.
#' @return A tibble with one row per .rds file, holding the filename, the
#'   `cell_id` (the filename without its extension, which is how the cell
#'   runners construct it), and the modification time. Returns an empty tibble,
#'   not an error, when nothing has completed: a run that has only just started
#'   is a normal state for a progress check.
#' @details `mtime` is included because it is the cheapest way to see whether a
#'   long array is still making progress or has stalled.
list_completed_cells <- function(out_dir = "outputs/design_analysis") {
  # Cells are written with saveRDS (see R/06_simulate_design.R), so the
  # extension is .rds. This globbed .qs until 2026-07, which silently matched
  # nothing and reported every run as having produced no cells.
  files <- list.files(out_dir, pattern = "\\.rds$", full.names = FALSE)
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
#'
#' @param design_grid Tibble; the full design grid, one row per requested cell.
#' @param completed Tibble from list_completed_cells().
#' @return The design grid with `completed` (logical) and `status` ("done" or
#'   "pending") added. A left join from the grid, so every requested cell appears
#'   whether or not it ran; `coalesce(completed, FALSE)` turns the join's NAs for
#'   unmatched rows into an explicit FALSE.
#'
#' @section Known limitation:
#'   The expected ID is rebuilt from the seven fields below, mirroring the base
#'   naming scheme in run_design_cell() (R/06_simulate_design.R). That function
#'   additionally appends "_gender" or "_genderX" for the gender model variations
#'   and "_<k>lang" for cross-language cells with an explicit language count.
#'   Those suffixes are not reproduced here, so cells of those kinds never match
#'   and are reported as "pending" even once they have finished. Progress figures
#'   for the gender and cross-language grids are therefore understated; the
#'   baseline single-language grids, which is what this report was written for,
#'   are unaffected. Counting the .rds files in the output directory is the
#'   reliable check for the affected grids until the two naming schemes are
#'   unified.
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

#' Query SLURM's accounting database for one job ID, via sacct.
#'
#' @param job_id A SLURM job ID, with or without an array suffix.
#' @return A tibble of one row per task with JobID, State, ExitCode, Elapsed and
#'   MaxRSS, or NULL when sacct is unavailable or returns nothing. NULL rather
#'   than an error, so the same scripts run unchanged on a laptop with no SLURM
#'   installed; callers must therefore handle NULL.
#' @details MaxRSS and Elapsed are requested because they are what the walltime
#'   and memory requests in hpc/ are tuned against: a task whose MaxRSS approaches
#'   the requested memory is the one that will be killed when the model grows.
#'   `--parsable2` gives pipe-delimited output with no padding and no trailing
#'   delimiter, and `--noheader` suppresses the header, so the column names below
#'   are supplied explicitly and stay correct regardless of the site's default
#'   sacct format.
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

#' Write the per-cell progress table to CSV and print a one-line tally.
#'
#' @param grid_status Output of check_grid_completion().
#' @param out_path Destination CSV; its directory is created if absent.
#' @return Invisibly, `out_path`.
#' @details The console tally exists so the numbers land in the SLURM job log,
#'   where they can be read without transferring the CSV off the cluster.
write_status_report <- function(grid_status, out_path = "outputs/job_status_report.csv") {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(grid_status, out_path)
  n_done    <- sum(grid_status$status == "done",    na.rm = TRUE)
  n_pending <- sum(grid_status$status == "pending", na.rm = TRUE)
  message("[status] Done: ", n_done, " | Pending: ", n_pending,
          " | Total: ", nrow(grid_status))
  invisible(out_path)
}

#' Record the provenance of a run: git SHA, timestamp and software versions.
#'
#' @param out_path Destination CSV.
#' @param additional_cols Named list of extra columns to append, for run-specific
#'   context such as the grid file or SLURM job ID.
#' @return Invisibly, the manifest tibble.
#' @details Written next to the outputs so that a result can be traced back to the
#'   exact code and package versions that produced it, which is the minimum needed
#'   for a computational result to be reproducible (Sandve et al., 2013,
#'   doi:10.1371/journal.pcbi.1003285). brms and cmdstanr are recorded by name
#'   because they determine the sampler's behaviour, and a change in either can
#'   move a Bayes factor without any change to this repository's code.
#'
#'   Each lookup is wrapped so that a missing tool degrades to "unknown" rather
#'   than aborting the run: an incomplete manifest is more useful than a job that
#'   died after computing its results. The timestamp keeps its UTC offset (%z), so
#'   runs on a cluster in one timezone remain comparable with local runs.
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
