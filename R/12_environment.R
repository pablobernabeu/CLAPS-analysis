# R/12_environment.R
#
# Purpose
#   One authoritative statement of what the analysis needs from its environment,
#   and the machinery to check a running session against it. Loaded by the test
#   suite, by scripts/check_environment.R and by scripts/record_environment.R, so
#   that "are we in a usable environment?" has a single answer rather than a
#   different implicit one in each place.
#
# Why this exists
#   Before it, a machine missing brms ran the test suite and reported success: the
#   three test files that need it failed to LOAD, and testthat reports a load
#   failure as a skip. Twenty-nine of fifty-three checks silently did not run. A
#   green suite therefore carried almost no information unless the reader also
#   inspected the skip count, which is not a reasonable thing to require.
#
# Deliberately dependency-free
#   Base R only. This file has to be loadable in exactly the degraded environments
#   it is meant to diagnose, so it cannot rely on the packages it is checking for.

# ---------------------------------------------------------------------------
# What the analysis requires
# ---------------------------------------------------------------------------
# Split by what breaks without them, so a partial environment can be described
# precisely rather than as a single pass/fail.
#
#   core      Needed to load and unit-test the deterministic logic: data
#             validation, factor coding, priors, the Bayes-factor arithmetic,
#             diagnostics classification.
#   modelling Needed to LOAD the analysis modules and therefore to run the test
#             suite at all. Kept deliberately narrow — brms and posterior — because
#             this is the group the test gate enforces, and demanding more would
#             block environments that can perfectly well run every test.
#   fitting   Needed to actually fit a model, over and above loading the code.
#             Separate from `modelling` because cmdstanr is not on CRAN (it comes
#             from the Stan repository at https://mc-stan.org/r-packages/), so a
#             CI runner reasonably has brms without it, and the unit tests never
#             sample anything.
#   pipeline  Needed by the runnable scripts and the {targets} orchestration, but
#             not by the unit tests.
#   reporting Needed to render the report and audit the bibliography.
#
# The split matters: the test gate checks core + modelling, whereas an ARC node
# about to spend days fitting should check everything, with --groups all.
CLAPS_REQUIREMENTS <- list(
  core      = c("testthat", "here", "withr", "dplyr", "tibble", "forcats",
                "readr", "assertr", "yaml", "purrr", "stringr", "tidyr"),
  modelling = c("brms", "posterior"),
  fitting   = c("cmdstanr", "MASS", "ordinal", "bridgesampling"),
  pipeline  = c("optparse", "targets", "tarchetypes"),
  reporting = c("quarto", "bib2df", "httr2", "ggplot2")
)

#' Which required packages are missing, by group.
#'
#' @param groups Requirement groups to check; defaults to all of them.
#' @return A named list, one character vector of missing package names per group.
#'   Empty vectors mean the group is satisfied.
#' @details Uses requireNamespace() rather than installed.packages(), because a
#'   package can be listed as installed yet fail to load — a broken shared library
#'   after a cluster upgrade is the usual cause, and it is precisely the case that
#'   must not be reported as present.
claps_missing_packages <- function(groups = names(CLAPS_REQUIREMENTS)) {
  groups <- intersect(groups, names(CLAPS_REQUIREMENTS))
  out <- lapply(groups, function(g) {
    pk <- CLAPS_REQUIREMENTS[[g]]
    pk[!vapply(pk, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
  })
  names(out) <- groups
  out
}

#' Describe the running environment in one printable block.
#'
#' @return Invisibly, a named list of the facts printed.
#' @details Printed at the head of a check or a job so the log records the
#'   environment that produced whatever follows it. Every value is read from the
#'   session rather than assumed, including the CmdStan path, which is the item
#'   most often silently different between a workstation and a cluster node.
claps_environment_summary <- function() {
  cmdstan <- Sys.getenv("CMDSTAN", "")
  facts <- list(
    r_version   = as.character(getRversion()),
    platform    = R.version$platform,
    os          = paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
    r_libs      = paste(.libPaths(), collapse = .Platform$path.sep),
    cmdstan     = if (nzchar(cmdstan)) cmdstan else "(CMDSTAN unset)",
    stan_threads = Sys.getenv("STAN_NUM_THREADS", "(unset)"),
    slurm_job   = Sys.getenv("SLURM_JOB_ID", "(not a SLURM job)"),
    r_module    = Sys.getenv("ARC_R_MODULE", "(unset)")
  )
  for (n in names(facts)) cat(sprintf("  %-13s %s\n", n, facts[[n]]))
  invisible(facts)
}

#' Compare the running session against a recorded environment.
#'
#' @param record_path CSV of Package,Version as written by
#'   scripts/record_environment.R. Comment lines beginning "#" are skipped.
#' @param packages Restrict the comparison to these packages; NULL compares every
#'   package named in the record that is also installed here.
#' @return A data frame with one row per compared package and a `status` of
#'   "match", "differs" or "absent_here".
#' @details Reports rather than judges. Two environments legitimately differ — the
#'   workstation used for editing is not the cluster that produces results — so the
#'   value of this comparison is in making the difference visible and recorded,
#'   not in forcing them to be identical. scripts/check_environment.R decides what
#'   to do about the result.
claps_compare_to_record <- function(record_path, packages = NULL) {
  stopifnot(file.exists(record_path))
  rec <- utils::read.csv(record_path, comment.char = "#", stringsAsFactors = FALSE)
  if (!all(c("Package", "Version") %in% names(rec))) {
    stop("[env] ", record_path, " must have Package and Version columns.")
  }
  if (!is.null(packages)) rec <- rec[rec$Package %in% packages, , drop = FALSE]

  inst <- utils::installed.packages()
  here_v <- ifelse(rec$Package %in% rownames(inst),
                   unname(inst[match(rec$Package, rownames(inst)), "Version"]),
                   NA_character_)
  data.frame(
    package        = rec$Package,
    recorded       = rec$Version,
    here           = here_v,
    status         = ifelse(is.na(here_v), "absent_here",
                     ifelse(here_v == rec$Version, "match", "differs")),
    stringsAsFactors = FALSE
  )
}
