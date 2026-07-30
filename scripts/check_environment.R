#!/usr/bin/env Rscript
# scripts/check_environment.R
#
# Purpose
#   Decide whether the current session can run the analysis, and say so with a
#   non-zero exit code if it cannot. Intended for three places: a developer's
#   machine before starting, CI before the tests, and the head of a SLURM job
#   before hours of compute are spent on a broken environment.
#
# Usage
#   Rscript scripts/check_environment.R                  # core + modelling
#   Rscript scripts/check_environment.R --groups all     # everything
#   Rscript scripts/check_environment.R --groups core
#   Rscript scripts/check_environment.R --compare config/arc_environment_recorded.csv
#   Rscript scripts/check_environment.R --lockfile renv.lock
#
# Exit codes
#   0  every requested group satisfied
#   1  something required is missing
#
# Deliberately base R only, with no optparse: this script must run in exactly the
# incomplete environments it exists to diagnose, so it cannot depend on a package
# that might itself be missing. That is also why the arguments are parsed by hand.

source("R/12_environment.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

groups_arg <- arg_value("--groups", "core,modelling")
# Braced so the `else` cannot start a new top-level expression, which is a syntax
# error in R outside a block.
groups <- if (identical(groups_arg, "all")) {
  names(CLAPS_REQUIREMENTS)
} else {
  trimws(strsplit(groups_arg, ",")[[1]])
}
compare_path  <- arg_value("--compare")
lockfile_path <- arg_value("--lockfile")

cat("CLAPS environment check\n")
cat("=======================\n")
claps_environment_summary()

# --- requirements ------------------------------------------------------------
cat("\nRequirements\n")
missing <- claps_missing_packages(groups)
for (g in names(missing)) {
  m <- missing[[g]]
  cat(sprintf("  %-10s %s\n", g,
              if (length(m)) paste0("MISSING: ", paste(m, collapse = ", ")) else "ok"))
}
n_missing <- sum(lengths(missing))

# --- optional: drift against a recorded environment --------------------------
# Reported, never fatal. The workstation used for editing is legitimately not the
# cluster that produces results; the point is that the difference is visible and
# ends up in the job log, not that the two are forced to match.
if (!is.null(compare_path) && file.exists(compare_path)) {
  cmp <- claps_compare_to_record(compare_path)
  cat("\nDrift against ", compare_path, "\n", sep = "")
  cat(sprintf("  compared %d packages: %d match, %d differ, %d absent here\n",
              nrow(cmp), sum(cmp$status == "match"),
              sum(cmp$status == "differs"), sum(cmp$status == "absent_here")))
  d <- cmp[cmp$status == "differs", , drop = FALSE]
  if (nrow(d)) {
    show <- utils::head(d[order(d$package), ], 15)
    for (i in seq_len(nrow(show))) {
      cat(sprintf("    %-16s recorded %-12s here %s\n",
                  show$package[i], show$recorded[i], show$here[i]))
    }
    if (nrow(d) > nrow(show)) cat(sprintf("    ... and %d more\n", nrow(d) - nrow(show)))
  }
}

# --- optional: agreement with renv.lock --------------------------------------
# Also reported rather than fatal, because the lockfile currently describes fewer
# packages than the analysis actually uses; see docs/reproducibility.md.
if (!is.null(lockfile_path) && file.exists(lockfile_path) &&
    requireNamespace("jsonlite", quietly = TRUE)) {
  lock <- jsonlite::fromJSON(lockfile_path, simplifyVector = FALSE)
  cat("\nLockfile ", lockfile_path, "\n", sep = "")
  cat(sprintf("  records R %s; this session is R %s%s\n",
              lock$R$Version, as.character(getRversion()),
              if (identical(lock$R$Version, as.character(getRversion()))) "" else "  <- DIFFERS"))
  cat(sprintf("  locks %d packages\n", length(lock$Packages)))
  needed <- unlist(CLAPS_REQUIREMENTS[intersect(c("core", "modelling"), groups)])
  unlocked <- setdiff(needed, names(lock$Packages))
  if (length(unlocked)) {
    cat("  required but NOT locked: ", paste(sort(unlocked), collapse = ", "), "\n", sep = "")
  }
}

cat("\n")
if (n_missing > 0) {
  cat("RESULT: environment incomplete —", n_missing, "package(s) missing.\n")
  cat("See docs/reproducibility.md for the expected environment and how to build it.\n")
  quit(status = 1)
}
cat("RESULT: environment satisfies the requested groups (",
    paste(groups, collapse = ", "), ").\n", sep = "")
quit(status = 0)
