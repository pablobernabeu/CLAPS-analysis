#!/usr/bin/env Rscript
# scripts/record_environment.R
#
# Purpose
#   Write down exactly what produced a result: R version, platform, OS, compiler,
#   CmdStan version, every installed package version, the SLURM context and the git
#   commit. Provenance of this kind is the minimum a computational result needs to
#   be reproducible (Sandve et al., 2013, doi:10.1371/journal.pcbi.1003285), and it
#   is far cheaper to capture at run time than to reconstruct afterwards.
#
# Usage
#   Rscript scripts/record_environment.R                       # to outputs/
#   Rscript scripts/record_environment.R --out config/arc_environment_recorded.csv
#   Rscript scripts/record_environment.R --format json
#
# Two outputs, deliberately
#   The CSV is a package table, diffable in review and readable by
#   scripts/check_environment.R. The JSON carries the same package table plus the
#   scalar facts (versions, toolchain, job identifiers) that do not fit a two-column
#   table. Neither is derived from the other at read time, so both are written.
#
# Where to run it
#   On ARC, inside a job or after loading the same modules a job would, so the
#   record describes the environment that actually fits the models. Running it on a
#   workstation records the workstation, which is useful for diagnosing drift but is
#   not the environment of record.
#
# Base R only, for the same reason as scripts/check_environment.R: it must run in a
# degraded environment, since that is often exactly when a record is wanted.

source("R/12_environment.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

out_path <- arg_value("--out", file.path("outputs", "environment_record.csv"))
format   <- arg_value("--format", "both")

# --- scalar facts ------------------------------------------------------------
git_sha <- tryCatch(
  paste(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
        collapse = ""),
  error = function(e) "unknown"
)
cmdstan_dir <- Sys.getenv("CMDSTAN", "")
cmdstan_ver <- if (nzchar(cmdstan_dir)) basename(cmdstan_dir) else "(CMDSTAN unset)"
cc <- tryCatch(
  paste(system2("gcc", "--version", stdout = TRUE, stderr = FALSE)[1], collapse = ""),
  error = function(e) "unknown"
)

facts <- list(
  recorded_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  git_sha        = git_sha,
  r_version      = as.character(getRversion()),
  platform       = R.version$platform,
  os             = paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
  compiler       = cc,
  cmdstan        = cmdstan_ver,
  r_module       = Sys.getenv("ARC_R_MODULE", "(unset)"),
  r_libs         = paste(.libPaths(), collapse = .Platform$path.sep),
  stan_threads   = Sys.getenv("STAN_NUM_THREADS", "(unset)"),
  slurm_job_id   = Sys.getenv("SLURM_JOB_ID", "(not a SLURM job)"),
  slurm_array_id = Sys.getenv("SLURM_ARRAY_TASK_ID", "(not an array task)"),
  slurm_account  = Sys.getenv("SLURM_JOB_ACCOUNT", "(unset)")
)
# n_packages and n_shadowed are added below, once the table has been deduplicated,
# so the printed count is the number of packages actually recorded rather than the
# raw row count of installed.packages(), which double-counts shadowed copies.

# --- package table -----------------------------------------------------------
# installed.packages() returns one row per package PER LIBRARY, so a package present
# in both the project library and the R module's site library appears twice with
# different versions. On ARC seven packages did (httr2, knitr, yaml, cpp11, ps,
# rlang, systemfonts). A record listing both is ambiguous about which version the
# analysis actually used, which defeats its purpose.
#
# R resolves that ambiguity by search order: the first .libPaths() entry providing a
# package wins. The record therefore keeps exactly that copy, one row per package,
# and names the library it came from. The shadowed copies are not silently dropped:
# they are listed in the preamble, because a shadowed newer version is a common
# reason for a result differing between two accounts on the same cluster.
inst <- utils::installed.packages()
ord  <- order(match(inst[, "LibPath"], .libPaths()))
inst <- inst[ord, , drop = FALSE]
dup  <- duplicated(inst[, "Package"])

pkgs <- data.frame(
  Package = unname(inst[!dup, "Package"]),
  Version = unname(inst[!dup, "Version"]),
  Library = unname(inst[!dup, "LibPath"]),
  stringsAsFactors = FALSE
)
pkgs <- pkgs[order(tolower(pkgs$Package)), ]

shadowed <- inst[dup, c("Package", "Version", "LibPath"), drop = FALSE]
facts$n_packages <- nrow(pkgs)
facts$n_shadowed <- nrow(shadowed)

cat("Recording environment\n=====================\n")
for (n in names(facts)) cat(sprintf("  %-15s %s\n", n, facts[[n]]))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

if (format %in% c("csv", "both")) {
  # A commented preamble so the file explains itself when read on its own, and is
  # still parseable: read.csv(comment.char = "#") skips these lines.
  con <- file(out_path, open = "wt")
  writeLines(c(
    paste0("# ", out_path),
    "# Environment record written by scripts/record_environment.R.",
    "# A RECORD of what was present, not a specification of what is required;",
    "# the requirements live in R/12_environment.R. See docs/reproducibility.md.",
    "#",
    paste0("#   recorded_at ", facts$recorded_at),
    paste0("#   git_sha     ", facts$git_sha),
    paste0("#   R           ", facts$r_version, " (", facts$platform, ")"),
    paste0("#   os          ", facts$os),
    paste0("#   compiler    ", facts$compiler),
    paste0("#   cmdstan     ", facts$cmdstan),
    paste0("#   r_module    ", facts$r_module),
    paste0("#   slurm_job   ", facts$slurm_job_id),
    "#",
    paste0("# Libraries searched, in the order R resolves them:"),
    paste0("#   ", seq_along(.libPaths()), ". ", .libPaths()),
    "#",
    if (nrow(shadowed)) {
      paste0("# Shadowed copies (present in a later library, so NOT the version used):")
    } else "# No package was shadowed by a later library.",
    if (nrow(shadowed)) {
      paste0("#   ", shadowed[, "Package"], " ", shadowed[, "Version"],
             "  in ", shadowed[, "LibPath"])
    } else character(0)
  ), con)
  utils::write.csv(pkgs, con, row.names = FALSE)
  close(con)
  cat("\nwritten ->", out_path, sprintf("(%d packages)\n", nrow(pkgs)))
}

if (format %in% c("json", "both")) {
  json_path <- sub("[.]csv$", ".json", out_path)
  if (identical(json_path, out_path)) json_path <- paste0(out_path, ".json")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      c(facts, list(packages = pkgs)),
      json_path, auto_unbox = TRUE, pretty = TRUE
    )
    cat("written ->", json_path, "\n")
  } else {
    # Not fatal: the CSV already carries the package table and the scalar facts as
    # its preamble, so a missing jsonlite costs machine-readability, not the record.
    cat("note: jsonlite absent, JSON record skipped (CSV written)\n")
  }
}
