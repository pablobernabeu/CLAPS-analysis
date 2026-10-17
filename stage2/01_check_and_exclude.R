#!/usr/bin/env Rscript
# stage2/01_check_and_exclude.R
# ---------------------------------------------------------------------------
# Stage 2, step 1 of 6: check the teams' data and apply the registered exclusions.
#
# Reads every *acceptability*.csv and *ratings*.csv file in the raw directory (the
# canonical format described in stage2/README.md), checks them, applies the three
# exclusion criteria of the Stage 1 report to both tasks, and keeps the registered
# number of analysed participants per variety (80 and 20). A variety that has not
# reached both numbers, or that has fewer than the minimum number of verbs, has not
# completed data collection and goes no further.
#
# Writes, under the processed and results directories of config/analysis_config.yaml:
#   processed/acceptability_analysed.csv   trials of the analysed participants
#   processed/ratings_analysed.csv         trials of the analysed raters
#   results/exclusions_acceptability.csv   every participant, with each criterion
#   results/exclusions_ratings.csv
#   results/samples.csv                    per variety: recruited, excluded, analysed
#
# Usage (from design_analysis/):
#   Rscript stage2/01_check_and_exclude.R [--dry-run]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/02_preprocess_factors.R")
source("R/03_define_priors.R")
source("R/14_stage2.R")

args <- commandArgs(trailingOnly = TRUE)
cfg <- stage2_settings(stage2_config(), dry_run = "--dry-run" %in% args)
P <- cfg$paths
dir.create(P$processed, recursive = TRUE, showWarnings = FALSE)
dir.create(P$results, recursive = TRUE, showWarnings = FALSE)
lang <- read_language_table(cfg$languages_table)

read_task <- function(pattern) {
  files <- list.files(P$raw, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("[stage2] no files matching ", pattern, " in ", P$raw)
  dplyr::bind_rows(lapply(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }))
}
acc <- validate_task_data(read_task("acceptability.*[.]csv$"), "acceptability")
rat <- validate_task_data(read_task("ratings.*[.]csv$"), "ratings")
acc$Attention_Target <- as.integer(acc$Attention_Target)
rat$Attention_Target <- as.integer(rat$Attention_Target)

unknown <- setdiff(unique(c(acc$Variety, rat$Variety)), lang$variety)
if (length(unknown)) {
  stop("[stage2] varieties not in ", cfg$languages_table, ": ", paste(unknown, collapse = ", "),
       ". Add them, with their Glottolog language and family, before running.")
}

# The order in which participants completed the task, needed only if a team has more
# analysed participants than the registered number.
order_file <- file.path(P$raw, "sessions.csv")
sessions <- if (file.exists(order_file)) readr::read_csv(order_file, show_col_types = FALSE) else NULL

X <- cfg$exclusions
excl_acc <- participant_exclusions(acc, "acceptability", X$attention_trials$acceptability,
                                   X$attention_min_correct, X$straightline_share) |>
  cap_sample(cfg$samples$acceptability, sessions)
excl_rat <- participant_exclusions(rat, "ratings", X$attention_trials$ratings,
                                   X$attention_min_correct, X$straightline_share) |>
  cap_sample(cfg$samples$ratings, sessions)

verbs <- acc |>
  dplyr::filter(Trial_Type == "test") |>
  dplyr::group_by(Variety) |>
  dplyr::summarise(n_verbs = dplyr::n_distinct(Verb), .groups = "drop")
summ <- function(e, what) {
  e |>
    dplyr::group_by(Variety) |>
    dplyr::summarise(recruited = dplyr::n(), incomplete = sum(incomplete),
                     failed_attention = sum(failed_attention), straightlined = sum(straightlined),
                     excluded = sum(excluded), analysed = sum(analysed), .groups = "drop") |>
    dplyr::rename_with(~ paste0(what, "_", .x), -Variety)
}
samples <- dplyr::full_join(summ(excl_acc, "acc"), summ(excl_rat, "rat"), by = "Variety") |>
  dplyr::left_join(verbs, by = "Variety") |>
  dplyr::mutate(
    enough_verbs = n_verbs >= cfg$samples$min_verbs & n_verbs <= cfg$samples$max_verbs,
    complete = Variety %in% attr(excl_acc, "complete") & Variety %in% attr(excl_rat, "complete") &
      enough_verbs
  )

keep <- samples$Variety[samples$complete]
ids_acc <- excl_acc |> dplyr::filter(analysed, Variety %in% keep) |> dplyr::select(Variety, Participant)
ids_rat <- excl_rat |> dplyr::filter(analysed, Variety %in% keep) |> dplyr::select(Variety, Participant)
readr::write_csv(dplyr::semi_join(acc, ids_acc, by = c("Variety", "Participant")),
                 file.path(P$processed, "acceptability_analysed.csv"))
readr::write_csv(dplyr::semi_join(rat, ids_rat, by = c("Variety", "Participant")),
                 file.path(P$processed, "ratings_analysed.csv"))
readr::write_csv(excl_acc, file.path(P$results, "exclusions_acceptability.csv"))
readr::write_csv(excl_rat, file.path(P$results, "exclusions_ratings.csv"))
readr::write_csv(samples, file.path(P$results, "samples.csv"))

message(sprintf("[stage2:01] %d of %d varieties complete: %s", length(keep), nrow(samples),
                paste(keep, collapse = ", ")))
if (any(!samples$complete)) {
  message("[stage2:01] not complete, so not analysed: ",
          paste(samples$Variety[!samples$complete], collapse = ", "))
}
