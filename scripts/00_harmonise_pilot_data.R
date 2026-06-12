#!/usr/bin/env Rscript
# scripts/00_harmonise_pilot_data.R
# ---------------------------------------------------------------------------
# Harmonise the raw per-language pilot CSVs into a single CLAPS-standard file.
#
# PREREGISTERED DECISIONS:
#
#   Semantics predictor  : affectedness_scores_all  (verb-level, not agent-specific)
#                          Matched to OSF study (V6_English_Passives.R: Data$Semantics=Data$FACA,
#                          where FACA = PCA1 of whole-event affectedness questions, Semantics.csv).
#                          affectedness_scores_agent varies per trial (agent-specific) → not used.
#
#   Hypothesis 1 (H1)    : S_TypeActive:Semantics < 0  — ONE-TAILED, negative direction.
#
#   Hypothesis 2 (H2)    : S_TypePseudo_Passive:Semantics — TWO-TAILED.
#                          Secondary prediction; Turkish pilot showed opposite direction.
#                          Implemented in R/05_hypothesis_tests.R as H2a (positive) + H2b (negative).
#
#   Languages with        : English, Turkish
#   Pseudo_Passive        :
#
#   Norwegian             : Synthetic_Passive rows excluded before analysis.
#                          English and Turkish have Pseudo_Passive; Norwegian does not.
#
#   Verb labels           : Standardised cross-language identifiers in the form
#                          'Language_VERB'  (e.g. "Turkish_JUMP").
#                          Verb in English for all languages — confirmed from pilot files.
#                          Verb_ID is used as the random-effects 'Verb' grouping factor
#                          so that Norwegian_ADMIRE and English_ADMIRE are distinct in the
#                          multi-language models (R/04_model_formulas.R L5_cross_maximal).
#
#   Affectedness scale    : English/Turkish: raw ~1–7.
#                          Norwegian: raw 0–100 (converted from continuous slider).
#                          Both are z-scored per language by scale_semantics() in
#                          R/02_preprocess_factors.R, so raw-scale differences are harmless.
#                          No manual rescaling applied here.
#
#   Balinese pilot        : NOT YET AVAILABLE in data/pilot/.
#                          Must be fetched from Google Drive:
#                          https://drive.google.com/drive/folders/1dkL60DCPIws-bvZksZ20fov5Xc2ba1OB
#                          (Data Resources and Analysis > Pilot).
#                          The old OSF Balinese_Passives.csv (10-point scale, no Pseudo_Passive)
#                          must NOT be used — it is from the previous study with a different design.
#                          When available, add an entry to PILOT_FILES below.
#
# Output: data/pilot/claps_pilot_harmonised.csv
#         CLAPS-standard columns:
#           Participant, Language, Verb, Verb_ID, Item,
#           S_Type, Semantics, Response, Is_Pilot
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(assertr)
  library(here)
})

source("R/00_check_data_consistency.R")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PILOT_FILES <- list(
  English = list(
    path               = here::here("data/pilot/dataEng_Final.csv"),
    has_pseudo_passive = TRUE,
    # Norwegian has an extra 'sentence' col before 'sentence_type', shifting column positions.
    # English/Turkish layout: id, trial, verb, agent_theme, verb_agent,
    #                         sentence_type, verb_type, notes, rating,
    #                         affectedness_scores_all, affectedness_scores_agent
    has_sentence_col   = FALSE
  ),
  Turkish = list(
    path               = here::here("data/pilot/dataTur_Final.csv"),
    has_pseudo_passive = TRUE,
    has_sentence_col   = FALSE
  ),
  Norwegian = list(
    path               = here::here("data/pilot/dataNor_Final.csv"),
    has_pseudo_passive = FALSE,
    has_sentence_col   = TRUE
    # Norwegian layout adds 'sentence' at col 5:
    # id, trial, verb, agent_theme, sentence, sentence_type, verb_type,
    # rating_100, verb_agent, affectedness_scores_all, affectedness_scores_agent, rating
  )
  # Balinese = list(
  #   path               = here::here("data/pilot/dataBal_Final.csv"),
  #   has_pseudo_passive = TRUE,   # confirm once file is available
  #   has_sentence_col   = FALSE   # confirm once file is available
  # )
)

# Canonical S_Type value mapping (raw pilot → CLAPS standard)
S_TYPE_MAP <- c(
  "Active"                    = "Active",
  "Passive"                   = "Passive",
  "PseudoPassive_Topicalization" = "Pseudo_Passive",
  "passive - analytical"      = "Passive",
  "passive - syntetic"        = "Synthetic_Passive",  # typo in source retained for tracing
  "active"                    = "Active"
)

OUT_PATH <- here::here("data/pilot/claps_pilot_harmonised.csv")

# ---------------------------------------------------------------------------
# Helper: read and harmonise one language file
# ---------------------------------------------------------------------------

harmonise_one <- function(language, cfg) {
  stopifnot(file.exists(cfg$path))
  message("[harmonise] Reading: ", cfg$path)

  # read.csv handles quoted fields with embedded commas (notes column) correctly
  raw <- read.csv(cfg$path, header = TRUE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8", comment.char = "")

  # ------------------------------------------------------------------
  # Apply robust UTF-8 fix (encoding safety: coerce all character cols to valid UTF-8)
  # ------------------------------------------------------------------
  robust_utf8 <- function(x) {
    if (is.character(x)) {
      x <- iconv(x, to = "UTF-8", sub = "byte")
      Encoding(x) <- "UTF-8"
    }
    x
  }
  raw <- as.data.frame(lapply(raw, robust_utf8))

  # ------------------------------------------------------------------
  # Rename columns to CLAPS schema
  # ------------------------------------------------------------------
  # Norwegian has 'rating' at col 12 and 'rating_100' at col 8.
  # English/Turkish have 'rating' at col 9.
  # Both schemas have 'affectedness_scores_all' and 'affectedness_scores_agent'.

  # Validate presence of required raw columns
  required_raw <- c("id", "verb", "agent_theme", "sentence_type",
                     "rating", "affectedness_scores_all", "affectedness_scores_agent")
  missing_cols <- setdiff(required_raw, names(raw))
  if (length(missing_cols) > 0) {
    stop("[harmonise] ", language, ": missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  out <- raw |>
    dplyr::rename(
      Participant              = id,
      Response                 = rating,
      Semantics                = affectedness_scores_all
      # affectedness_scores_agent retained as a separate column for reference
    ) |>
    dplyr::mutate(
      Language = language,
      Is_Pilot = TRUE,

      # Standardised S_Type: remap raw values to CLAPS vocabulary
      S_Type = dplyr::recode(sentence_type, !!!S_TYPE_MAP, .default = NA_character_),

      # Cross-language verb identifier: Language_VERB (verb in English, uppercase)
      # Used as the 'Verb' grouping factor in random effects so that
      # Norwegian_ADMIRE and English_ADMIRE are separate random-effect units
      # in multi-language models (R/04_model_formulas.R L5_cross_maximal).
      Verb_ID = paste0(language, "_", toupper(trimws(verb))),

      # Preserve original lower-case verb name for readability
      Verb = trimws(verb),

      # Item: verb + agent_theme combination (identifies the specific trial stimulus)
      Item = paste(trimws(verb), trimws(agent_theme), sep = "_")
    )

  # ------------------------------------------------------------------
  # Flag unmapped S_Type values
  # ------------------------------------------------------------------
  unmapped <- unique(raw$sentence_type[is.na(out$S_Type)])
  if (length(unmapped) > 0) {
    warning("[harmonise] ", language, ": unmapped sentence_type values: ",
            paste(unmapped, collapse = ", "))
  }

  # ------------------------------------------------------------------
  # Norwegian: exclude Synthetic_Passive (preregistered design decision)
  # ------------------------------------------------------------------
  if (!cfg$has_pseudo_passive) {
    n_before <- nrow(out)
    out <- dplyr::filter(out, S_Type != "Synthetic_Passive")
    message("[harmonise] Norwegian: excluded ", n_before - nrow(out),
            " Synthetic_Passive rows (", n_before, " → ", nrow(out), " rows)")
  }

  # ------------------------------------------------------------------
  # Select and order final columns
  # ------------------------------------------------------------------
  out <- out |>
    dplyr::select(
      Participant, Language, Verb, Verb_ID, Item,
      S_Type, Semantics, Response, Is_Pilot,
      # Retain for diagnostics; not used in models
      trial, verb_type, affectedness_scores_agent
    )

  out
}

# ---------------------------------------------------------------------------
# Harmonise all languages
# ---------------------------------------------------------------------------

parts <- mapply(harmonise_one,
                language = names(PILOT_FILES),
                cfg      = PILOT_FILES,
                SIMPLIFY = FALSE)

harmonised <- dplyr::bind_rows(parts)

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

message("[validate] Running post-harmonisation checks...")

harmonised |>
  assertr::verify(
    !is.na(Participant), error_fun = assertr::just_warn
  ) |>
  assertr::verify(
    !is.na(S_Type),
    description = "No unmapped S_Type values",
    error_fun   = assertr::just_warn
  ) |>
  assertr::verify(
    Response %in% 1:7,
    description = "Response in 1–7",
    error_fun   = assertr::just_warn
  ) |>
  assertr::verify(
    is.numeric(Semantics) & !is.na(Semantics),
    description = "Semantics is numeric, no NAs",
    error_fun   = assertr::just_warn
  ) |>
  assertr::verify(
    grepl("^[A-Za-z]+_[A-Z]+", Verb_ID),
    description = "Verb_ID has Language_VERB format",
    error_fun   = assertr::just_warn
  )

# Summary
message("\n[harmonise] Combined dataset: ", nrow(harmonised), " rows")
message("[harmonise] Languages       : ", paste(sort(unique(harmonised$Language)), collapse = ", "))
message("[harmonise] S_Type levels   : ", paste(sort(unique(harmonised$S_Type)), collapse = ", "))
message("[harmonise] Participants    : ", dplyr::n_distinct(harmonised$Participant, harmonised$Language))

cat("\n--- Rows per language × S_Type ---\n")
print(as.data.frame(table(harmonised$Language, harmonised$S_Type)))

cat("\n--- Verbs per language ---\n")
print(harmonised |>
  dplyr::group_by(Language) |>
  dplyr::summarise(n_verbs = dplyr::n_distinct(Verb), .groups = "drop"))

cat("\n--- Verb_ID sample (first 5 per language) ---\n")
print(harmonised |>
  dplyr::group_by(Language) |>
  dplyr::slice(1) |>
  dplyr::select(Language, Verb, Verb_ID))

cat("\n--- Semantics scale per language (raw, before z-scoring) ---\n")
print(harmonised |>
  dplyr::group_by(Language) |>
  dplyr::summarise(
    min    = round(min(Semantics, na.rm = TRUE), 3),
    max    = round(max(Semantics, na.rm = TRUE), 3),
    mean   = round(mean(Semantics, na.rm = TRUE), 3),
    sd     = round(sd(Semantics, na.rm = TRUE), 3),
    .groups = "drop"
  ))

# ---------------------------------------------------------------------------
# Write atomic output
# ---------------------------------------------------------------------------

tmp_path <- paste0(OUT_PATH, ".tmp")
readr::write_csv(harmonised, tmp_path)
file.rename(tmp_path, OUT_PATH)

message("\n[harmonise] Written: ", OUT_PATH)

# ---------------------------------------------------------------------------
# Cross-language consistency checks
# ---------------------------------------------------------------------------
message("[harmonise] Running cross-language consistency checks...")
check_crosslanguage_consistency(harmonised)

message("[harmonise] REMINDER: Balinese pilot data not yet included.")
message("            Fetch from Google Drive (Data Resources and Analysis > Pilot)")
message("            then uncomment the Balinese entry in PILOT_FILES above.")
