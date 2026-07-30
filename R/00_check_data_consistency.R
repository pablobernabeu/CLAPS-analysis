# R/00_check_data_consistency.R
#
# Purpose
#   Checks that only make sense once the languages are combined. Whereas
#   R/01_read_validate_data.R validates each file against the schema, this file
#   asks whether the files fit together: whether identifiers collide, whether the
#   sentence types present match what each language is supposed to have, and
#   whether the encoding survived the merge. These are the failures that produce a
#   model which fits without complaint and answers the wrong question.
#
# Checks performed
#   Per language
#     1. Required columns present                     -> error
#     2. Response integer within [1, 7]               -> error
#     3. S_Type levels expected for that language      -> error, except the
#        unexpected-pseudo-passive case, which only messages (see below)
#     4. Character columns are valid UTF-8            -> warning
#     5. Semantics numeric, complete and non-constant -> error
#   Across languages
#     6. Verb_ID unique to one language                -> error
#     7. Participant IDs unique to one language        -> warning
#
# Severity policy
#   Most checks stop, because the condition cannot be true of usable data. Three
#   deliberately do not:
#     - Participant ID collisions warn rather than stop, because a genuinely
#       shared bare ID such as "001" is possible if a site numbered participants
#       from 1 in each language. It is far more often a merge error, so it is
#       reported loudly and left to a human.
#     - Encoding problems warn, because robust_utf8_df() below can repair them and
#       the caller may intend to.
#     - Pseudo-passive rows in a language not listed as having them only message,
#       because the likelier cause is that LANGUAGES_WITH_PSEUDO_PASSIVE has not
#       been updated for a newly added language.
#   Nothing here modifies the data. Repair is always a separate, explicit call.
#
# Usage
#   source("R/00_check_data_consistency.R")
#   check_crosslanguage_consistency(df_all)   # df_all is the combined data frame

suppressPackageStartupMessages({
  library(dplyr)
})

# These four constants intentionally repeat the definitions in
# R/01_read_validate_data.R and config/analysis_config.yaml, so that this file can
# be sourced on its own as a standalone data check. The cost is that a change must
# be made in more than one place: adding a language means updating
# config/analysis_config.yaml AND the two language constants below.
REQUIRED_COLUMNS <- c(
  "Participant", "Language", "Verb", "Verb_ID", "Item",
  "S_Type", "Semantics", "Response"
)

VALID_S_TYPES   <- c("Active", "Passive", "Pseudo_Passive",
                     "Synthetic_Passive")  # Synthetic_Passive excluded upstream
RESPONSE_RANGE  <- c(1L, 7L)

# Languages that have Pseudo_Passive. Mirrors the per-language has_pseudo_passive
# flags in config/analysis_config.yaml, where Norwegian and Balinese are false.
LANGUAGES_WITH_PSEUDO_PASSIVE <- c("English", "Turkish")

# Languages where Synthetic_Passive must be absent by the time data reach a model.
# Norwegian pilot data contain synthetic passives, and their exclusion is
# preregistered; this constant is what makes a failure to apply that exclusion an
# error rather than an oversight.
LANGUAGES_EXCLUDE_SYNTHETIC  <- c("Norwegian")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.check_columns <- function(df, lang) {
  missing <- setdiff(REQUIRED_COLUMNS, names(df))
  if (length(missing) > 0) {
    stop("[consistency:", lang, "] Missing required columns: ",
         paste(missing, collapse = ", "))
  }
  invisible(NULL)
}

.check_response_range <- function(df, lang) {
  bad <- df$Response[!is.na(df$Response) &
                     (df$Response < RESPONSE_RANGE[1] |
                      df$Response > RESPONSE_RANGE[2])]
  if (length(bad) > 0) {
    stop("[consistency:", lang, "] Response out of [",
         RESPONSE_RANGE[1], ",", RESPONSE_RANGE[2], "]: ",
         paste(head(bad, 5), collapse = ", "))
  }
  if (any(df$Response != as.integer(df$Response), na.rm = TRUE)) {
    stop("[consistency:", lang, "] Response is not integer-valued.")
  }
  invisible(NULL)
}

.check_s_type_levels <- function(df, lang) {
  observed <- unique(as.character(df$S_Type))
  # The explicit "Synthetic_Passive" is redundant, as VALID_S_TYPES already
  # contains it; harmless, and kept so the intent stays legible if that constant
  # is ever narrowed.
  invalid   <- setdiff(observed, c(VALID_S_TYPES, "Synthetic_Passive"))
  if (length(invalid) > 0) {
    stop("[consistency:", lang, "] Unexpected S_Type values: ",
         paste(invalid, collapse = ", "))
  }
  # Norwegian must not have Synthetic_Passive (preregistered exclusion)
  if (lang %in% LANGUAGES_EXCLUDE_SYNTHETIC &&
      "Synthetic_Passive" %in% observed) {
    stop("[consistency:", lang, "] Synthetic_Passive rows found; ",
         "these must be excluded upstream via exclude_norwegian_synthetic_passive().")
  }
  # Languages without pseudo-passive should not have Pseudo_Passive rows
  if (!lang %in% LANGUAGES_WITH_PSEUDO_PASSIVE &&
      "Pseudo_Passive" %in% observed) {
    message("[consistency:", lang, "] WARNING: Pseudo_Passive rows present ",
            "but this language is not listed in LANGUAGES_WITH_PSEUDO_PASSIVE. ",
            "Update that constant or check the data.")
  }
  invisible(NULL)
}

.check_verb_labels <- function(df_list) {
  # Verb_ID must be globally unique per language: Language_VERB format.
  # Each Verb_ID should appear in exactly one language.
  # Bare Verb strings WILL repeat across languages (same English concept, different words)
  # — that is expected and correct. Only Verb_ID must be non-colliding.
  all_verb_id_lang <- lapply(names(df_list), function(lang) {
    data.frame(
      verb_id  = unique(as.character(df_list[[lang]]$Verb_ID)),
      language = lang,
      stringsAsFactors = FALSE
    )
  })
  vid_tbl <- do.call(rbind, all_verb_id_lang)

  # Verb_ID format check: must match Language_VERB. The pattern is unanchored at
  # the end, so a trailing suffix is tolerated, and requires the verb part to be
  # upper case, which is the repository's convention for distinguishing a Verb_ID
  # from the surface Verb form.
  bad_format <- vid_tbl$verb_id[!grepl("^[A-Za-z]+_[A-Z]+", vid_tbl$verb_id)]
  if (length(bad_format) > 0) {
    stop("[consistency] Verb_ID values do not match 'Language_VERB' format: ",
         paste(head(bad_format, 5), collapse = ", "))
  }

  # Verb_ID must not appear in more than one language
  duplicates <- vid_tbl |>
    dplyr::group_by(verb_id) |>
    dplyr::filter(dplyr::n_distinct(language) > 1) |>
    dplyr::ungroup()

  if (nrow(duplicates) > 0) {
    dup_ids <- unique(duplicates$verb_id)
    stop("[consistency] Verb_ID values shared across languages (should never happen): ",
         paste(head(dup_ids, 10), collapse = ", "))
  }
  invisible(NULL)
}

.check_participant_ids <- function(df_list) {
  # Participant IDs must be language-scoped; a bare numeric ID appearing
  # in multiple languages is almost certainly a data-merge error.
  all_ppt_lang <- lapply(names(df_list), function(lang) {
    data.frame(
      participant = unique(as.character(df_list[[lang]]$Participant)),
      language    = lang,
      stringsAsFactors = FALSE
    )
  })
  ppt_tbl <- do.call(rbind, all_ppt_lang)

  duplicates <- ppt_tbl |>
    dplyr::group_by(participant) |>
    dplyr::filter(dplyr::n_distinct(language) > 1) |>
    dplyr::ungroup()

  if (nrow(duplicates) > 0) {
    warning("[consistency] ", dplyr::n_distinct(duplicates$participant),
            " participant ID(s) appear in multiple languages. ",
            "Verify that IDs are language-scoped ",
            "(e.g., 'EN_001' not just '001').")
  }
  invisible(NULL)
}

.check_utf8 <- function(df, lang) {
  char_cols <- names(df)[sapply(df, is.character)]
  for (col in char_cols) {
    # Round-tripping UTF-8 to UTF-8 looks pointless but is the standard idiom for
    # detecting invalid byte sequences: with sub = NA, iconv() returns NA for any
    # value it cannot interpret as UTF-8 and leaves valid values untouched.
    raw_bytes <- iconv(df[[col]], from = "UTF-8", to = "UTF-8", sub = NA)
    # Compare against the original NAs so that genuinely missing values are not
    # counted as encoding failures.
    n_invalid <- sum(is.na(raw_bytes) & !is.na(df[[col]]))
    if (n_invalid > 0) {
      warning("[consistency:", lang, "] Column '", col, "' contains ",
              n_invalid, " non-UTF-8 value(s). ",
              "Apply robust_utf8() to fix encoding.")
    }
  }
  invisible(NULL)
}

.check_semantics <- function(df, lang) {
  if (!is.numeric(df$Semantics)) {
    stop("[consistency:", lang, "] Semantics column is not numeric.")
  }
  if (any(is.na(df$Semantics))) {
    stop("[consistency:", lang, "] NA(s) in Semantics column.")
  }
  # Zero variance is fatal, not merely odd. scale_semantics() divides by twice the
  # standard deviation, so a constant Semantics column would yield Inf or NaN for
  # every row of the focal predictor, and the model would fail later with an error
  # that says nothing about its cause.
  if (stats::sd(df$Semantics, na.rm = TRUE) == 0) {
    stop("[consistency:", lang, "] Semantics is constant (zero variance); ",
         "check data preparation.")
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Run all cross-language consistency checks on a combined data frame.
#'
#' @param df A combined data frame with a Language column. Every language present
#'   in the data is checked; the language constants above determine what is
#'   *expected* of each, not which are examined.
#' @return Invisibly TRUE if no check raised an error. Note that this is not a
#'   clean bill of health: the warning-level checks (encoding, participant-ID
#'   collisions) can have fired and the function still returns TRUE, so the
#'   console output must be read, not just the return value.
#'
#' @examples
#' \dontrun{
#'   df_all <- dplyr::bind_rows(df_english, df_turkish, df_norwegian)
#'   check_crosslanguage_consistency(df_all)
#' }
check_crosslanguage_consistency <- function(df) {
  if (!"Language" %in% names(df)) {
    stop("[consistency] 'Language' column not found in combined data frame.")
  }
  languages <- unique(as.character(df$Language))
  message("[consistency] Checking ", length(languages), " language(s): ",
          paste(languages, collapse = ", "))

  # split() keys the list by language, which the two cross-language checks below
  # rely on: they read names(df_list) to label each language's identifiers.
  df_list <- split(df, df$Language)

  # Per-language checks first. Each stops on the language that fails, so the error
  # message names the language rather than reporting an aggregate failure.
  for (lang in languages) {
    d <- df_list[[lang]]
    .check_columns(d, lang)
    .check_response_range(d, lang)
    .check_s_type_levels(d, lang)
    .check_utf8(d, lang)
    .check_semantics(d, lang)
  }

  # Cross-language checks
  .check_verb_labels(df_list)
  .check_participant_ids(df_list)

  message("[consistency] All checks passed for: ", paste(languages, collapse = ", "))
  invisible(TRUE)
}

#' Apply a robust UTF-8 encoding fix to every character column.
#'
#' @param df A data frame.
#' @return A base data frame with character columns coerced to valid UTF-8.
#' @details Call immediately after reading a file, before any other processing.
#'   Relevant to CLAPS because the stimuli include Turkish and Norwegian
#'   orthography, and a file written on one platform and read on another can arrive
#'   with mis-decoded bytes in exactly the columns used as grouping factors, where
#'   two spellings of one verb would become two verbs.
#'
#'   `sub = "byte"` replaces an undecodable byte with a visible \\xNN escape rather
#'   than dropping it or returning NA. That is deliberate: the mangled value
#'   survives into the data where it can be seen and traced to its source file,
#'   whereas silent deletion would leave two subtly different labels that still
#'   look plausible.
#'
#'   Two consequences worth knowing. The return value is a base data.frame, so a
#'   tibble passed in comes back as a data.frame; and as.data.frame() applies R's
#'   name repair, so a column name that is not syntactically valid will be altered.
#'   Neither affects the CLAPS schema, whose column names are all plain.
robust_utf8_df <- function(df) {
  robust_utf8 <- function(x) {
    if (is.character(x)) {
      x <- iconv(x, to = "UTF-8", sub = "byte")
      # iconv() does not always tag the result, and an untagged string is read
      # using the session's locale, which is what reintroduces the problem on a
      # non-UTF-8 machine such as a Windows workstation.
      Encoding(x) <- "UTF-8"
    }
    x
  }
  as.data.frame(lapply(df, robust_utf8), stringsAsFactors = FALSE)
}
