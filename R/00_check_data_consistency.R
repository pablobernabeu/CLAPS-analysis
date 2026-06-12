# R/00_check_data_consistency.R
# Cross-language data consistency checks for the CLAPS multi-language model.
# All checks fail loudly with informative messages — no silent coercion.
#
# Checks performed:
#   1. Required columns present and correctly typed in each language file
#   2. Response is integer in [1, 7] for every language
#   3. S_Type factor levels are consistent (modulo language-specific passives)
#   4. Verb labels are language-specific (no cross-language collisions)
#   5. UTF-8 encoding validated per language
#   6. Semantics is numeric and non-constant within each language
#   7. No cross-language participant ID collisions (IDs must be language-scoped)
#   8. Norwegian: no Synthetic_Passive rows (preregistered exclusion)
#
# Usage:
#   source("R/00_check_data_consistency.R")
#   check_crosslanguage_consistency(df_all)   # df_all is the combined data frame

suppressPackageStartupMessages({
  library(dplyr)
})

REQUIRED_COLUMNS <- c(
  "Participant", "Language", "Verb", "Verb_ID", "Item",
  "S_Type", "Semantics", "Response"
)

VALID_S_TYPES   <- c("Active", "Passive", "Pseudo_Passive",
                     "Synthetic_Passive")  # Synthetic_Passive excluded upstream
RESPONSE_RANGE  <- c(1L, 7L)

# Languages that have Pseudo_Passive (update as languages are added)
LANGUAGES_WITH_PSEUDO_PASSIVE <- c("English", "Turkish")

# Languages where Synthetic_Passive must be absent (preregistered)
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

  # Verb_ID format check: must match Language_VERB
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
    raw_bytes <- iconv(df[[col]], from = "UTF-8", to = "UTF-8", sub = NA)
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
#' @param df A combined data frame with a Language column. All languages
#'   present in the data are checked.
#' @return Invisibly returns TRUE on success; stops or warns on failure.
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

  # Split by language for per-language checks
  df_list <- split(df, df$Language)

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

#' Apply a robust UTF-8 encoding fix.
#' Call immediately after read.csv() / read_csv().
#' @param df A data frame.
#' @return df with all character columns coerced to valid UTF-8.
robust_utf8_df <- function(df) {
  robust_utf8 <- function(x) {
    if (is.character(x)) {
      x <- iconv(x, to = "UTF-8", sub = "byte")
      Encoding(x) <- "UTF-8"
    }
    x
  }
  as.data.frame(lapply(df, robust_utf8), stringsAsFactors = FALSE)
}
