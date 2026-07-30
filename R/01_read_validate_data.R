# R/01_read_validate_data.R
#
# Purpose
#   The gate every raw CLAPS data file passes through before any analysis sees
#   it. Checks the schema, the response range and the sentence-type levels, and
#   applies the one preregistered row exclusion.
#
# Design principle: fail loudly, never coerce
#   Each check below stops with a message naming the file and the offending
#   values. Nothing is dropped, repaired or silently converted. A response of 8
#   on a 7-point scale is a data-entry error whose cause must be found in the
#   source file, and coercing or discarding it here would hide a problem that
#   affects how the rest of that participant's data should be read.
#
# Entry points
#   read_raw_data(path)            One file, no validation.
#   validate_raw_data(df, label)   Schema checks; stops on violation.
#   read_all_raw_data(glob)        Read and validate everything matching a glob.
#   exclude_norwegian_synthetic_passive(df)  The preregistered exclusion.
#   split_pilot_confirmatory(df)   Partition by the Is_Pilot flag.
#
# Related
#   R/00_check_data_consistency.R  Cross-file consistency, run after this.
#   R/02_preprocess_factors.R      Contrasts and scaling, run after this.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(assertr)
  library(here)
})

# The schema. Verb and Verb_ID are both required and are not redundant: Verb is
# the surface form, which differs across languages, while Verb_ID identifies the
# concept that translation-equivalent verbs share and is what the by-verb random
# effects group on.
REQUIRED_COLUMNS <- c(
  "Participant", "Language", "Verb", "Verb_ID", "Item",
  "S_Type", "Semantics", "Response"
)

# Every sentence type that may legitimately appear in a raw file. Synthetic_Passive
# is accepted at this stage and removed further down the pipeline by
# exclude_norwegian_synthetic_passive(); validation and exclusion are kept apart so
# that an unexpected value is still an error rather than being quietly filtered.
VALID_S_TYPES  <- c("Active", "Passive", "Pseudo_Passive",
                    "Synthetic_Passive")

# Acceptability was collected on a 7-point Likert scale, which is why the ordinal
# model in R/04_model_formulas.R estimates six thresholds (categories minus one).
RESPONSE_RANGE <- c(1L, 7L)

#' Read a single raw data file (CSV or TSV) and return a tibble.
#'
#' @param path Path to a .csv or .tsv file.
#' @return A tibble, unvalidated. Validation is a separate step so that a file can
#'   be inspected before it is judged.
#' @details An unrecognised extension is an error rather than a guess, because a
#'   delimiter inferred wrongly produces a single-column data frame that would
#'   then fail the schema check with a confusing message about missing columns.
read_raw_data <- function(path) {
  stopifnot(file.exists(path))
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
    csv = readr::read_csv(path, show_col_types = FALSE),
    tsv = readr::read_tsv(path, show_col_types = FALSE),
    stop("Unsupported file extension: ", ext)
  )
  df
}

#' Validate that a raw data frame conforms to the CLAPS schema.
#'
#' @param df A data frame from read_raw_data().
#' @param source_label Label used in error messages, normally the file path. Worth
#'   passing when validating several files in a loop, since otherwise a failure
#'   does not say which file caused it.
#' @return Invisibly, `df` unchanged. This function is called for its side effect
#'   of stopping on bad data; it never returns a modified data frame.
#' @details The five checks run in a deliberate order, from the most structural to
#'   the most specific, so the first error a user sees is the most informative one.
#'   Missing columns are reported before anything else, because every check after
#'   that would otherwise fail with an unhelpful "object not found".
validate_raw_data <- function(df, source_label = "data") {
  # 1. Required columns. Reported as a set, so one run names every missing column
  #    rather than stopping at the first.
  missing_cols <- setdiff(REQUIRED_COLUMNS, names(df))
  if (length(missing_cols) > 0) {
    stop("[validate] Missing required columns in ", source_label, ": ",
         paste(missing_cols, collapse = ", "))
  }

  # 2. Response must be a whole number within the scale. Both parts matter: the
  #    bounds check catches out-of-scale codes, and the integrality check catches
  #    a mean or interpolated value that has been written into a raw file by
  #    mistake. An ordinal model would accept 4.5 as a factor level and silently
  #    estimate an extra threshold for it.
  df |>
    assertr::assert(
      assertr::within_bounds(RESPONSE_RANGE[1], RESPONSE_RANGE[2]),
      Response,
      error_fun = assertr::error_stop
    ) |>
    assertr::assert(
      function(x) x == as.integer(x),
      Response,
      error_fun = assertr::error_stop
    )

  # 3. S_Type levels. An unexpected value is an error, not a new category: a typo
  #    such as "passive" would otherwise become a fourth sentence type with its
  #    own coefficient.
  invalid_s_type <- setdiff(unique(df$S_Type), VALID_S_TYPES)
  if (length(invalid_s_type) > 0) {
    stop("[validate] Unexpected S_Type values in ", source_label, ": ",
         paste(invalid_s_type, collapse = ", "))
  }

  # 4. No missing values in the columns the model cannot proceed without. These
  #    five are the grouping factors, the focal predictor and the outcome; brms
  #    would drop such rows silently through complete-case deletion, changing the
  #    effective sample size without any record of it. Language and Item are
  #    absent from this list because they are checked in
  #    R/00_check_data_consistency.R against the item inventory.
  key_cols <- c("Participant", "S_Type", "Semantics", "Response", "Verb")
  for (col in key_cols) {
    n_na <- sum(is.na(df[[col]]))
    if (n_na > 0) {
      stop("[validate] ", n_na, " NA(s) in column '", col, "' in ", source_label)
    }
  }

  # 5. Semantics must be numeric, or at least coercible without warning. Note that
  #    the coercion here is a test only: its result is discarded and `df` is
  #    returned unchanged, so a character Semantics column passes validation and
  #    is still character afterwards. Converting it is the caller's job, in
  #    R/02_preprocess_factors.R. A column that cannot be coerced warns, and the
  #    warning is promoted to an error, which is the case worth catching: it means
  #    the column holds something like a decimal comma or a stray label.
  if (!is.numeric(df$Semantics)) {
    tryCatch(
      as.numeric(df$Semantics),
      warning = function(w) stop("[validate] Semantics is not numeric in ", source_label)
    )
  }

  invisible(df)
}

#' Read all raw data files matching a glob pattern and bind them.
read_all_raw_data <- function(glob_pattern = "data/raw/*.csv") {
  paths <- Sys.glob(glob_pattern)
  if (length(paths) == 0) {
    stop("[read] No files matched pattern: ", glob_pattern)
  }
  message("[read] Found ", length(paths), " file(s): ", paste(paths, collapse = ", "))

  df_list <- purrr::map(paths, function(p) {
    df <- read_raw_data(p)
    validate_raw_data(df, source_label = p)
    df
  })

  dplyr::bind_rows(df_list)
}

#' Exclude Norwegian synthetic passive rows.
#' Preregistered decision: Norwegian pilot data include only Active and
#' analytical Passive; synthetic passive rows are excluded before any
#' analysis (including pilot).
#' @param df A data frame with Language and S_Type columns.
#' @return df with Norwegian Synthetic_Passive rows removed; S_Type levels dropped.
exclude_norwegian_synthetic_passive <- function(df) {
  if (!"Language" %in% names(df) || !"S_Type" %in% names(df)) {
    return(df)
  }
  if ("Norwegian" %in% df$Language) {
    n_before <- nrow(df)
    df <- df[!(df$Language == "Norwegian" & df$S_Type == "Synthetic_Passive"), ]
    n_removed <- n_before - nrow(df)
    if (n_removed > 0) {
      message("[exclude] Removed ", n_removed,
              " Norwegian Synthetic_Passive rows (preregistered exclusion).")
    }
    if (is.factor(df$S_Type)) {
      df$S_Type <- droplevels(df$S_Type)
    }
  }
  df
}

#' Separate pilot data from confirmatory data.
#' Pilot data are identified by a flag column 'Is_Pilot' or by explicit participant IDs.
#' Returns a list(pilot = ..., confirmatory = ...).
split_pilot_confirmatory <- function(df, pilot_col = "Is_Pilot") {
  if (pilot_col %in% names(df)) {
    list(
      pilot         = dplyr::filter(df, .data[[pilot_col]] == TRUE),
      confirmatory  = dplyr::filter(df, .data[[pilot_col]] == FALSE)
    )
  } else {
    message("[split] No '", pilot_col, "' column found; treating all data as pilot.")
    list(pilot = df, confirmatory = df[0, ])
  }
}
