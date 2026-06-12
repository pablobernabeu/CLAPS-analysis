# R/01_read_validate_data.R
# Read and validate raw CLAPS data files.
# Validates column presence, factor levels, response range, and participant IDs.
# All checks fail loudly — no silent coercion.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(assertr)
  library(here)
})

REQUIRED_COLUMNS <- c(
  "Participant", "Language", "Verb", "Verb_ID", "Item",
  "S_Type", "Semantics", "Response"
)

VALID_S_TYPES  <- c("Active", "Passive", "Pseudo_Passive",
                    "Synthetic_Passive")  # Synthetic_Passive excluded later by
                                          # exclude_norwegian_synthetic_passive()
RESPONSE_RANGE <- c(1L, 7L)   # 7-point acceptability scale

#' Read a single raw data file (CSV or TSV) and return a tibble.
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
#' Stops on any violation.
validate_raw_data <- function(df, source_label = "data") {
  # 1. Required columns
  missing_cols <- setdiff(REQUIRED_COLUMNS, names(df))
  if (length(missing_cols) > 0) {
    stop("[validate] Missing required columns in ", source_label, ": ",
         paste(missing_cols, collapse = ", "))
  }

  # 2. Response must be integer in [1, 7]
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

  # 3. S_Type levels
  invalid_s_type <- setdiff(unique(df$S_Type), VALID_S_TYPES)
  if (length(invalid_s_type) > 0) {
    stop("[validate] Unexpected S_Type values in ", source_label, ": ",
         paste(invalid_s_type, collapse = ", "))
  }

  # 4. No NA in key columns
  key_cols <- c("Participant", "S_Type", "Semantics", "Response", "Verb")
  for (col in key_cols) {
    n_na <- sum(is.na(df[[col]]))
    if (n_na > 0) {
      stop("[validate] ", n_na, " NA(s) in column '", col, "' in ", source_label)
    }
  }

  # 5. Semantics must be numeric or coercible
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
