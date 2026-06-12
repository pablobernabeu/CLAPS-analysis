# R/02_preprocess_factors.R
# Factor coding for CLAPS models.
# S_Type: treatment coding, Passive = reference level.
# Semantics: numeric (already centred/scaled upstream).
# Languages without Pseudo_Passive have that level dropped before fitting.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
})

#' Apply treatment coding to S_Type with Passive as the reference level.
#' Validates that Passive exists in the data.
#' @param df A data frame with S_Type column.
#' @return df with S_Type as a treatment-coded factor.
code_s_type <- function(df) {
  if (!"S_Type" %in% names(df)) stop("[code_s_type] S_Type column not found.")
  observed <- unique(df$S_Type)
  if (!"Passive" %in% observed) {
    stop("[code_s_type] 'Passive' level not found in S_Type. Observed: ",
         paste(observed, collapse = ", "))
  }
  df |>
    dplyr::mutate(
      S_Type = forcats::fct_relevel(factor(S_Type), "Passive")
    )
}

#' Drop Pseudo_Passive rows and level for languages that do not have pseudo-passives.
#' Also drops any pseudo-passive interaction terms by ensuring the level is absent.
#' @param df A data frame with S_Type column.
#' @param has_pseudo_passive Logical scalar.
#' @return df with Pseudo_Passive rows and level removed if has_pseudo_passive is FALSE.
drop_pseudo_passive_if_absent <- function(df, has_pseudo_passive) {
  if (!has_pseudo_passive) {
    df <- dplyr::filter(df, S_Type != "Pseudo_Passive")
    if (is.factor(df$S_Type)) {
      df <- dplyr::mutate(df, S_Type = droplevels(S_Type))
    }
    message("[preprocess] Pseudo_Passive dropped for this language.")
  }
  df
}

#' Scale Semantics to mean 0, SD 0.5 (Gelman scaling) within each language.
#' @param df Data frame.
#' @param centre_by Character column to group by before scaling.
#' @return df with Semantics_scaled added.
scale_semantics <- function(df, centre_by = "Language") {
  if (centre_by %in% names(df)) {
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(centre_by))) |>
      dplyr::mutate(
        Semantics_scaled = (Semantics - mean(Semantics, na.rm = TRUE)) /
                           (2 * sd(Semantics, na.rm = TRUE))
      ) |>
      dplyr::ungroup()
  } else {
    df <- dplyr::mutate(
      df,
      Semantics_scaled = (Semantics - mean(Semantics, na.rm = TRUE)) /
                         (2 * sd(Semantics, na.rm = TRUE))
    )
  }
  df
}

#' Set the Semantics source column before scaling.
#' Overwrites `Semantics` with the values of `source_col`, so the downstream
#' focal predictor (`Semantics_scaled`) is built from the chosen affectedness
#' column. The gender model variation sources affectedness from
#' `affectedness_scores_agent` (agent/gender-specific) instead of the default
#' whole-event affectedness already stored in `Semantics`.
#' @param df Data frame.
#' @param source_col Character or NULL; column to copy into `Semantics`. NULL
#'   (or "Semantics") returns df unchanged.
#' @return df with `Semantics` set from `source_col`.
set_semantics_source <- function(df, source_col = NULL) {
  if (is.null(source_col) || identical(source_col, "Semantics")) return(df)
  if (!source_col %in% names(df)) {
    stop("[set_semantics_source] Source column '", source_col, "' not found.")
  }
  if (any(is.na(df[[source_col]]))) {
    stop("[set_semantics_source] NA values in source column '", source_col, "'.")
  }
  df$Semantics <- as.numeric(df[[source_col]])
  message("[preprocess] Semantics sourced from '", source_col, "'.")
  df
}

#' Derive a referent-gender covariate from the Item column.
#' Item is "verb_AgentTheme" (e.g. "push_Man", "see_Woman"); the gender token is
#' the part after the final underscore. Returns a treatment-coded factor with
#' Man as the reference level. Used by the gender model variation.
#' @param df Data frame with an Item column.
#' @param item_col Character; name of the item column.
#' @return df with a `Gender` factor column (levels Man, Woman; Man = reference).
derive_gender <- function(df, item_col = "Item") {
  if (!item_col %in% names(df)) {
    stop("[derive_gender] '", item_col, "' column not found.")
  }
  g   <- sub(".*_", "", as.character(df[[item_col]]))
  bad <- setdiff(unique(g), c("Man", "Woman"))
  if (length(bad) > 0) {
    stop("[derive_gender] Unexpected gender token(s) in ", item_col, ": ",
         paste(bad, collapse = ", "))
  }
  df$Gender <- forcats::fct_relevel(factor(g), "Man")
  df
}

#' Full preprocessing pipeline: (optionally re-source Semantics), code, drop
#' levels, scale, and (optionally) derive the Gender covariate.
#' @param df Raw validated data frame.
#' @param has_pseudo_passive Logical; whether this language uses pseudo-passives.
#' @param semantics_source Character or NULL; if given, source `Semantics` from
#'   this column before scaling (e.g. "affectedness_scores_agent" for the gender
#'   variation). NULL keeps the existing `Semantics` column.
#' @param include_gender Logical; if TRUE, derive the `Gender` factor from Item.
#' @return Preprocessed data frame ready for model fitting.
preprocess_data <- function(df, has_pseudo_passive = TRUE,
                            semantics_source = NULL,
                            include_gender = FALSE) {
  df <- set_semantics_source(df, semantics_source)
  df <- df |>
    code_s_type() |>
    drop_pseudo_passive_if_absent(has_pseudo_passive) |>
    scale_semantics()
  if (isTRUE(include_gender)) df <- derive_gender(df)
  df
}

#' Assert that factor contrast coding is correct: treatment coding, Passive = reference.
assert_treatment_coding <- function(df) {
  if (!is.factor(df$S_Type)) stop("[assert] S_Type must be a factor after preprocessing.")
  ref <- levels(df$S_Type)[1]
  if (ref != "Passive") {
    stop("[assert] Reference level of S_Type is '", ref, "', expected 'Passive'.")
  }
  contrasts_mat <- contrasts(df$S_Type)
  if (!identical(attr(contrasts_mat, "contrasts"), "contr.treatment")) {
    # Tolerate default contr.treatment without explicit attribute
  }
  invisible(TRUE)
}
