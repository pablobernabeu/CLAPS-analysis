# R/02_preprocess_factors.R
#
# Purpose
#   Turn validated raw data into the exact variable coding the models assume.
#   Everything here is a decision about what the fitted coefficients will *mean*,
#   which is why it lives in one file rather than being scattered through the
#   fitting scripts.
#
# The decisions, and why they matter
#   S_Type is treatment-coded with Passive as the reference level. This is not
#   cosmetic. Under treatment coding the interaction coefficient
#   S_TypeActive:Semantics_scaled is the *difference* between the affectedness
#   slope in actives and that in passives, which is exactly what H1b predicts to
#   be negative. Change the reference level and the same coefficient answers a
#   different question, so the hypothesis tests in R/05_hypothesis_tests.R, the
#   priors in R/03_define_priors.R and the simulator in R/06_simulate_design.R all
#   depend on this choice; each of them names coefficients built on this
#   assumption.
#
#   Semantics is rescaled to Semantics_scaled by centring and dividing by two
#   standard deviations (Gelman, 2008, doi:10.1002/sim.3107). On that scale a unit
#   change corresponds to moving across most of the observed range of
#   affectedness, which puts a continuous predictor on a footing roughly
#   comparable with a binary one and makes the prior scales in
#   R/03_define_priors.R interpretable as effects of a substantial change in
#   affectedness rather than of an arbitrary raw unit.
#
#   Languages without pseudo-passives have that level dropped, rather than
#   retained and estimated from no data.
#
# Entry point
#   preprocess_data() runs the whole sequence in the required order. The
#   individual steps are exported mainly so the tests can exercise them.
#
# Order dependence
#   The steps are not commutative. The Semantics source must be chosen before
#   scaling, since scaling reads whichever column is then in place; and rows must
#   be dropped before scaling, or the mean and SD would be computed over rows that
#   the model never sees.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
})

#' Apply treatment coding to S_Type with Passive as the reference level.
#'
#' @param df A data frame with an S_Type column.
#' @return df with S_Type as a factor whose first level is Passive.
#' @details The absence of Passive is an error rather than a fallback to whatever
#'   level happens to sort first. Without it there is no baseline against which
#'   the H1b interaction is defined, so a silently relevelled factor would produce
#'   coefficients that carry the expected names while answering a different
#'   question. R's default would be alphabetical, making Active the reference and
#'   inverting the sign of the focal interaction.
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
#'
#' @param df A data frame with an S_Type column.
#' @param has_pseudo_passive Logical scalar, taken from the per-language settings
#'   in config/analysis_config.yaml. FALSE for Norwegian and Balinese.
#' @return df with Pseudo_Passive rows and the unused level removed when
#'   has_pseudo_passive is FALSE; unchanged otherwise.
#' @details Dropping the *level* as well as the rows is the point of this function.
#'   A factor level with no observations still generates a model term, which brms
#'   then tries to estimate from nothing, and it also makes the model formula
#'   reference a coefficient that the priors in R/03_define_priors.R would have to
#'   supply. Removing the level keeps the model matrix consistent with the data, so
#'   the same ladder and prior code serves every language.
drop_pseudo_passive_if_absent <- function(df, has_pseudo_passive) {
  if (!has_pseudo_passive) {
    df <- dplyr::filter(df, S_Type != "Pseudo_Passive")
    # droplevels() only applies once S_Type is a factor. Called on a character
    # column it would be a no-op, hence the guard rather than an unconditional call.
    if (is.factor(df$S_Type)) {
      df <- dplyr::mutate(df, S_Type = droplevels(S_Type))
    }
    message("[preprocess] Pseudo_Passive dropped for this language.")
  }
  df
}

#' Scale Semantics to mean 0 and SD 0.5, within each language by default.
#'
#' @param df Data frame with a numeric Semantics column.
#' @param centre_by Column to group by before scaling. When the column is absent
#'   from `df`, scaling is applied to the whole data set instead; this is what makes
#'   the function usable on both single-language and pooled data without a branch
#'   at every call site.
#' @return df with Semantics_scaled added. Semantics itself is left untouched, so
#'   the raw affectedness values remain available for descriptive plots.
#' @details Dividing by twice the standard deviation, rather than once, is the
#'   convention of Gelman (2008, doi:10.1002/sim.3107): it leaves a continuous
#'   predictor with SD 0.5, which is the SD of a balanced binary predictor coded
#'   -0.5/+0.5. A prior scale is then interpretable in the same units whether the
#'   predictor is continuous or categorical, which is what allows the single set of
#'   prior scales in R/03_define_priors.R to apply to both.
#'
#'   Scaling within language is deliberate. Affectedness ratings were collected
#'   separately per language and their raw spread differs, so a pooled scaling
#'   would let one language's rating variance set the units for another's
#'   coefficients. The consequence to keep in mind when reading cross-language
#'   models is that a one-unit change means "two SDs of affectedness *in that
#'   language*", so slopes are comparable in standardised terms and not in raw
#'   rating points.
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
  # Greedy ".*_" strips everything up to and including the LAST underscore, which
  # is what makes this correct for verbs that themselves contain an underscore.
  g   <- sub(".*_", "", as.character(df[[item_col]]))
  # Anything other than the two expected tokens is an error. Silently accepting a
  # third value would add a level to the covariate, and a mis-parsed item would
  # otherwise surface only as an unexplained extra coefficient at fit time.
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
  # Order matters here; see the note in the file header. Re-sourcing Semantics
  # must precede scaling, and dropping unused rows must precede it too, so that
  # the scaling constants are computed over exactly the rows that will be fitted.
  df <- set_semantics_source(df, semantics_source)
  df <- df |>
    code_s_type() |>
    drop_pseudo_passive_if_absent(has_pseudo_passive) |>
    scale_semantics()
  if (isTRUE(include_gender)) df <- derive_gender(df)
  df
}

#' Assert that S_Type is treatment-coded with Passive as its reference level.
#'
#' @param df A preprocessed data frame.
#' @return Invisibly TRUE; called for the error it raises on bad coding.
#' @details A guard for the assumption the whole file rests on, intended to be
#'   called after preprocessing and before fitting, so that a data frame assembled
#'   by some other route cannot reach a model with the wrong baseline.
#'
#'   Four things are checked, in order of how basic they are: that S_Type is a
#'   factor at all, that it has at least two levels, that Passive is the first of
#'   them, and that its contrast matrix is the treatment one. The last two are
#'   separate properties. Reordering the levels changes which level is the baseline;
#'   assigning contr.sum keeps Passive first but changes the coefficients from
#'   "Active minus Passive" to deviations from the grand mean. Either would leave the
#'   coefficient NAMES intact while altering what they estimate, which is precisely
#'   the failure this function exists to make loud.
#'
#'   Note that droplevels() discards a previously assigned contrast, so a factor
#'   passed through drop_pseudo_passive_if_absent() arrives here carrying R's
#'   default regardless of what it had before.
assert_treatment_coding <- function(df) {
  if (!is.factor(df$S_Type)) stop("[assert] S_Type must be a factor after preprocessing.")
  ref <- levels(df$S_Type)[1]
  if (ref != "Passive") {
    stop("[assert] Reference level of S_Type is '", ref, "', expected 'Passive'.")
  }
  # Contrasts are undefined below two levels, and both calls below would otherwise
  # fail with R's opaque "contrasts not defined for 0 degrees of freedom".
  if (nlevels(df$S_Type) < 2L) {
    stop("[assert] S_Type has ", nlevels(df$S_Type),
         " level(s); at least two are needed to define a contrast.")
  }

  # Compare the contrast MATRIX, not the factor's "contrasts" attribute. That
  # attribute is NULL both for a factor carrying R's default contr.treatment and
  # for one explicitly assigned contr.sum, so it cannot tell them apart; an earlier
  # version of this check tested it and was therefore inert.
  #
  # contr.treatment() is given levels(), not nlevels(). Passing the level names
  # reproduces the dimnames that contrasts() returns, so a correctly coded factor
  # compares equal. Passing the count instead yields dimnames "1", "2", "3", and
  # all.equal() would then report a dimnames mismatch for every input, including
  # correct ones.
  observed <- stats::contrasts(df$S_Type)
  expected <- stats::contr.treatment(levels(df$S_Type))
  cmp <- all.equal(observed, expected)
  if (!isTRUE(cmp)) {
    stop("[assert] S_Type does not carry treatment contrasts (levels: ",
         paste(levels(df$S_Type), collapse = ", "), "). ",
         "A non-default contrast changes what every S_Type coefficient means, so ",
         "the priors in R/03_define_priors.R and the hypothesis tests in ",
         "R/05_hypothesis_tests.R would no longer refer to the intended quantities. ",
         "Mismatch: ", paste(cmp, collapse = "; "))
  }
  invisible(TRUE)
}
