# tests/testthat/test-data-validation.R
# Unit tests for data validation (R/01_read_validate_data.R)
# and factor coding (R/02_preprocess_factors.R).

library(testthat)
source(here::here("R", "01_read_validate_data.R"))
source(here::here("R", "02_preprocess_factors.R"))

# Helper: build a minimal valid CLAPS data frame
make_valid_df <- function(has_pp = TRUE) {
  s_types <- if (has_pp) c("Passive","Active","Pseudo_Passive") else c("Passive","Active")
  n <- length(s_types) * 5
  data.frame(
    Participant    = rep(paste0("P", seq_len(5)), each = length(s_types)),
    Language       = "English",
    Verb           = paste0("V", seq_len(n)),
    # Verb_ID is the language-prefixed verb key; the harmonised pilot carries
    # both it and Verb, and validate_raw_data requires both.
    Verb_ID        = paste0("English_V", seq_len(n)),
    Item           = seq_len(n),
    S_Type         = rep(s_types, 5),
    Semantics      = runif(n, -0.5, 0.5),
    Response       = sample(1L:7L, n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# validate_raw_data
# ---------------------------------------------------------------------------

test_that("validate_raw_data passes on valid data", {
  df <- make_valid_df()
  expect_silent(validate_raw_data(df, "test"))
})

test_that("validate_raw_data fails on missing required column", {
  df <- make_valid_df()
  df$Response <- NULL
  expect_error(validate_raw_data(df, "test"), regexp = "Missing required columns")
})

test_that("validate_raw_data fails on out-of-range Response", {
  df <- make_valid_df()
  df$Response[1] <- 8L
  expect_error(validate_raw_data(df, "test"))
})

test_that("validate_raw_data fails on invalid S_Type", {
  df <- make_valid_df()
  df$S_Type[1] <- "Impersonal"
  expect_error(validate_raw_data(df, "test"), regexp = "Unexpected S_Type")
})

test_that("validate_raw_data fails on NA in key column", {
  df <- make_valid_df()
  df$Participant[1] <- NA
  expect_error(validate_raw_data(df, "test"), regexp = "NA")
})

# ---------------------------------------------------------------------------
# code_s_type
# ---------------------------------------------------------------------------

test_that("code_s_type sets Passive as reference level", {
  df <- make_valid_df()
  result <- code_s_type(df)
  expect_true(is.factor(result$S_Type))
  expect_equal(levels(result$S_Type)[1], "Passive")
})

test_that("code_s_type fails without Passive level", {
  df <- make_valid_df()
  df$S_Type[df$S_Type == "Passive"] <- "Active"
  expect_error(code_s_type(df), regexp = "Passive")
})

# ---------------------------------------------------------------------------
# drop_pseudo_passive_if_absent
# ---------------------------------------------------------------------------

test_that("drop_pseudo_passive_if_absent removes Pseudo_Passive rows", {
  df <- make_valid_df(has_pp = TRUE) |> code_s_type()
  result <- drop_pseudo_passive_if_absent(df, has_pseudo_passive = FALSE)
  expect_false("Pseudo_Passive" %in% levels(result$S_Type))
  expect_false(any(result$S_Type == "Pseudo_Passive"))
})

test_that("drop_pseudo_passive_if_absent retains all levels when TRUE", {
  df <- make_valid_df(has_pp = TRUE) |> code_s_type()
  result <- drop_pseudo_passive_if_absent(df, has_pseudo_passive = TRUE)
  expect_true("Pseudo_Passive" %in% levels(result$S_Type))
})

# ---------------------------------------------------------------------------
# scale_semantics
# ---------------------------------------------------------------------------

test_that("scale_semantics produces Semantics_scaled column", {
  df <- make_valid_df()
  result <- scale_semantics(df, centre_by = "Language")
  expect_true("Semantics_scaled" %in% names(result))
})

test_that("scale_semantics produces near-zero mean and ~0.5 SD", {
  df <- make_valid_df()
  result <- scale_semantics(df, centre_by = "Language")
  expect_lt(abs(mean(result$Semantics_scaled)), 0.05)
  # Gelman scaling divides the centred predictor by twice its standard
  # deviation, so the rescaled predictor has an SD of 0.5, not 0.25. The
  # analysis priors are set on this scale, so the constant is load-bearing.
  expect_lt(abs(sd(result$Semantics_scaled) - 0.5), 0.1)
})

# ---------------------------------------------------------------------------
# preprocess_data
# ---------------------------------------------------------------------------

test_that("preprocess_data produces treatment-coded S_Type with Passive reference", {
  df <- make_valid_df()
  result <- preprocess_data(df, has_pseudo_passive = TRUE)
  expect_true(assert_treatment_coding(result))
})

test_that("assert_treatment_coding accepts the two-level (no pseudo-passive) case", {
  # Norwegian and Balinese have only Passive and Active, so the assertion must hold
  # for a two-column contrast matrix as well as a three-column one.
  df <- make_valid_df(has_pp = FALSE)
  result <- preprocess_data(df, has_pseudo_passive = FALSE)
  expect_true(assert_treatment_coding(result))
})

# ---------------------------------------------------------------------------
# assert_treatment_coding: the negative cases.
#
# These are what make the assertion worth having. Each below leaves the
# coefficient NAMES that the priors and hypothesis tests refer to intact while
# changing what those coefficients estimate, so nothing downstream would fail
# loudly on its own.
# ---------------------------------------------------------------------------

test_that("assert_treatment_coding rejects a non-default contrast (contr.sum)", {
  df <- preprocess_data(make_valid_df(), has_pseudo_passive = TRUE)
  # Passive stays the first level, so the reference-level check still passes and
  # only the contrast-matrix check can catch this.
  stats::contrasts(df$S_Type) <- stats::contr.sum(nlevels(df$S_Type))
  expect_error(assert_treatment_coding(df),
               regexp = "does not carry treatment contrasts")
})

test_that("assert_treatment_coding rejects contr.helmert", {
  df <- preprocess_data(make_valid_df(), has_pseudo_passive = TRUE)
  stats::contrasts(df$S_Type) <- stats::contr.helmert(nlevels(df$S_Type))
  expect_error(assert_treatment_coding(df),
               regexp = "does not carry treatment contrasts")
})

test_that("assert_treatment_coding rejects a wrong reference level", {
  df <- preprocess_data(make_valid_df(), has_pseudo_passive = TRUE)
  # R's alphabetical default would make Active the baseline, which inverts the sign
  # of the focal H1b interaction.
  df$S_Type <- factor(as.character(df$S_Type),
                      levels = c("Active", "Passive", "Pseudo_Passive"))
  expect_error(assert_treatment_coding(df), regexp = "Reference level")
})

test_that("assert_treatment_coding rejects a non-factor S_Type", {
  df <- make_valid_df()   # S_Type is still character here
  expect_error(assert_treatment_coding(df), regexp = "must be a factor")
})

test_that("assert_treatment_coding rejects a single-level S_Type", {
  # Guarded explicitly, because contrasts() on a one-level factor otherwise fails
  # with R's opaque "contrasts not defined for 0 degrees of freedom".
  df <- make_valid_df()
  df <- df[df$S_Type == "Passive", ]
  df$S_Type <- factor(df$S_Type, levels = "Passive")
  expect_error(assert_treatment_coding(df), regexp = "at least two")
})

test_that("preprocess_data without Pseudo_Passive drops that level", {
  df <- make_valid_df(has_pp = TRUE)
  result <- preprocess_data(df, has_pseudo_passive = FALSE)
  expect_false("Pseudo_Passive" %in% levels(result$S_Type))
})
