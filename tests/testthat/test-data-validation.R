# tests/testthat/test-data-validation.R
# Unit tests for data validation (R/01_read_validate_data.R)
# and factor coding (R/02_preprocess_factors.R).

library(testthat)
source("R/01_read_validate_data.R")
source("R/02_preprocess_factors.R")

# Helper: build a minimal valid CLAPS data frame
make_valid_df <- function(has_pp = TRUE) {
  s_types <- if (has_pp) c("Passive","Active","Pseudo_Passive") else c("Passive","Active")
  n <- length(s_types) * 5
  data.frame(
    Participant    = rep(paste0("P", seq_len(5)), each = length(s_types)),
    Language       = "English",
    Verb           = paste0("V", seq_len(n)),
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

test_that("scale_semantics produces near-zero mean and ~0.25 SD", {
  df <- make_valid_df()
  result <- scale_semantics(df, centre_by = "Language")
  expect_lt(abs(mean(result$Semantics_scaled)), 0.05)
  # Gelman scaling: SD ≈ 0.25 (= 0.5 / 2)
  expect_lt(abs(sd(result$Semantics_scaled) - 0.25), 0.1)
})

# ---------------------------------------------------------------------------
# preprocess_data
# ---------------------------------------------------------------------------

test_that("preprocess_data produces treatment-coded S_Type with Passive reference", {
  df <- make_valid_df()
  result <- preprocess_data(df, has_pseudo_passive = TRUE)
  assert_treatment_coding(result)
})

test_that("preprocess_data without Pseudo_Passive drops that level", {
  df <- make_valid_df(has_pp = TRUE)
  result <- preprocess_data(df, has_pseudo_passive = FALSE)
  expect_false("Pseudo_Passive" %in% levels(result$S_Type))
})
