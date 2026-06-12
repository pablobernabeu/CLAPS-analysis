# tests/testthat/test-gender-variation.R
# Unit tests for the gender model variation:
#   R/02_preprocess_factors.R: set_semantics_source(), derive_gender(), preprocess_data()
#   R/04_model_formulas.R:     build_model_ladder(include_gender = TRUE)

library(testthat)
source("R/02_preprocess_factors.R")

mk_df <- function() {
  tibble::tibble(
    Participant = rep(1:2, each = 4),
    Language    = "English",
    Verb        = rep(c("push", "see"), 4),
    Item        = rep(c("push_Man", "see_Woman", "push_Woman", "see_Man"), 2),
    S_Type      = rep(c("Passive", "Active"), 4),
    Semantics   = as.numeric(1:8),
    affectedness_scores_agent = as.numeric(8:1),
    Response    = rep(c(3L, 5L), 4)
  )
}

test_that("derive_gender extracts Man/Woman from Item as a Man-reference factor", {
  d <- derive_gender(mk_df())
  expect_true("Gender" %in% names(d))
  expect_setequal(levels(d$Gender), c("Man", "Woman"))
  expect_identical(levels(d$Gender)[1], "Man")
})

test_that("derive_gender errors on unexpected gender tokens", {
  bad <- mk_df(); bad$Item[1] <- "push_Robot"
  expect_error(derive_gender(bad), "Unexpected gender token")
})

test_that("set_semantics_source switches the Semantics column; NULL is a no-op", {
  d0 <- mk_df()
  d1 <- set_semantics_source(d0, "affectedness_scores_agent")
  expect_equal(d1$Semantics, as.numeric(d0$affectedness_scores_agent))
  expect_identical(set_semantics_source(d0, NULL)$Semantics, d0$Semantics)
  expect_error(set_semantics_source(d0, "nope"), "not found")
})

test_that("preprocess_data gender variation adds Gender and re-sources Semantics", {
  d0 <- mk_df()
  d  <- preprocess_data(d0, has_pseudo_passive = TRUE,
                        semantics_source = "affectedness_scores_agent",
                        include_gender = TRUE)
  expect_true(all(c("Gender", "Semantics_scaled") %in% names(d)))
  expect_s3_class(d$Gender, "factor")
  expect_equal(d$Semantics, as.numeric(d0$affectedness_scores_agent))
})

test_that("preprocess_data default is unchanged (no Gender, Semantics kept)", {
  d0 <- mk_df()
  d  <- preprocess_data(d0, has_pseudo_passive = TRUE)
  expect_false("Gender" %in% names(d))
  expect_equal(d$Semantics, d0$Semantics)
})

test_that("build_model_ladder includes a Gender term only when requested", {
  skip_if_not_installed("brms")
  source("R/04_model_formulas.R")
  base <- build_model_ladder(has_pseudo_passive = TRUE, include_gender = FALSE)
  gend <- build_model_ladder(has_pseudo_passive = TRUE, include_gender = TRUE)
  base_chr <- paste(deparse(base[["L5_correlated_maximal"]]$formula), collapse = " ")
  gend_chr <- paste(deparse(gend[["L5_correlated_maximal"]]$formula), collapse = " ")
  expect_false(grepl("Gender", base_chr))
  expect_true(grepl("\\+ Gender", gend_chr))
  expect_identical(attr(gend, "include_gender"), TRUE)
})
