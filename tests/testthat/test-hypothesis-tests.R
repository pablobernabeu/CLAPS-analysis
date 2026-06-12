# tests/testthat/test-hypothesis-tests.R
# Unit tests for Savage-Dickey BF implementation (R/05_hypothesis_tests.R).

library(testthat)
source("R/05_hypothesis_tests.R")

# ---------------------------------------------------------------------------
# classify_bf
# ---------------------------------------------------------------------------

test_that("classify_bf returns extreme_H1 for BF > 100", {
  expect_equal(classify_bf(150), "extreme_H1")
})

test_that("classify_bf returns very_strong_H1 for BF 30-100", {
  expect_equal(classify_bf(50), "very_strong_H1")
})

test_that("classify_bf returns strong_H1 for BF 10-30", {
  expect_equal(classify_bf(15), "strong_H1")
})

test_that("classify_bf returns moderate_H1 for BF 3-10", {
  expect_equal(classify_bf(5), "moderate_H1")
})

test_that("classify_bf returns anecdotal_H1 for BF 1-3", {
  expect_equal(classify_bf(2), "anecdotal_H1")
})

test_that("classify_bf returns anecdotal_H0 for BF 1/3 to 1", {
  expect_equal(classify_bf(0.5), "anecdotal_H0")
})

test_that("classify_bf returns moderate_H0 for BF 1/10 to 1/3", {
  expect_equal(classify_bf(0.2), "moderate_H0")
})

test_that("classify_bf returns strong_or_more_H0 for BF < 1/10", {
  expect_equal(classify_bf(0.05), "strong_or_more_H0")
})

test_that("classify_bf handles exactly 1 (anecdotal_H0)", {
  expect_equal(classify_bf(1.0), "anecdotal_H0")
})

# ---------------------------------------------------------------------------
# savage_dickey_directional_bf
# ---------------------------------------------------------------------------

test_that("savage_dickey_directional_bf errors if prior_samples is missing", {
  expect_error(
    savage_dickey_directional_bf(
      posterior_samples = rnorm(2000, 0.5, 0.1),
      prior_samples     = NULL,
      direction         = "positive"
    ),
    regexp = "prior_samples"
  )
})

test_that("savage_dickey_directional_bf returns list with bf and log_bf", {
  set.seed(1)
  result <- savage_dickey_directional_bf(
    posterior_samples = rnorm(4000, 0.5, 0.3),
    prior_samples     = rnorm(4000, 0.0, 1.0),
    direction         = "positive"
  )
  expect_true(is.list(result))
  expect_true(!is.null(result$bf))
  expect_true(!is.null(result$log_bf))
  expect_true(is.finite(result$bf))
  expect_true(is.finite(result$log_bf))
})

test_that("savage_dickey_directional_bf BF > 1 when posterior clearly positive", {
  set.seed(2)
  result <- savage_dickey_directional_bf(
    posterior_samples = rnorm(4000, 2.0, 0.2),
    prior_samples     = rnorm(4000, 0.0, 1.0),
    direction         = "positive"
  )
  expect_gt(result$bf, 1.0)
})

test_that("savage_dickey_directional_bf BF < 1 when posterior clearly negative (wrong direction)", {
  set.seed(3)
  result <- savage_dickey_directional_bf(
    posterior_samples = rnorm(4000, -2.0, 0.2),
    prior_samples     = rnorm(4000,  0.0, 1.0),
    direction         = "positive"
  )
  expect_lt(result$bf, 1.0)
})
