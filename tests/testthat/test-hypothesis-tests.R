# tests/testthat/test-hypothesis-tests.R
# Unit tests for Savage-Dickey BF implementation (R/05_hypothesis_tests.R).

library(testthat)
source(here::here("R", "05_hypothesis_tests.R"))

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

test_that("classify_bf puts the boundary value 1 on the H1 side", {
  # A Bayes factor of exactly 1 is no evidence either way. The bands are
  # half-open upwards, so 1 falls in anecdotal_H1 and anything below it in
  # anecdotal_H0.
  expect_equal(classify_bf(1.0), "anecdotal_H1")
  expect_equal(classify_bf(0.999), "anecdotal_H0")
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

# These exercise directional_bf_from_draws(), the arithmetic core of
# savage_dickey_directional_bf(). The wrapper itself takes a fitted brms model
# and only extracts the two sets of draws before handing them over, so testing
# the core on plain vectors covers the inferential content without fitting.

test_that("directional_bf_from_draws returns finite bf and log_bf", {
  set.seed(1)
  result <- directional_bf_from_draws(
    post_vals  = rnorm(4000, 0.5, 0.3),
    prior_vals = rnorm(4000, 0.0, 1.0),
    direction  = "positive"
  )
  expect_true(is.list(result))
  expect_true(is.finite(result$bf_10))
  expect_true(is.finite(result$log_bf))
  expect_equal(result$log_bf, log(result$bf_10))
})

test_that("directional_bf_from_draws gives BF > 1 when the posterior has the predicted sign", {
  set.seed(2)
  result <- directional_bf_from_draws(
    post_vals  = rnorm(4000, 2.0, 0.2),
    prior_vals = rnorm(4000, 0.0, 1.0),
    direction  = "positive"
  )
  expect_gt(result$bf_10, 1.0)
})

test_that("directional_bf_from_draws gives BF < 1 when the posterior has the opposite sign", {
  set.seed(3)
  result <- directional_bf_from_draws(
    post_vals  = rnorm(4000, -2.0, 0.2),
    prior_vals = rnorm(4000,  0.0, 1.0),
    direction  = "positive"
  )
  expect_lt(result$bf_10, 1.0)
})

test_that("directional_bf_from_draws is symmetric under reversing the direction", {
  set.seed(4)
  post  <- rnorm(4000, 1.0, 0.5)
  prior <- rnorm(4000, 0.0, 1.0)
  pos <- directional_bf_from_draws(post, prior, "positive")
  neg <- directional_bf_from_draws(post, prior, "negative")
  # With continuous draws no value is exactly zero, so the two directions
  # partition the probability. Both the posterior and the prior odds invert,
  # and the two Bayes factors are therefore reciprocals whatever the prior.
  expect_equal(pos$bf_10 * neg$bf_10, 1, tolerance = 1e-6)
})
