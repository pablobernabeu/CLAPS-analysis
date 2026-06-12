# tests/testthat/test-priors.R
# Unit tests for prior construction (R/03_define_priors.R).

library(testthat)
library(brms)
source("R/03_define_priors.R")

# ---------------------------------------------------------------------------
# build_brms_prior
# ---------------------------------------------------------------------------

test_that("build_brms_prior returns a brmsprior for each regime × threshold", {
  grid <- prior_sensitivity_grid()
  for (i in seq_len(nrow(grid))) {
    p <- build_brms_prior(
      regime_name    = grid$regime_name[i],
      threshold_mode = grid$threshold_mode[i]
    )
    expect_s3_class(p, "brmsprior")
    expect_gt(nrow(p), 0L)
  }
})

test_that("build_brms_prior with ceiling_calibrated uses threshold_params", {
  tp <- list(
    threshold_means = c(-2.5, -1.0, 0.0, 1.0, 2.0, 3.0),
    threshold_sds   = rep(0.8, 6)
  )
  p <- build_brms_prior("primary", "ceiling_calibrated", threshold_params = tp)
  expect_s3_class(p, "brmsprior")
  # Should contain Intercept class priors
  expect_true(any(p$class == "Intercept"))
})

test_that("build_brms_prior fails on unknown regime", {
  expect_error(build_brms_prior("nonexistent", "broad"), regexp = "regime_name")
})

# ---------------------------------------------------------------------------
# compute_ceiling_calibrated_thresholds
# ---------------------------------------------------------------------------

test_that("compute_ceiling_calibrated_thresholds returns 6 threshold means", {
  set.seed(1)
  pilot_df <- data.frame(Response = sample(1:7, 100, replace = TRUE,
                                            prob = c(0.01, 0.02, 0.05, 0.12, 0.20, 0.30, 0.30)))
  tp <- compute_ceiling_calibrated_thresholds(pilot_df, "test_lang")
  expect_length(tp$threshold_means, 6L)
  expect_length(tp$threshold_sds,   6L)
  expect_true(all(is.finite(tp$threshold_means)))
})

test_that("compute_ceiling_calibrated_thresholds thresholds are monotone increasing", {
  set.seed(42)
  pilot_df <- data.frame(Response = sample(1:7, 200, replace = TRUE))
  tp <- compute_ceiling_calibrated_thresholds(pilot_df, "test_lang")
  expect_true(all(diff(tp$threshold_means) > 0))
})

# ---------------------------------------------------------------------------
# prior_sensitivity_grid
# ---------------------------------------------------------------------------

test_that("prior_sensitivity_grid has 8 rows (4 regimes x 2 modes)", {
  grid <- prior_sensitivity_grid()
  expect_equal(nrow(grid), 8L)
  expect_setequal(grid$regime_name, names(PRIOR_REGIMES))
  expect_setequal(grid$threshold_mode, THRESHOLD_MODES)
})
