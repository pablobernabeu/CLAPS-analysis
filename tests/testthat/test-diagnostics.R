# tests/testthat/test-diagnostics.R
# Unit tests for convergence diagnostics (R/07_extract_diagnostics.R).
# Criteria are publication-grade: "converged" requires R-hat < 1.01, bulk & tail
# ESS >= 400, zero divergences and zero max-treedepth saturations.

library(testthat)
source("R/07_extract_diagnostics.R")

# Field names match the tibble returned by extract_convergence_diagnostics().
make_diag <- function(max_rhat = 1.00, min_ess_bulk = 2000, min_ess_tail = 1800,
                      n_divergent = 0, n_max_treedepth = 0) {
  tibble::tibble(
    max_rhat = max_rhat, min_ess_bulk = min_ess_bulk, min_ess_tail = min_ess_tail,
    n_divergent = n_divergent, n_max_treedepth = n_max_treedepth
  )
}

test_that("clean diagnostics are converged", {
  expect_equal(classify_convergence(make_diag()), "converged")
  expect_equal(classify_convergence(make_diag(max_rhat = 1.005)), "converged")
})

test_that("R-hat tiers", {
  expect_equal(classify_convergence(make_diag(max_rhat = 1.06)),  "failed_rhat")
  expect_equal(classify_convergence(make_diag(max_rhat = 1.03)),  "marginal_rhat")
  expect_equal(classify_convergence(make_diag(max_rhat = 1.011)), "marginal_rhat")
})

test_that("divergence tiers (publication-grade requires zero)", {
  expect_equal(classify_convergence(make_diag(n_divergent = 11)), "failed_divergences")
  expect_equal(classify_convergence(make_diag(n_divergent = 5)),  "marginal_divergences")
  expect_equal(classify_convergence(make_diag(n_divergent = 1)),  "marginal_divergences")
})

test_that("ESS tiers use both bulk and tail", {
  expect_equal(classify_convergence(make_diag(min_ess_bulk = 50)),   "failed_ess")
  expect_equal(classify_convergence(make_diag(min_ess_tail = 80)),   "failed_ess")
  expect_equal(classify_convergence(make_diag(min_ess_bulk = 300)),  "marginal_ess")
  expect_equal(classify_convergence(make_diag(min_ess_tail = 350)),  "marginal_ess")
})

test_that("treedepth saturation is not publication-grade", {
  expect_equal(classify_convergence(make_diag(n_max_treedepth = 5)), "marginal_treedepth")
})

test_that("R-hat failure takes precedence over divergences", {
  expect_equal(classify_convergence(make_diag(max_rhat = 1.10, n_divergent = 100)), "failed_rhat")
})

# ---------------------------------------------------------------------------
# select_highest_feasible_model (reads convergence_status; full ladder names)
# ---------------------------------------------------------------------------

test_that("returns the highest strictly-converged level", {
  ladder <- tibble::tibble(
    model_level = c("L5_correlated_maximal", "L4_uncorrelated_maximal",
                    "L3_no_participant_interaction_slope"),
    convergence_status = c("marginal_rhat", "converged", "converged")
  )
  expect_equal(select_highest_feasible_model(ladder), "L4_uncorrelated_maximal")
})

test_that("falls through to L0 when only L0 converges", {
  ladder <- tibble::tibble(
    model_level = c("L5_correlated_maximal", "L4_uncorrelated_maximal",
                    "L3_no_participant_interaction_slope",
                    "L2_sentence_type_slopes_only",
                    "L1_random_intercepts_plus_participant_semantics",
                    "L0_random_intercepts_only"),
    convergence_status = c("failed_rhat", "marginal_divergences", "failed_ess",
                           "marginal_rhat", "marginal_treedepth", "converged")
  )
  expect_equal(select_highest_feasible_model(ladder), "L0_random_intercepts_only")
})

test_that("returns NA when no level converges", {
  ladder <- tibble::tibble(
    model_level = c("L5_correlated_maximal", "L4_uncorrelated_maximal"),
    convergence_status = c("failed_rhat", "marginal_divergences")
  )
  expect_true(is.na(select_highest_feasible_model(ladder)))
})
