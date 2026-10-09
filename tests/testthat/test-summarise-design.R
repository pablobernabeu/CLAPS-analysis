# tests/testthat/test-summarise-design.R
# Unit tests for the design summaries in R/08_summarise_design.R: that the gender arm
# is never averaged with the baseline arm, and that the maximal feasible model and the
# recommended sample size are read from the baseline arm (code audit of September 2026,
# finding SEP-3).

library(testthat)

# R/08_summarise_design.R sources R/05_hypothesis_tests.R, which calls library(brms),
# and does so by a path relative to design_analysis/, so it is sourced from there.
skip_if_not_installed("brms")
withr::with_dir(here::here(), source("R/08_summarise_design.R"))

# Two arms at the same design points: the baseline (no gender_spec) and a gender arm
# whose H1b rates are much lower, so pooling them would be visible.
cells <- function() {
  grid <- expand.grid(n_participants = c(50L, 80L), seed = 1:10, arm = c("none", "main"),
                      hypothesis = c("H1a_semantics_positive", "H1b_active_interaction_negative"),
                      stringsAsFactors = FALSE)
  # Baseline: every study reaches BF 10 at N = 80, none at N = 50. Gender arm: H1b
  # reaches it in only 2 of 10 studies at N = 80.
  gender_h1b_miss <- grid$arm == "main" & grid$hypothesis == "H1b_active_interaction_negative" & grid$seed > 2L
  hit <- grid$n_participants == 80L & !gender_h1b_miss
  data.frame(cell_id = seq_len(nrow(grid)), status = "success", language = "English",
             model_level = "L5_correlated_maximal", n_verbs = 72L, prior_regime = "primary",
             n_participants = grid$n_participants, seed = grid$seed, hypothesis = grid$hypothesis,
             gender_spec = grid$arm, BF_10 = ifelse(hit, 30, 2), convergence_ok = TRUE)
}

test_that("the exceedance table keeps the two arms apart", {
  exc <- compute_bf_exceedance(cells())
  expect_true("gender_spec" %in% names(exc))
  expect_equal(nrow(exc), 2L * 2L * 2L)
  h1b80 <- exc[exc$hypothesis == "H1b_active_interaction_negative" & exc$n_participants == 80L, ]
  expect_equal(h1b80$p_bf_primary[h1b80$gender_spec == "none"], 1)
  expect_equal(h1b80$p_bf_primary[h1b80$gender_spec == "main"], 0.2)
  expect_true(all(exc$n_sims == 10L))
})

test_that("cells without the arm columns are summarised as before", {
  d <- cells()
  d <- d[d$gender_spec == "none", setdiff(names(d), "gender_spec")]
  exc <- compute_bf_exceedance(d)
  expect_false("gender_spec" %in% names(exc))
  expect_equal(nrow(exc), 2L * 2L)
})

test_that("the recommendation is read from the baseline arm", {
  d <- cells()
  exc <- compute_bf_exceedance(d)
  rec <- recommended_sample_size(exc, maximal_feasible_model(d))
  expect_equal(rec$recommended_n_participants, 80L)
  expect_true(rec$meets_target)
})

test_that("baseline_arm keeps rows with no gender term and drops the rest", {
  d <- data.frame(gender_spec = c(NA, "none", "main", "interaction"), include_gender = c(NA, FALSE, TRUE, NA))
  expect_equal(nrow(baseline_arm(d)), 2L)
  expect_equal(nrow(baseline_arm(data.frame(x = 1:3))), 3L)
})
