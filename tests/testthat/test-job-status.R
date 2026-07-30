# tests/testthat/test-job-status.R
#
# Tests for R/10_job_status.R, chiefly .design_cell_ids(), which rebuilds the output
# filenames a design grid is expected to produce so that progress can be measured
# from the files on disk.
#
# The expected names below are written out LITERALLY rather than computed. That is
# the point of these tests: .design_cell_ids() mirrors the cell_id construction in
# run_design_cell() (R/06_simulate_design.R), and the two cannot be shared because
# that module loads brms while the status scripts must run without a Stan toolchain.
# Deriving the expectations from either implementation would let both drift together
# unnoticed. Until 2026-07-30 the mirror reproduced only the seven base fields, so
# every gender and cross-language cell was reported as "pending" however long ago it
# had finished.

library(testthat)
source(here::here("R", "10_job_status.R"))

# A grid covering the four naming shapes. All rows carry all columns, with NA where
# a column does not apply, which is what a real grid looks like once the variants are
# crossed in.
mk_grid <- function() {
  tibble::tibble(
    language        = c("English", "English", "English", "AllLanguages"),
    model_level     = c(rep("L5_correlated_maximal", 3), "L4_cross_uncorrelated"),
    prior_regime    = "primary",
    threshold_mode  = "broad",
    n_participants  = c(30L, 30L, 30L, 30L),
    n_verbs         = c(72L, 72L, 72L, 20L),
    seed            = c(600001L, 300001L, 300002L, 900002L),
    include_gender  = c(NA, TRUE, NA, NA),
    gender_spec     = c(NA, NA, "interaction", NA),
    n_languages     = c(NA, NA, NA, 3L)
  )
}

EXPECTED <- c(
  "English_L5_correlated_maximal_primary_broad_30_72_600001",
  "English_L5_correlated_maximal_primary_broad_30_72_300001_gender",
  "English_L5_correlated_maximal_primary_broad_30_72_300002_genderX",
  "AllLanguages_L4_cross_uncorrelated_primary_broad_30_20_900002_3lang"
)

test_that(".design_cell_ids reproduces every naming shape run_design_cell writes", {
  expect_equal(.design_cell_ids(mk_grid()), EXPECTED)
})

test_that(".design_cell_ids treats a grid without the variant columns as baseline", {
  g <- mk_grid()[1, c("language", "model_level", "prior_regime", "threshold_mode",
                      "n_participants", "n_verbs", "seed")]
  expect_equal(.design_cell_ids(g), EXPECTED[1])
})

test_that(".design_cell_ids gives include_gender = TRUE the main-effect suffix", {
  g <- mk_grid()[2, ]
  g$gender_spec <- NA          # explicit spec absent, so include_gender decides
  expect_match(.design_cell_ids(g), "_gender$")
})

test_that(".design_cell_ids lets an explicit gender_spec override include_gender", {
  g <- mk_grid()[2, ]
  g$gender_spec <- "interaction"   # both set; the explicit spec must win
  expect_match(.design_cell_ids(g), "_genderX$")
})

test_that(".design_cell_ids appends the language count after any gender suffix", {
  g <- mk_grid()[4, ]
  g$gender_spec <- "main"
  expect_match(.design_cell_ids(g), "_gender_3lang$")
})

test_that(".design_cell_ids handles a zero-row grid", {
  expect_equal(.design_cell_ids(mk_grid()[0, ]), character(0))
})

test_that(".design_cell_ids preserves R's scientific rendering of a round seed", {
  # NOT a wart to tidy up. readr types the seed column as double, and paste() renders
  # a round double in scientific notation, so run_design_cell() has already written
  # files named "..._7e+05.rds". Reformatting the seed here would stop those cells
  # from ever being matched. Pinned so that changing it becomes a deliberate act.
  g <- mk_grid()[1, ]
  g$seed <- 700000            # double, exactly as readr would supply it
  expect_match(.design_cell_ids(g), "_7e\\+05$")
})

# ---------------------------------------------------------------------------
# check_grid_completion, end to end against files on disk
# ---------------------------------------------------------------------------

test_that("check_grid_completion marks gender and cross-language cells done", {
  dir <- withr::local_tempdir()
  # One file per expected cell, so every naming shape must match for this to pass.
  for (nm in EXPECTED) saveRDS(list(status = "success"), file.path(dir, paste0(nm, ".rds")))

  status <- check_grid_completion(mk_grid(), list_completed_cells(dir))
  expect_equal(status$status, rep("done", 4))
  expect_true(all(status$completed))
})

test_that("check_grid_completion reports a cell with no output as pending", {
  dir <- withr::local_tempdir()
  # Only the baseline cell has been written.
  saveRDS(list(status = "success"), file.path(dir, paste0(EXPECTED[1], ".rds")))

  status <- check_grid_completion(mk_grid(), list_completed_cells(dir))
  expect_equal(status$status, c("done", "pending", "pending", "pending"))
})

test_that("check_grid_completion returns one row per grid row, none dropped", {
  dir <- withr::local_tempdir()
  saveRDS(list(status = "success"), file.path(dir, paste0(EXPECTED[3], ".rds")))

  status <- check_grid_completion(mk_grid(), list_completed_cells(dir))
  expect_equal(nrow(status), 4L)
})

test_that("list_completed_cells returns an empty tibble for an empty directory", {
  dir <- withr::local_tempdir()
  expect_equal(nrow(list_completed_cells(dir)), 0L)
})
