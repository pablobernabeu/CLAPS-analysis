# tests/testthat/test-environment.R
#
# The gate. Runs first (testthat takes files alphabetically, and "environment"
# sorts ahead of the rest) and fails the suite outright when the session cannot
# actually exercise it.
#
# The problem this solves. Three test files source modules that call
# library(brms). Without brms those files fail to LOAD, and testthat reports a load
# failure as a SKIP, not a failure. The suite then printed
#   [ FAIL 0 | WARN 0 | SKIP 4 | PASS 37 ]
# which reads as success while 29 of the 66 checks had not run. Nothing in that
# output distinguishes "everything passed" from "most of it never executed". A test
# suite that reports green in an environment where it cannot do its job is worse
# than no suite, because it is trusted.
#
# The rule now: a missing core or modelling package FAILS here, unless the reader
# explicitly says they know. That opt-out exists because a partial run is genuinely
# useful while editing on a machine with no Stan toolchain:
#
#   CLAPS_ALLOW_PARTIAL_TESTS=1 Rscript tests/testthat.R
#
# which downgrades the gate to a skip carrying the same message. CI sets nothing, so
# CI cannot pass while its environment is incomplete.

library(testthat)
source(here::here("R", "12_environment.R"))

.partial_allowed <- function() {
  v <- tolower(Sys.getenv("CLAPS_ALLOW_PARTIAL_TESTS", ""))
  v %in% c("1", "true", "yes")
}

.report <- function(missing) {
  paste0(
    "The session is missing packages this suite needs, so part of it cannot run.\n",
    paste(vapply(names(missing), function(g) {
      m <- missing[[g]]
      sprintf("    %-10s %s", g, if (length(m)) paste(m, collapse = ", ") else "(complete)")
    }, character(1)), collapse = "\n"),
    "\n  Install them, or set CLAPS_ALLOW_PARTIAL_TESTS=1 to accept a partial run.\n",
    "  See docs/reproducibility.md for the environment this analysis expects."
  )
}

test_that("the session can actually run this suite", {
  # core and modelling are the two groups the tests themselves depend on. pipeline
  # and reporting are not checked here: no test exercises them, so demanding them
  # would block a perfectly adequate test environment.
  missing <- claps_missing_packages(c("core", "modelling"))
  incomplete <- sum(lengths(missing)) > 0

  if (incomplete) {
    cat("\n--- environment ---\n")
    claps_environment_summary()
    cat("-------------------\n")
    if (.partial_allowed()) {
      # Deliberately visible: a partial run should be obvious in the log, not a
      # quiet line among the skips.
      cat("\n!! PARTIAL TEST RUN — CLAPS_ALLOW_PARTIAL_TESTS is set !!\n")
      cat(.report(missing), "\n\n")
      skip("partial run accepted via CLAPS_ALLOW_PARTIAL_TESTS")
    }
    fail(.report(missing))
  }
  succeed()
})

test_that("the requirement groups are well formed", {
  # Guards the declaration itself: a typo that emptied a group would make the gate
  # above pass vacuously.
  expect_true(all(c("core", "modelling") %in% names(CLAPS_REQUIREMENTS)))
  expect_true(all(lengths(CLAPS_REQUIREMENTS) > 0))
  expect_false(anyDuplicated(unlist(CLAPS_REQUIREMENTS)) > 0)
  expect_true(all(vapply(unlist(CLAPS_REQUIREMENTS), nzchar, logical(1))))
})

test_that("the recorded ARC environment is readable and covers the analysis", {
  rec <- here::here("config", "arc_environment_recorded.csv")
  skip_if_not(file.exists(rec), "no recorded ARC environment")

  cmp <- claps_compare_to_record(rec)
  expect_gt(nrow(cmp), 100)                       # a real library, not a stub
  expect_true(all(cmp$status %in% c("match", "differs", "absent_here")))

  # The record must describe an environment in which the analysis could run, or it
  # is not a record of the environment that produces the results.
  need <- c(CLAPS_REQUIREMENTS$core, CLAPS_REQUIREMENTS$modelling)
  expect_equal(setdiff(need, cmp$package), character(0))
})
