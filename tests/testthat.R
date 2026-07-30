# tests/testthat.R
#
# Entry point for the unit-test suite, invoked by CI
# (.github/workflows/static-checks.yaml) and runnable by hand from
# design_analysis/ with:
#
#   Rscript tests/testthat.R
#
# What the suite does and does not cover
#   These are fast unit tests of the deterministic parts of the workflow: schema
#   validation, factor coding and scaling, prior construction, the Bayes-factor
#   arithmetic, and diagnostic classification. They run in seconds and need no
#   sampling.
#
#   They deliberately do NOT fit models, but three of the six files cannot even be
#   loaded without brms installed, because they source a module whose own
#   library(brms) call fails. testthat reports that load failure as a SKIP rather
#   than an error, so the whole file is silently passed over:
#
#     test-diagnostics.R        all  9 test_that blocks skipped
#     test-hypothesis-tests.R   all 14 test_that blocks skipped
#     test-priors.R             all  6 test_that blocks skipped
#     test-gender-variation.R    1 of 6 skipped (an explicit skip_if_not_installed)
#
#   So without brms roughly 30 checks never run, and only test-data-validation.R and
#   test-seed-disjointness.R are fully exercised. CI does install brms (the "Install
#   test dependencies" step in .github/workflows/static-checks.yaml lists it), so
#   those checks do run there.
#
#   READ THE SKIP COUNT, NOT JUST THE FAIL COUNT. A local run reporting no failures
#   has exercised far less than CI. Note also that only the gender-variation skip is
#   a deliberate one; the other three are load failures being reported leniently, so
#   a brms that is installed but broken would present the same way.
#
#   Correctness of the simulation engines is not unit-tested; it is checked by the
#   smoke tests in hpc/ and by the convergence diagnostics recorded for every fit.

library(testthat)
testthat::test_dir("tests/testthat")
