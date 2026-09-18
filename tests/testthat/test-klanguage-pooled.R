# tests/testthat/test-klanguage-pooled.R
#
# Tests for the K-language pooled design engine (R/13_simulate_klanguage_pooled.R).
#
# The engine's job is to stand in for a fit that costs 29 to 96 hours, so the
# things worth testing are the ones that would let it stand in WRONGLY: the
# quadrature not integrating to what it claims, the Bayes factor not meaning what
# R/05_hypothesis_tests.R means by it, and the shared-verb component not producing
# the correlation it is supposed to produce. Each test below has a closed form or
# an exact combinatorial answer, so none of them is a tolerance dressed up as a check.

source(testthat::test_path("..", "..", "R", "13_simulate_klanguage_pooled.R"))

test_that("with tau pinned at zero the posterior is the fixed-effect answer", {
  # A random-effects meta-analysis whose tau cannot move is a precision-weighted
  # mean with a closed form. Forcing the grid to the single point zero is the only
  # way to isolate the quadrature from the tau prior.
  y <- c(-0.44, -0.15, -0.27, 0.33)
  s <- c(0.10, 0.08, 0.09, 0.13)
  post <- ra_meta_posterior(y, s, tau_grid = 0)

  w <- 1 / s^2
  expect_equal(post$mu_mean, sum(w * y) / sum(w), tolerance = 1e-10)
  expect_equal(post$mu_sd,   sqrt(1 / sum(w)),    tolerance = 1e-10)
  expect_equal(post$p_mu_negative,
               pnorm(0, sum(w * y) / sum(w), sqrt(1 / sum(w))), tolerance = 1e-10)
})

test_that("a proper normal prior on mu shrinks the estimate by the closed-form factor", {
  y <- c(-0.5, -0.3); s <- c(0.2, 0.2); p <- 0.4
  post <- ra_meta_posterior(y, s, prior_mu_sd = p, tau_grid = 0)
  w    <- 1 / s^2
  prec <- sum(w) + 1 / p^2
  expect_equal(post$mu_mean, sum(w * y) / prec, tolerance = 1e-10)
  expect_equal(post$mu_sd,   sqrt(1 / prec),    tolerance = 1e-10)
})

test_that("the directional Bayes factor is the posterior odds on the sign", {
  # Same definition as R/05_hypothesis_tests.R::directional_bf_from_draws, so a
  # Bayes factor from this engine and one from a fitted brms model are the same
  # quantity. With a symmetric prior the prior odds are 1 and BF is P/(1-P).
  for (p in c(0.5, 0.75, 10 / 11, 0.95, 0.99)) {
    bf <- directional_bf_meta(list(p_mu_negative = p), direction = "negative")
    expect_equal(bf, p / (1 - p), tolerance = 1e-10)
  }
  # The threshold the whole design turns on: BF = 10 is exactly P = 10/11.
  expect_equal(directional_bf_meta(list(p_mu_negative = 10 / 11)), 10, tolerance = 1e-9)
  # And in z units, which is how the calibration against the brms runs is stated.
  expect_equal(qnorm(10 / 11), 1.33517773, tolerance = 1e-6)
})

test_that("the direction argument is not silently ignored", {
  post <- list(p_mu_negative = 0.9)
  expect_equal(directional_bf_meta(post, "negative"), 9,     tolerance = 1e-9)
  expect_equal(directional_bf_meta(post, "positive"), 1 / 9, tolerance = 1e-9)
})

test_that("the tau grid weights are a normalised distribution", {
  post <- ra_meta_posterior(c(-0.4, 0.3, -0.8), c(0.1, 0.13, 0.2))
  expect_equal(sum(post$weights), 1, tolerance = 1e-10)
  expect_true(all(post$weights >= 0))
  expect_true(post$tau_median > 0)
})

test_that("heterogeneous data widen the posterior relative to homogeneous data", {
  # Three languages that agree, against three that do not, at identical precision
  # and identical mean. The only thing that differs is the spread, so any
  # difference in mu_sd is the random-effects term doing its job.
  s <- rep(0.1, 3)
  agree    <- ra_meta_posterior(c(-0.30, -0.31, -0.29), s)
  disagree <- ra_meta_posterior(c(-0.80, -0.30,  0.20), s)
  expect_equal(agree$mu_mean, disagree$mu_mean, tolerance = 0.02)
  expect_gt(disagree$mu_sd, agree$mu_sd)
  expect_gt(disagree$tau_median, agree$tau_median)
})

test_that("the shared verb component induces the intended cross-language correlation", {
  # split_se as the driver computes it: a common draw of SD sqrt(rho)*se_verb plus
  # independent draws. The resulting correlation between two languages' errors must
  # be rho once the participant term is set to zero.
  set.seed(11)
  se_verb <- 0.074; rho <- 0.4
  shared <- se_verb * sqrt(rho); indep <- se_verb * sqrt(1 - rho)
  n <- 200000L
  common <- rnorm(n, 0, shared)
  e1 <- common + rnorm(n, 0, indep)
  e2 <- common + rnorm(n, 0, indep)
  expect_equal(cor(e1, e2), rho, tolerance = 0.01)
  # And the total per-language error is unchanged by how it is split.
  expect_equal(sd(e1), se_verb, tolerance = 0.01)
})

test_that("simulate_klanguage_replicate honours fixed language effects", {
  set.seed(3)
  r <- simulate_klanguage_replicate(3L, mu = 0, tau = 99, se_lang = 1e-8,
                                    theta_fixed = c(-0.5, -0.4, -0.45))
  # With negligible sampling error the estimate must sit at the mean of the fixed
  # effects, whatever tau says, because tau is not used when theta is supplied.
  expect_equal(r$mu_hat, mean(c(-0.5, -0.4, -0.45)), tolerance = 0.02)
})

test_that("power rises with K when the languages are drawn from one population", {
  set.seed(5)
  small <- run_klanguage_cell(K = 4,  n_reps = 400L, se_lang = 0.09, mode = "point",
                              mu_point = -0.3, tau_point = 0.3, seed = 1L)
  big   <- run_klanguage_cell(K = 16, n_reps = 400L, se_lang = 0.09, mode = "point",
                              mu_point = -0.3, tau_point = 0.3, seed = 2L)
  expect_gt(big$power, small$power)
  expect_true(small$power >= 0 && big$power <= 1)
})

test_that("enumerate_language_subsets returns every subset exactly once", {
  pop <- data.frame(language = c("A", "B", "C", "D", "E"),
                    y = c(-0.4, -0.2, -0.8, 0.3, -0.1), s = rep(0.1, 5))
  out <- enumerate_language_subsets(pop, K = 3, se_indep = 0.08, se_shared = 0.03,
                                    n_reps = 40L, seed = 7L)
  expect_equal(nrow(out), choose(5, 3))
  expect_equal(length(unique(out$languages)), choose(5, 3))
  # Each row names the languages it left out, and the two lists must partition.
  expect_true(all(vapply(seq_len(nrow(out)), function(i) {
    inc <- strsplit(out$languages[i], "\\+")[[1]]
    omt <- strsplit(out$omitted[i],   "\\+")[[1]]
    setequal(c(inc, omt), pop$language) && length(inc) == 3L
  }, logical(1))))
})

test_that("the calibrated tau prior scale is the one the design analysis documents", {
  # Guards against the constant drifting away from the value fitted to the 135
  # three-language brms replicates without the calibration being redone.
  expect_equal(CALIBRATED_TAU_PRIOR_SCALE, 1.1)
  expect_equal(CALIBRATION_BRACKET, c(1.0, 1.2))
  expect_true(CALIBRATED_TAU_PRIOR_SCALE > CALIBRATION_BRACKET[1] &&
              CALIBRATED_TAU_PRIOR_SCALE < CALIBRATION_BRACKET[2])
})

test_that("run_partly_known_cell reduces to the exchangeable arm when nothing is known", {
  pop <- data.frame(language = c("A", "B", "C", "D"),
                    y = c(-0.44, -0.15, -0.27, 0.33), s = c(0.12, 0.09, 0.06, 0.13))
  a <- run_partly_known_cell(K = 6, known = pop[0, ], n_reps = 1500L, se_lang = 0.089,
                             pop = pop, seed = 21L)
  b <- run_klanguage_cell(K = 6, n_reps = 1500L, se_lang = 0.089, mode = "assurance",
                          pop = pop, seed = 21L)
  expect_equal(a$n_known, 0L)
  # Same generative process, different seeds inside, so agreement is to Monte-Carlo
  # error: three standard errors of the difference.
  tol <- 3 * sqrt(a$mcse^2 + b$mcse^2)
  expect_lt(abs(a$power - b$power), tol)
})

test_that("knowing a favourable language's effect in advance raises power", {
  # Two sets differing only in whether three strongly-negative languages are known
  # in advance or drawn from a population with the same mean. Fixing them removes
  # tau's contribution for those places, so power must rise.
  pop <- data.frame(language = c("A", "B", "C", "D"),
                    y = c(-0.44, -0.15, -0.27, 0.33), s = c(0.12, 0.09, 0.06, 0.13))
  known <- data.frame(language = c("A", "B", "C"),
                      y = c(-0.44, -0.40, -0.36), s = c(0.12, 0.09, 0.06))
  none <- run_partly_known_cell(K = 8, known = pop[0, ], n_reps = 2000L,
                                se_lang = 0.089, pop = pop, seed = 31L)
  some <- run_partly_known_cell(K = 8, known = known, n_reps = 2000L,
                                se_lang = 0.089, pop = pop, seed = 32L)
  expect_equal(some$n_known, 3L)
  expect_gt(some$power, none$power)
})

test_that("doubt about protocol transfer removes the benefit of knowing a language", {
  # transfer_sd widens the prior on an already-measured language. Pushed well above
  # tau it should return that language to being effectively an unknown draw, so the
  # power advantage over the nothing-known case must shrink towards zero.
  pop <- data.frame(language = c("A", "B", "C", "D"),
                    y = c(-0.44, -0.15, -0.27, 0.33), s = c(0.12, 0.09, 0.06, 0.13))
  known <- data.frame(language = c("A", "B", "C"),
                      y = c(-0.44, -0.40, -0.36), s = c(0.12, 0.09, 0.06))
  tight <- run_partly_known_cell(K = 8, known = known, n_reps = 2000L, se_lang = 0.089,
                                 pop = pop, transfer_sd = 0, seed = 41L)
  loose <- run_partly_known_cell(K = 8, known = known, n_reps = 2000L, se_lang = 0.089,
                                 pop = pop, transfer_sd = 1.5, seed = 42L)
  expect_gt(tight$power, loose$power)
  expect_equal(loose$transfer_sd, 1.5)
})

test_that("a fully known set does not depend on the population it is not drawn from", {
  # With K equal to the number of known languages, nothing is drawn, so changing the
  # population must leave the answer alone up to Monte-Carlo error. This is what
  # guards against the known languages being silently redrawn.
  known <- data.frame(language = c("A", "B", "C", "D", "E", "F"),
                      y = c(-0.44, -0.15, -0.27, -0.80, -0.82, -0.07),
                      s = rep(0.1, 6))
  pop1 <- data.frame(language = c("X", "Y"), y = c(-0.9, -0.8), s = c(0.1, 0.1))
  pop2 <- data.frame(language = c("X", "Y"), y = c(0.9, 0.8), s = c(0.1, 0.1))
  a <- run_partly_known_cell(K = 6, known = known, n_reps = 2000L, se_lang = 0.089,
                             pop = pop1, seed = 51L)
  b <- run_partly_known_cell(K = 6, known = known, n_reps = 2000L, se_lang = 0.089,
                             pop = pop2, seed = 52L)
  expect_lt(abs(a$power - b$power), 3 * sqrt(a$mcse^2 + b$mcse^2))
})
