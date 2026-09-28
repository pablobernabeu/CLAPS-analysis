# tests/testthat/test-klanguage-pooled.R
#
# Tests for the K-language pooled design engine (R/13_simulate_klanguage_pooled.R).
#
# The engine's job is to stand in for a fit whose median runtime is 34 hours at
# N = 50 and 170 hours at N = 150 (outputs/design_pooled_v2), so the things worth
# testing are the ones that would let it stand in WRONGLY: the quadrature not
# integrating to what it claims, the Bayes factor not meaning what
# R/05_hypothesis_tests.R means by it, the shared-verb component not producing the
# correlation it is supposed to produce, and a known language being treated as
# something it is not. Where a closed form or an exact combinatorial answer exists
# the test uses it. The rest either compare two routes that must agree exactly at a
# common seed, or are Monte Carlo comparisons at a stated tolerance, usually three
# Monte Carlo standard errors.
#
# A few tests pin the reported numbers themselves. They need the gitignored pilot
# posteriors or brms replicates, so they skip on a fresh clone and in CI and run
# wherever those files are present.

source(testthat::test_path("..", "..", "R", "13_simulate_klanguage_pooled.R"))

da_path <- function(...) testthat::test_path("..", "..", ...)
have_pilot <- function() {
  all(file.exists(da_path("outputs", "pilot_models",
                          paste0("pilot_dgp_v2_pilot_", c("English", "Turkish", "Norwegian"), ".rds"))))
}
harmonised_csv <- da_path("outputs", "design_summary_pilot", "cross_language_effects_harmonised.csv")

# The engine's inputs exactly as scripts/run_klanguage_design_analysis.R builds them.
load_design_inputs <- function() {
  H  <- read.csv(harmonised_csv)
  PP <- read.csv(da_path("outputs", "design_summary_pilot", "pilot_params_ceilings.csv"))
  langs <- c("English", "Turkish", "Norwegian")
  # Read from the pilot table, not restated: a literal here would go on passing after
  # the table changed (audit item klanguage-engine#12). The run script writes the
  # values out to seven digits and checks them against this column at 1e-6. The
  # difference, about 3e-8, moves none of the pinned powers below; the one
  # continuous pin, the calibration's mean_z, moves by about 1e-8 and is compared
  # at a tolerance that allows for it.
  sdx <- PP$affectedness_sd[match(langs, PP$language)]
  stopifnot(!anyNA(sdx))
  raw <- read_claps_interaction_draws(da_path("outputs", "pilot_models"), langs)
  claps_se <- setNames(vapply(seq_along(langs), function(i) sd(raw[[i]]) * sdx[i], numeric(1)), langs)
  se_verb <- mean(PP$se_floor_h1b * PP$affectedness_sd)
  se_ppt <- sqrt(max(mean(PP$slope_post_sd * PP$affectedness_sd)^2 - se_verb^2, 0))
  list(hc = prepare_language_table(H, claps_se), sp = split_se(se_verb, se_ppt, 0.164),
       pilot_draws = lapply(seq_along(langs), function(i) raw[[i]] * sdx[i]), PP = PP)
}

# ---------------------------------------------------------------------------
# The meta-analytic posterior
# ---------------------------------------------------------------------------

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

test_that("the grid posterior matches direct numerical integration at a real tau prior", {
  # The tests above pin tau; this one checks the half-normal prior term, the
  # trapezoidal weights and the mixture over tau > 0 against integrate(), on the
  # three CLAPS pilot languages at the design prior. The reference integrates over
  # the grid's own range, [0, 3], so the two differ only by quadrature error.
  y <- c(-0.4439, -0.1453, -0.2678); s <- c(0.12, 0.0888, 0.0645); scale <- 0.5
  marg <- function(tau) {
    vapply(tau, function(t) {
      v <- s^2 + t^2; w <- 1 / v; m <- sum(w * y) / sum(w)
      exp(-0.5 * sum(log(v)) - 0.5 * log(sum(w)) - 0.5 * sum(w * (y - m)^2)) *
        dnorm(t, 0, scale)
    }, numeric(1))
  }
  p_neg_given <- function(tau) {
    vapply(tau, function(t) {
      w <- 1 / (s^2 + t^2)
      pnorm(0, sum(w * y) / sum(w), sqrt(1 / sum(w)))
    }, numeric(1))
  }
  z <- integrate(marg, 0, 3, rel.tol = 1e-10)$value
  ref_p <- integrate(function(t) marg(t) * p_neg_given(t), 0, 3, rel.tol = 1e-10)$value / z
  cdf <- function(q) integrate(marg, 0, q, rel.tol = 1e-10)$value / z
  ref_median <- uniroot(function(q) cdf(q) - 0.5, c(1e-6, 3), tol = 1e-8)$root

  post <- ra_meta_posterior(y, s, prior_tau_scale = scale)
  expect_equal(post$p_mu_negative, ref_p, tolerance = 1e-4)
  # The median is read off the nearest grid point, 0.005 apart.
  expect_lt(abs(post$tau_median - ref_median), 0.005)
})

test_that("the batch posterior returns exactly what the one-at-a-time posterior does", {
  # Every production cell now scores its replicates in one batch. If the batch
  # arithmetic drifted from ra_meta_posterior(), every reported number would move
  # without any test noticing, so the comparison is at machine precision, with
  # unequal SEs, a proper prior on mu and a one-point grid.
  set.seed(8)
  s <- c(0.12, 0.09, 0.06, 0.13, 0.20)
  Y <- matrix(rnorm(300 * 5, -0.3, 0.4), 300, 5)
  one <- t(vapply(seq_len(nrow(Y)), function(i) {
    p <- ra_meta_posterior(Y[i, ], s)
    c(p$p_mu_negative, p$tau_median)
  }, numeric(2)))
  b <- ra_meta_batch(Y, s, want_tau = TRUE, chunk = 70L)
  expect_equal(b$p_mu_negative, one[, 1], tolerance = 1e-12)
  expect_equal(b$tau_median, one[, 2])

  pm <- vapply(1:50, function(i) ra_meta_posterior(Y[i, ], s, prior_mu_sd = 0.4)$p_mu_negative,
               numeric(1))
  expect_equal(ra_meta_batch(Y[1:50, ], s, prior_mu_sd = 0.4)$p_mu_negative, pm, tolerance = 1e-12)

  p0 <- vapply(1:20, function(i) ra_meta_posterior(Y[i, ], s, tau_grid = 0)$p_mu_negative, numeric(1))
  expect_equal(ra_meta_batch(Y[1:20, ], s, tau_grid = 0)$p_mu_negative, p0, tolerance = 1e-12)
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
  # The vectorised form agrees element by element, clamp included.
  ps <- c(0, 1e-15, 0.3, 10 / 11, 1 - 1e-15, 1)
  expect_identical(bf_from_p_negative(ps),
                   vapply(ps, function(p) directional_bf_meta(list(p_mu_negative = p)), numeric(1)))
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
  # split_se() as the run script calls it: a common draw of SD sqrt(rho) * se_verb
  # plus independent draws. With the participant term at zero, the correlation
  # between two languages' errors must be rho and the total must be unchanged.
  se_verb <- 0.074; rho <- 0.4
  sp <- split_se(se_verb, se_ppt = 0, rho = rho)
  expect_equal(sp$shared^2 + sp$indep^2, se_verb^2, tolerance = 1e-12)
  set.seed(11)
  n <- 200000L
  common <- rnorm(n, 0, sp$shared)
  e1 <- common + rnorm(n, 0, sp$indep)
  e2 <- common + rnorm(n, 0, sp$indep)
  expect_equal(cor(e1, e2), rho, tolerance = 0.01)
  expect_equal(sd(e1), se_verb, tolerance = 0.01)
  # The participant term is never shared, so it dilutes the correlation.
  sp2 <- split_se(se_verb, se_ppt = 0.05, rho = rho)
  expect_equal(sp2$shared, sp$shared)
  expect_equal(sp2$indep^2, sp$indep^2 + 0.05^2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Partly known sets
# ---------------------------------------------------------------------------

pop4 <- data.frame(language = c("A", "B", "C", "D"), source = "Glossa2023",
                   y = c(-0.44, -0.15, -0.27, 0.33), s = c(0.12, 0.09, 0.06, 0.13))
claps_like <- data.frame(language = c("En", "Tr", "No"), source = "CLAPS_pilot",
                         y = c(-0.44, -0.40, -0.36), s = c(0.12, 0.09, 0.06))

test_that("run_partly_known_cell reduces to the exchangeable arm when nothing is known", {
  a <- run_partly_known_cell(K = 6, known = pop4[0, ], n_reps = 1500L, se_lang = 0.089,
                             pop = pop4, seed = 21L)
  b <- run_klanguage_cell(K = 6, n_reps = 1500L, se_lang = 0.089, mode = "assurance",
                          pop = pop4, seed = 21L)
  expect_equal(a$n_known, 0L)
  # Without a shared verb offset the two arms draw the same numbers in the same order,
  # so they agree to the replicate; the tolerance only allows for a BLAS that sums a
  # matrix product in a different order from a vector one.
  expect_lt(abs(a$power - b$power), 3 / 1500)
})

test_that("knowing a favourable language's effect in advance raises power", {
  # Two sets differing only in whether three strongly-negative languages are known
  # in advance or drawn from a population with the same mean. Fixing them removes
  # tau's contribution for those places, so power must rise.
  none <- run_partly_known_cell(K = 8, known = pop4[0, ], n_reps = 2000L,
                                se_lang = 0.089, pop = pop4, seed = 31L)
  some <- run_partly_known_cell(K = 8, known = claps_like, n_reps = 2000L,
                                se_lang = 0.089, pop = pop4, seed = 32L)
  expect_equal(some$n_known, 3L)
  expect_gt(some$power, none$power)
})

test_that("doubt about protocol transfer removes the benefit of knowing a language", {
  # transfer_sd widens the prior on an already-measured language. Pushed well above
  # tau it should return that language to being effectively an unknown draw, so the
  # power advantage over the nothing-known case must shrink. All three synthetic
  # languages here stand for ones measured under another protocol, so each gets it.
  other <- transform(claps_like, source = "Glossa2023")
  tight <- run_partly_known_cell(K = 8, known = other, n_reps = 2000L, se_lang = 0.089,
                                 pop = pop4, transfer_sd = transfer_sd_by_source(other, 0), seed = 41L)
  loose <- run_partly_known_cell(K = 8, known = other, n_reps = 2000L, se_lang = 0.089,
                                 pop = pop4, transfer_sd = transfer_sd_by_source(other, 1.5), seed = 42L)
  expect_gt(tight$power, loose$power)
  expect_equal(loose$transfer_sd, 1.5)
  expect_equal(loose$n_widened, 3L)
})

test_that("a CLAPS pilot language is never widened for protocol transfer", {
  # The defect behind partly-known rows c to e: a single transfer SD was recycled
  # over every known language, CLAPS pilot ones included. A CLAPS-only set must give
  # the same answer, to the replicate, whatever transfer doubt is set for the others.
  runs <- lapply(c(0, 0.15, 0.40), function(t) {
    run_partly_known_cell(K = 6, known = claps_like, n_reps = 800L, se_lang = 0.089,
                          pop = pop4, transfer_sd = transfer_sd_by_source(claps_like, t),
                          seed = 5226L)
  })
  expect_identical(runs[[1]]$power, runs[[2]]$power)
  expect_identical(runs[[1]]$power, runs[[3]]$power)
  expect_identical(runs[[3]]$n_widened, 0L)

  # In a mixed set only the non-CLAPS rows are widened.
  mixed <- rbind(claps_like, pop4[c(1, 4), ])
  expect_equal(transfer_sd_by_source(mixed, 0.4), c(0, 0, 0, 0.4, 0.4))
  r <- run_partly_known_cell(K = 8, known = mixed, n_reps = 400L, se_lang = 0.089, pop = pop4,
                             transfer_sd = transfer_sd_by_source(mixed, 0.4), seed = 3L)
  r_explicit <- run_partly_known_cell(K = 8, known = mixed, n_reps = 400L, se_lang = 0.089,
                                      pop = pop4, transfer_sd = c(0, 0, 0, 0.4, 0.4), seed = 3L)
  expect_identical(r$power, r_explicit$power)
})

test_that("a per-row transfer SD returns one row, and a recycled one is refused", {
  mixed <- rbind(claps_like, pop4[c(1, 4), ])
  r <- run_partly_known_cell(K = 7, known = mixed, n_reps = 200L, se_lang = 0.089, pop = pop4,
                             transfer_sd = transfer_sd_by_source(mixed, 0.15), seed = 2L)
  # data.frame() would recycle a length-5 column into five identical rows.
  expect_equal(nrow(r), 1L)
  expect_equal(r$transfer_sd, 0.15)
  expect_equal(r$n_widened, 2L)
  # A single non-zero value would widen the CLAPS rows too, so it is an error.
  expect_error(run_partly_known_cell(K = 7, known = mixed, n_reps = 20L, se_lang = 0.089,
                                     pop = pop4, transfer_sd = 0.15, seed = 2L),
               "one value per known language")
  expect_error(run_partly_known_cell(K = 7, known = mixed, n_reps = 20L, se_lang = 0.089,
                                     pop = pop4, transfer_sd = c(0, 0.15), seed = 2L),
               "one value per known language")
  # A single zero means no widening anywhere, which is unambiguous.
  expect_equal(run_partly_known_cell(K = 7, known = mixed, n_reps = 20L, se_lang = 0.089,
                                     pop = pop4, transfer_sd = 0, seed = 2L)$n_widened, 0L)
  expect_error(transfer_sd_by_source(mixed[, c("y", "s")], 0.15), "source column")
})

test_that("a lone CLAPS pilot language is refused a transfer doubt too", {
  # With one known language a single value is already one per row, so the length
  # check cannot tell it from a recycled scalar. The rows' source has to.
  en <- claps_like[1, ]
  expect_error(run_partly_known_cell(K = 6, known = en, n_reps = 20L, se_lang = 0.089,
                                     pop = pop4, transfer_sd = 0.40, seed = 9L),
               "CLAPS pilot language")
  # The same holds for a per-row vector that puts a doubt on a CLAPS row by hand.
  mixed <- rbind(claps_like, pop4[c(1, 4), ])
  expect_error(run_partly_known_cell(K = 7, known = mixed, n_reps = 20L, se_lang = 0.089,
                                     pop = pop4, transfer_sd = c(0.4, 0, 0, 0.4, 0.4), seed = 9L),
               "CLAPS pilot language")
  # Built by source, the lone CLAPS row is left alone and a lone 2023 row is widened.
  expect_identical(run_partly_known_cell(K = 6, known = en, n_reps = 20L, se_lang = 0.089,
                                         pop = pop4, transfer_sd = transfer_sd_by_source(en, 0.40),
                                         seed = 9L)$n_widened, 0L)
  g1 <- pop4[1, ]
  expect_identical(run_partly_known_cell(K = 6, known = g1, n_reps = 20L, se_lang = 0.089,
                                         pop = pop4, transfer_sd = 0.40, seed = 9L)$n_widened, 1L)
})

test_that("a fully known set does not consult the population, in either mode", {
  # With K equal to the number of known languages nothing is drawn, so the population
  # and the mode must leave the answer exactly unchanged at a common seed. This is
  # what guards against the known languages being silently redrawn, and it is why
  # the PI's construction at K = 7 is the same under every tau regime.
  known <- data.frame(language = c("A", "B", "C", "D", "E", "F"), source = "CLAPS_pilot",
                      y = c(-0.44, -0.15, -0.27, -0.80, -0.82, -0.07), s = rep(0.1, 6))
  pop1 <- data.frame(language = c("X", "Y"), y = c(-0.9, -0.8), s = c(0.1, 0.1))
  pop2 <- data.frame(language = c("X", "Y"), y = c(0.9, 0.8), s = c(0.1, 0.1))
  a <- run_partly_known_cell(K = 6, known = known, n_reps = 1000L, se_lang = 0.089,
                             pop = pop1, seed = 51L)
  b <- run_partly_known_cell(K = 6, known = known, n_reps = 1000L, se_lang = 0.089,
                             pop = pop2, seed = 51L)
  c_pt <- run_partly_known_cell(K = 6, known = known, n_reps = 1000L, se_lang = 0.089,
                                mode = "point", mu_point = 0.5, tau_point = 2, seed = 51L)
  expect_identical(a$power, b$power)
  expect_identical(a$power, c_pt$power)
  # And it needs no population at all.
  expect_identical(run_partly_known_cell(K = 6, known = known, n_reps = 1000L, se_lang = 0.089,
                                         seed = 51L)$power, a$power)
})

test_that("known effects can be held at their estimates, and are then not widened", {
  known <- data.frame(language = LETTERS[1:5], source = "Glossa2023", y = rep(-0.2, 5), s = rep(0.3, 5))
  # Held exactly at five identical negative estimates and measured almost without
  # error, the set leaves nothing for tau to explain, so the pooled test must succeed
  # every time. Drawn from priors of SD 0.3 the five effects spread out, and it cannot.
  fixed <- run_partly_known_cell(K = 5, known = known, n_reps = 300L, se_lang = 0.01,
                                 known_truth = "fixed", seed = 4L)
  drawn <- run_partly_known_cell(K = 5, known = known, n_reps = 300L, se_lang = 0.01,
                                 known_truth = "drawn", seed = 4L)
  expect_equal(fixed$power, 1)
  expect_lt(drawn$power, 1)
  expect_equal(fixed$known_truth, "fixed")
  expect_error(run_partly_known_cell(K = 5, known = known, n_reps = 10L, se_lang = 0.1,
                                     known_truth = "fixed",
                                     transfer_sd = transfer_sd_by_source(known, 0.15)),
               "transfer doubt")
})

test_that("point mode draws the unknown places at the stated (mu, tau)", {
  # With tau at zero and a mean far below zero, every drawn place carries the same
  # strongly negative effect, so the pooled test must succeed every time.
  r <- run_partly_known_cell(K = 8, known = claps_like, n_reps = 300L, se_lang = 0.05,
                             mode = "point", mu_point = -1, tau_point = 0, seed = 9L)
  expect_equal(r$power, 1)
  expect_equal(r$mu_point, -1)
  expect_error(run_partly_known_cell(K = 8, known = claps_like, n_reps = 10L, se_lang = 0.05,
                                     mode = "point"))
})

# ---------------------------------------------------------------------------
# Where the outcome variance comes from
# ---------------------------------------------------------------------------

test_that("the variance decomposition adds up and runs under the design prior", {
  pop <- data.frame(y = c(-0.44, -0.15, -0.27, -0.10, -0.30), s = c(0.12, 0.09, 0.06, 0.10, 0.08))
  d <- decompose_outcome_variance(6, pop, se_indep = 0.083, se_shared = 0.03,
                                  n_par = 60L, n_sets = 5L, n_reps = 40L, n_boot = 20L, seed = 3L)
  expect_equal(nrow(d), 1L)
  shares <- c(d$share_population_params, d$share_composition, d$share_noise)
  expect_true(all(shares >= 0))
  # The three unbiased components sum to the total variance up to a term of order
  # 1 / n_par, and the clipping at zero only ever adds.
  expect_lt(abs(sum(shares) - 1), 0.08)
  expect_equal(d$share_true_effects, d$share_population_params + d$share_composition)

  # The population is drawn under the DESIGN prior unless told otherwise: the
  # default reproduces an explicit design-prior call, and the analysis prior differs.
  expect_identical(formals(decompose_outcome_variance)$design_tau_prior_scale,
                   quote(DESIGN_TAU_PRIOR_SCALE))
  small <- function(...) decompose_outcome_variance(6, pop, 0.083, 0.03, n_par = 20L, n_sets = 3L,
                                                    n_reps = 20L, n_boot = 5L, seed = 5L, ...)
  expect_identical(small(), small(design_tau_prior_scale = DESIGN_TAU_PRIOR_SCALE))
  cols <- c("mean_power", "share_population_params", "share_composition")
  expect_false(identical(small()[cols], small(design_tau_prior_scale = 1.1)[cols]))
})

test_that("a population with nothing left to learn has no population share", {
  # Seven identical, precisely measured languages pin mu and tau (at zero), so
  # neither the population nor the composition of a set can vary the outcome.
  pop <- data.frame(y = rep(-0.07, 7), s = rep(0.001, 7))
  d <- decompose_outcome_variance(6, pop, se_indep = 0.083, se_shared = 0.03,
                                  n_par = 40L, n_sets = 4L, n_reps = 60L, n_boot = 10L, seed = 7L)
  expect_gt(d$mean_power, 0.05)
  expect_lt(d$mean_power, 0.95)
  expect_lt(d$share_population_params, 0.03)
  expect_lt(d$share_composition, 0.03)
  expect_gt(d$share_noise, 0.9)
})

# ---------------------------------------------------------------------------
# The observed languages and the passive-reading fork
# ---------------------------------------------------------------------------

toy_table <- function() {
  data.frame(
    source   = c("CLAPS_pilot", "CLAPS_pilot", "Glossa2023", "Glossa2023", "Glossa2023", "Glossa2023", "Glossa2023"),
    language = c("English", "Turkish", "Balinese", "English", "Hebrew", "Indonesian", "Mandarin"),
    se_kind  = c("slope_sd", "slope_sd", rep("posterior_sd", 5)),
    raw_est  = c(-0.88, -0.29, 0.33, -0.54, -0.81, -0.07, -0.80),
    raw_se   = c(0.19, 0.16, 0.13, 0.18, 0.20, 0.08, 0.21),
    sd_x     = c(0.5035, 0.5035, 1.0078, 1.0069, 1.0090, 1.0070, 1.0089),
    y_shared = 0, se_shared = 0, stringsAsFactors = FALSE)
}
with_y <- function(h) transform(h, y = raw_est * sd_x, se = raw_se * sd_x)

test_that("one_per_language prefers the CLAPS estimate and keeps the order", {
  h <- with_y(toy_table())
  o <- one_per_language(h)
  expect_equal(o$language, c("English", "Turkish", "Balinese", "Hebrew", "Indonesian", "Mandarin"))
  expect_equal(o$source[o$language == "English"], "CLAPS_pilot")
})

test_that("the alternative reading changes exactly the three forked languages", {
  h <- with_y(toy_table())
  expect_identical(apply_passive_reading(h, "canonical"), h)
  a <- apply_passive_reading(h, "alternative", "independent")
  changed <- which(a$y != h$y)
  expect_equal(sort(a$language[changed]), c("Balinese", "Indonesian", "Mandarin"))
  expect_true(all(a$source[changed] == "Glossa2023"))
  # The contrast is Active minus the other construction, on each language's own scale.
  ix <- match(c("Balinese", "Indonesian", "Mandarin"), a$language)
  expect_equal(a$y[ix], c(0.86, 0.27, 0.29) * h$sd_x[ix], tolerance = 1e-12)
  expect_equal(a$se[ix], sqrt(c(0.13, 0.08, 0.21)^2 + c(0.13, 0.10, 0.18)^2) * h$sd_x[ix],
               tolerance = 1e-12)
  expect_equal(apply_passive_reading(h, "alternative", "canonical")$se, h$se)
  # A contrast table that no longer matches the canonical estimates is refused.
  h2 <- h
  h2$raw_est[h2$language == "Mandarin"] <- -0.75
  expect_error(apply_passive_reading(with_y(h2), "alternative"), "no longer matches")
})

test_that("the populations keep their source and drop Glossa English", {
  p <- build_populations(with_y(toy_table()))
  expect_equal(names(p), c("claps_protocol", "seven", "glossa_protocol", "drop_balinese", "drop_hebrew"))
  expect_equal(p$glossa_protocol$language, c("Balinese", "Hebrew", "Indonesian", "Mandarin"))
  expect_true(all(c("language", "source", "y", "s") %in% names(p$seven)))
  expect_false("Balinese" %in% p$drop_balinese$language)
})

test_that("on the real inputs, two of the seven change sign under the alternative reading", {
  skip_if_not(file.exists(harmonised_csv))
  H <- read.csv(harmonised_csv)
  # CLAPS SEs are irrelevant to the effects themselves, so the table's own stand in.
  hc <- prepare_language_table(H, claps_se = c(English = 0.12, Turkish = 0.089, Norwegian = 0.064))
  canon <- build_populations(hc, "canonical")$seven
  alt <- build_populations(hc, "alternative")$seven
  expect_equal(mean(canon$y), -0.317, tolerance = 0.001)
  expect_equal(mean(alt$y), -0.035, tolerance = 0.001)
  flipped <- canon$language[sign(canon$y) != sign(alt$y)]
  expect_equal(sort(flipped), c("Indonesian", "Mandarin"))
  expect_gt(alt$y[alt$language == "Balinese"], 0)
  expect_gt(canon$y[canon$language == "Balinese"], 0)

  # The same mean in closed form. The alternative reading replaces each forked
  # language's Active coefficient b by b minus the other construction's, so on the
  # per-SD scale the seven-language mean moves by minus the sum of b_other * sd_x,
  # over seven. The coefficients are published to two decimals, each good to 0.005,
  # which leaves the mean good to about 0.002, so the third decimal is as far as it
  # can honestly be pinned.
  alt_tab <- ALT_PASSIVE_CONTRASTS
  g <- hc[hc$source == "Glossa2023", ]
  sdx_alt <- g$sd_x[match(alt_tab$language, g$language)]
  closed <- mean(canon$y) - sum(alt_tab$b_other * sdx_alt) / nrow(canon)
  expect_equal(nrow(canon), 7L)
  expect_equal(mean(alt$y), closed, tolerance = 1e-12)
  expect_equal(round(closed, 3), -0.035)
})

# ---------------------------------------------------------------------------
# The reported numbers, where their inputs are present
# ---------------------------------------------------------------------------

test_that("the corrected partly-known rows reproduce the audit's recomputation", {
  skip_if_not(have_pilot())
  inp <- load_design_inputs()
  pops <- build_populations(inp$hc)
  claps3 <- pops$claps_protocol
  g4 <- pops$glossa_protocol
  cell <- function(known, tsd, j, K) {
    run_partly_known_cell(K = K, known = known, n_reps = 4000L, se_lang = inp$sp$indep,
                          pop = pops$seven, shared_verb_sd = inp$sp$shared,
                          transfer_sd = transfer_sd_by_source(known, tsd), seed = 5200L + 10L * j + K)$power
  }
  # Row c, Hebrew and Mandarin known at transfer 0.15, and row e, Balinese and
  # Indonesian at 0.40, at six languages, with the run script's seeds.
  expect_equal(cell(rbind(claps3, g4[g4$language %in% c("Hebrew", "Mandarin"), ]), 0.15, 3L, 6), 0.93275)
  expect_equal(cell(rbind(claps3, g4[match(c("Balinese", "Indonesian"), g4$language), ]), 0.40, 5L, 6),
               0.26425)
})

test_that("the PI's construction reproduces the scratch cross-check of 23 September", {
  skip_if_not(have_pilot())
  inp <- load_design_inputs()
  pops <- build_populations(inp$hc)
  k7 <- pops$seven
  run <- function(K) {
    run_partly_known_cell(K = K, known = k7, n_reps = 4000L, se_lang = inp$sp$indep,
                          pop = pops$seven, shared_verb_sd = inp$sp$shared,
                          transfer_sd = transfer_sd_by_source(k7, 0.15), seed = 7000L + K)$power
  }
  expect_equal(round(run(10), 3), 0.831)
  expect_equal(round(run(20), 3), 0.889)
})

test_that("the PI's construction under the alternative reading uses that reading throughout", {
  # The alternative reading changes both the known Balinese, Indonesian and Mandarin
  # effects and the population the new places are drawn from. A driver that paired
  # the alternative known set with the canonical population, or the reverse, would
  # still run, so one cell is pinned: the seven regime, drawn with 0.15 transfer
  # doubt, assurance, K = 20, seed 7020, as in klanguage_known7.csv.
  skip_if_not(have_pilot())
  inp <- load_design_inputs()
  pops <- build_populations(inp$hc, "alternative", "independent")
  k7 <- pops$seven
  expect_equal(round(mean(k7$y), 3), -0.035)
  r <- run_partly_known_cell(K = 20, known = k7, n_reps = 4000L, se_lang = inp$sp$indep,
                             pop = pops$seven, shared_verb_sd = inp$sp$shared,
                             transfer_sd = transfer_sd_by_source(k7, 0.15), seed = 7020L)
  expect_equal(r$n_widened, 4L)
  expect_equal(r$power, 0.2330)
})

test_that("the calibration reproduces the reported engine rows and computes its target", {
  skip_if_not(have_pilot())
  inp <- load_design_inputs()
  ix <- match(c("English", "Turkish", "Norwegian"), inp$PP$language)
  s_floor <- inp$PP$se_floor_h1b[ix] * inp$PP$affectedness_sd[ix]
  r <- calibrate_against_pooled_k3(inp$pilot_draws, s_floor, n_reps = 4000L, prior_tau_scale = 1.1)
  expect_equal(r$p_bf10, 0.58225)
  expect_equal(r$p_bf6, 0.86375)
  # The reported value was computed with the run script's seven-digit per-SD
  # scaling, and load_design_inputs() reads the full-precision column, about 3e-8
  # away. The pilot draws scale with it, so mean_z moves by about 1e-8; 1e-7 still
  # catches any real change.
  expect_equal(r$mean_z, 1.40776519516873, tolerance = 1e-7)

  pooled_dir <- da_path("outputs", "design_pooled_v2")
  skip_if_not(dir.exists(pooled_dir))
  tg <- pooled_k3_target(read_pooled_k3_replicates(pooled_dir))
  expect_equal(tg$n_success, 135L)
  expect_equal(unname(round(tg$target, 4)), c(1.3865, 0.3113, 0.5778, 0.8519, 0.9778))
})
