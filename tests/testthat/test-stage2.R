# tests/testthat/test-stage2.R
# Unit tests for the registered Stage 2 analysis (R/14_stage2.R), on constructed data
# whose correct answer is known. The model fits themselves are exercised by the dry
# run (stage2/dry_run/), not here.

library(testthat)
skip_if_not_installed("brms")
suppressPackageStartupMessages(library(dplyr))
source(here::here("R", "02_preprocess_factors.R"))
source(here::here("R", "03_define_priors.R"))
source(here::here("R", "14_stage2.R"))

# One participant's trials: each of `verbs` once per sentence type, then 10
# attention trials, the first `n_right` answered with their target.
make_participant <- function(id, verbs = sprintf("v%02d", 1:10),
                             types = c("Active", "Passive", "Pseudo_Passive"),
                             responses = NULL, n_right = 10, variety = "X") {
  d <- expand.grid(Verb = verbs, S_Type = types, stringsAsFactors = FALSE)
  d$Response <- if (is.null(responses)) rep_len(1:7, nrow(d)) else rep_len(responses, nrow(d))
  d$Agent <- "Man"
  att <- data.frame(Verb = NA, S_Type = NA, Agent = NA, Response = rep(c(1L, 7L), 5))
  att$Response[seq_len(10 - n_right)] <- 4L
  out <- rbind(cbind(d, Trial_Type = "test", Attention_Target = NA),
               cbind(att, Trial_Type = "attention", Attention_Target = rep(c(1L, 7L), 5)))
  out$Trial <- seq_len(nrow(out))
  out$Participant <- id
  out$Variety <- variety
  out
}

test_that("validate_task_data accepts the canonical format and a missing rating", {
  d <- make_participant("p1")
  d$Response[2] <- NA
  expect_silent(v <- validate_task_data(d, "acceptability"))
  expect_type(v$Response, "integer")
})

test_that("validate_task_data stops on malformed input", {
  d <- make_participant("p1")
  expect_error(validate_task_data(d[, setdiff(names(d), "Agent")], "acceptability"), "missing columns")
  bad <- d; bad$Response[1] <- 8L
  expect_error(validate_task_data(bad, "acceptability"), "outside 1-7")
  bad <- d; bad$S_Type[1] <- "Middle"
  expect_error(validate_task_data(bad, "acceptability"), "S_Type")
  bad <- d; bad$Attention_Target[nrow(bad)] <- 4L
  expect_error(validate_task_data(bad, "acceptability"), "target of 1 or 7")
  bad <- rbind(d, d[1, ])
  expect_error(validate_task_data(bad, "acceptability"), "duplicated")
})

test_that("each exclusion criterion is applied at its registered threshold", {
  d <- rbind(
    make_participant("ok"),
    make_participant("missing_rating"),
    make_participant("missing_trial"),
    make_participant("seven_right", n_right = 7),
    make_participant("eight_right", n_right = 8),
    make_participant("flat90", responses = c(rep(7L, 27), 1:3)),
    make_participant("flat87", responses = c(rep(7L, 26), 1:4))
  )
  d$Response[d$Participant == "missing_rating" & d$Trial == 5] <- NA
  d <- d[!(d$Participant == "missing_trial" & d$Trial == 6), ]
  e <- participant_exclusions(validate_task_data(d, "acceptability"), "acceptability",
                              n_attention = 10, attention_min_correct = 8, straightline_share = 0.9)
  ex <- setNames(e$excluded, e$Participant)
  expect_false(ex[["ok"]])
  expect_true(ex[["missing_rating"]])
  expect_true(ex[["missing_trial"]])
  expect_true(ex[["seven_right"]])
  expect_false(ex[["eight_right"]])
  expect_true(ex[["flat90"]])
  expect_false(ex[["flat87"]])
  expect_true(e$straightlined[e$Participant == "flat90"])
  expect_true(e$incomplete[e$Participant == "missing_trial"])
})

test_that("cap_sample keeps the registered number and needs an order to trim", {
  e <- data.frame(Variety = "X", Participant = sprintf("p%d", 1:5), excluded = c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(attr(cap_sample(e, 5), "complete"), character(0))
  four <- cap_sample(e, 4)
  expect_equal(sum(four$analysed), 4)
  expect_equal(attr(four, "complete"), "X")
  expect_error(cap_sample(e, 3), "Session_Order")
  ord <- data.frame(Variety = "X", Participant = sprintf("p%d", 1:5), Session_Order = 5:1)
  three <- cap_sample(e, 3, ord)
  expect_setequal(three$Participant[three$analysed], c("p2", "p4", "p5"))
})

test_that("the affectedness predictor is scaled over verbs, to an SD of 0.5", {
  r <- expand.grid(Variety = "X", Participant = sprintf("r%d", 1:4), Verb = sprintf("v%d", 1:6),
                   Agent = c("Man", "Woman"), stringsAsFactors = FALSE)
  r$Trial_Type <- "test"
  r$Response <- as.integer((match(r$Verb, sprintf("v%d", 1:6)) + seq_len(nrow(r)) %% 2) %% 7 + 1)
  m <- verb_affectedness(r)
  expect_equal(nrow(m), 6)
  expect_equal(m$n_ratings, rep(8L, 6))
  s <- scale_verb_affectedness(m)
  expect_equal(mean(s$Semantics_scaled), 0, tolerance = 1e-12)
  expect_equal(sd(s$Semantics_scaled), 0.5, tolerance = 1e-12)
})

test_that("the outcome-neutral check passes an SD at the threshold and fails one below", {
  m <- data.frame(Variety = rep(c("A", "B"), each = 4),
                  Semantics = c(c(4, 4, 4, 4) + c(-1, 1, -1, 1) * 0.5 * sqrt(3 / 4), 4 + c(-0.1, 0.1, 0, 0)))
  ch <- outcome_neutral_check(m, min_sd = 0.5)
  expect_equal(ch$sd_verb_means[ch$Variety == "A"], 0.5, tolerance = 1e-12)
  expect_true(ch$passes[ch$Variety == "A"])
  expect_false(ch$passes[ch$Variety == "B"])
})

test_that("extract_interactions halves the coefficients and copes without a pseudo-passive", {
  set.seed(1)
  draws <- data.frame(`b_S_TypeActive:Semantics_scaled` = rnorm(4000, -0.8, 0.2),
                      `b_S_TypePseudo_Passive:Semantics_scaled` = rnorm(4000, -0.6, 0.3),
                      check.names = FALSE)
  e <- extract_interactions(draws)
  expect_equal(e$hypothesis, c("H1", "H2"))
  expect_equal(e$per_sd_mean, e$coef_mean / 2)
  expect_equal(e$per_sd_sd, e$coef_sd / 2)
  e1 <- extract_interactions(draws[, 1, drop = FALSE])
  expect_equal(e1$hypothesis, "H1")
  expect_error(extract_interactions(draws[, 2, drop = FALSE]), "no b_S_TypeActive")
})

test_that("varieties are averaged within a language: mean of means, mean of SDs", {
  tab <- data.frame(variety = c("A1", "A2", "B"), language = c("A", "A", "B"),
                    family = c("F", "F", "G"))
  est <- data.frame(Variety = c("A1", "A2", "B", "A1", "B"),
                    hypothesis = c("H1", "H1", "H1", "H2", "H2"),
                    per_sd_mean = c(-0.4, -0.2, -0.1, -0.5, 0.1),
                    per_sd_sd = c(0.10, 0.14, 0.2, 0.12, 0.3))
  a <- average_varieties(est, tab)
  h1a <- a[a$language == "A" & a$hypothesis == "H1", ]
  expect_equal(h1a$y, -0.3)
  expect_equal(h1a$se, 0.12)
  expect_equal(h1a$n_varieties, 2L)
  h2a <- a[a$language == "A" & a$hypothesis == "H2", ]
  expect_equal(h2a$n_varieties, 1L)   # only A1 used a pseudo-passive
  expect_equal(h2a$y, -0.5)
  expect_error(average_varieties(rbind(est, data.frame(Variety = "C", hypothesis = "H1",
                                                       per_sd_mean = 0, per_sd_sd = 1)), tab),
               "missing from the language table")
})

test_that("Bayes factors are read against the registered threshold of 10", {
  expect_equal(bf_from_p_negative(10 / 11), 10)
  expect_equal(interpret_bf(c(10, 9.99, 0.1, 0.1001, 1)),
               c("for", "inconclusive", "against", "inconclusive", "inconclusive"))
  b <- bf_from_draws(rep(-1, 100))
  expect_equal(b$bf, Inf)
  expect_equal(b$bound, "BF > 99")
  expect_equal(bf_from_draws(c(-1, -1, -1, 1))$bf, 3)
})

test_that("the meta-analytic priors are the per-language priors halved: 0.25 and 0.3", {
  expect_equal(meta_prior_sd(), c(H1 = 0.25, H2 = 0.3))
})

test_that("the Stage 2 settings parse and match the Stage 1 report", {
  cfg <- stage2_config(here::here("config", "analysis_config.yaml"))
  expect_equal(cfg$sampler$chains, 16L)
  expect_equal(cfg$sampler$iter, 5000L)
  expect_equal(cfg$sampler$warmup, 2000L)
  expect_equal(cfg$sampler$adapt_delta, 0.99)
  expect_equal(cfg$sampler$max_treedepth, 12L)
  expect_equal(cfg$samples$acceptability, 80L)
  expect_equal(cfg$samples$ratings, 20L)
  expect_equal(cfg$exclusions$attention_min_correct, 8L)
  expect_equal(cfg$exclusions$straightline_share, 0.9)
  expect_equal(cfg$inference$bf_threshold, 10L)
  expect_equal(cfg$inference$tau_prior_halfnormal_scale, 0.5)
  tab <- read_language_table(here::here(cfg$languages_table))
  expect_true(all(c("variety", "language", "family") %in% names(tab)))
  expect_equal(length(unique(tab$language)), 40L)
  expect_equal(sum(tapply(tab$family, tab$language, `[`, 1) == "Indo-European"), 24L)
})

test_that("the single-model ladder starts at the registered maximal model", {
  lad <- build_single_model_ladder()
  expect_equal(names(lad)[1], "R0_maximal")
  f <- function(x) gsub("[[:space:]]+", "", paste(deparse(x$formula), collapse = ""))
  expect_equal(f(lad$R0_maximal),
               gsub(" ", "", paste("Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled | Language)",
                                   "+ (1 + S_Type * Semantics_scaled | Participant) + (1 + S_Type | Verb)")))
  expect_equal(f(lad$R3_sentence_type_slopes),
               gsub(" ", "", paste("Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled || Language)",
                                   "+ (1 + S_Type || Participant) + (1 + S_Type || Verb)")))
  # No level carries a by-verb affectedness slope.
  expect_false(any(vapply(lad, function(x) grepl("Semantics_scaled[^|]*\\|+Verb", f(x)), logical(1))))
})

test_that("pooled data nest participants and verbs within variety", {
  tab <- data.frame(variety = c("A1", "A2"), language = c("A", "A"), family = "F")
  mk <- function(v) {
    d <- data.frame(Variety = v, Participant = factor(c("p1", "p1")), Verb = factor(c("eat", "eat")),
                    Agent = "Man", S_Type = factor(c("Passive", "Active")), Response = c(5L, 6L),
                    Semantics = 4, Semantics_scaled = c(0.1, 0.1))
    d$S_Type <- forcats::fct_relevel(d$S_Type, "Passive")
    d
  }
  p <- build_pooled_data(list(mk("A1"), mk("A2")), tab)
  expect_equal(nlevels(p$Participant), 2L)
  expect_equal(nlevels(p$Verb), 2L)
  expect_equal(levels(p$Language), "A")
  expect_equal(levels(p$S_Type)[1], "Passive")
})

test_that("run_meta_test agrees with the design analysis engine", {
  skip_if_not_installed("bayesmeta")
  skip_if(Sys.getenv("CLAPS_SLOW_TESTS") != "true", "set CLAPS_SLOW_TESTS=true to run (about 35 s per data set)")
  source(here::here("R", "13_simulate_klanguage_pooled.R"))
  set.seed(2026)
  for (i in 1:5) {
    s <- runif(20, 0.06, 0.2)
    y <- rnorm(20, rnorm(1, -0.25, 0.15), 0.39) + rnorm(20, 0, s)
    m <- run_meta_test(y, s, prior_mu_sd = 0.25, tau_scale = 0.5)
    e <- ra_meta_batch(matrix(y, 1), s, prior_mu_sd = 0.25, prior_tau_scale = 0.5)
    expect_lt(abs(m$p_mu_negative - e$p_mu_negative), 1e-3)
  }
})
