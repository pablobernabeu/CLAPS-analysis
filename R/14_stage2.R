# R/14_stage2.R
#
# Purpose
#   The confirmatory (Stage 2) analysis registered in the Stage 1 Registered Report,
#   as functions. The numbered scripts in stage2/ call them in the registered order,
#   and tests/testthat/test-stage2.R exercises each one on constructed data.
#
# What the Stage 1 report fixes, and where it is implemented here
#   Exclusions (Exclusion criteria)           participant_exclusions()
#   Sample of 80 + 20 per language (Sampling)  cap_sample()
#   Outcome-neutral check (Exclusion criteria) outcome_neutral_check()
#   Affectedness predictor (Pre-processing)    verb_affectedness(), scale_verb_affectedness()
#   Individual-language models and ladder      fit_variety_model(), via R/09_model_ladder.R
#   Per-SD interaction estimates               extract_interactions()
#   Varieties averaged within a language       average_varieties()
#   Meta-analysis and directional BF           run_meta_test(), bf_from_p_negative()
#   Single model across languages (robustness) build_single_model_ladder(), build_pooled_data()
#   Families (robustness)                      fit_family_meta()
#
# Scale of the estimates
#   Semantics_scaled is centred and divided by two SDs of the VERBS' mean ratings,
#   so a unit is two SDs of affectedness and the verb-level SD of the predictor is
#   0.5 exactly. A model coefficient is therefore converted to an effect per SD of
#   affectedness by halving it, which is what the report states. The pilot models
#   scaled over trial rows instead, which in a balanced design gives the same centre
#   and an SD within a fraction of a percent of it.
#
# Inputs
#   Two long CSV files per variety, in the canonical format checked by
#   validate_task_data(): one for the acceptability task and one for the semantic
#   ratings. stage2/README.md describes the columns. The variety-to-language table
#   is config/stage2_languages.csv.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

#' Read the registered Stage 2 settings from config/analysis_config.yaml.
#'
#' @param path Path to the analysis configuration.
#' @return The `stage2` block as a list, with the meta-analytic prior SDs added
#'   (`meta_prior_sd`), derived from the per-language interaction priors so that
#'   the two cannot drift apart.
stage2_config <- function(path = "config/analysis_config.yaml") {
  cfg <- yaml::read_yaml(path)$stage2
  if (is.null(cfg)) stop("[stage2] no stage2 block in ", path)
  needed <- c("sampler", "exclusions", "samples", "outcome_neutral", "inference",
              "languages_table", "seed_base")
  missing <- setdiff(needed, names(cfg))
  if (length(missing)) stop("[stage2] stage2 block lacks: ", paste(missing, collapse = ", "))
  cfg$meta_prior_sd <- meta_prior_sd()
  cfg
}

#' Settings for a real run or for the dry run on simulated data.
#'
#' The dry run exercises every script end to end on data from
#' stage2/dry_run/simulate_data.R. It is the only route by which anything below
#' differs from the registered values: smaller samples and fewer verbs, so that the
#' simulated data stay small, and a shorter sampler, so that it finishes on a laptop.
#' The exclusion rules, the outcome-neutral threshold, the priors, the convergence
#' criteria and the ladders are the registered ones in both modes.
#'
#' @param cfg Output of stage2_config().
#' @param dry_run Logical.
#' @return `cfg` with `paths`, `languages_table`, `samples` and `sampler` resolved.
stage2_settings <- function(cfg, dry_run = FALSE) {
  if (!dry_run) return(cfg)
  base <- "outputs/stage2_dry_run"
  cfg$paths <- list(raw = file.path(base, "raw"), processed = file.path(base, "processed"),
                    fits = file.path(base, "fits"), results = file.path(base, "results"))
  cfg$languages_table <- file.path(base, "languages.csv")
  cfg$samples <- list(acceptability = 24L, ratings = 8L, min_verbs = 12L, max_verbs = 12L)
  cfg$sampler <- list(chains = 4L, iter = 1000L, warmup = 500L, cores = 4L,
                      adapt_delta = cfg$sampler$adapt_delta,
                      max_treedepth = cfg$sampler$max_treedepth)
  cfg$dry_run <- TRUE
  cfg
}

#' The meta-analytic prior SD of each interaction, per SD of affectedness.
#'
#' The per-language models give each interaction a zero-centred normal prior on the
#' coefficient of Semantics_scaled (R/03_define_priors.R, primary regime). Halving
#' a coefficient expresses it per SD of affectedness, so the same halving applies to
#' the prior SD. The report gives the results, 0.25 for H1 and 0.3 for H2.
#'
#' @param regime Prior regime in PRIOR_REGIMES.
#' @return Named numeric vector, H1 and H2.
meta_prior_sd <- function(regime = "primary") {
  if (!exists("PRIOR_REGIMES")) source("R/03_define_priors.R")
  r <- PRIOR_REGIMES[[regime]]
  sd_of <- function(spec) {
    x <- suppressWarnings(as.numeric(sub("^normal\\(0, *([0-9.]+)\\)$", "\\1", spec)))
    if (is.na(x)) stop("[stage2] expected a zero-centred normal prior, got ", spec)
    x
  }
  c(H1 = sd_of(r$b_active_int) / 2, H2 = sd_of(r$b_pseudo_int) / 2)
}

#' Read the variety-to-language table and check it.
#'
#' @param path CSV with columns variety, language and family (Glottolog names).
#' @return A data frame with one row per variety.
read_language_table <- function(path) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, comment.char = "#")
  need <- c("variety", "language", "family")
  if (!all(need %in% names(tab))) {
    stop("[stage2] ", path, " needs columns ", paste(need, collapse = ", "))
  }
  if (anyDuplicated(tab$variety)) {
    stop("[stage2] duplicated varieties in ", path, ": ",
         paste(unique(tab$variety[duplicated(tab$variety)]), collapse = ", "))
  }
  # A language in two families would make the family model ill-defined.
  fam <- tapply(tab$family, tab$language, function(x) length(unique(x)))
  if (any(fam > 1)) {
    stop("[stage2] languages listed under more than one family: ",
         paste(names(fam)[fam > 1], collapse = ", "))
  }
  tab
}

# ---------------------------------------------------------------------------
# Input format
# ---------------------------------------------------------------------------

STAGE2_COLUMNS <- list(
  acceptability = c("Variety", "Participant", "Trial", "Trial_Type", "Verb", "Agent",
                    "S_Type", "Attention_Target", "Response"),
  ratings       = c("Variety", "Participant", "Trial", "Trial_Type", "Verb", "Agent",
                    "Attention_Target", "Response")
)

#' Check one task's data against the canonical format, loudly.
#'
#' @param df Data frame for one task, one row per trial, attention checks included.
#' @param task "acceptability" or "ratings".
#' @return `df`, with Response as integer, or an error naming what is wrong.
#' @details A missing Response is allowed here, because it is what the completeness
#'   criterion excludes a participant for; every other malformed value stops the run,
#'   since none of the criteria below could be applied to it as registered.
validate_task_data <- function(df, task = c("acceptability", "ratings")) {
  task <- match.arg(task)
  need <- STAGE2_COLUMNS[[task]]
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("[stage2:", task, "] missing columns: ", paste(miss, collapse = ", "))
  bad <- function(what, rows) {
    stop("[stage2:", task, "] ", what, " in ", sum(rows), " row(s), e.g. row ", which(rows)[1])
  }
  if (any(!df$Trial_Type %in% c("test", "attention"))) {
    bad("Trial_Type other than test/attention", !df$Trial_Type %in% c("test", "attention"))
  }
  resp <- suppressWarnings(as.integer(df$Response))
  if (any(!is.na(df$Response) & (is.na(resp) | !resp %in% 1:7))) {
    bad("Response outside 1-7", !is.na(df$Response) & (is.na(resp) | !resp %in% 1:7))
  }
  df$Response <- resp
  test <- df$Trial_Type == "test"
  att  <- !test
  rows <- test & !df$Agent %in% c("Man", "Woman")
  if (any(rows)) bad("Agent other than Man/Woman", rows)
  rows <- test & (is.na(df$Verb) | df$Verb == "")
  if (any(rows)) bad("missing Verb", rows)
  if (task == "acceptability") {
    rows <- test & !df$S_Type %in% c("Active", "Passive", "Pseudo_Passive")
    if (any(rows)) bad("S_Type other than Active/Passive/Pseudo_Passive", rows)
  }
  rows <- att & !df$Attention_Target %in% c(1, 7)
  if (any(rows)) bad("attention trial without a target of 1 or 7", rows)
  key <- paste(df$Variety, df$Participant, df$Trial, sep = "\r")
  if (anyDuplicated(key)) bad("duplicated Variety-Participant-Trial", duplicated(key))
  # The same participant must not rate the same verb twice in the same version.
  if (task == "acceptability") {
    k2 <- paste(df$Variety, df$Participant, df$Verb, df$S_Type, sep = "\r")[test]
  } else {
    k2 <- paste(df$Variety, df$Participant, df$Verb, df$Agent, sep = "\r")[test]
  }
  if (anyDuplicated(k2)) stop("[stage2:", task, "] a participant rated the same stimulus twice")
  df
}

# ---------------------------------------------------------------------------
# Exclusions (Stage 1 report, Exclusion criteria)
# ---------------------------------------------------------------------------

#' Apply the registered participant exclusions to one task.
#'
#' A participant is excluded if any of three criteria holds:
#'   incomplete       a trial without a rating, or fewer trials than the design gives
#'                    (each verb once in every sentence type for the acceptability
#'                    task, and in both AGENT versions for the ratings);
#'   failed_attention fewer than `attention_min_correct` of the attention trials
#'                    answered with their target (1 or 7);
#'   straightlined    the same rating on `straightline_share` or more of the test
#'                    trials.
#'
#' @param df Validated data for one task.
#' @param task "acceptability" or "ratings".
#' @param n_attention Attention trials each participant should have.
#' @param attention_min_correct,straightline_share Registered thresholds.
#' @return One row per participant and variety, with each criterion and `excluded`.
participant_exclusions <- function(df, task = c("acceptability", "ratings"), n_attention,
                                   attention_min_correct = 8, straightline_share = 0.9) {
  task <- match.arg(task)
  versions <- if (task == "acceptability") "S_Type" else "Agent"
  # The design count: each verb in every version the variety uses.
  expected <- df |>
    dplyr::filter(Trial_Type == "test") |>
    dplyr::group_by(Variety) |>
    dplyr::summarise(n_test_expected = dplyr::n_distinct(Verb) * dplyr::n_distinct(.data[[versions]]),
                     .groups = "drop")
  df |>
    dplyr::group_by(Variety, Participant) |>
    dplyr::summarise(
      n_test            = sum(Trial_Type == "test"),
      n_missing         = sum(is.na(Response)),
      n_attention       = sum(Trial_Type == "attention"),
      n_attention_right = sum(Trial_Type == "attention" & Response == Attention_Target, na.rm = TRUE),
      modal_share       = {
        r <- Response[Trial_Type == "test" & !is.na(Response)]
        if (length(r)) max(tabulate(r, nbins = 7)) / length(r) else NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::left_join(expected, by = "Variety") |>
    dplyr::mutate(
      incomplete       = n_missing > 0 | n_test < n_test_expected | n_attention < n_attention,
      failed_attention = n_attention_right < attention_min_correct,
      straightlined    = !is.na(modal_share) & modal_share >= straightline_share,
      excluded         = incomplete | failed_attention | straightlined
    )
}

#' Keep the registered number of included participants per variety.
#'
#' Teams replace excluded participants until they have `n_target` usable ones, so a
#' variety normally arrives with exactly that many. If it has more, the first
#' `n_target` in order of data collection are kept, which needs `order_col`. If it
#' has fewer, the variety has not completed data collection and is not analysed.
#'
#' @param excl Output of participant_exclusions().
#' @param n_target Registered sample per variety (80 or 20).
#' @param order Optional data frame with Variety, Participant and Session_Order.
#' @return `excl` with a logical `analysed` column and a per-variety attribute
#'   "complete" naming the varieties that reached `n_target`.
cap_sample <- function(excl, n_target, order = NULL) {
  excl$analysed <- FALSE
  complete <- character(0)
  for (v in unique(excl$Variety)) {
    i <- which(excl$Variety == v & !excl$excluded)
    if (length(i) < n_target) next
    if (length(i) > n_target) {
      if (is.null(order)) {
        stop("[stage2] ", v, " has ", length(i), " included participants, more than ",
             n_target, ", and no Session_Order to keep the first ", n_target)
      }
      o <- order$Session_Order[match(paste(v, excl$Participant[i]),
                                     paste(order$Variety, order$Participant))]
      if (anyNA(o)) stop("[stage2] Session_Order missing for some participants of ", v)
      i <- i[order(o)][seq_len(n_target)]
    }
    excl$analysed[i] <- TRUE
    complete <- c(complete, v)
  }
  attr(excl, "complete") <- complete
  excl
}

# ---------------------------------------------------------------------------
# Affectedness (Stage 1 report, Pre-processing) and the outcome-neutral check
# ---------------------------------------------------------------------------

#' Mean affectedness rating of each verb, per variety.
#'
#' @param ratings Test trials of the ratings task from the analysed raters.
#' @return One row per variety and verb: Semantics (the mean over raters and both
#'   AGENT versions) and n_ratings.
verb_affectedness <- function(ratings) {
  ratings |>
    dplyr::filter(Trial_Type == "test") |>
    dplyr::group_by(Variety, Verb) |>
    dplyr::summarise(Semantics = mean(Response), n_ratings = dplyr::n(), .groups = "drop")
}

#' The outcome-neutral check: do the verbs differ enough in affectedness?
#'
#' The hypotheses concern the affectedness slope, which cannot be estimated in a
#' variety whose verbs are rated alike. The check reads only the ratings task, so it
#' is independent of the acceptability judgments it protects.
#'
#' @param verb_means Output of verb_affectedness().
#' @param min_sd Registered minimum SD of the verbs' mean ratings, in scale points.
#' @return One row per variety with the SD, the number of verbs and `passes`.
outcome_neutral_check <- function(verb_means, min_sd) {
  verb_means |>
    dplyr::group_by(Variety) |>
    dplyr::summarise(n_verbs = dplyr::n(), sd_verb_means = stats::sd(Semantics), .groups = "drop") |>
    dplyr::mutate(passes = sd_verb_means >= min_sd)
}

#' Centre the verbs' mean ratings and divide by two of their SDs, per variety.
#'
#' @param verb_means Output of verb_affectedness().
#' @return `verb_means` with Semantics_scaled, whose verb-level SD is 0.5.
scale_verb_affectedness <- function(verb_means) {
  verb_means |>
    dplyr::group_by(Variety) |>
    dplyr::mutate(Semantics_scaled = (Semantics - mean(Semantics)) / (2 * stats::sd(Semantics))) |>
    dplyr::ungroup()
}

#' Assemble one variety's analysis data.
#'
#' @param acc Test trials of the acceptability task for the analysed participants.
#' @param verb_scores Output of scale_verb_affectedness() for the same variety.
#' @return Data frame ready for the models: Response (integer 1-7), S_Type
#'   treatment-coded with Passive as reference, Semantics_scaled, Participant and
#'   Verb as factors. Attribute has_pseudo_passive records whether the variety used
#'   the construction.
build_variety_data <- function(acc, verb_scores) {
  if (!exists("code_s_type")) source("R/02_preprocess_factors.R")
  if (length(unique(acc$Variety)) != 1L) stop("[stage2] build_variety_data takes one variety")
  acc <- dplyr::filter(acc, Trial_Type == "test")
  miss <- setdiff(unique(acc$Verb), verb_scores$Verb)
  if (length(miss)) stop("[stage2] verbs without affectedness ratings: ", paste(miss, collapse = ", "))
  has_pp <- "Pseudo_Passive" %in% acc$S_Type
  d <- acc |>
    dplyr::select(Variety, Participant, Verb, Agent, S_Type, Response) |>
    dplyr::left_join(dplyr::select(verb_scores, Verb, Semantics, Semantics_scaled), by = "Verb") |>
    dplyr::mutate(Participant = factor(Participant), Verb = factor(Verb))
  d <- code_s_type(d)
  d <- drop_pseudo_passive_if_absent(d, has_pp)
  assert_treatment_coding(d)
  attr(d, "has_pseudo_passive") <- has_pp
  d
}

# ---------------------------------------------------------------------------
# Individual-language models (Stage 1 report, Individual-language models)
# ---------------------------------------------------------------------------

STAGE2_TERMS <- c(H1 = "b_S_TypeActive:Semantics_scaled",
                  H2 = "b_S_TypePseudo_Passive:Semantics_scaled")

#' Fit one variety's model by the registered ladder.
#'
#' @param d Output of build_variety_data().
#' @param cfg Output of stage2_config().
#' @param out_dir Directory for the fits and their diagnostics.
#' @param label File label for the variety.
#' @param seed Sampler seed for this variety.
#' @param sampler Optional list overriding cfg$sampler, for the dry run only.
#' @param start_level NULL to start at the maximal model. A lower level only to
#'   resume a ladder whose previous level hit the cluster's wall time, one of the
#'   permitted fallback triggers in R/09_model_ladder.R.
#' @return The list returned by run_model_ladder().
fit_variety_model <- function(d, cfg, out_dir, label, seed, sampler = NULL, start_level = NULL) {
  if (!exists("run_model_ladder")) source("R/09_model_ladder.R")
  s <- if (is.null(sampler)) cfg$sampler else sampler
  has_pp <- isTRUE(attr(d, "has_pseudo_passive"))
  run_model_ladder(
    df            = d,
    prior_obj     = build_brms_prior("primary", "broad", has_pseudo_passive = has_pp),
    sampling_args = list(iter = s$iter, warmup = s$warmup, chains = s$chains,
                         cores = s$cores, seed = seed),
    control_args  = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth),
    has_pseudo_passive = has_pp,
    start_level   = start_level,
    out_dir       = out_dir,
    language      = label,
    prior_label   = "primary"
  )
}

#' Posterior mean and SD of each focal interaction, per SD of affectedness.
#'
#' @param fit A brmsfit, or a matrix or data frame of posterior draws with the
#'   coefficient names as columns (which is what the tests pass).
#' @return One row per hypothesis present in the fit: the coefficient's posterior
#'   mean and SD, and both halved (per SD of affectedness).
extract_interactions <- function(fit) {
  if (inherits(fit, "brmsfit")) {
    terms <- unname(STAGE2_TERMS)
    dr <- as.data.frame(fit, variable = terms[terms %in% brms::variables(fit)])
  } else {
    dr <- as.data.frame(fit)
  }
  present <- STAGE2_TERMS[STAGE2_TERMS %in% names(dr)]
  if (!"H1" %in% names(present)) stop("[stage2] the fit has no ", STAGE2_TERMS[["H1"]])
  tibble::tibble(
    hypothesis  = names(present),
    term        = unname(present),
    coef_mean   = vapply(present, function(t) mean(dr[[t]]), numeric(1)),
    coef_sd     = vapply(present, function(t) stats::sd(dr[[t]]), numeric(1)),
    per_sd_mean = coef_mean / 2,
    per_sd_sd   = coef_sd / 2
  )
}

# ---------------------------------------------------------------------------
# Varieties and the meta-analysis (Stage 1 report, Meta-analytic synthesis)
# ---------------------------------------------------------------------------

#' Average a language's varieties into one estimate per hypothesis.
#'
#' The mean of the varieties' posterior means is the language's estimate, and the
#' mean of their posterior SDs its standard error, which treats the varieties'
#' errors as fully correlated (they share verbs). A variety contributes to H2 only if
#' it used a pseudo-passive.
#'
#' @param est One row per variety and hypothesis, from extract_interactions(), with
#'   a Variety column.
#' @param lang_table Output of read_language_table().
#' @return One row per language and hypothesis: language, family, y, se and the
#'   number of varieties averaged.
average_varieties <- function(est, lang_table) {
  unknown <- setdiff(unique(est$Variety), lang_table$variety)
  if (length(unknown)) stop("[stage2] varieties missing from the language table: ",
                            paste(unknown, collapse = ", "))
  est |>
    dplyr::left_join(dplyr::select(lang_table, variety, language, family),
                     by = c("Variety" = "variety")) |>
    dplyr::group_by(language, family, hypothesis) |>
    dplyr::summarise(y = mean(per_sd_mean), se = mean(per_sd_sd),
                     n_varieties = dplyr::n(), varieties = paste(sort(Variety), collapse = "; "),
                     .groups = "drop")
}

#' Directional Bayes factor from the posterior probability that the mean is negative.
#'
#' With a zero-centred, symmetric prior the prior odds of either sign are 1, so the
#' Bayes factor for the predicted negative sign is the posterior odds.
#'
#' @param p Posterior probability that the parameter is negative.
#' @return The Bayes factor for "negative" against "positive".
bf_from_p_negative <- function(p) p / (1 - p)

#' Label a Bayes factor by the registered thresholds.
#'
#' @param bf Directional Bayes factor for the predicted (negative) sign.
#' @param threshold Registered threshold, 10.
#' @return "for" (BF >= threshold), "against" (BF <= 1/threshold: strong evidence
#'   that the effect runs the other way) or "inconclusive".
interpret_bf <- function(bf, threshold = 10) {
  ifelse(bf >= threshold, "for", ifelse(bf <= 1 / threshold, "against", "inconclusive"))
}

#' The registered random-effects meta-analysis of one hypothesis.
#'
#' Normal prior on the mean across languages and half-normal prior on the
#' between-language SD, computed by bayesmeta (Röver, 2020), which integrates
#' numerically, so the Bayes factor carries no Monte Carlo error.
#'
#' @param y,se Per-language estimates and standard errors, per SD of affectedness.
#' @param prior_mu_sd Prior SD of the mean (meta_prior_sd()).
#' @param tau_scale Scale of the half-normal prior on the between-language SD.
#' @return One-row tibble: k, P(mu < 0), the Bayes factor, the posterior mean, SD and
#'   95% interval of mu, and the posterior median and 95% interval of tau.
run_meta_test <- function(y, se, prior_mu_sd, tau_scale = 0.5) {
  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    stop("[stage2] the bayesmeta package is required for the confirmatory meta-analysis")
  }
  stopifnot(length(y) == length(se), length(y) >= 2, all(se > 0))
  b <- bayesmeta::bayesmeta(y = y, sigma = se, mu.prior.mean = 0, mu.prior.sd = prior_mu_sd,
                            tau.prior = function(t) bayesmeta::dhalfnormal(t, scale = tau_scale))
  p_neg <- b$pposterior(mu = 0)
  mu_ci <- b$post.interval(mu.level = 0.95)
  tau_ci <- b$post.interval(tau.level = 0.95)
  tibble::tibble(
    k = length(y), p_mu_negative = p_neg, bf = bf_from_p_negative(p_neg),
    mu_mean = b$summary["mean", "mu"], mu_sd = b$summary["sd", "mu"],
    mu_lo = mu_ci[1], mu_hi = mu_ci[2],
    tau_median = b$summary["median", "tau"], tau_lo = tau_ci[1], tau_hi = tau_ci[2]
  )
}

# ---------------------------------------------------------------------------
# Robustness 1: the single model across languages
# ---------------------------------------------------------------------------

#' The registered ladder for the single model across languages.
#'
#' R0 is the report's maximal model. Participant and Verb are the identifiers
#' made unique within variety by build_pooled_data(), which is how the report's
#' Language:Participant and Language:Verb nesting is written, and also how
#' participants and verbs come to be nested within variety. No level has a by-verb
#' affectedness slope, which the report rules out as not meaningful.
#'   R0  maximal, correlated
#'   R1  as R0, uncorrelated
#'   R2  as R1 without the by-participant interaction slope
#'   R3  sentence-type slopes only for participants and verbs (the report's
#'       anticipated model)
#'   R4  by-language sentence-type slopes, random intercepts otherwise
#'   R5  random intercepts only
#'
#' @param s Name of the affectedness predictor.
#' @return Named list of brmsformula objects, R0 first.
build_single_model_ladder <- function(s = "Semantics_scaled") {
  fam <- brms::cumulative(link = "logit", threshold = "flexible")
  fe <- paste0("Response ~ S_Type * ", s)
  f <- function(re) brms::bf(stats::as.formula(paste(fe, re)), family = fam)
  list(
    R0_maximal = f(paste0("+ (1 + S_Type * ", s, " | Language)",
                          " + (1 + S_Type * ", s, " | Participant)",
                          " + (1 + S_Type | Verb)")),
    R1_uncorrelated = f(paste0("+ (1 + S_Type * ", s, " || Language)",
                               " + (1 + S_Type * ", s, " || Participant)",
                               " + (1 + S_Type || Verb)")),
    R2_no_participant_interaction = f(paste0("+ (1 + S_Type * ", s, " || Language)",
                                             " + (1 + S_Type + ", s, " || Participant)",
                                             " + (1 + S_Type || Verb)")),
    R3_sentence_type_slopes = f(paste0("+ (1 + S_Type * ", s, " || Language)",
                                       " + (1 + S_Type || Participant)",
                                       " + (1 + S_Type || Verb)")),
    R4_language_sentence_type = f(paste0("+ (1 + S_Type || Language)",
                                         " + (1 | Participant)",
                                         " + (1 | Verb)")),
    R5_intercepts = f(paste0("+ (1 | Language)",
                             " + (1 | Participant)",
                             " + (1 | Verb)"))
  )
}

#' Stack the varieties' analysis data for the single model.
#'
#' @param variety_data List of build_variety_data() outputs.
#' @param lang_table Output of read_language_table().
#' @return One data frame with Language, and Participant and Verb identifiers that
#'   are unique within variety.
build_pooled_data <- function(variety_data, lang_table) {
  d <- dplyr::bind_rows(lapply(variety_data, function(x) {
    x$Participant <- as.character(x$Participant)
    x$Verb <- as.character(x$Verb)
    x$S_Type <- as.character(x$S_Type)
    x
  }))
  d$Language <- lang_table$language[match(d$Variety, lang_table$variety)]
  if (anyNA(d$Language)) stop("[stage2] varieties missing from the language table")
  d |>
    dplyr::mutate(Participant = factor(paste(Variety, Participant, sep = " / ")),
                  Verb = factor(paste(Variety, Verb, sep = " / ")),
                  Language = factor(Language)) |>
    code_s_type()
}

#' Directional Bayes factor from posterior draws, under a symmetric prior.
#'
#' @param draws Posterior draws of the coefficient.
#' @return List with p_negative, bf and, where every draw has the same sign, the
#'   bound the draws allow in place of a point value.
bf_from_draws <- function(draws) {
  n <- length(draws)
  k <- sum(draws < 0)
  if (k == n) return(list(p_negative = 1, bf = Inf, bound = paste0("BF > ", n - 1)))
  if (k == 0) return(list(p_negative = 0, bf = 0, bound = paste0("BF < 1/", n - 1)))
  p <- k / n
  list(p_negative = p, bf = bf_from_p_negative(p), bound = NA_character_)
}

# ---------------------------------------------------------------------------
# Robustness 2: languages nested within families
# ---------------------------------------------------------------------------

#' The meta-analysis with languages nested within their Glottolog families.
#'
#' Same likelihood and priors as the confirmatory meta-analysis, with a second
#' half-normal(0.5) prior on the between-family SD, fitted in brms.
#'
#' @param lang_est Rows of average_varieties() for one hypothesis.
#' @param prior_mu_sd,tau_scale As for run_meta_test().
#' @param sampler List with chains, iter, warmup, cores, adapt_delta, max_treedepth.
#' @param seed Sampler seed.
#' @return List with the fit, its convergence status and the directional Bayes factor.
fit_family_meta <- function(lang_est, prior_mu_sd, tau_scale, sampler, seed) {
  if (!exists("extract_convergence_diagnostics")) source("R/07_extract_diagnostics.R")
  d <- data.frame(y = lang_est$y, se = lang_est$se,
                  family = factor(lang_est$family), language = factor(lang_est$language))
  pr <- c(brms::prior_string(sprintf("normal(0, %s)", prior_mu_sd), class = "Intercept"),
          brms::prior_string(sprintf("normal(0, %s)", tau_scale), class = "sd", group = "family"),
          brms::prior_string(sprintf("normal(0, %s)", tau_scale), class = "sd", group = "family:language"))
  fit <- brms::brm(y | se(se) ~ 1 + (1 | family / language), data = d, prior = pr,
                   backend = getOption("brms.backend", "cmdstanr"),
                   chains = sampler$chains, iter = sampler$iter,
                   warmup = sampler$warmup, cores = sampler$cores, seed = seed,
                   control = list(adapt_delta = sampler$adapt_delta,
                                  max_treedepth = sampler$max_treedepth),
                   refresh = 0, silent = 2)
  diag <- extract_convergence_diagnostics(fit, max_treedepth = sampler$max_treedepth)
  mu <- as.data.frame(fit, variable = "b_Intercept")$b_Intercept
  c(list(fit = fit, diagnostics = diag, status = classify_convergence(diag),
         mu_mean = mean(mu), mu_sd = stats::sd(mu)),
    bf_from_draws(mu))
}
