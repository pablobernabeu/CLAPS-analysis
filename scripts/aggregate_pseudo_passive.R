#!/usr/bin/env Rscript
# scripts/aggregate_pseudo_passive.R
# ---------------------------------------------------------------------------
# Summarise the pseudo-passive prediction (H2) from the pilot-grounded decision
# arm, for the report (reports/preliminary_sample_size_analysis.qmd).
#
# Every decision-arm replicate stores directional Bayes factors for both signs
# of the pseudo-passive-by-affectedness interaction, H2a (positive) and H2b
# (negative), but scripts/aggregate_pilot.R summarises only H1a and H1b. The
# literature-anchored cross-check cannot stand in for them. Its grid generators
# set that interaction positive, where the meta-analysis (Ambridge, Arnon &
# Bekman, 2023) and the English pilot put it below zero, so its pseudo-passive
# rates describe an effect the design does not expect. The decision arm draws
# every effect from each language's own pilot posterior, so its rates carry the
# sign each pilot supports.
#
# The summary is kept apart from aggregate_pilot.R because
# joint_power_pilot.csv pools the decision arm with replicates of an earlier arm
# at 70 and 100 participants. Reading the decision arm alone keeps every H2 rate
# on the arm's own replicates. Its H1b counts are checked against
# joint_power_pilot.csv wherever the two files hold the same replicates, and that
# check ties the two files together.
#
# Outputs (to --outdir):
#   pseudo_passive_decision_arm.csv   One row per language x N: the replicate
#     and converged counts, and n_h1b, the H1b count at a Bayes factor of 10
#     that is checked against joint_power_pilot.csv. For each sign of H2, the
#     count reaching 10 (n_h2b, n_h2a) with its proportion, exact two-sided 95%
#     interval (lo_, hi_) and exact one-sided 95% lower bound (lb_), the count at
#     6 (_bf6) and the count over the converged replicates (_converged). Then
#     n_h2_inconclusive, where neither sign reaches 10, and n_h1b_and_h2b, where
#     H1b and H2b both reach 10 in the same simulated study. Counts are integers,
#     so a reader needs no round(p * reps). A language without a pseudo-passive
#     keeps its row for the H1b check, with every H2 column NA, so that "no such
#     prediction" cannot be read as "never detected".
#   pseudo_passive_pilot_posterior.csv   One row per language with a
#     pseudo-passive: the pilot posterior of the interaction per SD of
#     affectedness, the scale pilot_params_ceilings.csv uses for the focal
#     effects (mean, 95% interval, count and share of draws below zero). Written
#     only when --dgpdir holds the pilot DGP files, which are not tracked.
#
# Run from design_analysis/, with the replicates copied from the cluster:
#   Rscript scripts/aggregate_pseudo_passive.R \
#     --cells  <data>/outputs/design_databased_v2_decision \
#     --dgpdir <data>/outputs/pilot_models
# The grid, the reference and the output folder default to their places in this
# repository.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  # No default. The replicates are not tracked, so a default path would only
  # ever find nothing, and a summary of nothing must not be written.
  optparse::make_option("--cells"),
  # The grid defines the arm, both which seeds belong to it and what each
  # replicate should be. Taking the seeds from it keeps the completeness check
  # below in step with the grid if the grid ever changes.
  optparse::make_option("--grid",      default = "config/design_grid_databased_v2_decision.csv"),
  optparse::make_option("--reference", default = "outputs/design_summary_pilot/joint_power_pilot.csv"),
  optparse::make_option("--dgpdir",    default = "outputs/pilot_models"),
  optparse::make_option("--outdir",    default = "outputs/design_summary_pilot")
)))
tag <- "[aggregate pseudo-passive]"
if (is.null(opt$cells) || !dir.exists(opt$cells)) {
  stop(tag, " --cells must name the directory holding the decision-arm replicates")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Labels as R/05_hypothesis_tests.R writes them into each replicate's bf_results.
HYP <- c(h1b = "H1b_active_interaction_negative",
         h2a = "H2a_pseudo_positive",
         h2b = "H2b_pseudo_negative")
# What each label must have been computed for. A label only names a hypothesis, so a
# row is trusted only when its direction and coefficient agree with the name.
EXPECT <- list(
  h1b = c(direction = "negative", param = "b_S_TypeActive:Semantics_scaled"),
  h2a = c(direction = "positive", param = "b_S_TypePseudo_Passive:Semantics_scaled"),
  h2b = c(direction = "negative", param = "b_S_TypePseudo_Passive:Semantics_scaled"))

grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
if (anyDuplicated(grid$seed)) {
  stop(tag, " the grid repeats a seed, so a seed cannot identify a replicate")
}
# The flag is a property of the language. A grid that set it differently for two
# rows of one language would leave the H2 check below with no single answer.
if (any(tapply(grid$has_pseudo_passive, grid$language, function(x) length(unique(x))) != 1L)) {
  stop(tag, " the grid sets has_pseudo_passive differently within a language")
}

# ---- 1. Read every replicate --------------------------------------------------
# The first value a replicate recorded under a name, or NA when it recorded none,
# since a replicate whose fit failed carries a shorter summary than one that ran.
field <- function(x, name) {
  v <- x[[name]]
  if (is.null(v) || length(v) == 0L) NA else v[[1]]
}

files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0L) stop(tag, " no .rds files in ", opt$cells)

# Files that yield no summary at all are kept by name. A truncated .rds (a job
# killed mid-save) fails to deserialise, and it cannot report its own seed.
unreadable <- character(0)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  s <- if (is.list(d)) d$summary else NULL
  if (is.null(s) || is.null(s$seed)) {
    unreadable <- c(unreadable, basename(files[[i]]))
    next
  }
  bf <- d$bf_results
  has_bf <- is.data.frame(bf) && all(c("hypothesis", "BF_10") %in% names(bf))
  # How often each label appears. Zero is expected for H2 in a language without a
  # pseudo-passive. Two or more cannot be resolved, and the run stops on it below.
  n_lab <- vapply(HYP, function(h) if (has_bf) sum(bf$hypothesis == h) else 0L, integer(1))
  bf_val <- function(key) {
    if (n_lab[[key]] != 1L) return(NA_real_)
    suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == HYP[[key]]]))
  }
  # Absent and repeated labels are caught below. Here a present row must carry the
  # direction and coefficient its label promises.
  lab_ok <- has_bf && all(c("direction", "param") %in% names(bf)) &&
    all(vapply(names(HYP), function(key) {
      if (n_lab[[key]] != 1L) return(TRUE)
      r <- bf[bf$hypothesis == HYP[[key]], ]
      identical(as.character(r$direction), EXPECT[[key]][["direction"]]) &&
        identical(as.character(r$param), EXPECT[[key]][["param"]])
    }, logical(1)))
  g <- d$diagnostics
  rows[[i]] <- data.frame(
    file           = basename(files[[i]]),
    seed           = as.numeric(field(s, "seed")),
    status         = as.character(field(s, "status")),
    language       = as.character(field(s, "language")),
    n_participants = as.numeric(field(s, "n_participants")),
    mode           = as.character(field(s, "mode")),
    prior_source   = as.character(field(s, "prior_source")),
    prior_regime   = as.character(field(s, "prior_regime")),
    model_level    = as.character(field(s, "model_level")),
    threshold_mode = as.character(field(s, "threshold_mode")),
    draw_index     = as.numeric(field(s, "draw_index")),
    n_verbs        = as.numeric(field(s, "n_verbs")),
    chains         = as.numeric(field(s, "chains")),
    iter           = as.numeric(field(s, "iter")),
    warmup         = as.numeric(field(s, "warmup")),
    has_bf         = has_bf,
    dup_label      = any(n_lab > 1L),
    lab_ok         = lab_ok,
    # Read as aggregate_pilot.R reads it, so that n_converged here can be compared
    # with the n_converged column of joint_power_pilot.csv.
    converged      = if (is.null(g) || is.null(g$convergence_ok)) NA else isTRUE(g$convergence_ok[[1]]),
    bf_h1b         = bf_val("h1b"),
    bf_h2a         = bf_val("h2a"),
    bf_h2b         = bf_val("h2b"),
    stringsAsFactors = FALSE
  )
}
reps_all <- dplyr::bind_rows(rows)

# ---- 2. Keep the arm, and refuse anything short of all of it ------------------
# The directory may hold other arms' replicates too, as the merged directory on the
# cluster does, so only the seeds the grid lists are kept. Every committed grid has
# a seed set disjoint from every other (tests/testthat/test-seed-disjointness.R),
# so a seed identifies a decision-arm replicate.
r <- reps_all[reps_all$seed %in% grid$seed, ]
message(sprintf("%s %d file(s) read, %d in the arm, %d outside it and ignored, %d unreadable",
                tag, nrow(reps_all), nrow(r), nrow(reps_all) - nrow(r), length(unreadable)))

# The runner names each file after its seed, which is the only way to place a file
# that could not be read. One that belongs to the arm stops the run, because a
# summary over the rest would shrink a denominator without saying so.
name_seed <- suppressWarnings(as.numeric(sub("^.*_([0-9]+)\\.rds$", "\\1", unreadable)))
if (any(name_seed %in% grid$seed)) {
  stop(tag, " unreadable replicate(s) belonging to the arm: ",
       paste(unreadable[name_seed %in% grid$seed], collapse = ", "))
}

# A fit that fails is caught inside the replicate and recorded in its status, so the
# array task still exits 0 and the scheduler reports success. This is where such a
# failure first becomes visible.
failed <- r$seed[is.na(r$status) | r$status != "success"]
if (length(failed)) {
  stop(tag, " ", length(failed), " replicate(s) in the arm did not succeed, seeds ",
       paste(head(sort(failed), 20), collapse = ", "))
}

dup_seed <- unique(r$seed[duplicated(r$seed)])
missing_seed <- setdiff(grid$seed, r$seed)
if (length(dup_seed) || length(missing_seed)) {
  stop(tag, " the arm is not complete and unique: ", length(missing_seed), " grid seed(s) missing",
       if (length(missing_seed)) paste0(" (", paste(head(sort(missing_seed), 20), collapse = ", "), ")"),
       ", ", length(dup_seed), " seed(s) found more than once")
}

# Each replicate must be the one its grid row describes. The prior regime is left out
# of this comparison because its name differs between the cluster and the published
# copy of the grid. It is checked below for a single value across the arm instead.
chk <- dplyr::inner_join(r, grid, by = "seed", suffix = c("", ".grid"))
for (col in c("language", "n_participants", "mode", "prior_source", "model_level",
              "threshold_mode", "draw_index", "iter", "warmup", "chains")) {
  off <- is.na(chk[[col]]) | chk[[col]] != chk[[paste0(col, ".grid")]]
  if (any(off)) {
    stop(tag, " ", sum(off), " replicate(s) disagree with their grid row on ", col,
         ", for example seed ", chk$seed[which(off)[1]])
  }
}

# The replicate summary identifies a cell by language, N and draw, not by sampler
# settings, so two runs of one grid under different samplers would pass every check
# above. The output records one sampler and one model per arm, and says nothing true
# if the arm mixes them. n_verbs is fixed within a language, since Norwegian has fewer.
one_value <- function(x) length(unique(x)) == 1L
if (!one_value(paste(chk$chains, chk$iter, chk$warmup))) stop(tag, " the arm mixes sampler settings")
if (!one_value(chk$prior_regime)) stop(tag, " the arm mixes prior regimes")
if (!one_value(chk$model_level)) stop(tag, " the arm mixes model levels")
if (any(tapply(chk$n_verbs, chk$language, function(x) length(unique(x))) != 1L)) {
  stop(tag, " the verb count varies within a language")
}

# ---- 3. The Bayes factors themselves -------------------------------------------
if (!all(chk$has_bf) || any(chk$dup_label)) {
  stop(tag, " a replicate lacks a Bayes-factor table or repeats a hypothesis label")
}
if (!all(is.finite(chk$bf_h1b))) stop(tag, " a replicate has no finite H1b Bayes factor")
# Reciprocity cannot catch a swap of H2a and H2b, since swapped rows stay reciprocal.
# The direction and coefficient each row records can.
if (!all(chk$lab_ok)) {
  stop(tag, " a Bayes-factor row carries the wrong direction or coefficient for its label, ",
       "for example seed ", chk$seed[which(!chk$lab_ok)[1]])
}

# H2 must be present exactly where the grid says the language has a pseudo-passive.
# A language that had H2 in some replicates and not in others would give a count
# over an unknown subset.
chk$has_pp <- as.logical(chk$has_pseudo_passive)
h2_present <- !is.na(chk$bf_h2a) & !is.na(chk$bf_h2b)
h2_absent <- is.na(chk$bf_h2a) & is.na(chk$bf_h2b)
wrong <- (chk$has_pp & !h2_present) | (!chk$has_pp & !h2_absent)
if (any(wrong)) {
  stop(tag, " H2 is missing or unexpected in ", sum(wrong), " replicate(s), for example seed ",
       chk$seed[which(wrong)[1]])
}

# H2a and H2b are the two signs of one coefficient under a zero-centred prior, so
# each Bayes factor is the reciprocal of the other. Checking it makes "neither sign
# reaches 10" the same event as a Bayes factor for H2b between 1/10 and 10. It does not
# detect swapped labels, which stay reciprocal; the direction check above does that.
recip <- abs(log(chk$bf_h2a[chk$has_pp]) + log(chk$bf_h2b[chk$has_pp]))
if (!all(is.finite(recip)) || any(recip > 1e-6)) {
  stop(tag, " H2a and H2b are not reciprocal Bayes factors in every replicate")
}

# ---- 4. One row per language x N ------------------------------------------------
exact_ci <- function(k, n) stats::binom.test(k, n)$conf.int[1:2]
# The rule lower95() applies in the report, which judges adequacy on the exact
# one-sided lower bound.
lower95 <- function(k, n) stats::binom.test(k, n, alternative = "greater")$conf.int[1]

cell_summary <- function(d) {
  n <- nrow(d)
  pp <- d$has_pp[1]
  conv <- d$converged %in% TRUE
  # A count, or NA where the language has no pseudo-passive and so no H2 at all.
  k <- function(x) if (pp) as.integer(sum(x)) else NA_integer_
  ci <- function(kk) if (pp) exact_ci(kk, n) else c(NA_real_, NA_real_)
  lb <- function(kk) if (pp) lower95(kk, n) else NA_real_
  n_h2b <- k(d$bf_h2b >= 10)
  n_h2a <- k(d$bf_h2a >= 10)
  n_h2b_bf6 <- k(d$bf_h2b >= 6)
  n_h2a_bf6 <- k(d$bf_h2a >= 6)
  ci_b <- ci(n_h2b)
  ci_a <- ci(n_h2a)
  data.frame(
    language           = d$language[1],
    n_participants     = as.integer(d$n_participants[1]),
    reps               = n,
    n_converged        = as.integer(sum(conv)),
    has_pseudo_passive = pp,
    n_h1b              = as.integer(sum(d$bf_h1b >= 10)),
    n_h2b = n_h2b, p_h2b = n_h2b / n, lo_h2b = ci_b[1], hi_h2b = ci_b[2], lb_h2b = lb(n_h2b),
    n_h2a = n_h2a, p_h2a = n_h2a / n, lo_h2a = ci_a[1], hi_h2a = ci_a[2], lb_h2a = lb(n_h2a),
    n_h2_inconclusive  = k(d$bf_h2b < 10 & d$bf_h2a < 10),
    n_h2b_bf6 = n_h2b_bf6, p_h2b_bf6 = n_h2b_bf6 / n,
    n_h2a_bf6 = n_h2a_bf6, p_h2a_bf6 = n_h2a_bf6 / n,
    n_h2b_converged    = k(conv & d$bf_h2b >= 10),
    n_h2a_converged    = k(conv & d$bf_h2a >= 10),
    n_h1b_and_h2b      = k(d$bf_h1b >= 10 & d$bf_h2b >= 10),
    seed_min           = as.integer(min(d$seed)),
    seed_max           = as.integer(max(d$seed)),
    chains             = as.integer(d$chains[1]),
    iter               = as.integer(d$iter[1]),
    warmup             = as.integer(d$warmup[1]),
    stringsAsFactors   = FALSE
  )
}
summ <- dplyr::bind_rows(lapply(split(chk, list(chk$language, chk$n_participants), drop = TRUE),
                                cell_summary)) |>
  dplyr::arrange(language, n_participants)

# ---- 5. Tie the file to joint_power_pilot.csv ------------------------------------
# A cell can be compared only where the reference holds exactly these replicates,
# which its replicate count shows. Where it does, the H1b count and the converged
# count must agree to the replicate. Every language must have at least one such
# cell, so the check can never pass by comparing nothing.
ref <- readr::read_csv(opt$reference, show_col_types = FALSE) |>
  dplyr::filter(mode %in% unique(chk$mode)) |>
  dplyr::transmute(language, n_participants = as.integer(n_participants),
                   reference_reps = as.integer(reps),
                   reference_n_h1b = as.integer(round(p_h1b * reps)),
                   reference_n_converged = as.integer(n_converged))
summ <- dplyr::left_join(summ, ref, by = c("language", "n_participants"))
comparable <- !is.na(summ$reference_reps) & summ$reference_reps == summ$reps
if (!all(tapply(comparable, summ$language, any))) {
  stop(tag, " no cell of ", paste(names(which(!tapply(comparable, summ$language, any))), collapse = ", "),
       " can be compared with ", basename(opt$reference))
}
mism <- comparable & (summ$n_h1b != summ$reference_n_h1b |
                        summ$n_converged != summ$reference_n_converged)
if (any(mism)) {
  stop(tag, " H1b or convergence counts disagree with ", basename(opt$reference), " at ",
       paste(sprintf("%s N=%d", summ$language[mism], summ$n_participants[mism]), collapse = "; "))
}
ref_name <- tools::file_path_sans_ext(basename(opt$reference))
summ$validation <- ifelse(comparable, paste0("matches_", ref_name), "not_compared_reps_differ")
summ$source_grid <- basename(opt$grid)
summ$reference_n_h1b <- NULL
summ$reference_n_converged <- NULL
readr::write_csv(summ, file.path(opt$outdir, "pseudo_passive_decision_arm.csv"))
message(sprintf("%s %d replicates -> pseudo_passive_decision_arm.csv (%d rows, %d checked against %s)",
                tag, nrow(chk), nrow(summ), sum(comparable), basename(opt$reference)))

# ---- 6. The pilot posterior of the interaction ----------------------------------
# Per SD of affectedness, as pilot_params_ceilings.csv expresses the focal effects,
# so the report can state the direction each pilot supports without writing it out.
PP_TERM <- "S_TypePseudo_Passive:Semantics_scaled"
langs <- sort(unique(grid$language))
dgp_files <- file.path(opt$dgpdir, paste0("pilot_dgp_v2_pilot_", langs, ".rds"))
if (!all(file.exists(dgp_files))) {
  message(tag, " pilot DGP files not found in ", opt$dgpdir,
          ", so pseudo_passive_pilot_posterior.csv is left as it is")
} else {
  post <- dplyr::bind_rows(lapply(seq_along(langs), function(i) {
    d <- readRDS(dgp_files[[i]])
    # unclass() so that column extraction does not need the posterior package.
    draws <- unclass(d$fixef_draws)
    if (!PP_TERM %in% colnames(draws)) return(NULL)
    sd_x <- stats::sd(d$verb_affectedness)
    b <- as.numeric(draws[, PP_TERM])
    x <- b * sd_x
    data.frame(language = langs[[i]], affectedness_sd = sd_x, ndraws = length(x),
               fixef_mean = mean(b), per_sd_mean = mean(x),
               per_sd_lo95 = unname(stats::quantile(x, 0.025)),
               per_sd_hi95 = unname(stats::quantile(x, 0.975)),
               n_negative = sum(x < 0), p_negative = mean(x < 0),
               stringsAsFactors = FALSE)
  }))
  # The pilot models and the grid must agree on which languages have a
  # pseudo-passive, or the two output files would describe different languages.
  if (!setequal(post$language, unique(grid$language[grid$has_pseudo_passive]))) {
    stop(tag, " the pilot models and the grid disagree on which languages have a pseudo-passive")
  }
  readr::write_csv(post, file.path(opt$outdir, "pseudo_passive_pilot_posterior.csv"))
  message(sprintf("%s pseudo_passive_pilot_posterior.csv written (%d languages)", tag, nrow(post)))
}
