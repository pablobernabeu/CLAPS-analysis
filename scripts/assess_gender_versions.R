#!/usr/bin/env Rscript
# scripts/assess_gender_versions.R
# ---------------------------------------------------------------------------
# How much precision does the reversed agent-gender version of each stimulus
# buy? The PI asked on 16 September 2026, just before the Stage 1 Registered
# Report is finalised: dropping one of the two versions ("The woman pushed the
# man" / "The man pushed the woman", and likewise for the passive and
# pseudo-passive) would halve the session from 40 to 20 minutes, so what does it
# cost the per-language tests?
#
# WHY THE PILOT CAN ANSWER THIS DIRECTLY
#   Every pilot participant rated BOTH versions of every verb x sentence type
#   cell. Deleting one of the two rows per cell yields the data a one-version
#   session would have produced from the same people, so the information loss
#   can be measured by fitting one model to both datasets, with nothing assumed.
#
# WHAT A SECOND RATING CAN AND CANNOT BUY
#   The focal terms are verb-level contrasts over a fixed verb set: H1a is the
#   affectedness slope in passives and H1b (the primary per-language test) is its
#   difference in actives. A second rating of a cell averages away only the part
#   of per-rating noise that the two versions do not share. Between-verb
#   variation in the sentence-type effect and between-participant variation in
#   the slopes are untouched. So SE_full / SE_half must lie between 1/sqrt(2)
#   (pure rating noise) and 1 (pure verb and participant heterogeneity); this
#   script locates it.
#
# WHY IT MATTERS FOR THE NUMBERS ALREADY SENT
#   The decision-arm power figures (R/10_simulate_from_pilot.R) simulate ONE
#   rating per participant x verb x sentence type. They already describe the
#   halved design, so the question is how much the two-version design would add
#   to them, and the last output translates the measured ratio into that.
#
# WHAT IT DOES, PER LANGUAGE
#   1. Structure: both versions per cell, randomisation of trial order, the lag
#      between the two versions, and the Norwegian label anomaly (watch x Active,
#      both rows labelled woman-agent), resolved from the raw sentence text.
#   2. Agreement between the two versions for the same participant, against two
#      different participants rating the identical stimulus.
#   3. Variance decomposition: a participant x verb x sentence type cell term
#      added to the model, on the rating scale and on a latent logit scale.
#   4. Agent-gender effects, including the gender x S_Type x affectedness term
#      that would move H1b, and which verbs differ most between versions.
#   5. Order and fatigue: first against second exposure, first against second
#      half of the session.
#   6. Information loss: the same model fitted to the full data and to
#      one-version datasets built three ways (A: verbs split 50/50 into
#      woman-agent and man-agent sets per participant, complementary lists;
#      B: version drawn independently per cell; C: first exposure of each cell),
#      with an ordinal (clmm) check, a projection to N = 80 by participant
#      bootstrap, and a design-based bootstrap that does not rely on
#      variance-component estimates.
#   7. The SE ratio each design should have, computed from the variance
#      components on the rating scale and on the latent scale of the pilot brms
#      model, which is less noisy than any single refit.
#
# RANDOM-EFFECTS STRUCTURE
#   The preregistered structure is (1 + S_Type * Semantics_scaled | Participant)
#   + (1 + S_Type | Verb). A structure is only usable here if it is non-singular
#   on the full data AND on a one-version dataset, because a boundary estimate
#   in one of the two fits would move the SE ratio for reasons unrelated to the
#   design. The script walks a fixed ladder (maximal, then a zero-correlation
#   participant block with participant slopes estimated at zero pruned, then
#   the same with a zero-correlation verb block) and uses the first rung that
#   passes for EVERY fit in that language. The by-verb sentence-type slopes are
#   never pruned, because they set the H1b precision floor. Where the selected
#   structure is not the maximal one, the maximal structure is also fitted to
#   the full data and the Scheme A and C datasets as a sensitivity check.
#
# OUTPUTS (outputs/design_summary_pilot/gender_versions_*.csv)
#   structure, lag, structure_selection, agreement, variance,
#   variance_summary, expected_ratio, gender_effects, verb_differences, order,
#   fatigue, subset_fits, order_interactions, information_loss_fits,
#   information_loss_summary,
#   design_bootstrap, power_translation
#   Log: outputs/logs/gender_versions.log
#
# DATA PROTECTION
#   Reads participant-level pilot data and writes only aggregates. Nothing leaves
#   the machine.
#
# COST
#   About 90 minutes of wall time with three workers on a laptop-class machine.
#   The clmm fits dominate: 15-25 minutes each for the full English and Turkish
#   data. --cache=<dir> saves every finished fit so an interrupted run resumes.
#
# USAGE (from design_analysis/)
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" scripts/assess_gender_versions.R
#   Options: --workers=3 --reps-a=10 --reps-b=5 --boot-fits=3 --boot-design=2000
#            --clmm-reps=2 --latent-n=35 --languages=English,Turkish,Norwegian
#            --skip-ordinal --quick --cache=<dir> --progress=<file> --out=<dir>
#            --log=<file> --power=<csv> --dgp-dir=<dir>
#            --subsample-participants=<k> (smoke tests only)
# ---------------------------------------------------------------------------

options(stringsAsFactors = FALSE, width = 200, warn = 1)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
aval <- function(f, d) {
  h <- grep(paste0("^", f, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^", f, "="), "", h[1])
}
QUICK       <- "--quick" %in% args
SKIP_ORD    <- "--skip-ordinal" %in% args
WORKERS     <- as.integer(aval("--workers", 3L))
REPS_A      <- as.integer(aval("--reps-a", if (QUICK) 2L else 10L))
REPS_B      <- as.integer(aval("--reps-b", if (QUICK) 1L else 5L))
BOOT_FITS   <- as.integer(aval("--boot-fits", if (QUICK) 1L else 3L))
BOOT_DESIGN <- as.integer(aval("--boot-design", if (QUICK) 100L else 2000L))
CLMM_REPS   <- as.integer(aval("--clmm-reps", if (QUICK) 1L else 2L))
LATENT_N    <- as.integer(aval("--latent-n", if (QUICK) 12L else 35L))
# Smoke-test option only: keep this many participants per language. Results from
# a subsampled run are not the pilot results.
SUBSAMPLE   <- as.integer(aval("--subsample-participants", NA_integer_))
LANGS       <- strsplit(aval("--languages", "English,Turkish,Norwegian"), ",")[[1]]
OUT         <- aval("--out", "outputs/design_summary_pilot")
LOG         <- aval("--log", "outputs/logs/gender_versions.log")
DATA_PATH   <- aval("--data", "data/pilot/claps_pilot_harmonised.csv")
NOR_RAW     <- aval("--nor-raw", "data/pilot/dataNor_Final.csv")
POWER_PATH  <- aval("--power", "outputs/design_summary_pilot/joint_power_pilot.csv")
DGP_DIR     <- aval("--dgp-dir", "outputs/pilot_models")
CACHE       <- aval("--cache", NA_character_)
PROGRESS    <- aval("--progress", file.path(tempdir(), "gender_versions_progress.txt"))
if (is.na(CACHE)) CACHE <- NULL else dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

N_TARGET <- 80L
SEED     <- 20260916L
# Directional BF 10 with a symmetric prior is posterior P(correct sign) >= 10/11,
# i.e. z >= 1.3352 under a normal posterior; used only for the power translation.
Z_BF10   <- qnorm(10 / 11)
# A variance-component SD below this, on the 1-7 scale, is treated as estimated
# at zero when pruning. It is a hundredth of a scale point, far below anything
# that could move a fixed-effect SE.
SD_ZERO  <- 0.01
FOCAL    <- c(h1a = "Semantics_scaled",
              h1b = "S_TypeActive:Semantics_scaled",
              h2  = "S_TypePseudo_Passive:Semantics_scaled")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
log_con <- file(LOG, open = "wt")
sink(log_con, split = TRUE)

say  <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
show <- function(df, n = 60) print(utils::head(as.data.frame(df), n), row.names = FALSE)

write_out <- function(df, name) {
  df  <- as.data.frame(df)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = 6)
  path <- file.path(OUT, paste0("gender_versions_", name, ".csv"))
  write.csv(df, path, row.names = FALSE)
  say("wrote ", path, " (", nrow(df), " rows)")
  invisible(df)
}

T_START <- Sys.time()
rule("CLAPS: what does the reversed agent-gender version buy?")
say("R ", R.version$major, ".", R.version$minor, "; lme4 ", as.character(packageVersion("lme4")),
    if (requireNamespace("ordinal", quietly = TRUE))
      paste0("; ordinal ", as.character(packageVersion("ordinal"))) else "; ordinal NOT installed")
say("workers = ", WORKERS, "; Scheme A reps = ", REPS_A, "; Scheme B reps = ", REPS_B,
    "; bootstrap fits = ", BOOT_FITS, "; design bootstrap = ", BOOT_DESIGN,
    "; clmm reps = ", CLMM_REPS, "; seed = ", SEED)

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------

raw_all <- read.csv(DATA_PATH)
need <- c("Participant", "Language", "Verb", "Item", "S_Type", "Semantics", "Response", "trial")
if (length(setdiff(need, names(raw_all)))) {
  stop("[data] missing columns: ", paste(setdiff(need, names(raw_all)), collapse = ", "))
}
LANGS <- intersect(LANGS, unique(raw_all$Language))

# The agent is the sentence-initial noun in a Norwegian active and the noun after
# "av" (the last word) in an analytical passive. Only these two sentence types
# survive harmonisation, so the rule covers every row that reaches the models.
agent_from_sentence <- function(sentence, s_type) {
  w     <- strsplit(tolower(trimws(sentence)), "\\s+")
  first <- vapply(w, function(x) x[1], "")
  last  <- vapply(w, function(x) gsub("[[:punct:]]", "", x[length(x)]), "")
  noun  <- ifelse(s_type == "Active", first, last)
  unname(c(mannen = "Man", kvinnen = "Woman")[noun])
}

LABEL_NOTES <- list()

prepare_language <- function(L) {
  d <- raw_all[raw_all$Language == L, ]
  if (!is.na(SUBSAMPLE)) {
    set.seed(SEED)
    d <- d[d$Participant %in% sample(unique(d$Participant), min(SUBSAMPLE, n_distinct(d$Participant))), ]
    say("SMOKE TEST: ", L, " subsampled to ", n_distinct(d$Participant), " participants")
  }
  d$Gender_label <- sub(".*_", "", d$Item)
  d$Gender <- d$Gender_label
  d$Cell <- paste(d$Participant, d$Verb, d$S_Type, sep = ":")

  lab_ok   <- function(g, cell) ave(g, cell, FUN = function(x)
    if (length(x) == 2 && all(sort(x) == c("Man", "Woman"))) "ok" else "bad") == "ok"
  anomalous <- !lab_ok(d$Gender_label, d$Cell)
  note <- list(language = L, anomalous_rows = sum(anomalous),
               anomalous_verb_stype = paste(unique(paste(d$Verb[anomalous], d$S_Type[anomalous], sep = " x ")),
                                            collapse = "; "),
               rows_relabelled = 0L, resolution = "none needed", parse_mismatch_outside = NA,
               parse_unresolved = NA)

  # The raw Norwegian export keeps the sentence text, which names the agent
  # unambiguously. Relabelling from it is preferred to exclusion because it keeps
  # the anomalous cell in the balanced design; it is applied only if the parse
  # agrees with the label on every other row, so a wrong parsing rule cannot
  # silently rewrite good labels.
  if (L == "Norwegian" && any(anomalous) && file.exists(NOR_RAW)) {
    nr <- read.csv(NOR_RAW, fileEncoding = "UTF-8")
    m  <- match(paste(d$Participant, d$trial), paste(nr$id, nr$trial))
    ag <- agent_from_sentence(nr$sentence[m], d$S_Type)
    mismatch <- !is.na(ag) & ag != d$Gender_label
    note$parse_unresolved       <- sum(is.na(ag))
    note$parse_mismatch_outside <- sum(mismatch & !anomalous)
    if (note$parse_unresolved == 0 && note$parse_mismatch_outside == 0) {
      cand <- d$Gender_label
      cand[anomalous] <- ag[anomalous]
      if (all(lab_ok(cand, d$Cell))) {
        note$rows_relabelled <- sum(cand != d$Gender_label)
        note$resolution <- paste0("relabelled from raw sentence text (",
                                  paste(unique(paste(d$Verb, d$S_Type, sep = " x ")[cand != d$Gender_label]),
                                        collapse = "; "), ", agent 'Mannen')")
        d$Gender <- cand
      }
    }
  }
  if (any(anomalous) && note$rows_relabelled == 0) {
    note$resolution <- "anomalous cells excluded from gender-specific analyses"
  }
  d$gender_ok <- lab_ok(d$Gender, d$Cell)
  LABEL_NOTES[[L]] <<- note

  lv <- intersect(c("Passive", "Active", "Pseudo_Passive"), unique(d$S_Type))
  d$S_Type <- factor(d$S_Type, levels = lv)

  # Gelman scaling at VERB level, computed once on the full verb set and then
  # held fixed for every subset. Rescaling within a subset would change the units
  # of the coefficients between the full and the one-version fits and contaminate
  # the SE ratio.
  vs <- d[!duplicated(d$Verb), c("Verb", "Semantics")]
  if (any(tapply(d$Semantics, d$Verb, function(x) length(unique(x))) != 1)) {
    stop("[data] ", L, ": Semantics is not constant within verb")
  }
  d$Semantics_scaled <- (d$Semantics - mean(vs$Semantics)) / (2 * sd(vs$Semantics))

  # Numeric dummies let the random part be written with uncorrelated terms while
  # the fixed part keeps the preregistered factor coding and coefficient names.
  d$Act    <- as.numeric(d$S_Type == "Active")
  d$PP     <- as.numeric(d$S_Type == "Pseudo_Passive")
  d$Sem    <- d$Semantics_scaled
  d$ActSem <- d$Act * d$Sem
  d$PPSem  <- d$PP * d$Sem
  d$Gender_c <- ifelse(d$Gender == "Woman", 0.5, -0.5)
  d$PV <- paste(d$Participant, d$Verb, sep = ":")
  d$PS <- paste(d$Participant, d$S_Type, sep = ":")
  d$VS <- paste(d$Verb, d$S_Type, sep = ":")

  d <- d[order(d$Participant, d$trial), ]
  d$pos       <- ave(d$trial, d$Participant, FUN = seq_along)
  d$exposure  <- ave(d$trial, d$Cell, FUN = function(t) rank(t, ties.method = "first"))
  d$lag       <- ave(d$trial, d$Cell, FUN = function(t) if (length(t) == 2) abs(diff(t)) else NA_real_)
  d$lag_pos   <- ave(d$pos, d$Cell, FUN = function(t) if (length(t) == 2) abs(diff(t)) else NA_real_)
  d$max_trial <- ave(d$trial, d$Participant, FUN = max)
  # Session halves use the RAW trial number, which counts fillers and (in
  # Norwegian) the excluded synthetic passives, because fatigue follows time in
  # the session and not position among the analysed rows.
  d$half      <- ifelse(d$trial <= d$max_trial / 2, "first", "second")
  # Centred codes, so that the focal coefficients in the order-interaction models
  # stay averages over halves and exposures.
  d$Half_c     <- ifelse(d$half == "second", 0.5, -0.5)
  d$Exposure_c <- ifelse(d$exposure == 2, 0.5, -0.5)
  rownames(d) <- NULL
  d
}

DATA <- lapply(setNames(LANGS, LANGS), prepare_language)

# ---------------------------------------------------------------------------
# 2. Structure: versions per cell, randomisation, lag, label anomaly
# ---------------------------------------------------------------------------

rule("2. STRUCTURE")

structure_rows <- list()
lag_rows <- list()
for (L in LANGS) {
  d <- DATA[[L]]
  set.seed(SEED + 100L * match(L, LANGS) + 1L)
  cells <- d |> group_by(Cell) |> summarise(n = n(), ok = all(gender_ok), .groups = "drop")
  per_p <- table(d$Participant)
  nt    <- tapply(d$trial, d$Participant, max) - as.vector(per_p)

  # Stimulus positions per participant. If every participant saw one fixed order
  # the mean between-participant correlation of positions would be 1; independent
  # randomisation gives about 0.
  d$stim <- paste(d$Verb, d$Gender, d$S_Type)
  pm <- tidyr::pivot_wider(d[d$gender_ok, c("stim", "Participant", "pos")],
                           names_from = Participant, values_from = pos)
  pmat <- as.matrix(pm[, -1])
  cm   <- suppressWarnings(cor(pmat, use = "pairwise.complete.obs", method = "spearman"))
  seqs <- tapply(d$stim, d$Participant, paste, collapse = "|")

  d$quart <- ave(d$pos, d$Participant, FUN = function(p) ceiling(4 * p / max(p)))
  chi_st  <- suppressWarnings(chisq.test(table(d$S_Type, d$quart)))
  first_rows <- d[d$exposure == 1 & d$gender_ok, ]
  wf <- binom.test(sum(first_rows$Gender == "Woman"), nrow(first_rows))
  note <- LABEL_NOTES[[L]]

  structure_rows[[L]] <- data.frame(
    language = L,
    check = c("participants", "verbs", "sentence_types", "rows", "rows_per_participant_min",
              "rows_per_participant_max", "cells", "cells_with_two_rows",
              "anomalous_label_rows", "rows_relabelled_from_sentence_text",
              "sentence_parse_mismatch_outside_anomaly", "cells_gender_ok",
              "raw_trials_per_participant_max", "unanalysed_trials_per_participant_median",
              "distinct_trial_orders", "mean_between_participant_position_spearman",
              "stype_by_session_quartile_chisq_p", "position_semantics_correlation",
              "prop_cells_woman_version_first", "woman_first_binomial_p"),
    value = c(n_distinct(d$Participant), n_distinct(d$Verb), nlevels(d$S_Type), nrow(d),
              min(per_p), max(per_p), nrow(cells), sum(cells$n == 2),
              note$anomalous_rows, note$rows_relabelled,
              ifelse(is.na(note$parse_mismatch_outside), NA, note$parse_mismatch_outside),
              sum(cells$ok), max(d$trial), median(nt),
              length(unique(seqs)), mean(cm[upper.tri(cm)], na.rm = TRUE),
              chi_st$p.value, cor(d$pos, d$Semantics),
              mean(first_rows$Gender == "Woman"), wf$p.value),
    detail = c("", "", paste(levels(d$S_Type), collapse = ", "), "", "", "", "participant x verb x S_Type",
               "", note$anomalous_verb_stype, note$resolution, "rows where parsed agent disagrees with label",
               "cells with exactly one man-agent and one woman-agent row", "raw trial numbers include fillers",
               if (L == "Norwegian") "fillers plus the excluded synthetic passives" else "fillers and other non-analysed trials",
               "out of participants", "0 = independent orders, 1 = one fixed order",
               "S_Type proportions constant across session quarters", "raw Semantics vs position",
               "", "")
  )

  # Lag between the two versions of a cell, against the lag the same trial slots
  # would give under a fresh random order. A lag much shorter than the baseline
  # would mean the versions were presented as near-adjacent pairs, which would
  # inflate agreement and make the second rating partly a memory of the first.
  lag_obs <- d[d$exposure == 2, c("Participant", "S_Type", "lag", "lag_pos")]
  perm <- replicate(20, {
    tr <- ave(d$trial, d$Participant, FUN = function(t) t[sample.int(length(t))])
    abs(tapply(tr, d$Cell, function(t) diff(t))[d$Cell[d$exposure == 2]])
  })
  perm_st <- rep(as.character(lag_obs$S_Type), 20)
  lag_summary <- function(x, src, st) {
    q <- quantile(x, c(0, .05, .1, .25, .5, .75, .9, .95, 1), names = FALSE)
    data.frame(language = L, S_Type = st, source = src, n = length(x), mean = mean(x), sd = sd(x),
               min = q[1], q05 = q[2], q10 = q[3], q25 = q[4], median = q[5], q75 = q[6],
               q90 = q[7], q95 = q[8], max = q[9], prop_adjacent = mean(x == 1),
               prop_le_5 = mean(x <= 5), prop_le_20 = mean(x <= 20))
  }
  lr <- list(lag_summary(lag_obs$lag, "observed_raw_trials", "All"),
             lag_summary(as.vector(perm), "random_order_baseline", "All"),
             lag_summary(lag_obs$lag_pos, "observed_analysed_positions", "All"))
  for (st in levels(d$S_Type)) {
    lr[[length(lr) + 1]] <- lag_summary(lag_obs$lag[lag_obs$S_Type == st], "observed_raw_trials", st)
    lr[[length(lr) + 1]] <- lag_summary(as.vector(perm)[perm_st == st], "random_order_baseline", st)
  }
  lag_rows[[L]] <- do.call(rbind, lr)
}
structure_df <- do.call(rbind, structure_rows)
lag_df <- do.call(rbind, lag_rows)
show(structure_df, 100)
show(lag_df[lag_df$S_Type == "All", ])
write_out(structure_df, "structure")
write_out(lag_df, "lag")

# ---------------------------------------------------------------------------
# 3a. Agreement between the two versions
# ---------------------------------------------------------------------------

rule("3a. AGREEMENT BETWEEN VERSIONS")

agree_metrics <- function(x, y, xc = NULL, yc = NULL) {
  data.frame(n_pairs = length(x), pearson = cor(x, y), spearman = cor(x, y, method = "spearman"),
             exact_agreement = mean(x == y), within_one_point = mean(abs(x - y) <= 1),
             mean_abs_diff = mean(abs(x - y)),
             pearson_double_centred = if (is.null(xc)) NA_real_ else cor(xc, yc))
}

agreement_rows <- list()
for (L in LANGS) {
  d <- DATA[[L]][DATA[[L]]$gender_ok, ]
  set.seed(SEED + 100L * match(L, LANGS) + 2L)
  # Double centring removes each participant's mean and each stimulus's mean,
  # leaving the participant-specific reaction to a particular item. That residual
  # is what a second version can confirm or average away; the raw correlation is
  # inflated by participant and verb main effects that a second rating adds nothing
  # to.
  d$rc <- d$Response - ave(d$Response, d$Participant, d$Gender) -
    ave(d$Response, d$Verb, d$S_Type, d$Gender) + ave(d$Response, d$Gender)
  w <- d |> select(Participant, Verb, S_Type, Gender, Response, rc) |>
    tidyr::pivot_wider(names_from = Gender, values_from = c(Response, rc))
  for (st in c("All", levels(d$S_Type))) {
    ww <- if (st == "All") w else w[w$S_Type == st, ]
    agreement_rows[[length(agreement_rows) + 1]] <- cbind(
      language = L, S_Type = st, comparison = "same_participant_other_version", replicates = 1,
      agree_metrics(ww$Response_Man, ww$Response_Woman, ww$rc_Man, ww$rc_Woman))
  }
  # Baseline: two different participants rating the identical stimulus. A cyclic
  # shift of a random permutation is a derangement, so no participant is paired
  # with themselves.
  parts <- unique(d$Participant)
  key   <- paste(d$Participant, d$Verb, d$S_Type, d$Gender)
  base  <- lapply(seq_len(20), function(b) {
    p <- sample(parts)
    partner <- setNames(c(p[-1], p[1]), p)
    m <- match(paste(partner[as.character(d$Participant)], d$Verb, d$S_Type, d$Gender), key)
    data.frame(S_Type = d$S_Type, x = d$Response, y = d$Response[m], xc = d$rc, yc = d$rc[m], b = b)
  })
  for (st in c("All", levels(d$S_Type))) {
    mets <- do.call(rbind, lapply(base, function(bb) {
      bb <- if (st == "All") bb else bb[bb$S_Type == st, ]
      agree_metrics(bb$x, bb$y, bb$xc, bb$yc)
    }))
    agreement_rows[[length(agreement_rows) + 1]] <- cbind(
      language = L, S_Type = st, comparison = "different_participants_same_version", replicates = 20,
      as.data.frame(lapply(mets, mean)))
  }
}
agreement_df <- do.call(rbind, agreement_rows)
show(agreement_df)
write_out(agreement_df, "agreement")

# ---------------------------------------------------------------------------
# 3b. Order and fatigue descriptives
# ---------------------------------------------------------------------------

rule("3b. ORDER AND FATIGUE (DESCRIPTIVE)")

order_rows <- list()
fatigue_rows <- list()
for (L in LANGS) {
  d <- DATA[[L]]
  pr <- d |> filter(gender_ok) |> group_by(Participant, Verb, S_Type, Cell) |>
    summarise(first = Response[exposure == 1], second = Response[exposure == 2],
              first_gender = Gender[exposure == 1], lag = lag[1], .groups = "drop") |>
    mutate(shift = second - first)
  lag_cut <- quantile(pr$lag, c(1 / 3, 2 / 3))
  pr$lag_tertile <- cut(pr$lag, c(-Inf, lag_cut, Inf), labels = c("short", "medium", "long"))
  # Standard errors are clustered on participant (means of participant means),
  # because the same person contributes every cell.
  clus <- function(x, g) { m <- tapply(x, g, mean); c(mean(m), sd(m) / sqrt(length(m))) }
  one <- function(p, st, subset) {
    s   <- clus(p$shift, p$Participant)
    wfp <- p[p$first_gender == "Woman", ]; mfp <- p[p$first_gender == "Man", ]
    swf <- if (nrow(wfp)) clus(wfp$shift, wfp$Participant) else c(NA, NA)
    smf <- if (nrow(mfp)) clus(mfp$shift, mfp$Participant) else c(NA, NA)
    data.frame(language = L, S_Type = st, subset = subset, n_cells = nrow(p),
               prop_woman_first = mean(p$first_gender == "Woman"),
               woman_first_binomial_p = binom.test(sum(p$first_gender == "Woman"), nrow(p))$p.value,
               mean_first = mean(p$first), mean_second = mean(p$second),
               mean_shift_second_minus_first = s[1], se_shift = s[2],
               shift_if_woman_first = swf[1], shift_if_man_first = smf[1],
               # With both orders present, half their sum is the pure repetition
               # effect and half their difference is the woman-minus-man version
               # difference, each free of the other.
               repetition_effect = (swf[1] + smf[1]) / 2,
               woman_minus_man = (smf[1] - swf[1]) / 2,
               prop_identical = mean(p$shift == 0),
               sd_first = sd(p$first), sd_second = sd(p$second),
               prop_7_first = mean(p$first == 7), prop_7_second = mean(p$second == 7))
  }
  for (st in c("All", levels(d$S_Type))) {
    p <- if (st == "All") pr else pr[pr$S_Type == st, ]
    order_rows[[length(order_rows) + 1]] <- one(p, st, "all_cells")
  }
  for (lt in levels(pr$lag_tertile)) {
    order_rows[[length(order_rows) + 1]] <- one(pr[pr$lag_tertile == lt, ], "All",
                                                 paste0("lag_", lt, "_", paste(range(pr$lag[pr$lag_tertile == lt]), collapse = "-")))
  }

  # Straightlining shows up as long runs of one response and near-zero spread
  # within a participant; comparing session halves and quarters tells whether it
  # grows with time on task, which is the cost a 40-minute session would add.
  d$quart <- ave(d$trial, d$Participant, FUN = function(t) ceiling(4 * t / max(t)))
  seg_stats <- function(dd, label, contiguous = TRUE) {
    by_p <- dd |> arrange(Participant, trial) |> group_by(Participant) |>
      summarise(sd = sd(Response), same_next = mean(diff(Response) == 0),
                longest_run = max(rle(Response)$lengths),
                modal_share = max(table(Response)) / n(), .groups = "drop")
    # Run-based measures only mean something for a stretch of consecutive trials.
    if (!contiguous) { by_p$same_next <- NA_real_; by_p$longest_run <- NA_real_ }
    data.frame(language = L, segment = label, n_rows = nrow(dd),
               trial_range = paste(range(dd$trial), collapse = "-"),
               mean_rating = mean(dd$Response), sd_rating = sd(dd$Response),
               prop_7 = mean(dd$Response == 7),
               mean_within_participant_sd = mean(by_p$sd),
               prop_consecutive_identical = mean(by_p$same_next),
               mean_longest_run = mean(by_p$longest_run), max_longest_run = max(by_p$longest_run),
               n_participants_sd_below_0.5 = sum(by_p$sd < 0.5),
               n_participants_modal_share_above_0.9 = sum(by_p$modal_share > 0.9),
               n_participants = nrow(by_p))
  }
  fatigue_rows[[length(fatigue_rows) + 1]] <- seg_stats(d, "whole_session")
  for (h in c("first", "second")) fatigue_rows[[length(fatigue_rows) + 1]] <- seg_stats(d[d$half == h, ], paste0("half_", h))
  for (q in 1:4) fatigue_rows[[length(fatigue_rows) + 1]] <- seg_stats(d[d$quart == q, ], paste0("quarter_", q))
  for (e in 1:2) fatigue_rows[[length(fatigue_rows) + 1]] <- seg_stats(d[d$exposure == e, ], paste0("exposure_", e),
                                                                       contiguous = FALSE)
}
order_df <- do.call(rbind, order_rows)
fatigue_df <- do.call(rbind, fatigue_rows)
show(order_df)
show(fatigue_df)
write_out(order_df, "order")
write_out(fatigue_df, "fatigue")

# ---------------------------------------------------------------------------
# 3c. Verb-level version differences (descriptive, paired)
# ---------------------------------------------------------------------------

rule("3c. VERB-LEVEL VERSION DIFFERENCES")

verb_rows <- list()
het_rows <- list()
for (L in LANGS) {
  d <- DATA[[L]][DATA[[L]]$gender_ok, ]
  w <- d |> select(Participant, Verb, S_Type, Gender, Response) |>
    tidyr::pivot_wider(names_from = Gender, values_from = Response) |>
    mutate(diff = Woman - Man)
  summ <- function(g) {
    g |> summarise(n = n(), mean_man = mean(Man), mean_woman = mean(Woman),
                   woman_minus_man = mean(diff), se = sd(diff) / sqrt(n()), .groups = "drop") |>
      mutate(t = woman_minus_man / se, p = 2 * pt(-abs(t), n - 1))
  }
  vs <- summ(group_by(w, Verb, S_Type)) |> mutate(S_Type = as.character(S_Type))
  vall <- summ(group_by(w, Verb)) |> mutate(S_Type = "All")
  # The version difference that could add noise to H1b is its Active-minus-Passive
  # contrast within a verb, because H1b is built from that contrast.
  dd <- w |> filter(S_Type %in% c("Active", "Passive")) |>
    select(Participant, Verb, S_Type, diff) |>
    tidyr::pivot_wider(names_from = S_Type, values_from = diff) |>
    mutate(did = Active - Passive) |> group_by(Verb) |>
    summarise(n = n(), woman_minus_man = mean(did), se = sd(did) / sqrt(n()), .groups = "drop") |>
    mutate(S_Type = "Active_minus_Passive", mean_man = NA_real_, mean_woman = NA_real_,
           t = woman_minus_man / se, p = 2 * pt(-abs(t), n - 1))
  vv <- bind_rows(vs, vall, dd) |> group_by(S_Type) |>
    mutate(p_bh = p.adjust(p, "BH"), rank_abs_diff = rank(-abs(woman_minus_man), ties.method = "first")) |>
    ungroup() |> mutate(language = L, .before = 1) |> arrange(S_Type, rank_abs_diff)
  verb_rows[[L]] <- vv
  # Heterogeneity of verb-level differences beyond sampling error, by method of
  # moments (DerSimonian-Laird). tau is directly comparable with the by-verb SD of
  # the sentence-type effects in the model, which sets the H1b precision floor.
  for (st in unique(vv$S_Type)) {
    x <- vv[vv$S_Type == st, ]
    # A verb on which every participant gave both versions the same rating has a
    # zero SE; flooring it at the smallest positive SE keeps its weight finite.
    wts <- 1 / pmax(x$se, min(x$se[x$se > 0]))^2
    mbar <- sum(wts * x$woman_minus_man) / sum(wts)
    Q <- sum(wts * (x$woman_minus_man - mbar)^2)
    k <- nrow(x)
    tau2 <- max(0, (Q - (k - 1)) / (sum(wts) - sum(wts^2) / sum(wts)))
    het_rows[[length(het_rows) + 1]] <- data.frame(
      language = L, S_Type = st, k_verbs = k, weighted_mean_diff = mbar, sd_of_verb_diffs = sd(x$woman_minus_man),
      mean_sampling_se = sqrt(mean(x$se^2)), Q = Q, Q_p = pchisq(Q, k - 1, lower.tail = FALSE),
      tau_verb_diff = sqrt(tau2), n_p_bh_below_05 = sum(x$p_bh < .05, na.rm = TRUE),
      n_p_below_05 = sum(x$p < .05, na.rm = TRUE),
      top_verbs = paste(head(sprintf("%s (%+.2f)", x$Verb, x$woman_minus_man), 5), collapse = "; "))
  }
}
verb_df <- bind_rows(verb_rows)
het_df <- bind_rows(het_rows)
show(het_df)
write_out(verb_df, "verb_differences")

# ---------------------------------------------------------------------------
# 4. Model machinery shared by the master and the workers
# ---------------------------------------------------------------------------

p_terms_all <- function(L) if (L == "Norwegian") c("1", "Act", "Sem", "ActSem") else
  c("1", "Act", "PP", "Sem", "ActSem", "PPSem")
v_terms_all <- function(L) if (L == "Norwegian") c("1", "Act") else c("1", "Act", "PP")

re_block <- function(terms, group, corr) {
  if (!length(terms)) return(character(0))
  if (corr && length(terms) > 1) {
    lhs <- if ("1" %in% terms) paste(terms, collapse = " + ") else paste(c("0", terms), collapse = " + ")
    return(sprintf("(%s | %s)", lhs, group))
  }
  unname(vapply(terms, function(t) if (t == "1") sprintf("(1 | %s)", group) else
    sprintf("(0 + %s | %s)", t, group), ""))
}

structure_label <- function(st) {
  paste0("P[", if (st$p_corr) "corr" else "zcp", ":", paste(st$p_terms, collapse = ","), "] ",
         "V[", if (st$v_corr) "corr" else "zcp", ":", paste(st$v_terms, collapse = ","), "]")
}

build_formula <- function(st, fixed = "S_Type * Semantics_scaled", extra = NULL) {
  re <- c(re_block(st$p_terms, "Participant", st$p_corr), re_block(st$v_terms, "Verb", st$v_corr), extra)
  as.formula(paste("Response ~", fixed, "+", paste(re, collapse = " + ")))
}

fit_lmer <- function(formula, data) {
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  # calc.derivs = FALSE skips the finite-difference Hessian of the variance
  # parameters, which would cost more than the fit itself for the maximal models.
  # The optimiser's own convergence code is still recorded.
  fit <- withCallingHandlers(
    lmer(formula, data = data, REML = FALSE, control = lmerControl(calc.derivs = FALSE)),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  cf <- summary(fit)$coefficients
  vc <- as.data.frame(VarCorr(fit))
  vc$grp <- sub("\\.[0-9]+$", "", vc$grp)
  list(coefs = data.frame(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], z = cf[, 3], row.names = NULL),
       vc = vc, sigma = sigma(fit), singular = isSingular(fit),
       opt_conv = fit@optinfo$conv$opt, warnings = paste(unique(warns), collapse = " | "),
       nobs = nobs(fit), n_participants = length(unique(data$Participant)),
       secs = proc.time()[["elapsed"]] - t0)
}

fit_clmm <- function(formula, data) {
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  data$R <- factor(data$Response, levels = sort(unique(data$Response)), ordered = TRUE)
  fit <- withCallingHandlers(
    ordinal::clmm(formula, data = data),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  cf <- summary(fit)$coefficients
  cf <- cf[!grepl("\\|", rownames(cf)), , drop = FALSE]
  vcl <- ordinal::VarCorr(fit)
  vc <- data.frame(grp = names(vcl), var1 = "(Intercept)", var2 = NA_character_,
                   vcov = vapply(vcl, function(m) m[1, 1], 0),
                   sdcor = vapply(vcl, function(m) sqrt(m[1, 1]), 0), row.names = NULL)
  list(coefs = data.frame(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], z = cf[, 3], row.names = NULL),
       vc = vc, sigma = pi / sqrt(3), singular = any(vc$sdcor < 1e-3),
       opt_conv = fit$optRes$convergence, warnings = paste(unique(warns), collapse = " | "),
       max_grad = max(abs(fit$gradient)),
       nobs = nrow(data), n_participants = length(unique(data$Participant)),
       secs = proc.time()[["elapsed"]] - t0)
}

progress <- function(path, ...) {
  if (!is.null(path)) try(cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "",
                              file = path, append = TRUE), silent = TRUE)
}

# Walk the ladder on the full data and on one Scheme A dataset together. A rung
# passes only if both fits are non-singular, so whichever structure is chosen can
# be applied unchanged to the full and the one-version data.
#
# The by-verb sentence-type slopes are never pruned. They carry the between-verb
# variation that sets the H1b precision floor (Barr et al., 2013,
# doi:10.1016/j.jml.2012.11.001), so dropping one because a noisier one-version
# fit put it at zero would shrink both SEs for a reason unrelated to the design.
# Participant terms estimated at zero are pruned, following Matuschek et al.
# (2017, doi:10.1016/j.jml.2017.01.001); removing a zero-variance term leaves that
# fit unchanged. If no rung passes on both datasets, the richest rung that is
# non-singular on the full data is used and the table says so.
select_structure <- function(L, probe_rows, progress_path = NULL) {
  d <- DATA[[L]]
  probe <- d[probe_rows, ]
  hist <- list()
  fallback <- NULL
  try_rung <- function(st) {
    step <- length(hist) + 1L
    f1 <- fit_lmer(build_formula(st), d)
    f2 <- if (!f1$singular) fit_lmer(build_formula(st), probe) else NULL
    if (!f1$singular && is.null(fallback)) fallback <<- st
    foc <- function(f, t) if (is.null(f)) NA_real_ else f$coefs$se[f$coefs$term == t][1]
    hist[[step]] <<- data.frame(
      language = L, step = step, structure = structure_label(st),
      formula = paste(deparse(build_formula(st), width.cutoff = 500), collapse = ""),
      full_singular = f1$singular, probe_singular = if (is.null(f2)) NA else f2$singular,
      se_h1b_full = foc(f1, FOCAL[["h1b"]]), se_h1b_probe = foc(f2, FOCAL[["h1b"]]),
      se_h1a_full = foc(f1, FOCAL[["h1a"]]), se_h1a_probe = foc(f2, FOCAL[["h1a"]]),
      secs = f1$secs + if (is.null(f2)) 0 else f2$secs, selected = FALSE, note = "")
    progress(progress_path, L, " selection step ", step, ": ", structure_label(st),
             " full_singular=", f1$singular, " probe_singular=", if (is.null(f2)) NA else f2$singular)
    list(ok = !f1$singular && !is.null(f2) && !f2$singular, fits = Filter(Negate(is.null), list(f1, f2)))
  }
  zero_participant_terms <- function(fits) {
    unique(unlist(lapply(fits, function(f) {
      v <- f$vc[is.na(f$vc$var2) & f$vc$grp == "Participant" & f$vc$var1 != "(Intercept)" &
                  f$vc$sdcor < SD_ZERO, ]
      v$var1
    })))
  }
  done <- function(st, note) {
    h <- do.call(rbind, hist)
    k <- max(which(h$structure == structure_label(st)))
    h$selected[k] <- TRUE
    h$note[k] <- note
    list(structure = st, history = h)
  }
  st <- list(p_terms = p_terms_all(L), p_corr = TRUE, v_terms = v_terms_all(L), v_corr = TRUE)
  if (try_rung(st)$ok) return(done(st, "non-singular on full and probe"))
  st$p_corr <- FALSE
  for (v_corr in c(TRUE, FALSE)) {
    st$v_corr <- v_corr
    repeat {
      r <- try_rung(st)
      if (r$ok) return(done(st, "non-singular on full and probe"))
      nz <- zero_participant_terms(r$fits)
      if (!length(nz)) break
      st$p_terms <- setdiff(st$p_terms, nz)
    }
  }
  if (!is.null(fallback)) return(done(fallback, "no rung non-singular on both; richest non-singular on full data"))
  done(st, "singular on the full data at every rung; last rung used")
}

task_key <- function(task) {
  paste(task$id, if (!is.null(task$structure)) structure_label(task$structure) else "",
        task$fixed, paste(task$extra, collapse = "+"), task$clmm_formula, length(task$rows),
        sum(as.numeric(task$rows)), sum(as.numeric(task$rows) * seq_along(task$rows)) %% 1e9,
        if (is.null(task$pid)) "" else sum(as.numeric(task$pid) * seq_along(task$pid)) %% 1e9, sep = "|")
}

# Each finished fit is saved under the cache directory with the key that produced
# it, so a run interrupted part-way (another process on a shared machine can kill
# R sessions) resumes without refitting. A cached result is reused only if its key
# matches exactly, so changing a structure, subset or seed forces a refit. The
# cache holds coefficients and variance components only, never data.
run_task <- function(task, progress_path = NULL, cache_dir = NULL) {
  key <- task_key(task)
  cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_]", "_", task$id), ".rds"))
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$key, key) && identical(cached$status, "ok")) {
      progress(progress_path, task$id, " cached")
      return(cached)
    }
  }
  out <- tryCatch({
    d <- DATA[[task$language]][task$rows, ]
    if (!is.null(task$pid)) {
      d$Participant <- task$pid
      d$Cell <- paste(d$Participant, d$Verb, d$S_Type, sep = ":")
      d$PV <- paste(d$Participant, d$Verb, sep = ":")
      d$PS <- paste(d$Participant, d$S_Type, sep = ":")
    }
    if (task$engine == "lmer") {
      res <- fit_lmer(build_formula(task$structure, task$fixed, task$extra), d)
    } else {
      res <- fit_clmm(as.formula(task$clmm_formula), d)
    }
    res$status <- "ok"
    res
  }, error = function(e) list(status = paste("error:", conditionMessage(e))))
  progress(progress_path, task$id, " ", out$status,
           if (!is.null(out$secs)) sprintf(" %.0fs", out$secs) else "")
  out$task <- task[setdiff(names(task), c("rows", "pid", "structure"))]
  out$structure_label <- if (!is.null(task$structure)) structure_label(task$structure) else task$clmm_formula
  out$key <- key
  if (!is.null(cache_file) && identical(out$status, "ok")) try(saveRDS(out, cache_file), silent = TRUE)
  out
}

# ---------------------------------------------------------------------------
# 5. One-version datasets
# ---------------------------------------------------------------------------

# Scheme A: two complementary lists. Verbs are split 50/50 into sets V1 and V2;
# list 1 sees V1 with a woman agent and V2 with a man agent, list 2 the reverse,
# and all sentence types of a verb share the version. Every verb x version is
# then rated by half the participants, which is the counterbalancing a 20-minute
# study would use.
scheme_a_keep <- function(participant, verb, gender, gender_ok, exposure, seed = NULL) {
  # seed = NULL draws from the current stream, which the bootstrap needs so that
  # it does not reset the generator on every replicate.
  if (!is.null(seed)) set.seed(seed)
  verbs <- sort(unique(verb))
  parts <- sort(unique(participant))
  v1 <- sample(verbs, floor(length(verbs) / 2))
  p1 <- sample(parts, floor(length(parts) / 2))
  woman_here <- ifelse(participant %in% p1, verb %in% v1, !(verb %in% v1))
  # A cell whose labels could not be resolved contributes its first exposure, so
  # that the one-version dataset still covers exactly the cells of the full data.
  ifelse(gender_ok, (gender == "Woman") == woman_here, exposure == 1)
}

scheme_b_rows <- function(d, seed) {
  set.seed(seed)
  u <- runif(nrow(d))
  o <- order(d$Cell, u)
  sort(o[!duplicated(d$Cell[o])])
}

# ---------------------------------------------------------------------------
# 6. Structure selection (parallel, one worker per language)
# ---------------------------------------------------------------------------

rule("6. RANDOM-EFFECTS STRUCTURE SELECTION")

say("worker progress is appended to ", PROGRESS, if (!is.null(CACHE)) paste0("; fit cache in ", CACHE) else "")
cl <- makeCluster(min(WORKERS, max(length(LANGS), 1L)))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(lme4) }); NULL }))
clusterExport(cl, c("DATA", "FOCAL", "SD_ZERO", "p_terms_all", "v_terms_all", "re_block", "structure_label",
                    "build_formula", "fit_lmer", "fit_clmm", "progress", "select_structure"))

probe_rows <- lapply(setNames(LANGS, LANGS), function(L) {
  d <- DATA[[L]]
  which(scheme_a_keep(d$Participant, d$Verb, d$Gender, d$gender_ok, d$exposure,
                      SEED + 100L * match(L, LANGS) + 1000L + 1L))
})
sel <- parLapply(cl, LANGS, function(L, pr, pp, cache) {
  key <- paste(L, "ladder-v2", length(pr[[L]]), sum(as.numeric(pr[[L]])), SD_ZERO, sep = "|")
  f <- if (!is.null(cache)) file.path(cache, paste0("selection_", L, ".rds"))
  if (!is.null(f) && file.exists(f)) {
    s <- readRDS(f)
    if (identical(s$key, key)) return(s)
  }
  s <- select_structure(L, pr[[L]], pp)
  s$key <- key
  if (!is.null(f)) saveRDS(s, f)
  s
}, pr = probe_rows, pp = PROGRESS, cache = CACHE)
names(sel) <- LANGS
selection_df <- do.call(rbind, lapply(sel, `[[`, "history"))
show(selection_df)
write_out(selection_df, "structure_selection")
STRUCT <- lapply(sel, `[[`, "structure")
MAXIMAL <- lapply(setNames(LANGS, LANGS), function(L)
  list(p_terms = p_terms_all(L), p_corr = TRUE, v_terms = v_terms_all(L), v_corr = TRUE))
for (L in LANGS) say(L, ": selected ", structure_label(STRUCT[[L]]))
stopCluster(cl)

# ---------------------------------------------------------------------------
# 7. Task list
# ---------------------------------------------------------------------------

rule("7. FITS")

tasks <- list()
add_task <- function(...) {
  t <- list(...)
  t$id <- paste(t$language, t$analysis, t$scheme, t$replicate, t$variant, sep = "/")
  tasks[[length(tasks) + 1]] <<- t
}
same <- function(a, b) identical(structure_label(a), structure_label(b))

design_rows <- list()
for (L in LANGS) {
  d <- DATA[[L]]
  li <- match(L, LANGS)
  S <- STRUCT[[L]]
  M <- MAXIMAL[[L]]
  all_rows <- seq_len(nrow(d))
  variants <- if (same(S, M)) list(selected = S) else list(selected = S, maximal = M)

  for (v in names(variants)) {
    add_task(language = L, analysis = "information_loss", scheme = "full", replicate = 0L, variant = v,
             engine = "lmer", rows = all_rows, structure = variants[[v]], fixed = "S_Type * Semantics_scaled",
             cost = nrow(d) * 1.5)
    add_task(language = L, analysis = "information_loss", scheme = "C_first_exposure", replicate = 0L,
             variant = v, engine = "lmer", rows = which(d$exposure == 1), structure = variants[[v]],
             fixed = "S_Type * Semantics_scaled", cost = nrow(d))
    for (r in seq_len(REPS_A)) {
      rows <- which(scheme_a_keep(d$Participant, d$Verb, d$Gender, d$gender_ok, d$exposure,
                                  SEED + 100L * li + 1000L + r))
      add_task(language = L, analysis = "information_loss", scheme = "A_verb_split_lists", replicate = r,
               variant = v, engine = "lmer", rows = rows, structure = variants[[v]],
               fixed = "S_Type * Semantics_scaled", cost = nrow(d))
    }
  }
  for (r in seq_len(REPS_B)) {
    add_task(language = L, analysis = "information_loss", scheme = "B_random_per_cell", replicate = r,
             variant = "selected", engine = "lmer", rows = scheme_b_rows(d, SEED + 100L * li + 2000L + r),
             structure = S, fixed = "S_Type * Semantics_scaled", cost = nrow(d))
  }

  # Cell-term models: the variance decomposition, and an honest full-data SE that
  # does not treat two ratings of one cell as independent.
  add_task(language = L, analysis = "variance", scheme = "full_cell_term", replicate = 0L, variant = "selected",
           engine = "lmer", rows = all_rows, structure = S, fixed = "S_Type * Semantics_scaled",
           extra = "(1 | Cell)", cost = nrow(d) * 2.5)
  add_task(language = L, analysis = "variance", scheme = "full_cell_and_participant_verb", replicate = 0L,
           variant = "selected", engine = "lmer", rows = all_rows, structure = S,
           fixed = "S_Type * Semantics_scaled", extra = c("(1 | PV)", "(1 | Cell)"), cost = nrow(d) * 3)

  # Order and fatigue subsets, each fitted with the selected structure.
  subsets <- list(second_exposure = which(d$exposure == 2), session_first_half = which(d$half == "first"),
                  session_second_half = which(d$half == "second"),
                  first_half_first_exposure = which(d$half == "first" & d$exposure == 1))
  for (s in names(subsets)) {
    add_task(language = L, analysis = "subset", scheme = s, replicate = 0L, variant = "selected", engine = "lmer",
             rows = subsets[[s]], structure = S, fixed = "S_Type * Semantics_scaled", cost = nrow(d))
  }
  # The subset fits give an effect per half or per exposure, but the two subsets
  # share participants and verbs, so their SEs cannot be combined into a test of
  # the difference. These models test it directly. Trial order was randomised per
  # participant, so both contrasts are randomised within participant.
  order_codes <- c(session_half = "Half_c", exposure = "Exposure_c")
  for (oc in names(order_codes)) {
    add_task(language = L, analysis = "order_interaction", scheme = oc, replicate = 0L, variant = "selected",
             engine = "lmer", rows = all_rows, structure = S,
             fixed = paste("S_Type * Semantics_scaled *", order_codes[[oc]]),
             extra = sprintf("(0 + %s | Participant)", order_codes[[oc]]), cost = nrow(d) * 3)
  }

  # Gender model. Participants and verbs may each prefer one version, so both get
  # a gender slope; the Verb:S_Type slope lets the version difference vary by
  # sentence type within a verb, which is the variation that would add noise to H1b
  # in a one-version design.
  add_task(language = L, analysis = "gender", scheme = "gender_model", replicate = 0L, variant = "selected",
           engine = "lmer", rows = which(d$gender_ok), structure = S,
           fixed = "S_Type * Semantics_scaled * Gender_c",
           extra = c("(0 + Gender_c | Participant)", "(0 + Gender_c | Verb)", "(0 + Gender_c | VS)"),
           cost = nrow(d) * 3)

  # Projection to N = 80: resample participants with replacement, relabel them, and
  # fit the full data and a Scheme A dataset built on the relabelled sample.
  rows_by_p <- split(all_rows, d$Participant)
  for (b in seq_len(BOOT_FITS)) {
    set.seed(SEED + 100L * li + 3000L + b)
    pick <- sample(names(rows_by_p), N_TARGET, replace = TRUE)
    brow <- unlist(rows_by_p[pick], use.names = FALSE)
    bpid <- rep(seq_len(N_TARGET), times = lengths(rows_by_p[pick]))
    keep <- scheme_a_keep(bpid, d$Verb[brow], d$Gender[brow], d$gender_ok[brow], d$exposure[brow],
                          SEED + 100L * li + 3500L + b)
    add_task(language = L, analysis = "information_loss", scheme = "boot80_full", replicate = b,
             variant = "selected", engine = "lmer", rows = brow, pid = bpid, structure = S,
             fixed = "S_Type * Semantics_scaled", cost = length(brow) * 1.5)
    add_task(language = L, analysis = "information_loss", scheme = "boot80_A_verb_split_lists", replicate = b,
             variant = "selected", engine = "lmer", rows = brow[keep], pid = bpid[keep], structure = S,
             fixed = "S_Type * Semantics_scaled", cost = length(brow))
  }

  # Ordinal checks on the latent scale. Timed on this machine, clmm needs 4-5
  # minutes for random intercepts for participant and for verb x sentence type on
  # 15,000 rows, and a correlated slope structure did not finish, so the intercept
  # structure is the richest that fits in reasonable time. The verb x sentence type
  # intercept stands in for the by-verb sentence-type slopes that set the H1b floor.
  # The same structure is used for full and one-version data.
  if (!SKIP_ORD) {
    ORD_FORMULA <- "R ~ S_Type * Semantics_scaled + (1 | Participant) + (1 | VS)"
    add_task(language = L, analysis = "ordinal", scheme = "full", replicate = 0L, variant = "clmm",
             engine = "clmm", rows = all_rows, clmm_formula = ORD_FORMULA, cost = nrow(d) * 20)
    for (r in seq_len(CLMM_REPS)) {
      rows <- which(scheme_a_keep(d$Participant, d$Verb, d$Gender, d$gender_ok, d$exposure,
                                  SEED + 100L * li + 1000L + r))
      add_task(language = L, analysis = "ordinal", scheme = "A_verb_split_lists", replicate = r,
               variant = "clmm", engine = "clmm", rows = rows, clmm_formula = ORD_FORMULA, cost = nrow(d) * 10)
    }

    # Latent-scale cell variance. A participant x verb x sentence type term means
    # 15,000 random effects for English and Turkish; the Norwegian fit on all 57
    # participants took about 10 minutes and the English one would take several
    # times that. Variance components do not depend on the number of participants
    # in expectation, so the decomposition uses a seeded subset of LATENT_N
    # participants in every language. The participant x sentence type intercept is
    # there because participants differ widely in how much they accept each sentence
    # type, and without it that variation would be counted as cell variance. The
    # same structure is fitted with lmer on the same subset, so the latent and
    # rating-scale shares are like for like.
    set.seed(SEED + 100L * li + 4000L)
    sub_p <- sample(unique(d$Participant), min(LATENT_N, n_distinct(d$Participant)))
    sub_rows <- which(d$Participant %in% sub_p)
    add_task(language = L, analysis = "variance_latent", scheme = "clmm_cell_term_subset", replicate = 0L,
             variant = "clmm", engine = "clmm", rows = sub_rows,
             clmm_formula = "R ~ S_Type * Semantics_scaled + (1 | Participant) + (1 | PS) + (1 | VS) + (1 | Cell)",
             cost = length(sub_rows) * 40)
    add_task(language = L, analysis = "variance_latent", scheme = "lmer_cell_term_subset", replicate = 0L,
             variant = "intercepts", engine = "lmer", rows = sub_rows,
             structure = list(p_terms = "1", p_corr = FALSE, v_terms = character(0), v_corr = FALSE),
             fixed = "S_Type * Semantics_scaled", extra = c("(1 | PS)", "(1 | VS)", "(1 | Cell)"),
             cost = length(sub_rows))
  }
}

say(length(tasks), " fits queued")
tasks <- tasks[order(-vapply(tasks, function(t) t$cost, 0))]

cl <- makeCluster(WORKERS)
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(lme4); library(ordinal) }); NULL }))
clusterExport(cl, c("DATA", "FOCAL", "SD_ZERO", "p_terms_all", "v_terms_all", "re_block", "structure_label",
                    "build_formula", "fit_lmer", "fit_clmm", "progress", "run_task", "task_key"))
t_fit <- Sys.time()
# chunk.size = 1 matters. By default parLapplyLB cuts the list into 2 x workers
# contiguous chunks, and because the tasks are sorted by cost every clmm fit would
# land in the first chunk and queue on a single worker.
results <- parLapplyLB(cl, tasks, run_task, progress_path = PROGRESS, cache_dir = CACHE, chunk.size = 1)
stopCluster(cl)
say("fits finished in ", format(round(difftime(Sys.time(), t_fit, units = "mins"), 1)))
failed <- Filter(function(r) r$status != "ok", results)
if (length(failed)) {
  say(length(failed), " fits failed:")
  for (f in failed) say("  ", f$task$id, ": ", f$status)
}

# ---------------------------------------------------------------------------
# 8. Collect fits
# ---------------------------------------------------------------------------

rule("8. RESULTS")

fits_df <- bind_rows(lapply(results, function(r) {
  if (r$status != "ok") {
    return(data.frame(language = r$task$language, analysis = r$task$analysis, scheme = r$task$scheme,
                      replicate = r$task$replicate, variant = r$task$variant, engine = r$task$engine,
                      structure = r$structure_label, status = r$status))
  }
  cbind(data.frame(language = r$task$language, analysis = r$task$analysis, scheme = r$task$scheme,
                   replicate = r$task$replicate, variant = r$task$variant, engine = r$task$engine,
                   structure = r$structure_label, status = r$status, n_participants = r$n_participants,
                   n_rows = r$nobs, singular = r$singular, opt_conv = r$opt_conv,
                   max_grad = if (is.null(r$max_grad)) NA_real_ else r$max_grad,
                   warnings = r$warnings, secs = r$secs), r$coefs)
}))
vc_df <- bind_rows(lapply(results, function(r) {
  if (r$status != "ok") return(NULL)
  cbind(data.frame(language = r$task$language, analysis = r$task$analysis, scheme = r$task$scheme,
                   replicate = r$task$replicate, variant = r$task$variant, engine = r$task$engine),
        r$vc[, c("grp", "var1", "var2", "vcov", "sdcor")],
        residual_sd = r$sigma)
}))

focal_fits <- fits_df |> filter(analysis %in% c("information_loss", "ordinal", "variance"),
                                term %in% FOCAL) |>
  mutate(hypothesis = names(FOCAL)[match(term, FOCAL)])
write_out(focal_fits |> select(language, analysis, scheme, replicate, variant, engine, structure, hypothesis,
                               term, estimate, se, z, n_participants, n_rows, singular, opt_conv, max_grad,
                               warnings, secs), "information_loss_fits")

# SE ratio r = SE_full / SE_half and z ratio z_half / z_full, per replicate, then
# summarised. Replicates are paired with the full fit of the same engine and
# structure; bootstrap replicates are paired with the full fit on the same
# bootstrap sample.
ratio_block <- function(half, full, label_full) {
  j <- merge(half, full, by = c("language", "engine", "variant", "term", "hypothesis", "pair"),
             suffixes = c("_half", "_full"))
  # Ratios are formed per pair BEFORE summarising, because a bootstrap full fit
  # differs from pair to pair and a ratio of means would not be a mean ratio.
  j |> mutate(r = se_full / se_half, zr = z_half / z_full) |>
    group_by(language, engine, variant, hypothesis, term, scheme = scheme_half) |>
    summarise(compared_with = label_full, n_reps = n(), structure = first(structure_half),
              r_mean = mean(r), r_min = min(r), r_max = max(r),
              z_ratio_mean = mean(zr), z_ratio_min = min(zr), z_ratio_max = max(zr),
              est_full_mean = mean(estimate_full), se_full_mean = mean(se_full), z_full_mean = mean(z_full),
              est_half_mean = mean(estimate_half), se_half_mean = mean(se_half), z_half_mean = mean(z_half),
              n_half_singular = sum(singular_half), full_singular = any(singular_full), .groups = "drop")
}
ff <- focal_fits
ff$pair <- ifelse(grepl("^boot80", ff$scheme), ff$replicate, 0L)
full_std  <- ff |> filter(scheme %in% c("full", "boot80_full"), analysis %in% c("information_loss", "ordinal"))
half_sets <- ff |> filter(analysis %in% c("information_loss", "ordinal"), !scheme %in% c("full", "boot80_full"))
full_cell <- ff |> filter(scheme == "full_cell_term") |> mutate(variant = "selected")
summary_df <- bind_rows(
  ratio_block(half_sets, full_std, "full data, same model"),
  ratio_block(half_sets |> filter(engine == "lmer", !grepl("^boot80", scheme)), full_cell,
              "full data with participant x verb x S_Type cell term")
) |> arrange(language, hypothesis, engine, variant, compared_with, scheme)

# Analytic projection to N = 80. SE^2 = a + b / N, where a is the verb-level floor
# (by-verb variance of the relevant sentence-type effect over Sxx) and b collects
# everything that averages over participants. a is common to both designs, so b is
# recovered from each design's pilot SE and rescaled to N = 80.
verb_floor <- function(L) {
  v <- vc_df |> filter(language == L, analysis == "information_loss", scheme == "full", variant == "selected",
                       engine == "lmer", grp == "Verb", is.na(var2))
  x <- unique(DATA[[L]][, c("Verb", "Semantics_scaled")])$Semantics_scaled
  sxx <- sum((x - mean(x))^2)
  get <- function(nm) { k <- v$vcov[v$var1 == nm]; if (length(k)) k[1] else 0 }
  c(h1a = get("(Intercept)") / sxx, h1b = get("Act") / sxx, h2 = get("PP") / sxx, sxx = sxx)
}
FLOOR <- lapply(setNames(LANGS, LANGS), verb_floor)
proj <- summary_df |> filter(engine == "lmer", variant == "selected", compared_with == "full data, same model",
                             scheme == "A_verb_split_lists") |>
  rowwise() |>
  mutate(n_pilot = n_distinct(DATA[[language]]$Participant),
         verb_floor_var = FLOOR[[language]][[hypothesis]],
         b_full = n_pilot * max(se_full_mean^2 - verb_floor_var, 0),
         b_half = n_pilot * max(se_half_mean^2 - verb_floor_var, 0),
         r_analytic_N80 = sqrt((verb_floor_var + b_full / N_TARGET) / (verb_floor_var + b_half / N_TARGET)),
         verb_floor_share_full = verb_floor_var / se_full_mean^2) |>
  ungroup() |>
  select(language, hypothesis, term, n_pilot, verb_floor_var, verb_floor_share_full, r_analytic_N80)
summary_df <- summary_df |> left_join(proj |> select(language, term, r_analytic_N80, verb_floor_share_full) |>
                                        mutate(scheme = "A_verb_split_lists", engine = "lmer", variant = "selected",
                                               compared_with = "full data, same model"),
                                      by = c("language", "term", "scheme", "engine", "variant", "compared_with"))
show(summary_df |> filter(hypothesis %in% c("h1a", "h1b")) |>
       select(language, hypothesis, engine, variant, scheme, compared_with, n_reps, se_full_mean, se_half_mean,
              r_mean, r_min, r_max, z_ratio_mean, n_half_singular, r_analytic_N80), 80)
write_out(summary_df, "information_loss_summary")

# Subset fits (order and fatigue), with the first-exposure fit alongside for
# reference.
subset_df <- fits_df |>
  filter((analysis == "subset") |
         (analysis == "information_loss" & variant == "selected" & scheme %in% c("full", "C_first_exposure")),
         term %in% FOCAL) |>
  mutate(hypothesis = names(FOCAL)[match(term, FOCAL)]) |>
  select(language, scheme, structure, hypothesis, term, estimate, se, z, n_participants, n_rows, singular, secs) |>
  arrange(language, hypothesis, scheme)
show(subset_df |> filter(hypothesis == "h1b"))
write_out(subset_df, "subset_fits")

# Each order term is the second-minus-first difference in the coefficient it
# modifies. Whether that difference weakens or strengthens an effect depends on
# the sign of the effect itself, so the modified coefficient is carried alongside.
oi <- fits_df |> filter(analysis == "order_interaction")
order_int_df <- oi |>
  filter(grepl("Half_c|Exposure_c", term)) |>
  mutate(p_two_sided = 2 * pnorm(-abs(z)),
         modifies = sub(":?(Half_c|Exposure_c)$", "", term),
         modifies = ifelse(modifies == "", "(Intercept)", modifies)) |>
  left_join(oi |> select(language, scheme, modifies = term, estimate_of_modified_term = estimate),
            by = c("language", "scheme", "modifies")) |>
  select(language, scheme, structure, term, modifies, estimate, se, z, p_two_sided,
         estimate_of_modified_term, n_participants, n_rows, singular) |>
  arrange(language, scheme, term)
show(order_int_df)
write_out(order_int_df, "order_interactions")

# ---------------------------------------------------------------------------
# 9. Variance decomposition
# ---------------------------------------------------------------------------

rule("9. VARIANCE DECOMPOSITION")

var_out <- vc_df |> filter(analysis %in% c("variance", "variance_latent", "gender", "ordinal") |
                             (analysis == "information_loss" & scheme == "full"))
write_out(var_out, "variance")

vsum <- list()
for (L in LANGS) {
  for (sc in c("full_cell_term", "full_cell_and_participant_verb", "lmer_cell_term_subset",
               "clmm_cell_term_subset")) {
    v <- vc_df |> filter(language == L, scheme == sc, is.na(var2))
    if (!nrow(v)) next
    latent <- v$engine[1] == "clmm"
    cellv <- v$vcov[v$grp == "Cell"]
    pvv   <- if (any(v$grp == "PV")) v$vcov[v$grp == "PV"] else NA_real_
    # On the latent logit scale of a cumulative model the residual is the standard
    # logistic, variance pi^2 / 3, so the cell share is comparable across languages.
    resv  <- if (latent) pi^2 / 3 else v$vcov[v$grp == "Residual"]
    n_p   <- fits_df$n_participants[fits_df$language == L & fits_df$scheme == sc][1]
    vsum[[length(vsum) + 1]] <- data.frame(
      language = L, scale = if (latent) "latent_logit" else "response_1_to_7", model = sc,
      engine = v$engine[1], n_participants = n_p,
      random_structure = fits_df$structure[fits_df$language == L & fits_df$scheme == sc][1],
      extra_terms = if (sc == "full_cell_and_participant_verb") "(1 | PV) + (1 | Cell)" else
        if (sc == "full_cell_term") "(1 | Cell)" else "(1 | Participant) + (1 | PS) + (1 | VS) + (1 | Cell)",
      cell_variance = cellv, participant_verb_variance = pvv, residual_variance = resv,
      cell_variance_share = cellv / (cellv + resv),
      # Two ratings of a cell share the cell term and differ in the residual, so
      # the correlation of the two versions after all participant and verb effects
      # is the cell share, and averaging them removes half the residual: this is the
      # share of per-rating noise (beyond participant and verb effects) that a
      # second version averages away.
      noise_removed_by_second_rating = (resv / 2) / (cellv + resv))
  }
}
vsum_df <- bind_rows(vsum)
show(vsum_df)
write_out(vsum_df, "variance_summary")

# ---------------------------------------------------------------------------
# 9b. Expected SE ratio from the variance components
# ---------------------------------------------------------------------------

rule("9b. EXPECTED SE RATIO FROM VARIANCE COMPONENTS")

# The replicate fits give an empirical ratio, but every one-version fit re-estimates
# the variance components from half the data, which moves its SE for reasons that
# have nothing to do with the design. With a balanced design and a centred
# verb-level predictor, the sampling variance of each focal term follows from the
# components directly, with k ratings per cell:
#
#   H1a  tau2_verb_intercept / Sxx + sigma2_participant_Sem / N
#          + (sigma2_PV + sigma2_cell + sigma2_residual / k) / (N * Sxx)
#   H1b  tau2_verb_Active / Sxx + sigma2_participant_ActSem / N
#          + 2 * (sigma2_cell + sigma2_residual / k) / (N * Sxx)
#
# (H2 as H1b with the pseudo-passive terms.) The participant x verb term drops out
# of H1b because the active and the passive of a verb share it, and stays in H1a,
# a passive-only slope. Only the residual is divided by k, so the ratio of the
# k = 2 and k = 1 variances is what the second version buys.
#
# Components come from two places. On the rating scale, from the lmer fit that adds
# participant x verb and cell terms to the selected structure. On the latent logit
# scale, from the participant and verb covariance matrices of the pilot brms model
# behind the decision-arm power figures, with the residual fixed at pi^2 / 3 and the
# cell variance from the latent clmm. That clmm has no participant x verb term, so
# its cell variance also contains it; the latent ratio is therefore also given with
# the cell variance set to zero, the most a second version could possibly help.
expected_var <- function(h, N, k, cm) {
  switch(h,
         h1a = cm$tv0 / cm$sxx + cm$pS / N + (cm$pv + cm$cell + cm$res / k) / (N * cm$sxx),
         h1b = cm$tvA / cm$sxx + cm$pA / N + 2 * (cm$cell + cm$res / k) / (N * cm$sxx),
         h2  = cm$tvP / cm$sxx + cm$pP / N + 2 * (cm$cell + cm$res / k) / (N * cm$sxx))
}

expected_rows <- function(L, cm, scale, source) {
  hs <- if (cm$tvP > 0 || cm$pP > 0) c("h1a", "h1b", "h2") else c("h1a", "h1b")
  n_pilot <- n_distinct(DATA[[L]]$Participant)
  bind_rows(lapply(hs, function(h) bind_rows(lapply(unique(c(n_pilot, N_TARGET)), function(N) {
    v2 <- expected_var(h, N, 2, cm)
    v1 <- expected_var(h, N, 1, cm)
    floor_part <- switch(h, h1a = cm$tv0, h1b = cm$tvA, h2 = cm$tvP) / cm$sxx
    part_part  <- switch(h, h1a = cm$pS, h1b = cm$pA, h2 = cm$pP) / N
    data.frame(language = L, scale = scale, source = source, hypothesis = h, term = FOCAL[[h]],
               n_participants = N, se_two_versions = sqrt(v2), se_one_version = sqrt(v1),
               r_expected = sqrt(v2 / v1),
               share_verb_floor_one_version = floor_part / v1,
               share_participant_slope_one_version = part_part / v1,
               share_rating_noise_one_version = (v1 - floor_part - part_part) / v1)
  }))))
}

exp_rows <- list()
for (L in LANGS) {
  x <- unique(DATA[[L]][, c("Verb", "Semantics_scaled")])$Semantics_scaled
  sxx <- sum((x - mean(x))^2)
  v <- vc_df |> filter(language == L, scheme == "full_cell_and_participant_verb", is.na(var2))
  if (nrow(v)) {
    g <- function(grp, term = "(Intercept)") {
      k <- if (grp == "Residual") v$vcov[v$grp == grp] else v$vcov[v$grp == grp & v$var1 == term]
      if (length(k)) k[1] else 0
    }
    cm <- list(tv0 = g("Verb"), tvA = g("Verb", "Act"), tvP = g("Verb", "PP"),
               pS = g("Participant", "Sem"), pA = g("Participant", "ActSem"), pP = g("Participant", "PPSem"),
               pv = g("PV"), cell = g("Cell"), res = g("Residual"), sxx = sxx)
    exp_rows[[length(exp_rows) + 1]] <- expected_rows(L, cm, "response_1_to_7", "lmer selected + PV + cell")
  }
  dgp_path <- file.path(DGP_DIR, paste0("pilot_dgp_v2_pilot_", L, ".rds"))
  lat <- vc_df |> filter(language == L, scheme == "clmm_cell_term_subset", grp == "Cell")
  if (file.exists(dgp_path)) {
    dgp <- readRDS(dgp_path)
    dg  <- function(S, nm) if (nm %in% rownames(S)) S[nm, nm] else 0
    xa  <- dgp$verb_affectedness
    base <- list(tv0 = dg(dgp$Sigma_verb, "Intercept"), tvA = dg(dgp$Sigma_verb, "S_TypeActive"),
                 tvP = dg(dgp$Sigma_verb, "S_TypePseudo_Passive"),
                 pS = dg(dgp$Sigma_part, "Semantics_scaled"), pA = dg(dgp$Sigma_part, "S_TypeActive:Semantics_scaled"),
                 pP = dg(dgp$Sigma_part, "S_TypePseudo_Passive:Semantics_scaled"),
                 pv = 0, res = pi^2 / 3, sxx = sum((xa - mean(xa))^2))
    if (nrow(lat)) {
      exp_rows[[length(exp_rows) + 1]] <- expected_rows(L, c(base, cell = lat$vcov[1]), "latent_logit",
                                                        "pilot brms DGP + clmm cell variance")
    }
    exp_rows[[length(exp_rows) + 1]] <- expected_rows(L, c(base, cell = 0), "latent_logit",
                                                      "pilot brms DGP, cell variance 0 (upper bound on gain)")
  } else {
    say("pilot DGP not found for ", L, " (", dgp_path, "); latent expected ratio skipped")
  }
}
expected_df <- bind_rows(exp_rows)
show(expected_df)
write_out(expected_df, "expected_ratio")

# ---------------------------------------------------------------------------
# 10. Gender effects
# ---------------------------------------------------------------------------

gender_terms <- fits_df |> filter(analysis == "gender", grepl("Gender_c", term)) |>
  select(language, structure, term, estimate, se, z, n_participants, n_rows, singular) |>
  mutate(p_two_sided = 2 * pnorm(-abs(z)))
# Put the three-way term next to the H1b estimate it would perturb: in a design
# that showed only one version to everyone, H1b would move by half the three-way
# term. Under Scheme A counterbalancing the shift cancels in expectation.
h1b_full <- fits_df |> filter(analysis == "information_loss", scheme == "full", variant == "selected",
                              term == FOCAL[["h1b"]]) |> select(language, h1b_estimate = estimate, h1b_se = se)
h1a_full <- fits_df |> filter(analysis == "information_loss", scheme == "full", variant == "selected",
                              term == FOCAL[["h1a"]]) |> select(language, h1a_estimate = estimate, h1a_se = se)
gender_terms <- gender_terms |> left_join(h1b_full, by = "language") |> left_join(h1a_full, by = "language") |>
  mutate(shift_if_one_version_only = estimate / 2,
         shift_relative_to_matching_focal = case_when(
           term == "S_TypeActive:Semantics_scaled:Gender_c" ~ (estimate / 2) / h1b_estimate,
           term == "Semantics_scaled:Gender_c" ~ (estimate / 2) / h1a_estimate,
           TRUE ~ NA_real_))
gvc <- vc_df |> filter(analysis == "gender", grepl("Gender_c", var1), is.na(var2)) |>
  transmute(language, structure = NA_character_, term = paste0("random SD: Gender_c | ", grp), estimate = sdcor)
sel_vc <- vc_df |> filter(analysis == "information_loss", scheme == "full", variant == "selected",
                          grp == "Verb", is.na(var2)) |>
  transmute(language, structure = NA_character_, term = paste0("reference random SD: ", var1, " | Verb"),
            estimate = sdcor)
het_long <- het_df |> transmute(language, structure = NA_character_,
                                term = paste0("descriptive verb heterogeneity tau: ", S_Type),
                                estimate = tau_verb_diff, se = NA_real_, z = NA_real_,
                                p_two_sided = Q_p)
gender_df <- bind_rows(gender_terms, gvc, sel_vc, het_long)
show(gender_df |> select(language, term, estimate, se, z, p_two_sided, shift_relative_to_matching_focal), 80)
write_out(gender_df, "gender_effects")

# ---------------------------------------------------------------------------
# 11. Design-based participant bootstrap
# ---------------------------------------------------------------------------

rule("11. DESIGN-BASED PARTICIPANT BOOTSTRAP")

# The model-based SE ratio depends on how well each fit estimates the variance
# components, and one-version fits estimate participant slopes less well. This
# check uses no variance components for the participant part: the ordinary
# least-squares estimate of each focal term is recomputed on participant
# bootstrap samples, for the full data and for a freshly drawn Scheme A dataset,
# and its spread is the participant-driven sampling SD with verbs held fixed. The
# verb-level floor from the full model is then added to both, because the verb
# sample is the same in either design.
boot_rows <- list()
for (L in LANGS) {
  d  <- DATA[[L]]
  li <- match(L, LANGS)
  X  <- model.matrix(~ S_Type * Semantics_scaled, d)
  y  <- d$Response
  cols <- intersect(FOCAL, colnames(X))
  rows_by_p <- split(seq_len(nrow(d)), d$Participant)
  n_pilot <- length(rows_by_p)
  # Integer cell index within participant, so a bootstrap cell id is a single
  # arithmetic step instead of a string paste on every replicate.
  cidx  <- ave(seq_len(nrow(d)), d$Participant, FUN = function(i) match(paste(d$Verb[i], d$S_Type[i]),
                                                                        unique(paste(d$Verb[i], d$S_Type[i]))))
  max_c <- max(cidx)
  set.seed(SEED + 100L * li + 5000L)
  for (N in unique(c(n_pilot, N_TARGET))) {
    est <- replicate(BOOT_DESIGN, {
      pick <- sample.int(n_pilot, N, replace = TRUE)
      brow <- unlist(rows_by_p[pick], use.names = FALSE)
      bpid <- rep(seq_len(N), times = lengths(rows_by_p[pick]))
      keep <- scheme_a_keep(bpid, d$Verb[brow], d$Gender[brow], d$gender_ok[brow], d$exposure[brow])
      # Scheme B: with exactly two rows per cell, keeping the first or second
      # exposure on an independent coin flip per cell is a random version per cell.
      coin <- runif(N * max_c) < 0.5
      rb   <- (d$exposure[brow] == 1) == coin[(bpid - 1L) * max_c + cidx[brow]]
      full <- .lm.fit(X[brow, , drop = FALSE], y[brow])$coefficients
      half <- .lm.fit(X[brow[keep], , drop = FALSE], y[brow[keep]])$coefficients
      hb   <- .lm.fit(X[brow[rb], , drop = FALSE], y[brow[rb]])$coefficients
      names(full) <- names(half) <- names(hb) <- colnames(X)
      c(full[cols], half[cols], hb[cols])
    })
    est <- matrix(est, nrow = 3 * length(cols))
    k <- length(cols)
    sd_full <- apply(est[1:k, , drop = FALSE], 1, sd)
    sd_half <- apply(est[(k + 1):(2 * k), , drop = FALSE], 1, sd)
    sd_hb   <- apply(est[(2 * k + 1):(3 * k), , drop = FALSE], 1, sd)
    fl <- FLOOR[[L]][names(FOCAL)[match(cols, FOCAL)]]
    boot_rows[[length(boot_rows) + 1]] <- data.frame(
      language = L, n_participants = N, bootstrap_samples = BOOT_DESIGN,
      hypothesis = names(FOCAL)[match(cols, FOCAL)], term = cols,
      sd_participant_full = sd_full, sd_participant_schemeA = sd_half, sd_participant_schemeB = sd_hb,
      r_participant_part_schemeA = sd_full / sd_half, r_participant_part_schemeB = sd_full / sd_hb,
      verb_floor_var = unname(fl),
      r_total_schemeA = sqrt((fl + sd_full^2) / (fl + sd_half^2)),
      r_total_schemeB = sqrt((fl + sd_full^2) / (fl + sd_hb^2)),
      verb_floor_share_full = unname(fl) / (unname(fl) + sd_full^2))
  }
}
boot_df <- bind_rows(boot_rows)
show(boot_df)
write_out(boot_df, "design_bootstrap")

# ---------------------------------------------------------------------------
# 12. Power translation
# ---------------------------------------------------------------------------

rule("12. WHAT THE SECOND VERSION WOULD ADD TO THE REPORTED POWER")

# The decision-arm figures already describe one rating per cell. If a design's
# test statistic is roughly normal with mean mu and unit SD, power at BF 10 is
# Phi(mu - 1.3352); the two-version design multiplies mu by 1 / r. Assurance
# averages over effect-size uncertainty, so this is an approximation that is
# accurate when the power is driven by sampling noise and optimistic about how
# much a precision gain can help when it is driven by uncertainty in the effect.
#
# The primary translation uses the latent-scale expected ratio, because it comes
# from the same pilot model as the power figures. The upper translation uses the
# smallest ratio any method produced at N = 80, i.e. the largest gain the evidence
# allows for the second version.
jp_path <- POWER_PATH
if (file.exists(jp_path)) {
  jp <- read.csv(jp_path) |> filter(mode == "assurance", n_participants == N_TARGET)
  rA <- summary_df |> filter(engine == "lmer", variant == "selected", scheme == "boot80_A_verb_split_lists",
                             compared_with == "full data, same model") |>
    select(language, hypothesis, r_boot80_model = r_mean)
  rP <- summary_df |> filter(engine == "lmer", variant == "selected", scheme == "A_verb_split_lists",
                             compared_with == "full data, same model") |>
    select(language, hypothesis, r_pilot_model = r_mean, r_analytic_N80)
  rD <- boot_df |> filter(n_participants == N_TARGET) |> select(language, hypothesis, r_design_N80 = r_total_schemeA)
  rE <- expected_df |> filter(n_participants == N_TARGET) |>
    mutate(col = case_when(scale == "response_1_to_7" ~ "r_expected_rating_N80",
                           grepl("clmm cell", source) ~ "r_expected_latent_N80",
                           TRUE ~ "r_expected_latent_nocell_N80")) |>
    select(language, hypothesis, col, r_expected) |>
    tidyr::pivot_wider(names_from = col, values_from = r_expected)
  pt <- bind_rows(lapply(c("h1a", "h1b"), function(h) {
    jp |> transmute(language, hypothesis = h, n_participants, reps,
                    power_one_version_decision_arm = .data[[paste0("p_", h)]])
  })) |> left_join(rE, by = c("language", "hypothesis")) |> left_join(rP, by = c("language", "hypothesis")) |>
    left_join(rA, by = c("language", "hypothesis")) |> left_join(rD, by = c("language", "hypothesis"))
  r_cols <- intersect(c("r_expected_latent_N80", "r_expected_latent_nocell_N80", "r_expected_rating_N80",
                        "r_analytic_N80", "r_boot80_model", "r_design_N80"), names(pt))
  pt$r_primary <- if ("r_expected_latent_N80" %in% names(pt)) pt$r_expected_latent_N80 else NA_real_
  pt$r_smallest <- do.call(pmin, c(unname(as.list(pt[r_cols])), na.rm = TRUE))
  pt <- pt |>
    mutate(mu_one_version = qnorm(pmin(pmax(power_one_version_decision_arm, 0.001), 0.999)) + Z_BF10,
           power_two_versions_primary = pnorm(mu_one_version / r_primary - Z_BF10),
           power_two_versions_upper = pnorm(mu_one_version / r_smallest - Z_BF10),
           gain_primary = power_two_versions_primary - power_one_version_decision_arm,
           gain_upper = power_two_versions_upper - power_one_version_decision_arm)
  show(pt)
  write_out(pt, "power_translation")
} else {
  say("joint_power_pilot.csv not found; power translation skipped")
}

rule("DONE")
say("total time ", format(round(difftime(Sys.time(), T_START, units = "mins"), 1)))
sink()
close(log_con)
