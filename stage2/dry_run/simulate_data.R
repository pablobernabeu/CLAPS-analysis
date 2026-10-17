#!/usr/bin/env Rscript
# stage2/dry_run/simulate_data.R
# ---------------------------------------------------------------------------
# Simulated team data for the Stage 2 dry run, in the canonical input format.
#
# Not part of the registered analysis. It writes small data sets that exercise every
# rule of stage2/01-06 on a laptop:
#   - six varieties of five languages in three families, one language in two
#     varieties and one language without a pseudo-passive;
#   - in every variety three participants and two raters who meet an exclusion
#     criterion (one incomplete, one failing the attention checks, one giving the same
#     rating throughout), so that exactly the target number remain;
#   - one variety that has not finished data collection (too few participants);
#   - one variety whose verbs are rated alike, which the outcome-neutral check stops.
# The acceptability ratings follow a cumulative-logit model in which the affectedness
# slope is smaller for actives and pseudo-passives than for passives, so H1 and H2
# hold in the simulated population.
#
# Writes outputs/stage2_dry_run/raw/<variety>_acceptability.csv and _ratings.csv, and
# outputs/stage2_dry_run/languages.csv.
#
# Usage (from design_analysis/):
#   Rscript stage2/dry_run/simulate_data.R
# ---------------------------------------------------------------------------
set.seed(2027)
base <- "outputs/stage2_dry_run"
raw <- file.path(base, "raw")
unlink(base, recursive = TRUE)
dir.create(raw, recursive = TRUE)

N_ACC <- 24L      # analysed participants per variety (the dry-run target)
N_RAT <- 8L       # analysed raters per variety
N_VERB <- 12L
verbs <- sprintf("verb%02d", seq_len(N_VERB))

varieties <- data.frame(
  variety  = c("Simlang A", "Simlang B North", "Simlang B South", "Simlang C",
               "Simlang D", "Simlang E", "Simlang F", "Simlang G"),
  language = c("Simlang A", "Simlang B", "Simlang B", "Simlang C",
               "Simlang D", "Simlang E", "Simlang F", "Simlang G"),
  family   = c("Family 1", "Family 1", "Family 1", "Family 2",
               "Family 3", "Family 2", "Family 3", "Family 1"),
  pseudo   = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  n_acc    = c(N_ACC, N_ACC, N_ACC, N_ACC, N_ACC, N_ACC, N_ACC - 5L, N_ACC),
  flat     = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
utils::write.csv(varieties[, c("variety", "language", "family")],
                 file.path(base, "languages.csv"), row.names = FALSE)

thresholds <- c(-3, -2, -1.2, -0.4, 0.4, 1.4)
rate <- function(eta) {
  # Cumulative logit: P(Y <= k) = logistic(threshold_k - eta)
  p <- stats::plogis(outer(eta, thresholds, function(e, t) t - e))
  probs <- cbind(p[, 1], p[, -1, drop = FALSE] - p[, -ncol(p), drop = FALSE], 1 - p[, ncol(p)])
  apply(probs, 1, function(pr) sample.int(7, 1, prob = pmax(pr, 0)))
}
attention_rows <- function(id, start) {
  data.frame(Participant = id, Trial = start + seq_len(10), Trial_Type = "attention",
             Verb = NA, Agent = NA, S_Type = NA, Attention_Target = rep(c(1L, 7L), 5),
             Response = rep(c(1L, 7L), 5))
}

for (i in seq_len(nrow(varieties))) {
  V <- varieties[i, ]
  true_aff <- if (V$flat) stats::rnorm(N_VERB, 4, 0.1) else stats::rnorm(N_VERB, 4.2, 1.4)
  true_aff <- pmin(pmax(true_aff, 1.3), 6.7)
  z <- (true_aff - mean(true_aff)) / (2 * stats::sd(true_aff))
  types <- if (V$pseudo) c("Active", "Passive", "Pseudo_Passive") else c("Active", "Passive")
  # Language-specific interactions around the population values (per unit of z).
  b_act <- stats::rnorm(1, -0.6, 0.15)
  b_pp <- stats::rnorm(1, -0.5, 0.15)
  verb_int <- stats::rnorm(N_VERB, 0, 0.4)

  # Acceptability: the analysed participants plus three who will be excluded.
  n_rec <- V$n_acc + 3L
  acc <- do.call(rbind, lapply(seq_len(n_rec), function(p) {
    id <- sprintf("P%03d", p)
    agent <- sample(rep(c("Man", "Woman"), length.out = N_VERB))
    d <- expand.grid(Verb = verbs, S_Type = types, stringsAsFactors = FALSE)
    d$Agent <- agent[match(d$Verb, verbs)]
    zz <- z[match(d$Verb, verbs)]
    # Kept off the ceiling, so that no simulated participant meets the 90% criterion
    # by chance and every variety ends with exactly the target number.
    eta <- 0.4 + 0.9 * zz + stats::rnorm(1, 0, 0.4) + verb_int[match(d$Verb, verbs)] +
      ifelse(d$S_Type == "Active", 0.8 + b_act * zz, 0) +
      ifelse(d$S_Type == "Pseudo_Passive", -0.6 + b_pp * zz, 0)
    d$Response <- rate(eta)
    d <- d[sample(nrow(d)), ]
    d$Participant <- id
    d$Trial <- seq_len(nrow(d))
    d$Trial_Type <- "test"
    d$Attention_Target <- NA
    rbind(d[, c("Participant", "Trial", "Trial_Type", "Verb", "Agent", "S_Type",
                "Attention_Target", "Response")], attention_rows(id, nrow(d)))
  }))
  # The three to be excluded: the last three recruited.
  ids <- sprintf("P%03d", n_rec - 2:0)
  acc$Response[acc$Participant == ids[1] & acc$Trial == 3] <- NA                     # incomplete
  acc$Response[acc$Participant == ids[2] & acc$Trial_Type == "attention"][1:3] <- 4L # 7 of 10 right
  acc$Response[acc$Participant == ids[3] & acc$Trial_Type == "test"] <- 6L           # same rating
  acc$Variety <- V$variety
  utils::write.csv(acc[, c("Variety", setdiff(names(acc), "Variety"))],
                   file.path(raw, paste0(gsub(" ", "_", V$variety), "_acceptability.csv")),
                   row.names = FALSE, na = "")

  # Ratings: the analysed raters plus two who will be excluded.
  rat <- do.call(rbind, lapply(seq_len(N_RAT + 2L), function(r) {
    id <- sprintf("R%03d", r)
    d <- expand.grid(Verb = verbs, Agent = c("Man", "Woman"), stringsAsFactors = FALSE)
    d$Response <- pmin(pmax(round(true_aff[match(d$Verb, verbs)] + stats::rnorm(nrow(d), 0, 0.9)), 1), 7)
    d <- d[sample(nrow(d)), ]
    d$Participant <- id
    d$Trial <- seq_len(nrow(d))
    d$Trial_Type <- "test"
    d$Attention_Target <- NA
    d$S_Type <- NA
    rbind(d[, c("Participant", "Trial", "Trial_Type", "Verb", "Agent", "S_Type",
                "Attention_Target", "Response")], attention_rows(id, nrow(d)))
  }))
  rid <- sprintf("R%03d", N_RAT + 1:2)
  rat$Response[rat$Participant == rid[1] & rat$Trial_Type == "attention"][1:5] <- 4L
  rat$Response[rat$Participant == rid[2] & rat$Trial_Type == "test"] <- 5L
  rat$Variety <- V$variety
  rat$S_Type <- NULL
  utils::write.csv(rat[, c("Variety", setdiff(names(rat), "Variety"))],
                   file.path(raw, paste0(gsub(" ", "_", V$variety), "_ratings.csv")),
                   row.names = FALSE, na = "")
}
message("[dry run] simulated ", nrow(varieties), " varieties in ", raw)
