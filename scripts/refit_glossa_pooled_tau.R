#!/usr/bin/env Rscript
# scripts/refit_glossa_pooled_tau.R
# ---------------------------------------------------------------------------
# Re-estimate the between-language SD of the affectedness-by-sentence-type
# interaction from the Ambridge, Arnon & Bekman (2023) five-language data, under
# a correctly specified model.
#
# WHY THIS FIT EXISTS
#   The K-language recruitment rule the PI wants to preregister turns on one
#   number: tau, the SD of the interaction across languages. Every figure in
#   scripts/run_klanguage_design_analysis.R is a function of it. The published
#   value is
#
#       sd(S_TypeActive:Semantics | Language) = 0.45, 95% CI [0.19, 1.12]
#       (All_Sem_Only_Bayes.txt, doi:10.5070/G6011177)
#
#   and it comes from a model with two misspecifications, each of which pushes
#   variance INTO the by-Language block and so inflates exactly this term.
#
#   1. VERBS ARE FORCED TO BE SHARED. The published model fits
#      (1 + S_Type*Semantics | verb) with a single 72-level verb factor spanning
#      all five languages, so a verb's deviation is identical in every language.
#      That structure is rejected by the paper's own data: against a
#      verb-nested-only model with the same number of parameters, nesting wins by
#      dBIC = 2310 on the 37 verbs common to all five languages and 4300 on the
#      full set. The estimated cross-language correlation of verb-level effects is
#      0.22, 95% CI [0.08, 0.46], not 1. Whatever a verb does idiosyncratically in
#      Mandarin but not in Hebrew has nowhere to go in that model except the
#      by-Language slopes.
#
#   2. LANGUAGES SHARE ONE THRESHOLD VECTOR. Languages differ sharply in how they
#      use the response scale: the proportion of ratings in the top category runs
#      from 0.274 (Indonesian) to 0.593 (Hebrew). A cumulative model with one
#      global threshold vector must absorb that as differences in the latent
#      predictor, and the by-Language block is where it lands.
#
#   The published model additionally carries a by-verb slope on Semantics. That
#   slope is not identified: affectedness is a property of the verb and is constant
#   within it, so (Semantics | verb) aliases the by-verb intercept and
#   (S_Type:Semantics | verb) aliases (S_Type | verb). This is the same point
#   R/04_model_formulas.R makes about L5_cross_verbblock_aligned. The corrected
#   arms drop it.
#
# WHAT IS FITTED (one arm per SLURM array task)
#   M0_published   The published structure exactly, refitted under this project's
#                  backend and sampler so that the comparison is like for like.
#   M1_verbsplit   M0 with the verb block corrected: a SHARED by-verb term plus a
#                  LANGUAGE-SPECIFIC one, and no unidentified by-verb affectedness
#                  slope.
#   M2_thresholds  M1 with per-language thresholds, so scale use is modelled where
#                  it belongs instead of inside the by-Language slopes.
#
# WHAT TO DO WITH THE ANSWER
#   Read sd(S_TypeActive:Semantics | Language) off the preferred arm and feed it
#   back into scripts/run_klanguage_design_analysis.R as the Glossa-protocol tau.
#   The arms are ordered a priori: M2 is the only one that models both the verb
#   nesting and the per-language scale use, and the verb nesting is already
#   established on these data by dBIC 2310 and 4300 above. LOO is confirmatory
#   rather than load-bearing here, and it is now computed by a separate step from
#   the persisted draws (--loo, off by default) so that it cannot take the fit
#   down with it.
#   The design analysis is cheap and rerunning it is minutes; this fit is the
#   expensive input. If the corrected tau comes back near 0.25 rather than 0.45,
#   the minimum number of languages reaching 90% falls from the low twenties to
#   about eight, which is the difference between a rule the project can meet and
#   one it cannot.
#
# COST
#   The published fit took 490,250 s elapsed (136 h) with 16 chains under rstan
#   (V6_All_Semantics_Only.Rout). These arms use 4 chains under cmdstanr on
#   76,841 observations, which is the same order as the CLAPS pooled fits at
#   75,528 observations and 29-96 h. Walltime is requested at 10 days.
#
# USAGE
#   Rscript scripts/refit_glossa_pooled_tau.R --arm M1_verbsplit
#   sbatch hpc/submit_glossa_tau_refit.sh          (array 1-3, one arm each)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse); library(brms)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--arm",      default = NULL, type = "character"),
  optparse::make_option("--datadir",  default = "data/glossa2023"),
  optparse::make_option("--outdir",   default = "outputs/glossa_tau_refit"),
  optparse::make_option("--iter",     default = 3000L, type = "integer"),
  optparse::make_option("--warmup",   default = 1500L, type = "integer"),
  optparse::make_option("--chains",   default = 4L,    type = "integer"),
  optparse::make_option("--seed",     default = 20260910L, type = "integer"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE),
  # Off by default. The first run of these arms (job 13106297, Sept 2026) lost
  # 92-160 h of finished sampling per arm because brms::loo() was evaluated
  # inline while building the result list, exceeded the cgroup limit, and was
  # SIGKILLed before anything reached disk. A tryCatch cannot catch that signal.
  # LOO is now a separate step run from the persisted draws, so a failure there
  # costs minutes rather than the fit.
  optparse::make_option("--loo", action = "store_true", default = FALSE)
)))

# cmdstanr takes its default output directory from this R option, NOT from the
# environment. Every submit script in this project exports CMDSTANR_OUTPUT_DIR,
# but nothing ever read it, so the chain CSVs went to tempdir() -- which on ARC
# is a per-job private /tmp created by the auto_tmpdir SPANK plugin and destroyed
# when the job ends. That is why the killed arms left nothing to recover. Wiring
# the option here makes the existing convention do what it always claimed to.
csv_dir <- Sys.getenv("CMDSTANR_OUTPUT_DIR", unset = "")
if (nzchar(csv_dir)) {
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  options(cmdstanr_output_dir = csv_dir)
  message("[glossa-tau] chain CSVs -> ", csv_dir)
} else {
  warning("[glossa-tau] CMDSTANR_OUTPUT_DIR unset; draws will not survive the job.")
}

ARMS <- c("M0_published", "M1_verbsplit", "M2_thresholds")
arm <- opt$arm
if (is.null(arm)) {
  sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (!nzchar(sid)) stop("[glossa-tau] --arm or SLURM_ARRAY_TASK_ID required.")
  arm <- ARMS[as.integer(sid)]
}
if (!arm %in% ARMS) stop("[glossa-tau] unknown arm: ", arm)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(opt$outdir, paste0("glossa_tau_", arm, ".rds"))
if (file.exists(out_file) && !opt$overwrite) {
  message("[glossa-tau] exists, skipping: ", out_file); quit(save = "no")
}

# ---------------------------------------------------------------------------
# Data, assembled exactly as the published script does
# ---------------------------------------------------------------------------
# Every step below mirrors V6_All_Semantics_Only.R so that M0 is a genuine
# reproduction and the corrected arms differ only in the random-effects structure.
# The one thing NOT reproduced is the Balinese 10-point-to-5-point collapse being
# silently rounded: it is kept, because changing it would confound the comparison,
# but it is recorded in the output so the report can say the arms inherit it.

read_lang <- function(f) {
  d <- read.csv(file.path(opt$datadir, f), header = TRUE)
  d$X <- NULL; d$X.1 <- NULL
  d
}
D1 <- read_lang("Balinese_Passives.csv")
D2 <- read_lang("English_Passives.csv")
D3 <- read_lang("Hebrew_Passives.csv")
D4 <- read_lang("Indonesian_Passives.csv")
D5 <- read_lang("Mandarin_Passives.csv")

D3$Passive_or_not <- gsub("Active", "Not_Passive", D3$S_Type)
# Balinese was collected on a 10-point scale; the published analysis maps it onto
# 5 points and rounds. Reproduced, not corrected.
D1$Response <- (5 - 1) * (D1$Response - 1) / (10 - 1) + 1

Data <- rbind(D1, D2, D3, D4, D5)
Data$Response <- round(as.numeric(Data$Response), 0)
Data$S_Type   <- gsub("Non_Canonical_Passive", "Pseudo_Passive", Data$S_Type)
Data$S_Type   <- gsub("Passive_basic",         "Pseudo_Passive", Data$S_Type)
Data$S_Type   <- gsub("Notional_Passive",      "Pseudo_Passive", Data$S_Type)
Data$S_Type   <- gsub("Passive_ka", "Passive", Data$S_Type)
Data$S_Type   <- gsub("Passive_ma", "Passive", Data$S_Type)
Data$S_Type   <- gsub("Passive_a",  "Passive", Data$S_Type)
Data <- subset(Data, S_Type != "BA_Active")

Data$verb        <- as.factor(Data$verb)
Data$Participant <- as.factor(Data$Participant)
Data$Language    <- as.factor(Data$Language)
Data$Semantics   <- as.numeric(Data$Semantics)
Data$S_Type      <- relevel(as.factor(Data$S_Type), ref = "Passive")

# The language-specific verb term. Written explicitly rather than relying on
# brms's `verb:Language` interaction syntax, so that the grouping factor's levels
# are exactly the verb-by-language combinations that occur.
Data$VerbLang <- factor(paste(Data$Language, Data$verb, sep = "_"))

message("[glossa-tau] arm ", arm, " | obs ", nrow(Data),
        " | participants ", nlevels(droplevels(Data$Participant)),
        " | verbs ", nlevels(droplevels(Data$verb)),
        " | verb-by-language ", nlevels(droplevels(Data$VerbLang)),
        " | languages ", nlevels(droplevels(Data$Language)))

# ---------------------------------------------------------------------------
# The three arms
# ---------------------------------------------------------------------------
FORMULAS <- list(
  # As published: one shared verb factor carrying the full fixed structure.
  M0_published = brms::bf(
    Response ~ S_Type * Semantics +
      (1 + S_Type * Semantics | Participant) +
      (1 + S_Type * Semantics | verb) +
      (1 + S_Type * Semantics | Language),
    family = brms::cumulative()),

  # Verb block corrected. The by-verb affectedness slope is dropped as
  # unidentified, and verb deviations are split into a component shared across
  # languages and a component specific to a language.
  M1_verbsplit = brms::bf(
    Response ~ S_Type * Semantics +
      (1 + S_Type * Semantics | Participant) +
      (1 + S_Type | verb) +
      (1 + S_Type | VerbLang) +
      (1 + S_Type * Semantics | Language),
    family = brms::cumulative()),

  # M1 plus per-language thresholds, so that a language which uses the top of the
  # scale differently is modelled as such rather than as a different effect.
  M2_thresholds = brms::bf(
    Response | thres(gr = Language) ~ S_Type * Semantics +
      (1 + S_Type * Semantics | Participant) +
      (1 + S_Type | verb) +
      (1 + S_Type | VerbLang) +
      (1 + S_Type * Semantics | Language),
    family = brms::cumulative())
)

# The published prior, kept for every arm so the comparison isolates structure.
PRIOR <- brms::prior(normal(0, 5), class = "b") +
         brms::prior(normal(0, 5), class = "Intercept")

t0 <- proc.time()[["elapsed"]]
fit <- brms::brm(
  FORMULAS[[arm]], data = Data, prior = PRIOR,
  family = brms::cumulative(),
  backend = "cmdstanr", sample_prior = "yes",
  iter = opt$iter, warmup = opt$warmup, chains = opt$chains,
  cores = opt$chains, seed = opt$seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  silent = 1
)
runtime <- proc.time()[["elapsed"]] - t0
message("[glossa-tau] fitted in ", round(runtime / 3600, 2), " h")

# ---------------------------------------------------------------------------
# The extract. A full brms object for this model is several gigabytes, so only
# what the design analysis and the report need is written.
# ---------------------------------------------------------------------------
draws <- as.data.frame(fit, variable = c("^b_", "^sd_", "^cor_"), regex = TRUE)
tau_col <- grep("^sd_Language__S_TypeActive:Semantics$", colnames(draws), value = TRUE)
beta_col <- grep("^b_S_TypeActive:Semantics$", colnames(draws), value = TRUE)

qs <- function(x) c(median = median(x), mean = mean(x), sd = sd(x),
                    q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))

tau_est  <- if (length(tau_col))  qs(draws[[tau_col]])  else NA
beta_est <- if (length(beta_col)) qs(draws[[beta_col]]) else NA

# The number this fit exists to produce, written before anything else touches the
# fitted object. Everything below is summary and diagnostics, and each of those
# calls allocates on the scale of the draws array. If one of them is killed, the
# arm has still left its answer on disk instead of nothing at all.
tau_file <- file.path(opt$outdir, paste0("glossa_tau_", arm, "_tau.rds"))
tmp <- paste0(tau_file, ".tmp")
saveRDS(list(arm = arm, n_obs = nrow(Data), runtime_h = runtime / 3600,
             tau_language_interaction = tau_est, beta_interaction = beta_est,
             csv_dir = csv_dir, seed = opt$seed), tmp)
file.rename(tmp, tau_file)
message("[glossa-tau] wrote ", tau_file)
if (length(tau_col)) {
  message("[glossa-tau] tau(Language, interaction) = ",
          paste(sprintf("%.3f", tau_est), collapse = " "))
}

result <- list(
  arm      = arm,
  formula  = FORMULAS[[arm]],
  n_obs    = nrow(Data),
  runtime_h = runtime / 3600,
  # THE number this whole fit exists to produce.
  tau_language_interaction = tau_est,
  beta_interaction         = beta_est,
  # Everything else, for the report's table and for anyone checking the arm did
  # what it claims: all group-level SDs and the fixed effects.
  summary_fixed  = brms::fixef(fit),
  summary_random = summary(fit)$random,
  diagnostics = list(
    rhat_max     = max(brms::rhat(fit), na.rm = TRUE),
    ess_bulk_min = min(posterior::summarise_draws(fit)$ess_bulk, na.rm = TRUE),
    divergent    = sum(brms::nuts_params(fit)$Value[
                        brms::nuts_params(fit)$Parameter == "divergent__"]),
    treedepth_sat = mean(brms::nuts_params(fit)$Value[
                        brms::nuts_params(fit)$Parameter == "treedepth__"] >= 12)
  ),
  loo = if (opt$loo) brms::loo(fit) else "not computed; run scripts/loo_glossa_arms.R on the persisted draws",
  balinese_rescale_note = "Balinese collapsed 10 -> 5 points and rounded, as published.",
  sessionInfo = utils::sessionInfo()
)

tmp <- paste0(out_file, ".tmp")
saveRDS(result, tmp); file.rename(tmp, out_file)
message("[glossa-tau] wrote ", out_file)
