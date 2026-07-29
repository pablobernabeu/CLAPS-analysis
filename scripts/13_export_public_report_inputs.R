#!/usr/bin/env Rscript
# Derive the two disclosure-reviewed public report inputs from the private pilot:
# the already-published Figure 1 raster and language-level composition counts.
#
# Called by scrub_to_public.sh rather than by the analysis pipeline. It takes
# three arguments, the private harmonised pilot CSV and the two public output
# paths, and writes nothing else.
#
# Disclosure scope. Neither output carries participant identifiers, and neither
# is a participant-level table. The composition CSV holds language-level counts
# only. The figure is a raster copy of one already published in the report PDF.
# It plots individual judgements as unlabelled, jittered points, so the
# distribution of the ratings is visible while no rating can be traced to a
# participant.
#
# The ggplot specification below deliberately mirrors the fig-pilot chunk of
# reports/preliminary_sample_size_analysis.qmd, so that the raster the public
# report shows is the same figure the private report renders. Nothing checks that
# automatically. The awk guard in scrub_to_public.sh verifies the chunk labels it
# rewrites, not the plot definition behind them, so a change to fig-pilot has to be
# made here as well or the two reports will silently show different figures.
# scripts/12_plot_raw_gender_comparison.R holds a third, older copy of the same
# plot, which feeds only the retired archive report.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript scripts/13_export_public_report_inputs.R ",
    "<private-pilot.csv> <public-figure.png> <public-composition.csv>"
  )
}

data_path <- args[[1]]
figure_path <- args[[2]]
composition_path <- args[[3]]

pilot <- read_csv(data_path, show_col_types = FALSE)
required <- c(
  "Participant", "Language", "Verb_ID", "Item", "S_Type",
  "Response", "affectedness_scores_agent"
)
missing <- setdiff(required, names(pilot))
if (length(missing)) {
  stop("Pilot data missing required column(s): ", paste(missing, collapse = ", "))
}

dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(composition_path), recursive = TRUE, showWarnings = FALSE)

composition <- pilot |>
  group_by(Language) |>
  summarise(
    Participants = n_distinct(Participant),
    Verbs = n_distinct(Verb_ID),
    `Sentence types` = n_distinct(S_Type),
    Items = n_distinct(Item),
    Judgements = n(),
    .groups = "drop"
  ) |>
  arrange(match(Language, c("English", "Turkish", "Norwegian")))

# The appendix prose in the report, and the facet order of the figure, assume
# these three pilot languages and no others. A fourth language, or a renamed
# one, would otherwise be published in a table that the surrounding text no
# longer describes.
if (!identical(composition$Language, c("English", "Turkish", "Norwegian"))) {
  stop("Expected exactly English, Turkish and Norwegian composition rows")
}
write_csv(composition, composition_path)

# What follows reproduces the fig-pilot chunk of
# reports/preliminary_sample_size_analysis.qmd. The two are separate copies of
# one specification, so a change to either has to be made in the other, or the
# raster published to the mirror will stop matching the figure in the PDF.
pal_stype <- c(
  "Passive" = "#0072B2",
  "Active" = "#D55E00",
  "Pseudo-passive" = "#009E73"
)

d_pilot <- pilot |>
  mutate(Gender = sub(".*_", "", Item)) |>
  filter(
    S_Type %in% c("Passive", "Active", "Pseudo_Passive"),
    Gender %in% c("Man", "Woman")
  ) |>
  group_by(Language) |>
  mutate(affect_z = as.numeric(scale(affectedness_scores_agent))) |>
  ungroup() |>
  mutate(
    S_Type = factor(
      S_Type,
      levels = c("Passive", "Active", "Pseudo_Passive"),
      labels = c("Passive", "Active", "Pseudo-passive")
    ),
    Gender = factor(Gender, levels = c("Man", "Woman")),
    Language = factor(Language, levels = c("English", "Turkish", "Norwegian"))
  )

figure <- ggplot(
  d_pilot,
  aes(affect_z, Response, colour = S_Type, fill = S_Type)
) +
  geom_jitter(
    width = 0.04, height = 0.20, alpha = 0.05, size = 0.45,
    show.legend = FALSE
  ) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  facet_grid(Language ~ Gender) +
  scale_colour_manual(values = pal_stype) +
  scale_fill_manual(values = pal_stype) +
  scale_y_continuous(breaks = 1:7, limits = c(0.5, 7.5)) +
  labs(
    x = "Agent semantic affectedness (z-scored within language)",
    y = "Acceptability rating (1 to 7)",
    colour = "Sentence type",
    fill = "Sentence type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# geom_jitter draws its offsets when the plot is built. Fixing the RNG state
# first makes repeated private-to-public scrubs produce the same raster.
set.seed(20260727)
ggsave(figure_path, figure, width = 6.5, height = 2.9, dpi = 200)

message("[public report inputs] figure: ", figure_path)
message("[public report inputs] aggregate composition: ", composition_path)
