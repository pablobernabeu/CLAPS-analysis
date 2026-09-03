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
  # key_glyph = "path" draws the key as a bare line. The default glyph for a
  # ribbon-bearing smooth paints the confidence band into the key too, in the
  # band's own grey, which showed as a grey swatch behind each coloured line.
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE, linewidth = 1,
    key_glyph = "path"
  ) +
  facet_grid(Language ~ Gender) +
  scale_colour_manual(values = pal_stype) +
  scale_fill_manual(values = pal_stype) +
  scale_y_continuous(breaks = 1:7, limits = c(0.5, 7.5)) +
  labs(
    x = "Agent semantic affectedness (z-scored within language)",
    y = "Acceptability rating (1 to 7)",
    colour = "Sentence type"
  ) +
  # Passive leads the legend because it is the model's reference level, set in
  # R/02_preprocess_factors.R, so the reading order matches the contrast the
  # coefficients answer. Do not reorder it to put Active first.
  # The key carries the trend line alone. Mapping fill as well drew a confidence
  # band behind each key, which at this size read as a colour patch, and the note
  # in the report already explains what the bands are.
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.8, alpha = 1)),
    fill = "none"
  ) +
  theme_bw(base_size = 10) +
  # Kept deliberately identical to the fig-pilot chunk in
  # reports/preliminary_sample_size_analysis.qmd, so the raster published for
  # the public report matches the figure the private report draws. Change both
  # together or the two reports will show different figures.
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5),
    strip.text.y = element_text(margin = margin(l = 3, r = 3)),
    axis.text.y = element_text(size = 8),
    panel.spacing.y = grid::unit(5, "pt"),
    # Wider keys, on no background, so each one reads as a line. The grey key
    # panel was legible while a confidence band tinted it, but with the band
    # dropped from the key it became a swatch behind the line.
    legend.key = element_blank(),
    legend.key.width = grid::unit(22, "pt"),
    # The entries sat close enough that each label was as near the next key as to
    # its own, so the pairing had to be guessed. A modest gap settles it, and is
    # preferred to stacking, which would take three rows of height from a figure
    # that has little to spare. The title takes a gap of its own so it does not
    # read as the first entry's label.
    legend.key.spacing.x = grid::unit(10, "pt"),
    legend.title = element_text(margin = margin(r = 8))
  )

# geom_jitter draws its offsets when the plot is built. Fixing the RNG state
# first makes repeated private-to-public scrubs produce the same raster.
set.seed(20260727)
ggsave(figure_path, figure, width = 6.5, height = 4.4, dpi = 200)

message("[public report inputs] figure: ", figure_path)
message("[public report inputs] aggregate composition: ", composition_path)
