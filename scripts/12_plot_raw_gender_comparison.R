#!/usr/bin/env Rscript
# scripts/12_plot_raw_gender_comparison.R
# Raw-data visualisation complementing the gender model variation.
# Shows the study's core comparison — acceptability as a function of the
# agent/gender-specific semantic affectedness, by sentence type — with referent
# gender (Man/Woman) side by side and one row per language.
# Uses affectedness_scores_agent (the predictor the gender models use), z-scored
# within language so panels are comparable. Pilot data = calibration only.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

cfg_path  <- "config/analysis_config.yaml"
data_path <- if (file.exists(cfg_path)) yaml::read_yaml(cfg_path)$pilot_data_path else
             "data/pilot/claps_pilot_harmonised.csv"
out_dir <- "outputs/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(data_path, show_col_types = FALSE) |>
  mutate(Gender = sub(".*_", "", Item)) |>
  filter(S_Type %in% c("Passive", "Active", "Pseudo_Passive"),
         Gender %in% c("Man", "Woman")) |>
  group_by(Language) |>
  mutate(affect_z = as.numeric(scale(affectedness_scores_agent))) |>
  ungroup() |>
  mutate(
    S_Type   = factor(S_Type, levels = c("Passive", "Active", "Pseudo_Passive"),
                       labels = c("Passive", "Active", "Pseudo-Passive")),
    Gender   = factor(Gender, levels = c("Man", "Woman")),
    Language = factor(Language, levels = c("English", "Turkish", "Norwegian"))
  )

pal <- c("Passive" = "#0072B2", "Active" = "#D55E00", "Pseudo-Passive" = "#009E73")

p <- ggplot(d, aes(affect_z, Response, colour = S_Type, fill = S_Type)) +
  geom_jitter(width = 0.04, height = 0.20, alpha = 0.05, size = 0.45, show.legend = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  facet_grid(Language ~ Gender) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(breaks = 1:7, limits = c(0.5, 7.5)) +
  # Title/subtitle/caption are intentionally omitted here and supplied in the
  # report's figure caption instead. Margins are balanced top and bottom.
  labs(
    x      = "Agent Semantic Affectedness (z-Scored Within Language)",
    y      = "Acceptability Rating",
    colour = "Sentence Type",
    fill   = "Sentence Type"
  ) +
  guides(colour = guide_legend(title.position = "top"),
         fill   = guide_legend(title.position = "top")) +
  theme_bw(base_size = 10) +
  theme(
    legend.position    = "top",
    legend.direction   = "horizontal",
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(hjust = 0.5, face = "bold", size = 11),
    strip.text.x       = element_text(face = "bold", size = 11),
    strip.text.y       = element_text(face = "bold", size = 10,
                                      margin = margin(t = 12, r = 6, b = 12, l = 6, unit = "pt")),
    legend.margin      = margin(0, 0, 2, 0, unit = "pt"),
    legend.box.spacing = grid::unit(3, "pt"),
    axis.title.x       = element_text(margin = margin(t = 4, b = 8, unit = "pt")),
    axis.title.y       = element_text(margin = margin(r = 6, unit = "pt")),
    plot.margin        = margin(t = 2, r = 8, b = 5, l = 8, unit = "pt")
  )

out_png <- file.path(out_dir, "raw_gender_comparison.png")
ggsave(out_png, p, width = 6.5, height = 3.8, dpi = 200)
cat("Wrote", out_png, "\n")
cat("Rows plotted:", nrow(d), "\n")
print(d |> count(Language, Gender, S_Type) |> tidyr::pivot_wider(names_from = S_Type, values_from = n))
