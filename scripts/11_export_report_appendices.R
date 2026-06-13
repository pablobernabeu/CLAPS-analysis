# scripts/11_export_report_appendices.R
# Export appendix content for reports/claps_report.qmd.
#
# The report includes these generated markdown files via Quarto include tags,
# so model/prior definitions stay synchronized with R/03_define_priors.R and
# R/04_model_formulas.R without embedding large code blocks in the main text.

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(repo_root, "reports", "appendices", "generated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_env <- new.env(parent = globalenv())
source(file.path(repo_root, "R", "03_define_priors.R"), local = load_env)
source(file.path(repo_root, "R", "04_model_formulas.R"), local = load_env)

write_lines <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

fence_block <- function(lines, lang = "r") {
  c(paste0("```", lang), lines, "```")
}

# -----------------------------------------------------------------------------
# Appendix A: Prior specification
# -----------------------------------------------------------------------------
anchors <- load_env$EMPIRICAL_ANCHORS
regimes <- load_env$PRIOR_REGIMES

anchor_lines <- c(
  "@tbl-anchors lists the empirical effect-size anchors that calibrate the focal-slope",
  "priors.",
  "",
  "| Anchor | Value |",
  "|---|---:|"
)
for (nm in names(anchors)) {
  anchor_lines <- c(anchor_lines, paste0("| `", nm, "` | ", anchors[[nm]], " |"))
}
anchor_lines <- c(anchor_lines, "",
                  ": Empirical Effect-Size Anchors for the Focal-Slope Priors {#tbl-anchors}",
                  "",
                  "```{=latex}", "\\vspace{-0.9em}", "```",
                  "",
                  "*Note.* Anchors are sourced from R/03_define_priors.R.")

prior_rows <- lapply(names(regimes), function(reg_name) {
  list(
    regime = reg_name,
    b_default = regimes[[reg_name]]$b_default,
    b_semantics = regimes[[reg_name]]$b_semantics,
    b_active_int = regimes[[reg_name]]$b_active_int,
    b_pseudo_int = regimes[[reg_name]]$b_pseudo_int,
    Intercept = regimes[[reg_name]]$Intercept,
    sd = regimes[[reg_name]]$sd,
    cor = regimes[[reg_name]]$cor
  )
})

# Compact distribution notation keeps the eight-column table within the page width:
# normal(...) -> N(...), student_t(...) -> t(...).
compact_prior <- function(x) {
  x <- gsub("normal\\(", "N(", x)
  x <- gsub("student_t\\(", "t(", x)
  x
}

prior_lines <- c(
  "",
  "@tbl-priorregimes gives the four prior regimes by parameter.",
  "",
  "| Regime | default | semantics | active | pseudo | Intercept | sd | cor |",
  "|---|---|---|---|---|---|---|---|"
)

for (r in prior_rows) {
  reg_label <- if (identical(r$regime, "literature_centred")) "LC" else r$regime
  prior_lines <- c(
    prior_lines,
    paste0(
      "| ", reg_label, " | ", compact_prior(r$b_default), " | ",
      compact_prior(r$b_semantics), " | ", compact_prior(r$b_active_int), " | ",
      compact_prior(r$b_pseudo_int), " | ", compact_prior(r$Intercept), " | ",
      compact_prior(r$sd), " | ", compact_prior(r$cor), " |"
    )
  )
}
prior_lines <- c(prior_lines, "",
                 ": Prior Regimes by Parameter {#tbl-priorregimes}",
                 "",
                 "```{=latex}", "\\vspace{-0.9em}", "```",
                 "",
                 "*Note.* N denotes a normal prior, t a Student-t prior and LC the literature-centred regime. The regimes are sourced from R/03_define_priors.R.")

write_lines(
  file.path(out_dir, "prior_specifications.md"),
  c(anchor_lines, prior_lines)
)

# -----------------------------------------------------------------------------
# Appendix B: Single-language model ladder
# -----------------------------------------------------------------------------
single_ladder <- load_env$build_model_ladder(has_pseudo_passive = TRUE)
single_names <- load_env$ladder_names()

single_lines <- c(
  "The following ladder definitions are sourced from R/04_model_formulas.R ",
  "via build_model_ladder(has_pseudo_passive = TRUE).",
  ""
)

for (nm in single_names) {
  f <- single_ladder[[nm]]$formula
  single_lines <- c(
    single_lines,
    paste0("### `", nm, "`"),
    "",
    fence_block(deparse(f), lang = "r"),
    ""
  )
}

write_lines(
  file.path(out_dir, "single_language_ladder.md"),
  single_lines
)

# -----------------------------------------------------------------------------
# Appendix C: Cross-language model ladder
# -----------------------------------------------------------------------------
cross_ladder <- load_env$build_multilanguage_ladder()
cross_names <- load_env$multilanguage_ladder_names()

cross_lines <- c(
  "The following ladder definitions are sourced from R/04_model_formulas.R ",
  "via build_multilanguage_ladder().",
  ""
)

for (nm in cross_names) {
  f <- cross_ladder[[nm]]$formula
  cross_lines <- c(
    cross_lines,
    paste0("### `", nm, "`"),
    "",
    fence_block(deparse(f), lang = "r"),
    ""
  )
}

write_lines(
  file.path(out_dir, "cross_language_ladder.md"),
  cross_lines
)

message("Appendix markdown files written to: ", out_dir)
