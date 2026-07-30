# .lintr.R — lint configuration for the CLAPS analysis
#
# Single source of truth for what "clean" means here. Read automatically by
# lintr::lint_dir("R") and by the Lint step in .github/workflows/static-checks.yaml.
# The CI step previously passed its linters inline, which put the standard in a YAML
# file rather than with the code and meant `lintr::lint_dir("R")` run locally applied
# different rules from CI. Reproduce CI exactly with:
#
#   Rscript -e 'lintr::lint_dir("R")'
#
# Written in R rather than the older DCF `.lintr` format (lintr >= 3.2), because DCF
# requires every continuation line to be indented and silently misparses otherwise.
#
# Every deviation from lintr's defaults is justified below. None is a blanket
# silencing of inconvenient output: the mechanical linters (indentation, braces,
# semicolons, infix spacing, commas, line length) are all left on, and the codebase
# is clean against them.

linters <- linters_with_defaults(

  # --- disabled, with reasons ------------------------------------------------

  # 188 findings, none of them real. Two architectural facts defeat this linter:
  #
  #   1. Non-standard evaluation. dplyr refers to columns by bare name, so BF_10,
  #      language, status, model_level and Semantics are reported as undefined
  #      global variables. Silencing that honestly would mean .data$ pronouns or
  #      utils::globalVariables() throughout, which buys nothing and obscures the
  #      pipelines.
  #   2. Cross-file source(). Modules under R/ call source("R/0X_....R") at run
  #      time, so compute_all_bf(), align_prior_to_model(),
  #      extract_convergence_diagnostics() and production_sampling() are defined in
  #      a file lintr analyses separately, and are reported as undefined functions.
  #
  # It also emits "Could not find exported symbols for package" whenever brms and
  # the other heavy dependencies are absent, which makes its output depend on the
  # machine rather than on the code.
  object_usage_linter = NULL,

  # Nine findings, every one prose rather than dead code: usage examples such as
  # "#   check_crosslanguage_consistency(df_all)", the OSF reference model written
  # out as a formula in R/04_model_formulas.R, and citation lists whose semicolons
  # make them parse as expressions ("Bürkner & Vuorre (2019); Liddell & Kruschke
  # (2018); ..."). Removing any of them would delete documentation, which is the
  # opposite of the intent. Genuinely dead code is caught in review.
  commented_code_linter = NULL,

  # Names are checked for length below, but not for style. The project mixes
  # snake_case functions with the column names of the source data (S_Type, Verb_ID,
  # Semantics_scaled, Response), which are fixed by the data and cannot be renamed
  # without breaking the model formulas and the priors that reference them.
  object_name_linter = NULL,

  # --- relaxed, with reasons -------------------------------------------------

  # The default of 30 flags seven deliberately descriptive names, among them
  # compute_ceiling_calibrated_thresholds (37) and
  # exclude_norwegian_synthetic_passive (35). The model-ladder levels are longer
  # still: L1_random_intercepts_plus_participant_semantics is 47 characters and is a
  # naming contract, parsed for its L-number and matched across three files. In a
  # research codebase an unambiguous name is worth more than a short one.
  object_length_linter = object_length_linter(50L),

  # 80 is too tight for roxygen prose and for brms formulas written on one line.
  # 120 is what the CI step used before this file existed.
  line_length_linter = line_length_linter(120L)
)

exclusions <- list(
  # renv manages these; they are not ours to style.
  "renv"
)
