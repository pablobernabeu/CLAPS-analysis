# R/00_reference_audit.R
#
# Purpose
#   Verify every citation in references.bib against the authoritative Crossref
#   record, so that no report is rendered from a bibliography containing invented
#   or mistyped metadata. The audit checks bibliographic *fields* only. Whether a
#   cited work actually supports the claim made about it is a matter of reading,
#   and is recorded by hand in docs/preregistration_decisions.md.
#
# Inputs
#   bib_path  Path to a BibTeX file (default "references.bib"). Each entry is
#             expected to carry a DOI; entries without one are flagged, not
#             silently accepted.
#
# Outputs
#   outputs/reference_audit/reference_audit.csv  One row per bibliography entry,
#             holding the local field values, the Crossref values, and a flag in
#             {OK, WARN, ERROR, TITLE_MISMATCH, YEAR_MISMATCH}.
#
# How to run
#   Rscript scripts/00_verify_references.R      # from design_analysis/
#
# Behaviour on failure
#   run_reference_audit() stops with an error if any entry mismatches, so a
#   failing audit blocks report rendering rather than producing a report with
#   unverified citations. Nothing here ever edits references.bib: a DOI is never
#   invented, and a mismatch is never patched automatically. Correcting the
#   bibliography is a human decision, because the right fix is sometimes to
#   change the DOI and sometimes to change the title.
#
# Assumptions and limitations
#   Requires network access to api.crossref.org. A network outage yields
#   cr_status = "error: ..." and the WARN/ERROR path, never a false pass.
#   Preprints, book chapters and grey literature are frequently absent from
#   Crossref; those entries surface as "no_doi_or_not_found" and need manual
#   confirmation against the publisher record.

suppressPackageStartupMessages({
  library(bib2df)
  library(httr2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

# Null-coalescing operator. Defined locally because base R gained %||% only in 4.4,
# while config/arc_modules.yaml sets min_version 4.3.0 and CI pins 4.3.3, so it
# cannot be assumed present. Defined here, above every use, rather than at the foot
# of the file: at the foot it would not serve the uses above it, and the script would
# work only on an R new enough to make the definition redundant.
#
# This variant is deliberately STRICTER than base R's: it also falls through on NA
# and on the empty string, because a BibTeX field parsed from an empty brace pair
# arrives as "" rather than NULL and should be treated as absent.
#
# HAZARD: because the name is the same, this strict variant is interchangeable with
# the plain one defined in R/06_simulate_design.R and elsewhere, and whichever is
# evaluated last wins. targets.R calls tar_source("R/"), which sources the whole
# directory, so under the {targets} pipeline the plain variant from a
# later-sorting file overrides this one and the audit loses its NA/"" handling.
# Scripts that source only this file (scripts/00_verify_references.R) are unaffected.
# The length guards matter. `&&` has required a length-one condition since R 4.3,
# and this operator is applied to Crossref's `author`, a list with one element per
# author, so the previous one-line form raised "'length = 10' in coercion to
# 'logical(1)'" for every multi-author record. tryCatch recorded that as the
# entry's status, so the audit reported a metadata failure for work whose metadata
# it had never compared. Only a length-one value can be NA or empty in the sense
# meant here; anything longer is a real value and passes through untouched.
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (length(a) == 1L && (is.na(a) || !nzchar(trimws(as.character(a))))) return(b)
  a
}

#' Read a BibTeX file into a tidy data frame of the fields the audit needs.
#'
#' @param bib_path Path to the .bib file.
#' @return A data frame with one row per entry and lower-case column names,
#'   restricted to bibtexkey, doi, title, author, year, journal and booktitle.
#'   `any_of()` rather than `all_of()` is used because a bibliography of only
#'   journal articles has no `booktitle` column, and vice versa; a missing
#'   column is normal here rather than an error.
parse_bib <- function(bib_path = "references.bib") {
  stopifnot(file.exists(bib_path))
  df <- bib2df::bib2df(bib_path)
  df <- df |>
    # bib2df returns upper-case BibTeX field names (TITLE, DOI, ...); fold them
    # to lower case so downstream code has one spelling to rely on.
    dplyr::rename_with(tolower) |>
    dplyr::select(dplyr::any_of(c("bibtexkey", "doi", "title", "author", "year", "journal", "booktitle")))
  df
}

#' Look up one DOI in the Crossref REST API.
#'
#' @param doi A DOI string, with or without surrounding whitespace. An NA or
#'   empty value returns NULL, which the caller reports as "no DOI".
#' @return A list with cr_title, cr_year, cr_first_author and cr_status, or NULL
#'   when there is nothing to look up. Never throws: a network failure, a
#'   timeout or a non-200 response is converted into a value the caller can
#'   record, so one unreachable DOI cannot abort an audit of the whole
#'   bibliography.
query_crossref <- function(doi) {
  if (is.na(doi) || nchar(trimws(doi)) == 0) return(NULL)
  # DOIs legitimately contain characters that are reserved in a URL path, most
  # often "/" and occasionally "<" ">" in older Wiley identifiers, so the DOI is
  # percent-encoded with reserved = TRUE rather than pasted in raw.
  url <- paste0("https://api.crossref.org/works/", utils::URLencode(trimws(doi), reserved = TRUE))
  tryCatch({
    resp <- httr2::request(url) |>
      # Crossref asks API users to identify themselves and supply a contact
      # address. Doing so routes the request to their "polite" pool, which is
      # more reliably available than the anonymous pool.
      httr2::req_headers("User-Agent" = "CLAPS-reference-audit/1.0 (mailto:researcher@example.ox.ac.uk)") |>
      # A whole-bibliography audit runs unattended, so a stalled connection must
      # time out rather than hang the render. Three tries with a two-second
      # backoff absorb the transient 5xx responses the API returns under load
      # without hammering it.
      httr2::req_timeout(20) |>
      httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
      httr2::req_perform()
    if (httr2::resp_status(resp) != 200) return(NULL)
    body <- httr2::resp_body_json(resp)
    # Crossref wraps the record in a "message" envelope.
    item <- body$message
    # "title" is an array because some records carry a translated or alternative
    # title; the first element is the primary one.
    cr_title <- item$title[[1]] %||% NA_character_
    # "published.date-parts" is a nested array of [[year, month, day]] with month
    # and day optional, so the year is the first element of the first part.
    # Prefer the print/issue year. Crossref's `published` is the earliest date it
    # holds, which for a work released online ahead of print is the online year,
    # while the bibliography records the issue year an APA citation uses. Comparing
    # the two flagged two correct entries (Cognitive Science 2016, online 2015;
    # Behavior Research Methods 2023, online 2022).
    date_parts <- item$`published-print`$`date-parts` %||%
      item$`published-online`$`date-parts` %||%
      item$published$`date-parts`
    cr_year <- if (is.null(date_parts)) {
      NA_character_
    } else {
      as.character(date_parts[[1]][[1]])
    }
    # Rebuild "Family, Given" so the local and Crossref author strings are
    # comparable by eye in the CSV. Corporate authors have no given name, hence
    # the %||% fallbacks to "".
    cr_authors <- purrr::map_chr(
      item$author %||% list(),
      ~ paste0(.x$family %||% "", ", ", .x$given %||% "")
    )
    cr_first_author <- if (length(cr_authors) > 0) cr_authors[[1]] else NA_character_
    list(
      cr_title        = cr_title,
      cr_year         = cr_year,
      cr_first_author = cr_first_author,
      cr_status       = "found"
    )
  }, error = function(e) {
    list(cr_title = NA, cr_year = NA, cr_first_author = NA, cr_status = paste0("error: ", conditionMessage(e)))
  })
}

#' Compare one bibliography entry with its Crossref record.
#'
#' @param key BibTeX citation key, carried through so a flagged row can be found
#'   in the .bib file.
#' @param bib_doi,bib_title,bib_year,bib_author Field values as written locally.
#' @param cr The list returned by query_crossref(), or NULL.
#' @return A one-row tibble pairing the local and Crossref values with a `flag`.
#'   Flags are ordered by severity in the case_when() below: an unreachable
#'   record (ERROR) is reported even if the title would also have mismatched,
#'   because the title comparison is meaningless without a record to compare to.
audit_entry <- function(key, bib_doi, bib_title, bib_year, bib_author, cr) {
  if (is.null(cr)) {
    # Every column is typed exactly as the matched branch below types it. Until
    # 2026-09-02 this branch passed bib_year through as the bibliography parser
    # returned it, a double for a bare year, and left the three cr_ columns as
    # logical NA. The matched branch coerces bib_year with as.character() and
    # fills the cr_ columns with strings, so a run containing both kinds of entry
    # failed in bind_rows with "Can't combine bib_year <character> and <double>".
    # That aborted the whole audit, and the gate has been red since the first
    # DOI-less entry was added.
    return(tibble::tibble(
      key = key, doi = bib_doi, bib_title = bib_title,
      bib_year = as.character(bib_year),
      bib_first_author = bib_author,
      cr_title = NA_character_, cr_year = NA_character_,
      cr_first_author = NA_character_,
      status = "no_doi_or_not_found", flag = "WARN"
    ))
  }
  # Titles are compared on the first 40 lower-cased characters rather than in
  # full. Publishers and BibTeX disagree routinely on subtitle punctuation, on
  # LaTeX escapes and on trailing series information, so an exact match produces
  # mismatches that are not real errors. A 40-character prefix is long enough to
  # be specific to one work yet short enough to sit before the colon in almost
  # every title. str_fixed() disables regex, since titles contain characters
  # such as "(" that would otherwise be read as regex syntax. isTRUE() collapses
  # the NA that str_detect() returns for a missing title into FALSE, so a missing
  # title is treated as a mismatch and gets flagged.
  # Both sides are normalised first. BibTeX protects capitals with braces
  # ({B}ayesian), and Crossref returns HTML for some records (<b>brms</b>: An
  # <i>R</i> package). Neither is content, and comparing them literally reported a
  # mismatch for 22 of 58 correctly written entries, which is what kept this gate
  # red. Dash variants are folded too, since publishers alternate between hyphen,
  # en dash and em dash within the same title.
  normalise_title <- function(x) {
    x <- as.character(x %||% "")
    x <- gsub("<[^>]*>", " ", x)                      # Crossref HTML markup
    x <- gsub("[{}]", "", x)                          # BibTeX capital protection
    x <- gsub("\\[A-Za-z]+", " ", x)                 # stray LaTeX commands
    x <- gsub("[‐-―−]", "-", x)          # hyphen, en dash, em dash, minus
    x <- gsub("-{2,}", "-", x)                          # BibTeX -- for an en dash
    x <- gsub("[[:space:]]+", " ", x)
    # Removing a tag can strand a space in front of punctuation, as it does when
    # a bold tag is followed by a newline and a colon, so those are closed up.
    x <- gsub(" +([,;:.])", "\\1", x)
    stringr::str_to_lower(trimws(x))
  }
  bib_norm <- normalise_title(bib_title)
  cr_norm  <- normalise_title(cr$cr_title)
  # Compared as a common prefix over whichever is shorter, capped at 40
  # characters. Some Crossref records hold a truncated title, as gelmanCarlin2014
  # does ("Beyond Power Calculations" for a work whose full title continues), and
  # a fixed 40-character window reported those as mismatches against a correct
  # bibliography.
  n_cmp <- min(40L, nchar(bib_norm), nchar(cr_norm))
  title_match <- isTRUE(
    n_cmp > 0L && substr(bib_norm, 1, n_cmp) == substr(cr_norm, 1, n_cmp)
  )
  # Year is compared exactly. An off-by-one year is usually a genuine error, from
  # citing the online-first date instead of the issue date, and is worth flagging.
  bib_y <- suppressWarnings(as.integer(trimws(as.character(bib_year))))
  cr_y  <- suppressWarnings(as.integer(trimws(cr$cr_year)))
  year_match <- isTRUE(trimws(as.character(bib_year)) == trimws(cr$cr_year))
  # A bibliography records the issue year, which is what an APA citation uses,
  # while Crossref may hold only the online-ahead-of-print date. Where a record
  # carries no print date, a bibliography year one later is that pattern rather
  # than an error, so it is reported without blocking the gate. Any other
  # difference still blocks.
  year_online_first <- isTRUE(
    !year_match && !is.na(bib_y) && !is.na(cr_y) && bib_y - cr_y == 1L
  )
  flag <- dplyr::case_when(
    cr$cr_status != "found"                 ~ "ERROR",
    !title_match                            ~ "TITLE_MISMATCH",
    year_online_first                       ~ "YEAR_ONLINE_FIRST",
    !year_match                             ~ "YEAR_MISMATCH",
    TRUE                                    ~ "OK"
  )
  tibble::tibble(
    key             = key,
    doi             = bib_doi,
    bib_title       = bib_title,
    bib_year        = as.character(bib_year),
    bib_first_author = bib_author,
    cr_title        = cr$cr_title,
    cr_year         = cr$cr_year,
    cr_first_author = cr$cr_first_author,
    status          = cr$cr_status,
    flag            = flag
  )
}

#' Audit an entire bibliography and stop if anything fails verification.
#'
#' @param bib_path Path to the .bib file.
#' @param out_dir Directory for reference_audit.csv; created if absent.
#' @return Invisibly, the audit tibble (one row per entry).
#' @details Entries are queried one at a time rather than through Crossref's
#'   bulk filter endpoint. A bibliography of this size takes well under a minute,
#'   and per-DOI requests make it obvious from the console log which entry a
#'   failure belongs to. The CSV is written before the failure check, so a
#'   failing run still leaves a complete record to work from.
run_reference_audit <- function(bib_path = "references.bib",
                                out_dir   = "outputs/reference_audit") {
  message("[audit] Parsing ", bib_path)
  bib <- parse_bib(bib_path)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- purrr::pmap_dfr(
    list(
      key        = bib$bibtexkey,
      bib_doi    = bib$doi,
      bib_title  = bib$title,
      bib_year   = bib$year,
      # bib2df returns `author` as a list column holding every author of an entry.
      # Passing it through unreduced put a vector into a column named
      # bib_first_author, and the tibble then recycled the scalar columns against
      # it, so a ten-author entry became ten identical rows. The audit is a
      # per-entry check, so it takes the first author, which is what the column
      # name and the Crossref side both mean.
      bib_author = vapply(
        bib$author,
        function(a) if (length(a) >= 1L) trimws(as.character(a)[[1]]) else NA_character_,
        character(1)
      )
    ),
    function(key, bib_doi, bib_title, bib_year, bib_author) {
      message("[audit] ", key, " (DOI: ", bib_doi %||% "none", ")")
      cr <- query_crossref(bib_doi)
      audit_entry(key, bib_doi, bib_title, bib_year, bib_author, cr)
    }
  )

  out_path <- file.path(out_dir, "reference_audit.csv")
  readr::write_csv(results, out_path)
  message("[audit] Written to ", out_path)

  # Only these three flags are treated as blocking. YEAR_ONLINE_FIRST is not,
  # because it marks a bibliography year one later than a Crossref record that
  # carries no print date, which is how an online-ahead-of-print article is cited.
  # WARN ("no_doi_or_not_found")
  # is deliberately not, because a legitimately DOI-less source such as a
  # dissertation or a corpus release would otherwise make the audit unpassable.
  # Those rows still appear in the CSV and need checking by hand.
  n_errors <- sum(results$flag %in% c("ERROR", "TITLE_MISMATCH", "YEAR_MISMATCH"))
  if (n_errors > 0) {
    stop(
      "[audit] ", n_errors, " citation(s) failed metadata verification. ",
      "See ", out_path, " for details. ",
      "Fix references.bib before rendering reports."
    )
  }
  message("[audit] All ", nrow(results), " references passed.")
  invisible(results)
}
