# R/00_reference_audit.R
# Reference audit: verify every citation in references.bib against Crossref.
# Fails loudly on metadata inconsistency. Never invents or silently patches a DOI.
# Output: outputs/reference_audit/reference_audit.csv

suppressPackageStartupMessages({
  library(bib2df)
  library(httr2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

#' Parse references.bib and return a data frame
parse_bib <- function(bib_path = "references.bib") {
  stopifnot(file.exists(bib_path))
  df <- bib2df::bib2df(bib_path)
  df <- df |>
    dplyr::rename_with(tolower) |>
    dplyr::select(dplyr::any_of(c("bibtexkey", "doi", "title", "author", "year", "journal", "booktitle")))
  df
}

#' Query Crossref for a single DOI. Returns a list with verified fields.
#' Returns NULL on network error or HTTP failure — never throws.
query_crossref <- function(doi) {
  if (is.na(doi) || nchar(trimws(doi)) == 0) return(NULL)
  url <- paste0("https://api.crossref.org/works/", utils::URLencode(trimws(doi), reserved = TRUE))
  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_headers("User-Agent" = "CLAPS-reference-audit/1.0 (mailto:researcher@example.ox.ac.uk)") |>
      httr2::req_timeout(20) |>
      httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
      httr2::req_perform()
    if (httr2::resp_status(resp) != 200) return(NULL)
    body <- httr2::resp_body_json(resp)
    item <- body$message
    cr_title <- item$title[[1]] %||% NA_character_
    cr_year  <- as.character(item$published$`date-parts`[[1]][[1]] %||% NA_character_)
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

#' Compare bib entry to Crossref record and return audit status
audit_entry <- function(key, bib_doi, bib_title, bib_year, bib_author, cr) {
  if (is.null(cr)) {
    return(tibble::tibble(
      key = key, doi = bib_doi, bib_title = bib_title, bib_year = bib_year,
      bib_first_author = bib_author, cr_title = NA, cr_year = NA,
      cr_first_author = NA, status = "no_doi_or_not_found", flag = "WARN"
    ))
  }
  title_match <- isTRUE(
    stringr::str_detect(
      stringr::str_to_lower(cr$cr_title),
      stringr::fixed(substr(stringr::str_to_lower(bib_title %||% ""), 1, 40))
    )
  )
  year_match <- isTRUE(trimws(as.character(bib_year)) == trimws(cr$cr_year))
  flag <- dplyr::case_when(
    cr$cr_status != "found"                 ~ "ERROR",
    !title_match                            ~ "TITLE_MISMATCH",
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

#' Run full audit
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
      bib_author = bib$author
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

# Operator for null coalescing (base R < 4.4 compatibility)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nchar(as.character(a)) > 0) a else b
