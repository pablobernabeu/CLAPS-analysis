# tests/testthat/test-seed-disjointness.R
#
# Guards the seed-allocation invariant documented under "Seed allocation across
# grids" in docs/design_power_analysis_pipeline.md: no two design grids may hand the
# same seed to different cells.
#
# Why this matters. A seed is the only thing distinguishing one replicate of a design
# point from another, and it is also part of every output filename. Two grids sharing
# a seed range draw from the same RNG stream, so their replicates are not independent,
# and two cells that also share every design parameter would produce byte-identical
# data while appearing in the aggregation as two separate replicates. Three
# generators shared base 7e5 and two shared 9e5 until 2026-07-30; the comments in
# those files asserted they were collision-free, which is precisely the kind of claim
# a test should be making instead of a comment.
#
# Two independent checks, because each catches what the other cannot:
#   1. The committed grid CSVs must have pairwise-disjoint seed sets. Empirical, and
#      decisive about the grids actually in use.
#   2. No seed base literal may appear in two different generators. Catches a
#      collision at source, before any grid is regenerated, including for grids whose
#      CSV is not committed.

library(testthat)

grid_dir <- here::here("config")
gen_dir  <- here::here("scripts")

# ---------------------------------------------------------------------------
# 1. Committed grid CSVs: pairwise-disjoint seed sets
# ---------------------------------------------------------------------------

test_that("committed design grids have pairwise-disjoint seed sets", {
  csvs <- list.files(grid_dir, pattern = "^design_grid.*[.]csv$", full.names = TRUE)
  skip_if(length(csvs) < 2, "fewer than two grid CSVs present")

  seeds <- list()
  for (f in csvs) {
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (!"seed" %in% names(d)) next
    s <- suppressWarnings(as.integer(d$seed))
    s <- s[!is.na(s)]
    if (length(s)) seeds[[basename(f)]] <- unique(s)
  }
  skip_if(length(seeds) < 2, "fewer than two grid CSVs carry a seed column")

  nms <- names(seeds)
  overlaps <- character(0)
  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      if (j <= i) next
      shared <- intersect(seeds[[nms[i]]], seeds[[nms[j]]])
      if (length(shared)) {
        overlaps <- c(overlaps, sprintf(
          "%s <-> %s: %d shared seed(s), e.g. %s",
          nms[i], nms[j], length(shared),
          paste(utils::head(sort(shared), 3), collapse = ", ")))
      }
    }
  }
  expect_equal(overlaps, character(0),
               info = paste0("Grids sharing seeds:\n  ",
                             paste(overlaps, collapse = "\n  ")))
})

# ---------------------------------------------------------------------------
# 2. Generator sources: no seed base shared between two generators
# ---------------------------------------------------------------------------

# Pull the seed base literals out of one generator. A base is an integer literal of
# at least five digits appearing on a line that mentions `seed`. set.seed() calls are
# excluded: their argument is an RNG seed for the generator's own sampling of draw
# indices, not a base for the grid's cell seeds, and it is a date literal here.
extract_seed_bases <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- sub("#.*$", "", ln)                       # ignore comments
  ln <- ln[grepl("seed", ln, ignore.case = TRUE)]
  ln <- ln[!grepl("set\\.seed", ln)]
  m  <- regmatches(ln, gregexpr("\\b[0-9]{5,}L?\\b", ln))
  v  <- suppressWarnings(as.numeric(sub("L$", "", unlist(m))))
  sort(unique(v[!is.na(v) & v >= 1e5]))
}

test_that("no seed base is shared between two grid generators", {
  gens <- list.files(gen_dir, pattern = "^generate_.*grid[.]R$", full.names = TRUE)
  skip_if(length(gens) < 2, "fewer than two grid generators present")

  bases <- lapply(gens, extract_seed_bases)
  names(bases) <- basename(gens)
  bases <- bases[lengths(bases) > 0]
  skip_if(length(bases) < 2, "fewer than two generators expose a seed base")

  # One generator may legitimately use several bases (generate_design_grid.R builds
  # four grids, and generate_pooled_v2_grid.R adds a smoke-test seed). The invariant
  # is only that a base is not used by two DIFFERENT files.
  flat  <- data.frame(
    file = rep(names(bases), lengths(bases)),
    base = unlist(bases, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  dup <- flat$base[duplicated(flat$base)]
  clashes <- character(0)
  for (b in unique(dup)) {
    who <- sort(unique(flat$file[flat$base == b]))
    if (length(who) > 1) {
      clashes <- c(clashes, sprintf("base %.0f used by: %s", b,
                                    paste(who, collapse = ", ")))
    }
  }
  expect_equal(clashes, character(0),
               info = paste0("Generators sharing a seed base:\n  ",
                             paste(clashes, collapse = "\n  ")))
})

test_that("extract_seed_bases finds the documented base of a known generator", {
  # A canary for check 2. If the extraction silently stopped matching, both files
  # would report zero bases and the disjointness test above would vacuously pass.
  f <- file.path(gen_dir, "generate_floor50_power_grid.R")
  skip_if_not(file.exists(f), "floor50 generator not present")
  expect_true(7e5 %in% extract_seed_bases(f))
})
