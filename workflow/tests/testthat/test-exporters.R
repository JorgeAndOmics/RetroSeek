# testthat tests for the erv_like GFF3 exporter in
# workflow/scripts/range_analysis/exporters.R.
#
# Focus: erv_like_track_exporter must (a) emit a header-only but valid GFF3 when
# there are no candidates, and (b) write parent candidates with ID= and their
# child members with Parent=<that ID>, with the `type` rendered as GFF3 column 3.
# A lightweight readLines+grepl check avoids a full rtracklayer re-import.
#
# Run with: make test-r.

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "range_analysis", "erv_assembly.R"))
source(file.path(.script_dir, "range_analysis", "exporters.R"))

.full_erv <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1", strand = "+",
    ranges = IRanges::IRanges(start = c(100, 300, 500), end = c(200, 400, 600))
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    probe = c("GAG", "POL", "ENV"), virus = "HIV", label = "Lab",
    species = "Test_sp", max_bitscore = 100, max_identity = 90,
    min_evalue = 1e-10, query_coverage = 0.5,
    ID = c("GAG_1", "POL_1", "ENV_1")
  )
  assemble_erv_like(gr, c("GAG", "POL", "ENV"), "virus", 1500, FALSE, 1.0,
                    list(agg_concat_separator = "; "))
}


test_that("empty input writes a header-only, valid GFF3", {
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp))
  erv_like_track_exporter(.empty_erv_parents(), .empty_erv_children(), tmp,
                          "RetroSeek/test")
  lines <- readLines(tmp)
  expect_equal(lines[1], "##gff-version 3")
  expect_true(any(grepl("^##source-version RetroSeek/test", lines)))
  # No feature rows.
  expect_false(any(grepl("erv_like", lines)))
})

test_that("parent + children are written with ID / Parent linkage and types", {
  res <- .full_erv()
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp))
  erv_like_track_exporter(res$parents, res$children, tmp, "RetroSeek/test")
  lines <- readLines(tmp)

  expect_equal(lines[1], "##gff-version 3")
  # One parent feature (type erv_like) carrying ID=ERV_like_1.
  parent_lines <- grep("\terv_like\t", lines, value = TRUE)
  expect_length(parent_lines, 1L)
  expect_match(parent_lines, "ID=ERV_like_1")
  # Three child members, each parented to the candidate.
  member_lines <- grep("\terv_like_member\t", lines, value = TRUE)
  expect_length(member_lines, 3L)
  expect_true(all(grepl("Parent=ERV_like_1", member_lines)))
})
