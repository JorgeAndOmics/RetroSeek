# testthat coverage for load_dropped_elements (ranges/io.R): reading the per-element
# verdict table the Dfam filter (ltr_arm_filter.py, ADR-026) writes.
#
# What can go wrong quietly: Python writes the booleans as True/False and a reader
# that expects TRUE/FALSE keeps nothing; a header-only table (no element overlaps a
# hit) read as an error; a missing file read as "nothing to drop", which would
# silently turn the filter off.

suppressMessages({
  library(testthat)
})

.script_dir <- file.path("..", "..", "scripts", "ranges")
source(file.path(.script_dir, "io.R"))

.verdict_csv <- function(lines) {
  path <- tempfile(fileext = ".csv")
  writeLines(c(paste("element", "seqname", "arm1_family", "arm1_class",
                     "arm2_family", "arm2_class", "drop", "reason", sep = ","),
               lines), path)
  path
}

test_that("load_dropped_elements returns no ids when no table is given (filter off)", {
  expect_identical(load_dropped_elements(NULL), character(0))
})

test_that("load_dropped_elements keeps only the elements marked True", {
  path <- .verdict_csv(c(
    "LTR_retrotransposon1,chr1,famA,LTR/ERVK,famA,LTR/ERVK,False,",
    "LTR_retrotransposon2,chr1,famB,LINE/L1,famB,LINE/L1,True,arm is LINE/L1",
    "LTR_retrotransposon3,chr2,,,famC,SINE/Alu,True,arm is SINE/Alu"
  ))
  expect_identical(load_dropped_elements(path),
                   c("LTR_retrotransposon2", "LTR_retrotransposon3"))
})

test_that("load_dropped_elements also reads R-style TRUE", {
  path <- .verdict_csv("LTR_retrotransposon9,chr1,f,DNA,f,DNA,TRUE,arm is DNA")
  expect_identical(load_dropped_elements(path), "LTR_retrotransposon9")
})

test_that("load_dropped_elements on a header-only table returns no ids", {
  expect_identical(load_dropped_elements(.verdict_csv(character(0))), character(0))
})

test_that("load_dropped_elements stops on a missing file instead of dropping nothing", {
  expect_error(load_dropped_elements(file.path(tempdir(), "no-such-file.csv")),
               "ltr_arm_filter")
})

test_that("load_dropped_elements stops when the table lacks the verdict columns", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("element,seqname", "LTR_retrotransposon1,chr1"), path)
  expect_error(load_dropped_elements(path), "drop")
})
