# testthat coverage for utils/chrom_names.R.
#
# The normaliser turns FASTA header lines into the chromosome names the rest of
# the pipeline uses. It imposes NO pattern: names come out in whatever shape the
# genome FASTA provides, which is what BLAST (subject id), LTRharvest (-seqids)
# and extract_region_fasta.R ("as bedtools does") all use. Those names are what
# reach the tracks and the catalog, so this has to agree with them.
#
# An earlier revision required an NCBI accession with a version suffix, which
# aborted hotspot detection on Bat1K and DNA Zoo assemblies (bare scaffold
# numbers like "02077418") and could also have returned a name DIFFERENT from
# the one BLAST wrote, mismatching silently.

suppressMessages({
  library(testthat)
})

source(file.path("..", "..", "scripts", "utils", "chrom_names.R"))


test_that("a descriptive NCBI header still yields just the accession", {
  # First-token gives the same answer the old accession regex did for every
  # well-formed NCBI header, so this behaviour does not regress.
  expect_equal(
    normalise_chrom_names("CM050508.1 Antrozous pallidus isolate x chromosome 1"),
    "CM050508.1"
  )
  expect_equal(
    normalise_chrom_names("NC_000001.11 Homo sapiens chromosome 1, GRCh38"),
    "NC_000001.11"
  )
})

test_that("bare scaffold numbers are accepted", {
  # The Bat1K / DNA Zoo case that aborted the hotspot stage. These are valid
  # chromosome names; they simply are not NCBI accessions.
  expect_equal(
    normalise_chrom_names(c("02077418", "02077419", "02077420")),
    c("02077418", "02077419", "02077420")
  )
})

test_that("arbitrary assembly names survive, keeping only the first token", {
  expect_equal(normalise_chrom_names("scaffold_1 length=12345"), "scaffold_1")
  expect_equal(normalise_chrom_names("HiC_scaffold_12"), "HiC_scaffold_12")
  expect_equal(normalise_chrom_names("chr1"), "chr1")
})

test_that("nothing returns NA for a non-empty header", {
  # The hard stop() in load_genome_for_hotspot() fires on any NA, so a header
  # that names something must always yield a name.
  headers <- c("CM050508.1 desc", "02077418", "scaffold_9", "chrX")
  expect_false(any(is.na(normalise_chrom_names(headers))))
})

test_that("an empty or whitespace-only header is still NA", {
  # There is genuinely no name here, and silently inventing one would put an
  # empty seqlevel into the genome.
  expect_true(is.na(normalise_chrom_names("")))
  expect_true(is.na(normalise_chrom_names("   ")))
})

test_that("mixed header styles in one genome each resolve correctly", {
  expect_equal(
    normalise_chrom_names(c("CM050508.1 chromosome 1", "02077418", "scaffold_3 len=9")),
    c("CM050508.1", "02077418", "scaffold_3")
  )
})

test_that("the returned vector is always the same length as the input", {
  # Callers assign this straight onto seqlevels(), so a length change would
  # silently misalign every chromosome.
  headers <- c("CM050508.1 a", "", "02077418")
  expect_length(normalise_chrom_names(headers), 3L)
  expect_length(normalise_chrom_names(character()), 0L)
})


test_that("a name that merely resembles an accession is not reinterpreted", {
  # The old pattern would have clipped this to "CM050508.1", disagreeing with
  # the seqname BLAST writes into the tracks. Whatever the FASTA says, wins.
  expect_equal(normalise_chrom_names("CM050508.1extra"), "CM050508.1extra")
})
