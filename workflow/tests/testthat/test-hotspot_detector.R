# End-to-end test of workflow/scripts/hotspot_detector.R, the orchestrator.
#
# The pure transforms are unit-tested in test-hotspot_*.R. This test runs the
# real script through Rscript on a toy genome, so what it checks is the wiring
# between the phases: every output is written, the windows cover the genome,
# and each lineage's dense cluster comes out as a hotspot with its composition.
# About 35 s (two NB fits and the stage PDF); set RETROSEEK_SKIP_SLOW_TESTS=1 to
# skip it on a quick run.

suppressMessages({
  library(testthat)
})

.detector <- normalizePath(file.path("..", "..", "scripts", "hotspot_detector.R"))

# Three 300 kb chromosomes of random bases (fixed seed, no N runs).
.write_toy_genome <- function(path) {
  withr::local_seed(7)
  chrom <- function(name) {
    bases <- paste(sample(c("A", "C", "G", "T"), 300000L, replace = TRUE), collapse = "")
    c(paste0(">", name), substring(bases, seq(1L, 300000L, 60L), seq(60L, 300000L, 60L)))
  }
  writeLines(c(chrom("chrA"), chrom("chrB"), chrom("chrC")), path)
}

# Two lineages, so the per-group split and the join across chromosomes both run.
# Each has one locus in each of 60 windows spread over the genome plus a dense
# cluster of 40 loci: Gammaretrovirus on chrB:100-110 kb, Betaretrovirus on
# chrC:200-210 kb.
.write_toy_catalog <- function(path) {
  lineage <- function(segment, chrom, cluster_start) {
    spread <- data.frame(
      seqname = rep(c("chrA", "chrB", "chrC"), each = 20L),
      start = rep(seq(1000L, 295000L, length.out = 20L), 3L)
    )
    dense <- data.frame(seqname = chrom,
                        start = seq(cluster_start, cluster_start + 9000L, length.out = 40L))
    cbind(rbind(spread, dense), segment = segment)
  }
  loci <- rbind(lineage("Gammaretrovirus", "chrB", 100000L),
                lineage("Betaretrovirus", "chrC", 200000L))
  loci$start <- as.integer(loci$start)
  loci$end <- loci$start + 499L
  loci$species <- "Toyus_toyus"
  loci$strand <- "+"
  loci$source <- "ltr-flanked"
  loci$structure_class <- "full"
  loci$confidence <- "0.9"
  utils::write.csv(loci, path, row.names = FALSE)
}

test_that("the detector writes every output and calls each lineage's cluster", {
  skip_if(identical(Sys.getenv("RETROSEEK_SKIP_SLOW_TESTS"), "1"), "slow end-to-end test")
  dir <- withr::local_tempdir()
  fasta <- file.path(dir, "Toyus_toyus.fa")
  catalog <- file.path(dir, "catalog.csv")
  config <- file.path(dir, "config.yaml")
  .write_toy_genome(fasta)
  .write_toy_catalog(catalog)
  writeLines(c("species:", "  Toyus_toyus: 'Toyus toyus'",
               "hotspot:", "  window_size: 5000", "  group_by: segment"), config)

  out <- function(sub) file.path(dir, sub)
  console <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(.detector), "--fasta", shQuote(fasta), "--hits", shQuote(catalog),
      "--config", shQuote(config), "--parquet_dir", shQuote(out("pq")),
      "--csv_dir", shQuote(out("csv")), "--track_output_dir", shQuote(out("tracks")),
      "--pdf_output_dir", shQuote(out("pdf")), "--log", shQuote(out("job.log"))),
    stdout = TRUE, stderr = TRUE
  ))  # a non-zero exit is a warning here; the status is checked below
  status <- attr(console, "status")  # set only on a non-zero exit
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(utils::tail(console, 15), collapse = "\n"))

  for (f in c("csv/Toyus_toyus.csv", "csv/Toyus_toyus.hotspots.csv",
              "pq/Toyus_toyus.parquet", "pq/Toyus_toyus.manifest.yaml",
              "tracks/Toyus_toyus.gff3", "tracks/Toyus_toyus.bed",
              "pdf/Toyus_toyus.hotspots.pdf")) {
    expect_true(file.exists(out(f)), info = f)
  }
  windows <- utils::read.csv(out("csv/Toyus_toyus.csv"))
  # 300 kb in 5 kb windows on three chromosomes, once per lineage
  expect_equal(nrow(windows), 2L * 3L * 60L)
  expect_setequal(unique(windows$label), c("Gammaretrovirus", "Betaretrovirus"))
  regions <- utils::read.csv(out("csv/Toyus_toyus.hotspots.csv"))
  cluster <- function(chrom, from) {
    regions$seqnames == chrom & regions$start <= from + 10000 & regions$end >= from
  }
  expect_equal(unique(regions$dominant_taxon[cluster("chrB", 100000)]), "Gammaretrovirus")
  expect_equal(unique(regions$dominant_taxon[cluster("chrC", 200000)]), "Betaretrovirus")
})
