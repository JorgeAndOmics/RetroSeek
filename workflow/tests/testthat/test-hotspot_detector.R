# End-to-end test of workflow/scripts/hotspot_detector.R, the orchestrator.
#
# The pure transforms are unit-tested in test-hotspot_*.R. This test runs the
# real script through Rscript on a toy genome, so what it checks is the wiring
# between the phases: every output is written, the windows cover the genome,
# and a dense cluster of loci comes out as a hotspot with its composition.
# About 25 s: it fits the NB model and draws the stage PDF.

suppressMessages({
  library(testthat)
})

.detector <- normalizePath(file.path("..", "..", "scripts", "hotspot_detector.R"))

# Three 300 kb chromosomes of random bases (fixed seed, no N runs).
.write_toy_genome <- function(path) {
  set.seed(7)
  chrom <- function(name) {
    bases <- paste(sample(c("A", "C", "G", "T"), 300000L, replace = TRUE), collapse = "")
    c(paste0(">", name), substring(bases, seq(1L, 300000L, 60L), seq(60L, 300000L, 60L)))
  }
  writeLines(c(chrom("chrA"), chrom("chrB"), chrom("chrC")), path)
}

# One locus in each of 120 windows spread over the genome, plus 40 loci packed
# into chrB:100,000-110,000: the cluster the detector must call.
.write_toy_catalog <- function(path) {
  spread <- data.frame(
    seqname = rep(c("chrA", "chrB", "chrC"), each = 40L),
    start = rep(seq(1000L, 295000L, length.out = 40L), 3L)
  )
  dense <- data.frame(seqname = "chrB", start = seq(100000L, 109000L, length.out = 40L))
  loci <- rbind(spread, dense)
  loci$start <- as.integer(loci$start)
  loci$end <- loci$start + 499L
  loci$species <- "Toyus_toyus"
  loci$strand <- "+"
  loci$source <- "ltr-flanked"
  loci$segment <- "Gammaretrovirus"
  loci$structure_class <- "full"
  loci$confidence <- "0.9"
  utils::write.csv(loci, path, row.names = FALSE)
}

test_that("the detector writes every output and calls the dense cluster", {
  dir <- withr::local_tempdir()
  fasta <- file.path(dir, "Toyus_toyus.fa")
  catalog <- file.path(dir, "catalog.csv")
  config <- file.path(dir, "config.yaml")
  .write_toy_genome(fasta)
  .write_toy_catalog(catalog)
  writeLines(c("species:", "  Toyus_toyus: 'Toyus toyus'",
               "hotspot:", "  window_size: 5000", "  group_by: segment"), config)

  out <- function(sub) file.path(dir, sub)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(.detector), "--fasta", shQuote(fasta), "--hits", shQuote(catalog),
      "--config", shQuote(config), "--parquet_dir", shQuote(out("pq")),
      "--csv_dir", shQuote(out("csv")), "--track_output_dir", shQuote(out("tracks")),
      "--pdf_output_dir", shQuote(out("pdf")), "--log", shQuote(out("job.log"))),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(status, 0L)

  for (f in c("csv/Toyus_toyus.csv", "csv/Toyus_toyus.hotspots.csv",
              "pq/Toyus_toyus.parquet", "pq/Toyus_toyus.manifest.yaml",
              "tracks/Toyus_toyus.gff3", "tracks/Toyus_toyus.bed",
              "pdf/Toyus_toyus.hotspots.pdf")) {
    expect_true(file.exists(out(f)), info = f)
  }
  windows <- utils::read.csv(out("csv/Toyus_toyus.csv"))
  expect_equal(nrow(windows), 3L * 60L)  # 300 kb in 5 kb windows, three chromosomes
  regions <- utils::read.csv(out("csv/Toyus_toyus.hotspots.csv"))
  on_cluster <- regions$seqnames == "chrB" & regions$start <= 110000 & regions$end >= 100000
  expect_true(any(on_cluster))
  expect_equal(unique(regions$dominant_taxon[on_cluster]), "Gammaretrovirus")
})
