# testthat coverage for the ranges run manifest (ranges/exporters.R).
#
# The manifest is how a past run can be audited. It used to omit `main_probes`,
# the list that decides which probes are main and the completeness count, so
# two runs with different lists were indistinguishable afterwards.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(yaml)
})

source("../../scripts/ranges/exporters.R")

test_that("the manifest records the main probes it ran with", {
  inputs <- vapply(c("fasta", "blast", "ltrdigest", "probes", "config"),
                   function(name) {
                     path <- tempfile(name)
                     writeLines(name, path)
                     path
                   }, character(1))
  on.exit(unlink(inputs))
  opts <- list(
    seed = 67, bitscore_threshold = 0, identity_threshold = 0,
    merge_option = "virus", probe_min_length = c(GAG = 200),
    main_probes = c("POL", "GAG", "ENV"),
    agg_virus = "list", agg_label = "list", agg_probe = "list",
    agg_species = "list", agg_best_tiebreaker = "first"
  )
  out <- tempfile(fileext = ".yaml")
  on.exit(unlink(out), add = TRUE)
  emit_manifest(as.list(inputs), "test", opts, out)
  expect_equal(yaml::read_yaml(out)$options$main_probes, c("POL", "GAG", "ENV"))
})
