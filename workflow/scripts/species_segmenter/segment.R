# =============================================================================
# Pure helpers for species_segmenter.R
# =============================================================================
# Extracted so the partition / backfill logic can be unit-tested without the
# surrounding Parquet/CSV I/O. Sourced by both species_segmenter.R and
# workflow/tests/testthat/test-species_segmenter.R. No side effects, base R only.

# Partition a hit table by probe class.
#
#   df               : data frame / tibble with a character `probe` column.
#   main_probe_names : probe names treated as "main" (config$parameters$main_probes).
#
# Returns list(main = <rows whose probe is in main_probe_names>,
#              accessory = <the remaining rows>). Row order is preserved within
# each part; NA probes fall into `accessory` (NA %in% x is FALSE).
segment_by_probe <- function(df, main_probe_names) {
  is_main <- df[["probe"]] %in% main_probe_names
  list(
    main = df[is_main, , drop = FALSE],
    accessory = df[!is_main, , drop = FALSE]
  )
}

# Parse the Snakefile-supplied comma-separated SPECIES list into a clean
# character vector: split on commas, trim whitespace, drop empty tokens.
parse_species_list <- function(species_arg) {
  parts <- strsplit(species_arg, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

# Configured species that produced no hits, so the orchestrator can write an
# empty placeholder table for each — keeping the B1 checkpoint DAG's declared
# per-genome outputs materialised even when a genome has zero hits.
species_to_backfill <- function(all_species, written) {
  setdiff(all_species, written)
}
