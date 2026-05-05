# -----------------------------------------------------------------------------
# range_analysis / io.R
# -----------------------------------------------------------------------------
# Loaders + config-readers for the ranges_analysis pipeline. Each function
# returns plain R structures (tibbles, named vectors, lists, GRanges) and
# does no analysis — that lives in the other modules.

suppressMessages({
  library(yaml)
  library(arrow)
  library(GenomicRanges)
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(readr)
})


# Read the project config YAML and return as a nested list.
read_config <- function(path) {
  yaml::read_yaml(path)
}


# Read named-integer chrom lengths from a FASTA file.
load_chrom_lengths <- function(fasta_path) {
  fa <- Biostrings::readDNAStringSet(fasta_path)
  setNames(BiocGenerics::width(fa), names(fa))
}


# Load BLAST hits as a tibble. Schema is fixed by upstream species_segmenter.R.
load_blast_parquet <- function(path) {
  arrow::read_parquet(path)
}


# Load LTRdigest GFF3 as a GRanges.
load_ltrdigest_gff3 <- function(path) {
  rtracklayer::import(path, format = "gff3")
}


# Load probes CSV. Returns a list with:
#   $df       — full tibble (one row per probe)
#   $df_sum   — distinct (Name, Label, Abbreviation) tibble
load_probes <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  df_sum <- dplyr::distinct(df, Name, Label, Abbreviation)
  list(df = df, df_sum = df_sum)
}


# Load probe protein lengths from the post-fetch enriched probe_dict parquet
# (the file produced by obj2dict.py after Entrez retrieval populates
# `genbank_seq`). Returns a named integer vector keyed by
# `paste(virus, probe, sep = "|")` — the same key shape used by
# build_blast_gr to attach per-hit query_coverage.
#
# When the same (virus, probe) pair has multiple proteins (e.g. HIV-1 VIF
# carries two accessions), the longest sequence is kept — the most
# conservative coverage estimate, deterministic, and consistent with the
# pre-fix code that arbitrarily picked one length per group.
load_probe_lengths <- function(probe_dict_path) {
  pdict <- arrow::read_parquet(probe_dict_path)
  agg <- pdict %>%
    dplyr::mutate(.seq_len = nchar(genbank_seq)) %>%
    dplyr::group_by(virus, probe) %>%
    dplyr::summarise(.probe_len = max(.seq_len, na.rm = TRUE), .groups = "drop")
  setNames(agg$.probe_len, paste(agg$virus, agg$probe, sep = "|"))
}


# Read the `parameters` block + `domains` block + aggregation block from config
# and return a flat opts list with safe defaults applied. Centralising the
# `%||%` defaults here keeps the orchestrator readable.
read_pipeline_options <- function(config) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  agg <- config$parameters$aggregation %||% list()
  list(
    probe_min_length    = unlist(config$parameters$probe_min_length),
    bitscore_threshold  = as.numeric(config$parameters$bitscore_threshold %||% 0),
    identity_threshold  = as.numeric(config$parameters$identity_threshold %||% 0),
    ltr_resize          = as.numeric(config$parameters$ltr_resize         %||% 0),
    merge_option        = config$parameters$merge_option %||% "virus",
    main_probes         = unique(config$parameters$main_probes),
    domains             = config$domains,
    agg_virus           = agg$virus            %||% "list",
    agg_label           = agg$label            %||% "list",
    agg_probe           = agg$probe            %||% "list",
    agg_species         = agg$species          %||% "first",
    agg_best_tiebreaker = agg$best_tiebreaker  %||% "bitscore",
    agg_concat_separator = agg$concat_separator %||% "; ",
    agg_strict_marker   = agg$strict_marker    %||% "ambiguous"
  )
}
