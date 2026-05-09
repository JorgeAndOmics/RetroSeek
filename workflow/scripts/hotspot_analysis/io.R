# -----------------------------------------------------------------------------
# hotspot_analysis / io.R
# -----------------------------------------------------------------------------
# Loaders + config-readers for hotspot_detector.R. Each function returns plain
# R structures (named integer vectors, GRanges, lists) and does no analysis —
# that lives in the other modules.
#
# Mirrors the shape of range_analysis/io.R: a thin yaml::read_yaml wrapper, a
# `read_hotspot_options()` helper that centralises the `%||%` defaults, plus
# loaders for the FASTA + GFF3 inputs.

suppressMessages({
  library(yaml)
  library(GenomicRanges)
  library(Biostrings)
  library(rtracklayer)
  library(S4Vectors)
})


# Read the project config YAML and return as a nested list.
read_config <- function(path) {
  yaml::read_yaml(path)
}


#' Load genome FASTA, normalise headers to NCBI-accession tokens, and rename
#' the DNAStringSet in place so downstream `vmatchPattern` -> mask GRanges
#' carries the normalised seqlevels and aligns with windows automatically.
#'
#' Aborts loudly if any header fails to match the chrom-name pattern. The
#' previous code path silently propagated NA seqlevels into the genome
#' GRanges; here we surface the misconfiguration instead.
#'
#' Returns a list with:
#'   $seqs        — DNAStringSet (renamed in place)
#'   $seqlengths  — named integer vector keyed by normalised chrom name
#'   $headers_raw — original FASTA headers (kept for the manifest)
load_genome_for_hotspot <- function(fasta_path) {
  seqs <- Biostrings::readDNAStringSet(fasta_path)
  headers_raw <- names(seqs)
  chrom_names <- normalise_chrom_names(headers_raw)
  if (any(is.na(chrom_names))) {
    stop(sprintf(
      "load_genome_for_hotspot(): %d of %d FASTA headers in %s did not match the chromosome-name pattern. Check FASTA contents or extend utils/chrom_names.R::normalise_chrom_names().",
      sum(is.na(chrom_names)), length(chrom_names), fasta_path
    ))
  }
  names(seqs) <- chrom_names
  seqlengths <- setNames(BiocGenerics::width(seqs), chrom_names)
  list(seqs = seqs, seqlengths = seqlengths, headers_raw = headers_raw)
}


#' Load a hotspot input GFF3 (raw or valid track).
#'
#' Validates that `mcols$label` is present (this is the per-genus tag attached
#' by ranges_analysis.R; without it we cannot per-label split). Aborts with a
#' clear error if missing — easier to debug than the cryptic NULL-subset error
#' the loop would otherwise hit.
load_hits_gff <- function(gff_path) {
  hits <- rtracklayer::import(gff_path, format = "gff3")
  if (!"label" %in% colnames(S4Vectors::mcols(hits))) {
    stop(sprintf(
      "Input GFF3 %s does not contain a 'label' metadata column. The hotspot detector splits hits by retrovirus genus via mcols$label; without it, hotspot_group_split=true cannot operate.",
      gff_path
    ))
  }
  hits
}


#' Centralise hotspot-pipeline option defaults.
#'
#' Reads from `config$parameters$*` and applies safe defaults so the
#' orchestrator stays readable. Emits one logical knob per design decision in
#' the plan; default values match `data/config/config.yaml`.
read_hotspot_options <- function(config) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  p <- config$parameters
  list(
    seed                       = as.integer(p$seed                       %||% 67L),
    input                      = p$hotspot_input                          %||% "valid",
    group_split                = isTRUE(p$hotspot_group_split             %||% FALSE),
    window_size                = as.integer(p$hotspot_window_size         %||% 10000L),
    mask_size                  = as.integer(p$hotspot_mask_size           %||% 20L),
    mask_mismatch              = as.integer(p$hotspot_mask_mismatch       %||% 2L),
    permutations               = as.integer(p$hotspot_permutations        %||% 100L),
    pvalue_threshold           = as.numeric(p$hotspot_pvalue_threshold    %||% 0.05),
    validate_permutation       = isTRUE(p$hotspot_validate_permutation    %||% FALSE),
    min_hits_per_window        = as.integer(p$hotspot_min_hits_per_window %||% 2L),
    merge_adjacent             = isTRUE(p$hotspot_merge_adjacent          %||% TRUE),
    merge_gap                  = as.integer(p$hotspot_merge_gap           %||% 0L),
    strata_by_chromosome       = isTRUE(p$hotspot_strata_by_chromosome    %||% TRUE),
    unplaced_min_factor        = as.integer(p$hotspot_unplaced_min_factor %||% 10L)
  )
}
