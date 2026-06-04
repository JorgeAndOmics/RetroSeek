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


#' Load a hotspot input GFF3 (original / valid / erv_like track).
#'
#' Validates that `mcols$label` is present (this is the per-genus tag attached
#' by ranges_analysis.R; without it we cannot per-label split). Aborts with a
#' clear error if missing — easier to debug than the cryptic NULL-subset error
#' the loop would otherwise hit.
load_hits_gff <- function(gff_path) {
  hits <- rtracklayer::import(gff_path, format = "gff3")
  if (!"label" %in% colnames(S4Vectors::mcols(hits))) {
    stop(sprintf(
      "Input GFF3 %s does not contain a 'label' metadata column. The hotspot detector splits hits by retrovirus genus via mcols$label; without it, group_split=true cannot operate.",
      gff_path
    ))
  }
  hits
}


#' Fail loud when the input track and the genome FASTA come from different
#' assemblies / accession namespaces (the classic GenBank `CM*`/`JA*` vs RefSeq
#' `NC_*`/`NW_*` mismatch). Without this, zero hits overlap the windows and the
#' run silently degrades to `insufficient_data` — indistinguishable from a
#' genome that genuinely has no ERV clusters.
#'
#' Aborts when fewer than `min_frac` of hits sit on a contig present in the
#' genome. `seqlengths` is the named integer vector from
#' `load_genome_for_hotspot()` (names = normalised genome contigs).
assert_hits_on_genome <- function(hits, seqlengths, min_frac = 0.01) {
  hit_chr <- as.character(GenomicRanges::seqnames(hits))
  if (length(hit_chr) == 0L) {
    return(invisible(0))
  }
  genome_chr <- names(seqlengths)
  n_on <- sum(hit_chr %in% genome_chr)
  frac <- n_on / length(hit_chr)
  if (frac < min_frac) {
    ex_hit <- utils::head(unique(hit_chr), 3L)
    ex_gen <- utils::head(genome_chr, 3L)
    stop(sprintf(
      paste0(
        "Only %d/%d (%.2f%%) input hits map to a genome contig. The input track and ",
        "--fasta almost certainly come from different assemblies or accession ",
        "namespaces (e.g. RefSeq NC_*/NW_* vs GenBank CM*/JA*). Track seqnames e.g.: ",
        "%s ; genome seqnames e.g.: %s. Re-run the track and the FASTA against the ",
        "same assembly."
      ),
      n_on, length(hit_chr), 100 * frac,
      paste(ex_hit, collapse = ", "), paste(ex_gen, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(frac)
}


#' Centralise hotspot-pipeline option defaults.
#'
#' Reads the operational knobs from the dedicated top-level `config$hotspot`
#' section and the global RNG seed from `config$parameters$seed`. Applies safe
#' defaults so the orchestrator stays readable. Defaults match
#' `data/config/config.yaml`.
read_hotspot_options <- function(config) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  h <- config$hotspot
  list(
    seed                 = as.integer(config$parameters$seed     %||% 67L),
    input                = h$input                                %||% "erv_like",
    group_split          = isTRUE(h$group_split                   %||% FALSE),
    window_size          = as.integer(h$window_size              %||% 500000L),
    mask_size            = as.integer(h$mask_size                %||% 20L),
    mask_mismatch        = as.integer(h$mask_mismatch            %||% 3L),
    pvalue_threshold     = as.numeric(h$pvalue_threshold         %||% 0.05),
    min_hits             = as.integer(h$min_hits                 %||% 2L),
    merge_gap            = as.integer(h$merge_gap                %||% 0L),
    strata_by_chromosome = isTRUE(h$strata_by_chromosome          %||% TRUE),
    unplaced_min_factor  = as.integer(h$unplaced_min_factor      %||% 10L)
  )
}
