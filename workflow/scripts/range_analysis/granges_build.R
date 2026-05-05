# -----------------------------------------------------------------------------
# range_analysis / granges_build.R
# -----------------------------------------------------------------------------
# Convert tabular / GFF3 inputs into well-typed GRanges objects with the
# metadata the downstream phases expect. Pure construction — no filtering or
# aggregation here.

suppressMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})


# Build a GRanges from the BLAST tibble produced upstream by species_segmenter.
# Columns picked into mcols are exactly those the rest of the pipeline reads;
# unused columns from the parquet (genbank_*, hsp_query, hsp_sbjct, etc.) are
# intentionally dropped to keep the GRanges lean.
#
# Trusts the upstream contract: `accession` is already a clean NCBI accession
# (e.g. "CM138268.1"), having been split from the BLAST hit_def by seq_utils
# / obj2dict. This matches the seqid LTRdigest emits, so findOverlaps against
# LTRdigest GRanges work without per-stage stripping.
#
# `probe_lengths` is a named integer vector keyed by paste(virus, probe, '|')
# (built by load_probe_lengths). When supplied, attaches per-hit
# `query_coverage = pmin(1, hsp_align_length / probe_total_length)` as an
# mcols column. The reduction stage downstream just sums the per-hit values
# across each merged range.
build_blast_gr <- function(blast_df, probe_lengths = NULL) {
  gr <- GenomicRanges::GRanges(
    seqnames = blast_df$accession,
    ranges   = IRanges::IRanges(
      start = blast_df$hsp_sbjct_start,
      end   = blast_df$hsp_sbjct_end
    ),
    strand   = blast_df$strand
  )
  S4Vectors::mcols(gr)$label        <- blast_df$label
  S4Vectors::mcols(gr)$virus        <- blast_df$virus
  S4Vectors::mcols(gr)$bitscore     <- blast_df$hsp_bits
  # Percent identity over the alignment length.
  S4Vectors::mcols(gr)$identity     <- (blast_df$hsp_identity / blast_df$hsp_align_length) * 100
  S4Vectors::mcols(gr)$species      <- blast_df$species
  S4Vectors::mcols(gr)$probe        <- blast_df$probe
  S4Vectors::mcols(gr)$evalue       <- blast_df$hsp_evalue
  S4Vectors::mcols(gr)$align_length <- blast_df$hsp_align_length
  S4Vectors::mcols(gr)$query_start  <- blast_df$hsp_query_start
  S4Vectors::mcols(gr)$query_end    <- blast_df$hsp_query_end

  if (is.null(probe_lengths)) {
    S4Vectors::mcols(gr)$query_coverage <- NA_real_
  } else {
    key <- paste(blast_df$virus, blast_df$probe, sep = "|")
    plen <- probe_lengths[key]
    S4Vectors::mcols(gr)$query_coverage <- pmin(1, blast_df$hsp_align_length / plen)
  }

  gr
}


# Pull the LTR_retrotransposon features (the parent ERVs) from the LTRdigest
# GFF3. Optionally pad each ERV by `resize_bp` on both sides so downstream
# overlap detection has tolerance.
extract_retrotransposons <- function(ltr_data, resize_bp = 0L) {
  retros <- ltr_data[ltr_data$type == "LTR_retrotransposon"]
  if (length(retros) == 0L || resize_bp == 0L) return(retros)
  # Resize both ends by resize_bp without going below seqlength bounds.
  GenomicRanges::resize(
    retros,
    width = BiocGenerics::width(retros) + 2L * resize_bp,
    fix = "center"
  )
}


# Pull repeat_region features (the outer-bound LTR retrotransposon containers,
# used for some overlap analyses).
extract_repeat_regions <- function(ltr_data) {
  ltr_data[ltr_data$type == "repeat_region"]
}


# Build a domain-name → probe regex map from the config$domains list.
# Each probe's pattern alternates its domain regexes via "|".
build_domain_map <- function(config_domains) {
  vapply(config_domains, function(domain_regexes) {
    paste(domain_regexes, collapse = "|")
  }, character(1))
}


# Given a single domain `name` (the GFF3 `name` attribute, typically a Pfam
# domain), return the probe whose regex matches first. Among multiple matches,
# prefer the most-specific (longest probe name).
assign_probe_to_domain_name <- function(domain_name, domain_map) {
  matched <- vapply(names(domain_map), function(probe) {
    grepl(domain_map[[probe]], domain_name, ignore.case = TRUE)
  }, logical(1))
  if (!any(matched)) return(NA_character_)
  candidates <- names(domain_map)[matched]
  candidates[order(nchar(candidates), decreasing = TRUE)[1]]
}


# Pull domain features (those with a non-NA `name` attribute) and assign each
# a probe label via the domain_map. Domains without a probe assignment are
# dropped — they cannot contribute to validation.
extract_domains_with_probes <- function(ltr_data, domain_map) {
  doms <- ltr_data[!is.na(ltr_data$name)]
  if (length(doms) == 0L) {
    S4Vectors::mcols(doms)$probe <- character(0)
    return(doms)
  }
  S4Vectors::mcols(doms)$probe <- vapply(
    doms$name,
    assign_probe_to_domain_name,
    domain_map = domain_map,
    FUN.VALUE = character(1)
  )
  doms[!is.na(doms$probe)]
}


# Build a GRanges of flanking LTRs (one or two per parent retrotransposon).
# LTRdigest emits LTR features as `long_terminal_repeat` rows with a Parent
# attribute pointing to their enclosing retrotransposon. We keep both arms,
# tag them as "L" (5') or "R" (3') in strand-aware orientation, and propagate
# the parent ERV's identifier so downstream analyses can reassemble pairs.
extract_flanking_ltrs <- function(ltr_data) {
  ltrs <- ltr_data[ltr_data$type == "long_terminal_repeat"]
  if (length(ltrs) == 0L) return(ltrs)
  # Sort by chromosome then start so that within each parent, the lower-start
  # element is "left" on the + strand and the higher-start is "left" on the
  # - strand (strand-aware labeling).
  ltrs <- BiocGenerics::sort(ltrs, ignore.strand = TRUE)
  parent <- as.character(ltrs$Parent)
  # First/second-occurrence flag within parent groups.
  occurrence <- ave(seq_along(parent), parent, FUN = seq_along)
  arm <- ifelse(
    BiocGenerics::strand(ltrs) == "-",
    ifelse(occurrence == 1L, "R", "L"),
    ifelse(occurrence == 1L, "L", "R")
  )
  S4Vectors::mcols(ltrs)$arm <- arm
  ltrs
}
