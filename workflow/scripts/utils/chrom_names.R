# -----------------------------------------------------------------------------
# utils / chrom_names.R
# -----------------------------------------------------------------------------
# FASTA-header normaliser, sourced by hotspot_detector.R. Returns the first
# whitespace-delimited token of each header, the same length as the input, so it
# can be assigned directly to `seqlevels(...)` or used in
# `Seqinfo(seqnames = ..., ...)`.
#
# No pattern is imposed on the name. The pipeline takes chromosome names in
# whatever shape the genome FASTA provides: BLAST uses the first header token as
# its subject id, LTRharvest does the same under `-seqids`, and
# extract_region_fasta.R keys on it "as bedtools does". Those names are what end
# up in the tracks and the catalog, so this must agree with them rather than
# reinterpret them.
#
# An earlier revision required an NCBI accession with a version suffix
# ("CM034567.1", "NC_000001.11"). That is true of NCBI-downloaded genomes and
# false of most others - Bat1K and DNA Zoo assemblies carry bare scaffold
# numbers like "02077418" - so the hotspot stage aborted on them. Worse, where a
# header began with an accession followed immediately by more characters, the
# pattern returned a DIFFERENT name from the one BLAST had written into the
# tracks, which would have mismatched silently rather than failed.

#' Chromosome names from FASTA headers.
#'
#' @param headers character vector of FASTA header lines
#' @return character vector the same length as `headers`; NA only where a header
#'   carries no name at all (empty or whitespace-only), which is a malformed
#'   FASTA the caller stops on rather than inventing an empty seqlevel for.
normalise_chrom_names <- function(headers) {
  names_only <- sub("\\s.*$", "", trimws(headers))
  names_only[!nzchar(names_only)] <- NA_character_

  unnamed_idx <- which(is.na(names_only))
  if (length(unnamed_idx) > 0L) {
    message(sprintf(
      "normalise_chrom_names(): %d of %d FASTA headers carry no name at all.",
      length(unnamed_idx), length(headers)
    ))
  }
  names_only
}
