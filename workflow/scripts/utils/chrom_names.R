# -----------------------------------------------------------------------------
# utils / chrom_names.R
# -----------------------------------------------------------------------------
# Shared FASTA-header normaliser used by hotspot_detector.R and
# circle_plot_generator.R. Pulls a leading NCBI-accession token (e.g.
# "CM034567.1", "NC_000001.11") out of a full header line, returning a vector
# the same length as the input so it can be assigned directly to
# `seqlevels(...)` or used in `Seqinfo(seqnames = ..., ...)`.
#
# Behaviour mirrors the previous duplicated regex in both call sites, with one
# difference: when any header fails to match, we emit an explicit `message()`
# listing the offenders so the silent-NA failure mode is visible to callers
# rather than surfacing later as a cryptic `Seqinfo` error.

suppressMessages({
  library(stringr)
})


#' Extract NCBI-style accessions from FASTA headers.
#'
#' @param headers character vector of FASTA header lines
#' @param pattern regex for the accession token; the default matches genus-style
#'   tokens like "CM######.#" and "NC_######.#"
#' @return character vector same length as `headers`; positions where the
#'   pattern did not match contain NA
normalise_chrom_names <- function(headers,
                                  pattern = "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}") {
  normalised <- stringr::str_extract(headers, pattern)
  unmatched_idx <- which(is.na(normalised))
  if (length(unmatched_idx) > 0L) {
    sample_unmatched <- utils::head(headers[unmatched_idx], 5L)
    message(sprintf(
      "normalise_chrom_names(): %d of %d FASTA headers did not match pattern '%s'. First %d: %s",
      length(unmatched_idx), length(headers), pattern,
      length(sample_unmatched), paste(sample_unmatched, collapse = " | ")
    ))
  }
  normalised
}
