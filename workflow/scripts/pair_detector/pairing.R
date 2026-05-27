# =============================================================================
# Pure pairing helpers for pair_detector.R
# =============================================================================
# Extracted so the probe-pairing logic can be unit-tested without GFF3 I/O or
# plotting. Sourced by both pair_detector.R and
# workflow/tests/testthat/test-pair_detector.R. Assumes GenomicRanges / IRanges
# / S4Vectors / dplyr / tibble are attached by the caller. No side effects.

# 5' coordinate of a feature given its strand (used for biological distance).
get_5prime <- function(start, end, strand) {
  dplyr::case_when(
    strand == "+" ~ start,
    strand == "-" ~ end,
    strand == "*" ~ start, # unknown strand treated as +
    TRUE ~ start
  )
}

# Find every (probe_to_pair, other-probe) pair within `max_gap` bp.
#
#   ranges        : GRanges with mcols probe/label/virus (name/ID optional).
#   probe_to_pair : the "self" probe paired against all other probes.
#   max_gap       : maximum gap (bp) passed to findOverlaps.
#
# Returns a tibble, one row per self/other overlap, carrying coordinate, label,
# virus and probe columns plus coord_distance / bio_distance / span_distance.
find_pairs <- function(ranges, probe_to_pair, max_gap) {
  self_ranges <- ranges[ranges$probe == probe_to_pair]
  other_ranges <- ranges[ranges$probe != probe_to_pair]

  ov <- findOverlaps(self_ranges, other_ranges, maxgap = max_gap)
  qh <- queryHits(ov)
  sh <- subjectHits(ov)

  pairs_df <- tibble::tibble(
    self.seqname  = as.character(seqnames(self_ranges)[qh]),
    other.seqname = as.character(seqnames(other_ranges)[sh]),

    self.start  = start(self_ranges)[qh],
    other.start = start(other_ranges)[sh],

    self.end  = end(self_ranges)[qh],
    other.end = end(other_ranges)[sh],

    self.strand  = as.character(strand(self_ranges)[qh]),
    other.strand = as.character(strand(other_ranges)[sh]),

    self.label  = mcols(self_ranges)$label[qh],
    other.label = mcols(other_ranges)$label[sh],

    self.virus  = mcols(self_ranges)$virus[qh],
    other.virus = mcols(other_ranges)$virus[sh],

    self.probe  = mcols(self_ranges)$probe[qh],
    other.probe = mcols(other_ranges)$probe[sh],
  )

  # Conditionally carry `name` (domain) and `ID` if the GFF3 provided them.
  if ("name" %in% names(mcols(self_ranges))) {
    pairs_df$self.domain <- mcols(self_ranges)$name[qh]
    pairs_df$other.domain <- mcols(other_ranges)$name[sh]
  }
  if ("ID" %in% names(mcols(self_ranges))) {
    pairs_df$self.ID <- mcols(self_ranges)$ID[qh]
    pairs_df$other.ID <- mcols(other_ranges)$ID[sh]
  }

  # Coordinate distance between the two start positions.
  pairs_df <- dplyr::mutate(pairs_df, coord_distance = abs(self.start - other.start))

  # Biological distance between strand-aware 5' ends.
  pairs_df$self.5prime <- get_5prime(pairs_df$self.start, pairs_df$self.end, pairs_df$self.strand)
  pairs_df$other.5prime <- get_5prime(pairs_df$other.start, pairs_df$other.end, pairs_df$other.strand)
  pairs_df$bio_distance <- abs(pairs_df$self.5prime - pairs_df$other.5prime)
  pairs_df <- dplyr::select(pairs_df, -self.5prime, -other.5prime)

  # Span distance ignoring strand (outermost extent of the pair).
  pairs_df <- dplyr::mutate(
    pairs_df,
    span_distance = abs(pmax(self.end, other.end) - pmin(self.start, other.start))
  )

  pairs_df
}
