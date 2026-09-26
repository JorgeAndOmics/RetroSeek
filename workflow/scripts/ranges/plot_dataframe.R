# -----------------------------------------------------------------------------
# ranges / plot_dataframe.R
# -----------------------------------------------------------------------------
# Build the per-row tibble consumed by plot2sort.R / circle_plot_generator.R /
# pair_detector.R. Replaces the pre-refactor three-file fanout (main /
# accessory / full) with a single tibble carrying a `probe_type` column.

suppressMessages({
  library(GenomicRanges)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
})


# Convert a GRanges (typically the per-(probe, virus) reduced gr_virus) into
# a tibble suitable for plot scripts. The `virus` column may be a `concatenate`
# string ("BLV; FeLV") or a single value depending on aggregation strategy;
# this function explodes such rows so each output row carries a single virus
# identity, then attaches Label / Abbreviation by joining against probe metadata.
build_plot_dataframe <- function(gr_virus, probe_df_sum, main_probes,
                                 agg_virus_strategy, agg_concat_separator = "; ") {
  if (length(gr_virus) == 0L) {
    return(tibble::tibble(
      seqnames = character(0), start = integer(0), end = integer(0),
      width = integer(0), strand = character(0),
      probe = character(0), virus = character(0), label = character(0),
      species = character(0), abbreviation = character(0),
      probe_type = character(0)
    ))
  }

  df <- as.data.frame(gr_virus, stringsAsFactors = FALSE)
  # `virus` arrives in one of three shapes, depending on how it was aggregated:
  #   * a single value            - `best` / `first`
  #   * a separator-joined string - `concatenate`, and ALSO `list` through the
  #     plyranges reduce path (see "IMPORTANT - why list is concatenated here"
  #     in range_aggregation_strategies.R)
  #   * a genuine list-column     - `list` when aggregate_values is called direct
  # Explode the last two so every row carries one virus identity. Miss this and
  # the match() below returns NA, so label and abbreviation go silently empty for
  # every multi-virus locus - and ranges/io.R defaults agg_virus to "list".
  if (is.list(df$virus)) {
    df <- tidyr::unnest(df, virus)
  } else if (agg_virus_strategy %in% c("concatenate", "list")) {
    df <- tidyr::separate_rows(df, virus, sep = agg_concat_separator)
  }

  df <- df %>%
    dplyr::mutate(
      label        = probe_df_sum$Label[match(virus, probe_df_sum$Name)],
      abbreviation = probe_df_sum$Abbreviation[match(virus, probe_df_sum$Name)],
      probe_type   = ifelse(probe %in% main_probes, "main", "accessory")
    )
  tibble::as_tibble(df)
}


# Tag a GRanges with a `probe_category` mcols column: "main" / "accessory" /
# "mixed" depending on whether the (possibly multi-value) probe column lies
# entirely in the main set, entirely outside, or spans both.
attach_probe_category <- function(gr, main_set, concat_separator = "; ") {
  if (length(gr) == 0L) {
    S4Vectors::mcols(gr)$probe_category <- character(0)
    return(gr)
  }
  probes <- .probe_lists(S4Vectors::mcols(gr)$probe, concat_separator)
  S4Vectors::mcols(gr)$probe_category <- .probe_categories(probes, main_set)
  gr
}

# The probe column as one CharacterList (a joined string is split on
# `concat_separator`), empty names dropped. Built once for the whole column:
# indexing a CharacterList row by row costs about a millisecond per row.
.probe_lists <- function(probe_col, concat_separator) {
  probes <- if (inherits(probe_col, "CharacterList")) {
    probe_col
  } else {
    IRanges::CharacterList(strsplit(as.character(probe_col), concat_separator, fixed = TRUE))
  }
  # nzchar() has no CharacterList method; nchar() does, but reads NA as NA
  # where nzchar(NA) is TRUE, so an NA probe is kept explicitly.
  probes[nchar(probes) != 0L | is.na(probes)]
}

# "main" when every probe of a range is in `main_set`, "accessory" when none is,
# "mixed" otherwise; NA for a range with no probe.
.probe_categories <- function(probes, main_set) {
  n_all <- lengths(probes)
  n_main <- sum(probes %in% main_set)
  category <- ifelse(n_main == n_all, "main", ifelse(n_main == 0L, "accessory", "mixed"))
  category[n_all == 0L] <- NA_character_
  category
}
