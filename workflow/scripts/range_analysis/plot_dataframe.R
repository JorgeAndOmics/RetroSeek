# -----------------------------------------------------------------------------
# range_analysis / plot_dataframe.R
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
  # When virus is a concatenated string, split into multiple rows so each row
  # carries a single virus identity.
  if (identical(agg_virus_strategy, "concatenate")) {
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
  probe_col <- S4Vectors::mcols(gr)$probe
  cat_chr <- vapply(seq_along(probe_col), function(i) {
    probes <- if (inherits(probe_col, "CharacterList")) {
      as.character(probe_col[[i]])
    } else {
      strsplit(as.character(probe_col[[i]]), concat_separator, fixed = TRUE)[[1]]
    }
    probes <- probes[nzchar(probes)]
    if (length(probes) == 0L) return(NA_character_)
    in_main <- probes %in% main_set
    if (all(in_main))   return("main")
    if (!any(in_main))  return("accessory")
    "mixed"
  }, character(1))
  S4Vectors::mcols(gr)$probe_category <- cat_chr
  gr
}
