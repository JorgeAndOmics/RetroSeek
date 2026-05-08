# =============================================================================
# plot2sort/io.R — config + parquet ingest + logged ggsave
# =============================================================================
# Anything that reads from / writes to disk lives here. Logging hooks
# (`log_section`) are defined in the orchestrator and resolved via R's lexical
# scope at call time — same convention as `ranges_analysis.R`.


# Discover every {genome}.parquet under `input_dir` and concatenate them into
# one tibble. Reports per-input row counts and the total so an empty input
# surfaces immediately in the verbose log instead of failing 300 lines later.
load_plot_dataframes <- function(input_dir) {
  parquet_files <- list.files(input_dir, pattern = "\\.parquet$",
                              full.names = TRUE)
  if (length(parquet_files) == 0L) {
    stop("plot2sort: no per-genome plot parquet files found in ", input_dir)
  }
  log_section(sprintf("Loading %d plot dataframes from %s",
                      length(parquet_files), input_dir))
  frames <- purrr::map(parquet_files, function(path) {
    df <- arrow::read_parquet(path)
    log_section(sprintf("  %-40s  %d rows",
                        basename(path), nrow(df)))
    df
  })
  out <- dplyr::bind_rows(frames)
  log_section(sprintf("Combined: %d rows across %d genomes",
                      nrow(out), length(parquet_files)))
  out
}


# Single guard that lists *all* missing columns at once instead of stopping
# at the first one. Saves a round-trip if multiple required columns are
# absent (typically the symptom of a stale parquet predating a schema change).
verify_required_columns <- function(df, required_cols, source_label = "input") {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0L) {
    stop(sprintf(
      "plot2sort: %s missing required column(s): %s. Rebuild ranges_analysis outputs — the plot dataframe contract may have changed.",
      source_label, paste(missing, collapse = ", ")
    ))
  }
  invisible(df)
}


# Replace the per-row `species` accession code with the readable species name
# from the YAML config map. The result keeps the column name `species`; the
# accession code is dropped.
attach_species_name <- function(df, species_map) {
  species_df <- tibble::tibble(
    species      = names(species_map),
    species_name = unname(unlist(species_map))
  )
  df %>%
    dplyr::left_join(species_df, by = "species") %>%
    dplyr::select(-species) %>%
    dplyr::rename(species = species_name)
}


# Write a ggplot to disk and log the action. `dims` (a list with $w + $h)
# overrides the config defaults per call; this is how the auto-scaled
# species-axis plots ask for a wider canvas at large N.
save_plot <- function(name, plot, output_dir,
                      dims = NULL, base_w, base_h, dpi) {
  w <- if (!is.null(dims) && !is.null(dims$w)) dims$w else base_w
  h <- if (!is.null(dims) && !is.null(dims$h)) dims$h else base_h
  ggplot2::ggsave(
    filename = file.path(output_dir, name),
    plot     = plot,
    width    = w,
    height   = h,
    dpi      = dpi
  )
  log_section(sprintf("  wrote %-44s  (%5.1f x %5.1f in, %d dpi)",
                      name, w, h, dpi))
  invisible(file.path(output_dir, name))
}
