# =============================================================================
# plot2sort/io.R - config + parquet ingest + logged ggsave
# =============================================================================
# Anything that reads from / writes to disk lives here. The logging functions
# (log_section, log_info) come from utils/log.R, which the entry script sources.


# Discover every {genome}.final_loci.parquet under `input_dir` and concatenate
# them into one tibble. The ranges-analysis tables directory holds several
# table types per genome; plot2sort consumes only the final-tier loci. Reports
# per-input row counts and the total so an empty input surfaces immediately in
# the verbose log instead of failing 300 lines later.
load_plot_dataframes <- function(input_dir) {
  parquet_files <- list.files(input_dir, pattern = "\\.final_loci\\.parquet$",
                              full.names = TRUE)
  if (length(parquet_files) == 0L) {
    stop("plot2sort: no per-genome *.final_loci.parquet files found in ", input_dir)
  }
  log_section(sprintf("Loading %d plot dataframes from %s",
                      length(parquet_files), input_dir))
  frames <- purrr::map(parquet_files, function(path) {
    df <- arrow::read_parquet(path)
    log_info("%-40s  %d rows", basename(path), nrow(df))
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
      paste(
        "plot2sort: %s missing required column(s): %s. Rebuild ranges_analysis",
        "outputs: the plot dataframe contract may have changed."
      ),
      source_label, paste(missing, collapse = ", ")
    ))
  }
  invisible(df)
}


# Write one ggplot to an image file and log the action. Only the README demo
# figures use it (PNG, for GitHub); every pipeline stage writes a multi-page PDF
# through save_stage_pdf() in style.R. `dims` (a list with $w + $h) overrides
# the base canvas per call.
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
  log_info("wrote %-44s  (%5.1f x %5.1f in, %d dpi)", name, w, h, dpi)
  invisible(file.path(output_dir, name))
}
