# =============================================================================
# stage_plot_generator/io.R — stage-parquet + manifest ingest
# =============================================================================
# Stage-specific disk I/O. `save_plot` / `verify_required_columns` are reused
# from plot2sort/io.R (sourced by the orchestrator). `log_section` is defined
# in the orchestrator and resolved via lexical scope at call time — same
# convention as plot2sort.R / ranges_analysis.R.


# Load the three stage-dataframe parquet types from one directory. Files are
# suffix-discriminated: {genome}_hits.parquet, {genome}_ltr.parquet,
# {genome}_reduced.parquet. Returns a named list of three bound tibbles, each
# carrying a `genome` column. A missing type yields an empty tibble so the
# downstream builders fall through to their empty_plot() guards.
load_stage_dataframes <- function(input_dir) {
  read_suffixed <- function(suffix) {
    pattern <- sprintf("_%s\\.parquet$", suffix)
    files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0L) {
      log_section(sprintf("  (no *_%s.parquet found in %s)", suffix, input_dir))
      return(tibble::tibble())
    }
    frames <- purrr::map(files, function(path) {
      df <- arrow::read_parquet(path)
      genome <- sub(sprintf("_%s$", suffix), "",
                    tools::file_path_sans_ext(basename(path)))
      df$genome <- genome
      log_section(sprintf("  %-44s  %d rows", basename(path), nrow(df)))
      df
    })
    dplyr::bind_rows(frames)
  }
  log_section(sprintf("Loading stage dataframes from %s", input_dir))
  list(
    hits    = read_suffixed("hits"),
    ltr     = read_suffixed("ltr"),
    reduced = read_suffixed("reduced")
  )
}


# Read every {genome}.yaml manifest and pull the per-stage counts into a tidy
# long tibble: (genome, stage, count). `stage` is a factor in pipeline order
# so the funnel plots render homology → valid left-to-right without re-sorting.
load_manifest_counts <- function(manifest_dir) {
  # Pipeline-ordered stages of the homology → valid refinement funnel. Names
  # are manifest `counts:` keys; values are the human-readable axis labels.
  stage_keys <- c(
    filtered_blast_hits  = "homology hits",
    first_reduced_ranges = "first-reduced",
    candidate_ranges     = "candidate",
    valid_ranges         = "valid"
  )
  files <- list.files(manifest_dir, pattern = "\\.yaml$", full.names = TRUE)
  if (length(files) == 0L) {
    log_section(sprintf("  (no manifests found in %s)", manifest_dir))
    return(tibble::tibble(
      genome = character(0),
      stage  = factor(character(0), levels = unname(stage_keys)),
      count  = integer(0)
    ))
  }
  log_section(sprintf("Loading %d manifests from %s",
                      length(files), manifest_dir))
  rows <- purrr::map(files, function(path) {
    y      <- yaml::read_yaml(path)
    genome <- tools::file_path_sans_ext(basename(path))
    counts <- y$counts %||% list()
    tibble::tibble(
      genome = genome,
      stage  = unname(stage_keys),
      count  = vapply(names(stage_keys),
                      function(k) as.integer(counts[[k]] %||% NA_integer_),
                      integer(1))
    )
  })
  dplyr::bind_rows(rows) %>%
    dplyr::mutate(stage = factor(stage, levels = unname(stage_keys)))
}
