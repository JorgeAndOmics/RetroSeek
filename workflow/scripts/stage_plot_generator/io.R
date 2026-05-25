# =============================================================================
# stage_plot_generator/io.R — stage-parquet + manifest ingest
# =============================================================================
# Stage-specific disk I/O. `save_plot` / `verify_required_columns` are reused
# from plot2sort/io.R (sourced by the orchestrator). `log_section` is defined
# in the orchestrator and resolved via lexical scope at call time — same
# convention as plot2sort.R / ranges_analysis.R.


# Load the three stage-dataframe parquet types from one directory. Files are
# name-discriminated: {genome}.homology_loci.parquet,
# {genome}.ltr_structure.parquet, {genome}.reduction_multiplicity.parquet.
# Returns a named list of three bound tibbles, each carrying a `genome` column.
# A missing type yields an empty tibble so the downstream builders fall through
# to their empty_plot() guards.
load_stage_dataframes <- function(input_dir) {
  read_table <- function(table) {
    pattern <- sprintf("\\.%s\\.parquet$", table)
    files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0L) {
      log_section(sprintf("  (no *.%s.parquet found in %s)", table, input_dir))
      return(tibble::tibble())
    }
    frames <- purrr::map(files, function(path) {
      df <- arrow::read_parquet(path)
      genome <- sub(sprintf("\\.%s$", table), "",
                    tools::file_path_sans_ext(basename(path)))
      df$genome <- genome
      log_section(sprintf("  %-44s  %d rows", basename(path), nrow(df)))
      df
    })
    dplyr::bind_rows(frames)
  }
  log_section(sprintf("Loading stage dataframes from %s", input_dir))
  list(
    hits    = read_table("homology_loci"),
    ltr     = read_table("ltr_structure"),
    reduced = read_table("reduction_multiplicity")
  )
}


# Read every {genome}.counts.parquet and pull the refinement-funnel stages into
# a tidy long tibble: (genome, stage, count). `stage` is a factor in pipeline
# order so the funnel plots render homology → valid left-to-right without
# re-sorting. The counts table is long-form with columns `metric` + `value`.
load_counts_table <- function(input_dir) {
  # Pipeline-ordered stages of the homology → valid refinement funnel. Names
  # are counts-table `metric` keys; values are the human-readable axis labels.
  stage_keys <- c(
    filtered_blast_hits  = "homology hits",
    first_reduced_ranges = "first-reduced",
    candidate_ranges     = "candidate",
    valid_ranges         = "valid"
  )
  files <- list.files(input_dir, pattern = "\\.counts\\.parquet$",
                      full.names = TRUE)
  if (length(files) == 0L) {
    log_section(sprintf("  (no *.counts.parquet found in %s)", input_dir))
    return(tibble::tibble(
      genome = character(0),
      stage  = factor(character(0), levels = unname(stage_keys)),
      count  = integer(0)
    ))
  }
  log_section(sprintf("Loading %d counts tables from %s",
                      length(files), input_dir))
  rows <- purrr::map(files, function(path) {
    df     <- arrow::read_parquet(path)
    genome <- sub("\\.counts$", "", tools::file_path_sans_ext(basename(path)))
    lookup <- setNames(df$value, df$metric)
    tibble::tibble(
      genome = genome,
      stage  = unname(stage_keys),
      count  = vapply(names(stage_keys),
                      function(k) as.integer(lookup[k]), integer(1))
    )
  })
  dplyr::bind_rows(rows) %>%
    dplyr::mutate(stage = factor(stage, levels = unname(stage_keys)))
}
