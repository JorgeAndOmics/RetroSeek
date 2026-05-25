# =============================================================================
# erv_like_plot_generator/io.R — erv_like parquet ingest
# =============================================================================
# Disk I/O for the erv-like panel. `save_plot` / `verify_required_columns` /
# `empty_plot` etc. are reused from plot2sort/io.R + helpers.R (sourced by the
# orchestrator). `log_section` is defined in the orchestrator and resolved via
# lexical scope at call time — same convention as the other plot generators.


# Read every {genome}.<table>.parquet under `input_dir`, bind into one tibble,
# and tag each row with its genome. Returns an empty tibble when none are found
# so downstream builders fall through to their empty_plot() guards.
.read_erv_table <- function(input_dir, table) {
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
    log_section(sprintf("  %-46s  %d rows", basename(path), nrow(df)))
    df
  })
  dplyr::bind_rows(frames)
}

load_erv_like_loci    <- function(input_dir) .read_erv_table(input_dir, "erv_like_loci")
load_erv_like_members <- function(input_dir) .read_erv_table(input_dir, "erv_like_members")


# Sum the `erv_like_dropped_noncanonical` metric across every {genome}.counts
# parquet. Returns 0L when the metric or the files are absent. Lets the
# canonical-vs-rearranged plot show dropped candidates even though they never
# reach the erv_like_loci table (they were filtered by require_canonical_order).
load_dropped_noncanonical <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.counts\\.parquet$",
                      full.names = TRUE)
  if (length(files) == 0L) return(0L)
  total <- 0L
  for (path in files) {
    df <- arrow::read_parquet(path)
    hit <- df$value[df$metric == "erv_like_dropped_noncanonical"]
    if (length(hit) > 0L) total <- total + sum(as.integer(hit), na.rm = TRUE)
  }
  as.integer(total)
}
