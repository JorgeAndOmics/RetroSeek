# -----------------------------------------------------------------------------
# range_analysis / exporters.R
# -----------------------------------------------------------------------------
# Output writers: GFF3 tracks, BED6, parquet plot dataframe, YAML manifest.
# All functions produce empty-but-valid files when given empty GRanges so
# Snakemake's output declarations are always satisfied.

suppressMessages({
  library(rtracklayer)
  library(arrow)
  library(yaml)
  library(GenomicRanges)
  library(S4Vectors)
})


# Write a GRanges as GFF3. For empty GRanges, emits a header-only file with
# a `##gff-version 3` and a `##source-version` pragma so the output exists
# and is traceable to a generator version.
track_exporter <- function(track, path, generator_version = "RetroSeek/unknown") {
  if (length(track) > 0L) {
    rtracklayer::export(track, path, format = "gff3")
    # Prepend a source-version pragma so downstream consumers can identify
    # which RetroSeek build emitted the file.
    lines <- readLines(path)
    if (length(lines) >= 1L) {
      lines <- c(lines[1], paste0("##source-version ", generator_version), lines[-1])
      writeLines(lines, path)
    }
  } else {
    writeLines(
      c("##gff-version 3", paste0("##source-version ", generator_version)),
      path
    )
  }
  invisible(path)
}


# Write a GRanges as 6-column BED. Score column uses max_bitscore where
# available, falling back to 0; name column uses ID where available.
bed_exporter <- function(track, path) {
  if (length(track) == 0L) {
    file.create(path)
    return(invisible(path))
  }
  m <- S4Vectors::mcols(track)
  name_col <- if ("ID" %in% names(m)) as.character(m$ID) else paste0("range_", seq_along(track))
  score_col <- if ("max_bitscore" %in% names(m)) {
    as.numeric(m$max_bitscore)
  } else if ("bitscore" %in% names(m)) {
    as.numeric(m$bitscore)
  } else {
    rep(0, length(track))
  }
  bed <- data.frame(
    chrom      = as.character(GenomicRanges::seqnames(track)),
    chromStart = BiocGenerics::start(track) - 1L,  # BED is 0-based half-open
    chromEnd   = BiocGenerics::end(track),
    name       = name_col,
    score      = round(score_col, 3),
    strand     = as.character(BiocGenerics::strand(track)),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    bed, path,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  invisible(path)
}


# Compute a simple per-probe overlap count matrix and write to CSV.
# Rows: probe; Columns: candidate / valid / total. Lightweight summary that
# replaces the pre-refactor full overlap-matrix machinery (which is not
# consumed by any current downstream rule).
overlap_matrix_exporter <- function(gr_virus, gr_candidates, gr_valid, path) {
  probes <- sort(unique(c(
    as.character(S4Vectors::mcols(gr_virus)$probe),
    as.character(S4Vectors::mcols(gr_candidates)$probe),
    as.character(S4Vectors::mcols(gr_valid)$probe)
  )))
  counts <- function(g, p) sum(as.character(S4Vectors::mcols(g)$probe) == p)
  out <- data.frame(
    probe     = probes,
    total     = vapply(probes, function(p) counts(gr_virus,      p), integer(1)),
    candidate = vapply(probes, function(p) counts(gr_candidates, p), integer(1)),
    valid     = vapply(probes, function(p) counts(gr_valid,      p), integer(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(out, path, row.names = FALSE)
  invisible(path)
}


# Resolve the generator version from `git rev-parse --short HEAD`, falling
# back to "unknown" when not in a git working tree (e.g. installed packages).
resolve_generator_version <- function() {
  ver <- tryCatch(
    suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  if (length(ver) == 1L && nzchar(ver)) paste0("RetroSeek/", ver)
  else "RetroSeek/unknown"
}


# Compute md5 of a file, returning NA_character_ if the file is missing.
file_md5 <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_character_)
  tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
}


# Emit a YAML manifest describing inputs, key counts, and outputs.
emit_manifest <- function(args, counts, generator_version, opts, path) {
  manifest <- list(
    generator       = generator_version,
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    inputs = list(
      fasta            = list(path = args$fasta,           md5 = file_md5(args$fasta)),
      blast_parquet    = list(path = args$blast,           md5 = file_md5(args$blast)),
      ltrdigest_gff3   = list(path = args$ltrdigest,       md5 = file_md5(args$ltrdigest)),
      probes_csv       = list(path = args$probes,          md5 = file_md5(args$probes)),
      config_yaml      = list(path = args$config,          md5 = file_md5(args$config))
    ),
    counts = counts,
    options = list(
      bitscore_threshold  = opts$bitscore_threshold,
      identity_threshold  = opts$identity_threshold,
      merge_option        = opts$merge_option,
      probe_min_length    = as.list(opts$probe_min_length),
      aggregation         = list(
        virus = opts$agg_virus, label = opts$agg_label,
        probe = opts$agg_probe, species = opts$agg_species,
        best_tiebreaker = opts$agg_best_tiebreaker
      )
    ),
    outputs = list(
      original_ranges             = args$original_ranges,
      original_ranges_reduced     = args$original_ranges_reduced,
      candidate_ranges            = args$candidate_ranges,
      candidate_ranges_reduced    = args$candidate_ranges_reduced,
      valid_ranges                = args$valid_ranges,
      valid_ranges_reduced        = args$valid_ranges_reduced,
      flanking_ltr_ranges         = args$flanking_ltr_ranges,
      overlap_matrix              = args$overlap_matrix,
      plot_dataframe              = args$plot_dataframe
    )
  )
  yaml::write_yaml(manifest, path)
  invisible(path)
}
