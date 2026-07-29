# =============================================================================
# RetroSeek hotspot detector - Negative-Binomial GLM.
#
# Per-genome detection of windows enriched for ERV integrations beyond chance.
# A single deterministic NB GLM models per-window counts with a mask-aware
# offset and an optional chromosome covariate; per-window upper-tail p-values
# are BH-adjusted, thresholded, and merged into hotspot regions. There is no
# RNG in the core, so results are reproducible by construction; the global
# `parameters.seed` is set and recorded in the manifest for provenance only.
#
# Orchestrator only - pure transforms live in `hotspot/*.R`.
#
# Outputs (per genome):
#   {csv_dir}/{species}.csv                         per-window summary
#   {parquet_dir}/{species}.parquet                 same content, parquet
#   {parquet_dir}/{species}.manifest.yaml           provenance manifest
#   {track_output_dir}/{species}.gff3              merged hotspot regions
#   {track_output_dir}/{species}.bed               same regions, BED6
#   {pdf_output_dir}/{species}_manhattan.pdf       per-label Manhattan
#   {pdf_output_dir}/{species}_karyotype.pdf       ideogram + hotspots
#   {pdf_output_dir}/{species}_qq.pdf              per-label Q-Q diagnostic
#   {pdf_output_dir}/{species}_summary.pdf         density + width panel
#   {pdf_output_dir}/{species}_composition.pdf     per-hotspot structural mix
#   {csv_dir}/{species}.hotspots.csv               called regions + composition
# =============================================================================
options(warn = 1)
suppressMessages({
  library(argparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(BiocGenerics)
  library(dplyr)
  library(tibble)
  library(readr)
  library(arrow)
  library(yaml)
  library(tools)
})


# -----------------------------------------------------------------------------
# Locate sibling scripts + source modules
# -----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "utils",            "chrom_names.R"))
source(file.path(.script_dir, "hotspot", "io.R"))
source(file.path(.script_dir, "hotspot", "masking.R"))
source(file.path(.script_dir, "hotspot", "windowing.R"))
source(file.path(.script_dir, "hotspot", "models.R"))
source(file.path(.script_dir, "hotspot", "postprocess.R"))
source(file.path(.script_dir, "hotspot", "plots.R"))
source(file.path(.script_dir, "ranges",           "exporters.R"))


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
parser <- ArgumentParser(
  description = "RetroSeek hotspot detector (Negative-Binomial GLM)"
)
parser$add_argument("--fasta",            required = TRUE, help = "Genome FASTA")
parser$add_argument("--hits",             required = TRUE,
                    help = paste("Input events: catalog.csv (hotspot.input=catalog)",
                                 "or a raw-hit GFF3 track (hotspot.input=original)"))
parser$add_argument("--config",           required = TRUE, help = "Project YAML config")
parser$add_argument("--parquet_dir",      required = TRUE,
                    help = "Output dir for pipeline-internal {species}.parquet / .manifest.yaml")
parser$add_argument("--csv_dir",          required = TRUE,
                    help = "Output dir for user-facing {species}.csv")
parser$add_argument("--track_output_dir", required = TRUE,
                    help = "Output dir for {species}.gff3 / .bed")
parser$add_argument("--pdf_output_dir",   required = TRUE,
                    help = "Output dir for hotspot PDFs")
args <- parser$parse_args()


# -----------------------------------------------------------------------------
# Pipeline instrumentation
# -----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x


# -----------------------------------------------------------------------------
# Phase 1. Load config, seed, FASTA, hits
# -----------------------------------------------------------------------------
log_section("Phase 1: loading config, FASTA, hits")
config       <- read_config(args$config)
opts         <- read_hotspot_options(config)
# Deterministic model; seed set + recorded for provenance only.
set.seed(opts$seed)

species      <- tools::file_path_sans_ext(basename(args$fasta))
species_name <- if (!is.null(config$species[[species]])) {
  as.character(config$species[[species]])
} else {
  species
}
message(sprintf("  species: %s (display: %s)", species, species_name))
message(sprintf("  seed: %d | input: %s | window: %d",
                opts$seed, opts$input, opts$window_size))

genome     <- load_genome_for_hotspot(args$fasta)
seqs       <- genome$seqs
seqlengths <- genome$seqlengths
# ADR-012: `catalog` counts one row per INTEGRATION EVENT (the non-overlapping
# per-locus assembly); `original` counts raw tBLASTn hits, where a multi-gene
# provirus contributes several features and window counts are therefore weighted
# by gene content. The former is the defensible unit for "integration hotspot";
# the latter is kept because its density is what gives the NB power on sparse
# assemblies.
hits <- if (identical(opts$input, "catalog")) {
  load_catalog_loci(args$hits, species, config$species, opts$source)
} else {
  load_hits_gff(args$hits)
}
message(sprintf("  loaded %d chromosomes, %d events (%s tier)",
                length(seqlengths), length(hits), opts$input))
# Fail loud on an assembly / accession-namespace mismatch (GenBank vs RefSeq):
# otherwise zero hits overlap the windows and the run silently degrades to
# insufficient_data.
overlap_frac <- assert_hits_on_genome(hits, seqlengths)
message(sprintf("  %.1f%% of hits map to genome contigs", 100 * overlap_frac))


# -----------------------------------------------------------------------------
# Phase 2. Build mask, windows, count matrix
# -----------------------------------------------------------------------------
log_section("Phase 2: building mask, windows, count matrix")
mask <- build_n_mask(seqs, opts$mask_size, opts$mask_mismatch)
message(sprintf("  N-mask: %d intervals, %d total bp masked",
                length(mask), sum(BiocGenerics::width(mask))))

windows <- tile_genome_for_hotspot(seqlengths, opts$window_size)
effective_bp <- effective_bp_per_window(windows, mask)
chrom_stratum <- pool_small_scaffolds(seqlengths, opts$window_size,
                                       opts$unplaced_min_factor)
message(sprintf("  tiled into %d windows; %d strata after pooling",
                length(windows), length(unique(chrom_stratum))))


# -----------------------------------------------------------------------------
# Phase 3. Per-label NB GLM
# -----------------------------------------------------------------------------
log_section("Phase 3: per-group NB GLM")

# Grouping axis (ADR-012): `segment` / `taxon_call` are the calibrated calls the
# classifier produces; the legacy `label` is probe provenance, which ADR-007/008
# superseded. A raw-hit input carries none of these, so an absent column pools
# rather than erroring - a missing covariate must not abort a detection run.
.group_col <- opts$group_by
if (!identical(.group_col, "none") &&
    !(.group_col %in% colnames(S4Vectors::mcols(hits)))) {
  message(sprintf(
    "  WARNING: group_by='%s' but the input carries no such column; pooling instead",
    .group_col))
  .group_col <- "none"
}
gff_groups <- if (!identical(.group_col, "none")) {
  split(hits, as.character(S4Vectors::mcols(hits)[[.group_col]]))
} else {
  list(Ungrouped = hits)
}
message(sprintf("  %d group(s): %s",
                length(gff_groups),
                paste(names(gff_groups), collapse = ", ")))

per_label_window_dfs <- list()
per_label_hotspots   <- list()
fit_diagnostics      <- list()

for (label in names(gff_groups)) {
  log_section(sprintf("  [%s] %d hits", label, length(gff_groups[[label]])))
  events <- gff_groups[[label]]
  counts <- count_hits_per_window(windows, events)
  win_df <- assemble_window_table(windows, counts, effective_bp,
                                   chrom_stratum, label)

  fit <- fit_nb_model(win_df, opts$window_size,
                      strata_by_chromosome = opts$strata_by_chromosome)
  fit_diagnostics[[label]] <- list(
    status = fit$status, family = fit$family,
    theta  = if (is.na(fit$theta)) NULL else as.numeric(fit$theta),
    n_fit_rows = nrow(fit$fit_data)
  )
  if (identical(fit$status, "insufficient_data")) {
    message(sprintf("    WARNING: %s has too few non-zero windows to fit; emitting NA p-values",
                    label))
  } else if (identical(fit$status, "failed")) {
    message(sprintf("    WARNING: %s NB GLM failed to converge; emitting NA p-values", label))
  }

  # Score windows
  scored <- score_windows_nb(win_df, fit)

  # Postprocess: select -> merge -> recompute -> min-hits
  significant <- select_significant_windows(scored, opts$pvalue_threshold)
  merged      <- merge_adjacent_hotspots(significant, gap = opts$merge_gap)
  merged      <- recompute_merged_pvalue(merged, fit)
  merged      <- apply_min_hits_filter(merged, opts$min_hits)

  message(sprintf("    %d significant windows -> %d hotspot regions after merge & filter",
                  nrow(significant), length(merged)))

  per_label_window_dfs[[label]] <- scored
  per_label_hotspots[[label]]   <- merged
}


# -----------------------------------------------------------------------------
# Phase 4. Concatenate, assign IDs, attach to per-window table
# -----------------------------------------------------------------------------
log_section("Phase 4: concatenating per-label outputs")
all_windows_df <- dplyr::bind_rows(per_label_window_dfs)
all_hotspots <- if (length(per_label_hotspots) == 0L) {
  .empty_merged_gr()
} else {
  do.call(c, unname(per_label_hotspots))
}
all_hotspots <- assign_hotspot_ids(all_hotspots, species)
all_windows_df <- attach_hotspot_id_to_windows(all_windows_df, all_hotspots)
# Describe each called region by the loci inside it (ADR-012): structural class,
# tier, dominant lineage, mean confidence. Annotation only - no region is added,
# removed or re-scored by this.
all_hotspots <- annotate_hotspot_composition(all_hotspots, hits, .group_col)
message(sprintf("  total windows: %d | total hotspots: %d",
                nrow(all_windows_df), length(all_hotspots)))
if (length(all_hotspots) > 0L) {
  .comp <- S4Vectors::mcols(all_hotspots)
  message(sprintf("  composition: %d loci (full %d / partial %d / gene %d)",
                  sum(.comp$n_loci), sum(.comp$n_full), sum(.comp$n_partial),
                  sum(.comp$n_gene)))
}


# -----------------------------------------------------------------------------
# Phase 5. Emit CSV / Parquet / GFF3 / BED / manifest
# -----------------------------------------------------------------------------
log_section("Phase 5: emitting tables, tracks, manifest")

dir.create(args$parquet_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(args$csv_dir,          showWarnings = FALSE, recursive = TRUE)
dir.create(args$track_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(args$pdf_output_dir,   showWarnings = FALSE, recursive = TRUE)

csv_path      <- file.path(args$csv_dir,          paste0(species, ".csv"))
regions_path  <- file.path(args$csv_dir,          paste0(species, ".hotspots.csv"))
parquet_path  <- file.path(args$parquet_dir,      paste0(species, ".parquet"))
manifest_path <- file.path(args$parquet_dir,      paste0(species, ".manifest.yaml"))
gff_path      <- file.path(args$track_output_dir, paste0(species, ".gff3"))
bed_path      <- file.path(args$track_output_dir, paste0(species, ".bed"))

readr::write_csv(all_windows_df, csv_path)
arrow::write_parquet(all_windows_df, parquet_path)
# Per-REGION table: the called hotspots with their composition (ADR-012). The
# per-window CSV above answers "where is enrichment"; this one answers "what is
# the enrichment made of", without parsing GFF3 attributes.
readr::write_csv(as.data.frame(all_hotspots, row.names = NULL), regions_path)

generator_version <- resolve_generator_version()
track_exporter(all_hotspots, gff_path, generator_version = generator_version)

# BED uses the standard `bed_exporter()` from ranges/exporters.R, which
# expects `mcols$ID` (name) and `mcols$max_bitscore` (score). Munge a copy of
# the GRanges with those names so we DRY the writer rather than duplicating.
hotspots_for_bed <- all_hotspots
if (length(hotspots_for_bed) > 0L) {
  scores_raw <- as.numeric(S4Vectors::mcols(hotspots_for_bed)$pval_nb_region)
  scores_raw[is.na(scores_raw)] <- 1
  scores <- pmin(round(-log10(pmax(scores_raw, .Machine$double.xmin)) * 100),
                 1000)
  S4Vectors::mcols(hotspots_for_bed)$ID <-
    as.character(S4Vectors::mcols(hotspots_for_bed)$hotspot_id)
  S4Vectors::mcols(hotspots_for_bed)$max_bitscore <- as.integer(scores)
}
bed_exporter(hotspots_for_bed, bed_path)

# Hotspot-specific manifest emitter (the ranges one is shaped for that
# pipeline's args). Keeps provenance for reruns.
emit_hotspot_manifest <- function(args, opts, species, species_name,
                                  fit_diagnostics, counts,
                                  generator_version, path) {
  manifest <- list(
    generator       = generator_version,
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    species         = species,
    species_display = species_name,
    inputs = list(
      fasta  = list(path = args$fasta,  md5 = file_md5(args$fasta)),
      gff    = list(path = args$gff,    md5 = file_md5(args$gff)),
      config = list(path = args$config, md5 = file_md5(args$config))
    ),
    options = opts,
    counts = counts,
    per_label_diagnostics = fit_diagnostics,
    outputs = list(
      csv      = csv_path,
      parquet  = parquet_path,
      gff3     = gff_path,
      bed      = bed_path
    )
  )
  yaml::write_yaml(manifest, path)
  invisible(path)
}

emit_hotspot_manifest(
  args, opts, species, species_name, fit_diagnostics,
  counts = list(
    total_hits      = length(hits),
    total_windows   = length(windows),
    total_hotspots  = length(all_hotspots),
    n_label_groups  = length(gff_groups)
  ),
  generator_version = generator_version,
  path = manifest_path
)


# -----------------------------------------------------------------------------
# Phase 6. Plots
# -----------------------------------------------------------------------------
log_section("Phase 6: plots")
plot_w <- as.integer(config$plots$width  %||% 15L)
plot_h <- as.integer(config$plots$height %||% 12L)

manhattan_pdf <- file.path(args$pdf_output_dir, paste0(species, "_manhattan.pdf"))
karyotype_pdf <- file.path(args$pdf_output_dir, paste0(species, "_karyotype.pdf"))
qq_pdf        <- file.path(args$pdf_output_dir, paste0(species, "_qq.pdf"))
summary_pdf   <- file.path(args$pdf_output_dir, paste0(species, "_summary.pdf"))
composition_pdf <- file.path(args$pdf_output_dir, paste0(species, "_composition.pdf"))

# Manhattan: one page per label
labels_present <- unique(all_windows_df$label)
manhattan_pages <- lapply(labels_present, function(lbl) {
  plot_manhattan(
    dplyr::filter(all_windows_df, .data$label == lbl),
    threshold = opts$pvalue_threshold,
    title     = sprintf("Hotspot Manhattan plot - %s", species_name),
    subtitle  = sprintf("Label: %s", lbl)
  )
})
save_plots_pdf_pages(manhattan_pages, manhattan_pdf, plot_w, plot_h)

# Karyotype: one page covering all hotspots
save_plots_pdf_pages(
  list(plot_karyotype(
    seqlengths, all_hotspots,
    title    = sprintf("Hotspot karyotype - %s", species_name),
    subtitle = sprintf("%d hotspots across %d chromosomes",
                       length(all_hotspots), length(seqlengths))
  )),
  karyotype_pdf, plot_w, plot_h
)

# Q-Q: one page per label
qq_pages <- lapply(labels_present, function(lbl) {
  plot_qq(
    dplyr::filter(all_windows_df, .data$label == lbl),
    title    = sprintf("Hotspot Q-Q diagnostic - %s", species_name),
    subtitle = sprintf("Label: %s", lbl)
  )
})
save_plots_pdf_pages(qq_pages, qq_pdf, plot_w, plot_h)

# Summary panel
save_plots_pdf_pages(
  list(plot_summary_panel(
    all_hotspots, seqlengths,
    title    = sprintf("Hotspot summary - %s", species_name),
    subtitle = sprintf("%d hotspots; threshold q < %.3g",
                       length(all_hotspots), opts$pvalue_threshold)
  )),
  summary_pdf, plot_w, plot_h
)

# Composition panel (ADR-012): what each called hotspot is actually made of.
save_plots_pdf_pages(
  list(plot_hotspot_composition(
    all_hotspots,
    title    = sprintf("Hotspot composition - %s", species_name),
    subtitle = sprintf("%s tier, grouped by %s", opts$input, .group_col)
  )),
  composition_pdf, plot_w, plot_h
)

log_section("Done")
