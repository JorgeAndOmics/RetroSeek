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
#   {pdf_output_dir}/{species}.hotspots.pdf        key page, then per-label
#                                                   Manhattan and Q-Q pages, the
#                                                   karyotype, the summary and
#                                                   the per-hotspot composition
#   {csv_dir}/{species}.hotspots.csv               called regions + composition
# =============================================================================
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
source(file.path(.script_dir, "utils", "log.R"))          # line contract, run_main (ADR-021)
source(file.path(.script_dir, "plot2sort", "style.R"))    # palette, theme, stage PDFs
source(file.path(.script_dir, "plot2sort", "helpers.R"))  # empty_plot, add_titles
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
parser$add_argument("--log", default = NULL,
                    help = "job log file; the Snakemake log: path")
args <- parser$parse_args()
log_job(args$log, "hotspot_detector")

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x


# -----------------------------------------------------------------------------
# Phase 1. Load config, seed, FASTA, hits
# -----------------------------------------------------------------------------
.load_inputs <- function(args) {
  log_section("Phase 1: loading config, FASTA, hits")
  config <- read_config(args$config)
  opts   <- read_hotspot_options(config)
  # Deterministic model; seed set + recorded for provenance only.
  set.seed(opts$seed)

  species      <- tools::file_path_sans_ext(basename(args$fasta))
  species_name <- .display_name(species, config$species)
  log_info("species: %s (display: %s)", species, species_name)
  log_info("seed: %d, input: %s, window: %d", opts$seed, opts$input, opts$window_size)

  genome <- load_genome_for_hotspot(args$fasta)
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
  log_info("loaded %d chromosomes, %d events (%s tier)",
           length(genome$seqlengths), length(hits), opts$input)
  # Fail loud on an assembly / accession-namespace mismatch (GenBank vs RefSeq):
  # otherwise zero hits overlap the windows and the run silently degrades to
  # insufficient_data.
  overlap_frac <- assert_hits_on_genome(hits, genome$seqlengths)
  log_info("%.1f%% of hits map to genome contigs", 100 * overlap_frac)
  list(config = config, opts = opts, species = species, species_name = species_name,
       seqs = genome$seqs, seqlengths = genome$seqlengths, hits = hits)
}

# -----------------------------------------------------------------------------
# Phase 2. Build mask, windows, count matrix
# -----------------------------------------------------------------------------
.build_windows <- function(inputs) {
  log_section("Phase 2: building mask, windows, count matrix")
  opts <- inputs$opts
  mask <- build_n_mask(inputs$seqs, opts$mask_size, opts$mask_mismatch)
  log_info("N-mask: %d intervals, %d total bp masked",
           length(mask), sum(BiocGenerics::width(mask)))

  windows <- tile_genome_for_hotspot(inputs$seqlengths, opts$window_size)
  effective_bp <- effective_bp_per_window(windows, mask)
  chrom_stratum <- pool_small_scaffolds(inputs$seqlengths, opts$window_size,
                                        opts$unplaced_min_factor)
  log_info("tiled into %d windows; %d strata after pooling",
           length(windows), length(unique(chrom_stratum)))
  list(windows = windows, effective_bp = effective_bp, chrom_stratum = chrom_stratum)
}

# -----------------------------------------------------------------------------
# Phase 3. Per-group NB GLM
# -----------------------------------------------------------------------------
# The column hits are grouped by, or "none" when the input lacks it. Grouping
# axis (ADR-012): `segment` / `taxon_call` are the calibrated calls the
# classifier produces; the legacy `label` is probe provenance, which ADR-007/008
# superseded. A raw-hit input carries none of these, so an absent column pools
# rather than erroring - a missing covariate must not abort a detection run.
.group_column <- function(hits, group_by) {
  if (identical(group_by, "none") ||
      group_by %in% colnames(S4Vectors::mcols(hits))) {
    return(group_by)
  }
  log_warn("hotspot.group_by is '%s' but the input has no such column; pooling all loci instead",
           group_by)
  "none"
}

.split_hits <- function(hits, group_col) {
  if (identical(group_col, "none")) return(list(Ungrouped = hits))
  split(hits, as.character(S4Vectors::mcols(hits)[[group_col]]))
}

# Log what a fit's status means for its windows. An insufficient fit is expected
# for rare lineages (a handful of loci), and recorded in the tables as NA
# p-values: information, not something to act on.
.log_fit_status <- function(label, status) {
  if (identical(status, "insufficient_data")) {
    log_info("%s: too few non-zero windows to fit the model; its windows get NA p-values",
             label)
  } else if (identical(status, "failed")) {
    log_warn(paste("%s: the NB GLM did not converge, so its windows get NA p-values;",
                   "a larger hotspot.window_size usually lets it fit"), label)
  }
}

# One group's scored windows, merged hotspots and fit diagnostics.
.scan_group <- function(label, events, layout, opts) {
  log_info("%s: %d hits", label, length(events))
  counts <- count_hits_per_window(layout$windows, events)
  win_df <- assemble_window_table(layout$windows, counts, layout$effective_bp,
                                  layout$chrom_stratum, label)
  fit <- fit_nb_model(win_df, opts$window_size,
                      strata_by_chromosome = opts$strata_by_chromosome)
  .log_fit_status(label, fit$status)

  scored <- score_windows_nb(win_df, fit)
  # Postprocess: select -> merge -> recompute -> min-hits
  significant <- select_significant_windows(scored, opts$pvalue_threshold)
  merged      <- merge_adjacent_hotspots(significant, gap = opts$merge_gap)
  merged      <- recompute_merged_pvalue(merged, fit)
  merged      <- apply_min_hits_filter(merged, opts$min_hits)
  log_info("%s: %d significant windows, %d hotspot regions after merge and filter",
           label, nrow(significant), length(merged))
  list(
    windows = scored,
    hotspots = merged,
    diagnostics = list(
      status = fit$status, family = fit$family,
      theta  = if (is.na(fit$theta)) NULL else as.numeric(fit$theta),
      n_fit_rows = nrow(fit$fit_data)
    )
  )
}

# -----------------------------------------------------------------------------
# Phase 4. Concatenate, assign IDs, attach to per-window table
# -----------------------------------------------------------------------------
# Lineages can have hotspots on different chromosomes; joining GRanges whose
# sequence levels differ warns ("no sequence levels in common"). Give every
# piece the full set first; the joined level order is the one c() made anyway.
.join_hotspots <- function(per_label_hotspots) {
  if (length(per_label_hotspots) == 0L) return(.empty_merged_gr())
  every_level <- unique(unlist(lapply(per_label_hotspots, GenomeInfoDb::seqlevels)))
  shared <- lapply(per_label_hotspots, function(gr) {
    GenomeInfoDb::seqlevels(gr) <- union(GenomeInfoDb::seqlevels(gr), every_level)
    gr
  })
  do.call(c, unname(shared))
}

.combine_groups <- function(scans, inputs, group_col) {
  log_section("Phase 4: concatenating per-label outputs")
  all_windows_df <- dplyr::bind_rows(lapply(scans, `[[`, "windows"))
  all_hotspots <- .join_hotspots(lapply(scans, `[[`, "hotspots"))
  all_hotspots <- assign_hotspot_ids(all_hotspots, inputs$species)
  all_windows_df <- attach_hotspot_id_to_windows(all_windows_df, all_hotspots)
  # Describe each called region by the loci inside it (ADR-012): structural class,
  # tier, dominant lineage, mean confidence. Annotation only - no region is added,
  # removed or re-scored by this.
  all_hotspots <- annotate_hotspot_composition(all_hotspots, inputs$hits, group_col)
  log_info("total windows: %d, total hotspots: %d",
           nrow(all_windows_df), length(all_hotspots))
  if (length(all_hotspots) > 0L) {
    .comp <- S4Vectors::mcols(all_hotspots)
    log_info("composition: %d loci (full %d, partial %d, gene %d)",
             sum(.comp$n_loci), sum(.comp$n_full), sum(.comp$n_partial), sum(.comp$n_gene))
  }
  list(windows = all_windows_df, hotspots = all_hotspots)
}

# -----------------------------------------------------------------------------
# Phase 5. Emit CSV / Parquet / GFF3 / BED / manifest
# -----------------------------------------------------------------------------
# BED uses the standard `bed_exporter()` from ranges/exporters.R, which expects
# `mcols$ID` (name) and `mcols$max_bitscore` (score). Munge a copy of the
# GRanges with those names so we DRY the writer rather than duplicating.
.hotspots_for_bed <- function(hotspots) {
  if (length(hotspots) == 0L) return(hotspots)
  scores_raw <- as.numeric(S4Vectors::mcols(hotspots)$pval_nb_region)
  scores_raw[is.na(scores_raw)] <- 1
  scores <- pmin(round(-log10(pmax(scores_raw, .Machine$double.xmin)) * 100), 1000)
  S4Vectors::mcols(hotspots)$ID <- as.character(S4Vectors::mcols(hotspots)$hotspot_id)
  S4Vectors::mcols(hotspots)$max_bitscore <- as.integer(scores)
  hotspots
}

.output_paths <- function(args, species) {
  list(
    csv      = file.path(args$csv_dir,          paste0(species, ".csv")),
    regions  = file.path(args$csv_dir,          paste0(species, ".hotspots.csv")),
    parquet  = file.path(args$parquet_dir,      paste0(species, ".parquet")),
    manifest = file.path(args$parquet_dir,      paste0(species, ".manifest.yaml")),
    gff      = file.path(args$track_output_dir, paste0(species, ".gff3")),
    bed      = file.path(args$track_output_dir, paste0(species, ".bed"))
  )
}

.write_outputs <- function(args, inputs, layout, result, scans) {
  log_section("Phase 5: emitting tables, tracks, manifest")
  # The PDF directory too: .write_plots() writes into it next.
  for (dir in c(args$parquet_dir, args$csv_dir, args$track_output_dir,
                args$pdf_output_dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  paths <- .output_paths(args, inputs$species)
  readr::write_csv(result$windows, paths$csv)
  arrow::write_parquet(result$windows, paths$parquet)
  # Per-REGION table: the called hotspots with their composition (ADR-012). The
  # per-window CSV above answers "where is enrichment"; this one answers "what is
  # the enrichment made of", without parsing GFF3 attributes.
  readr::write_csv(as.data.frame(result$hotspots, row.names = NULL), paths$regions)

  generator_version <- resolve_generator_version()
  track_exporter(result$hotspots, paths$gff, generator_version = generator_version)
  bed_exporter(.hotspots_for_bed(result$hotspots), paths$bed)

  emit_hotspot_manifest(
    inputs  = list(fasta = args$fasta, hits = args$hits, config = args$config),
    outputs = list(csv = paths$csv, parquet = paths$parquet,
                   gff3 = paths$gff, bed = paths$bed),
    opts = inputs$opts, species = inputs$species, species_name = inputs$species_name,
    fit_diagnostics = lapply(scans, `[[`, "diagnostics"),
    counts = list(
      total_hits      = length(inputs$hits),
      total_windows   = length(layout$windows),
      total_hotspots  = length(result$hotspots),
      n_label_groups  = length(scans)
    ),
    generator_version = generator_version,
    path = paths$manifest
  )
}

# -----------------------------------------------------------------------------
# Phase 6. Plots
# -----------------------------------------------------------------------------
# Per label: a Manhattan page, then its Q-Q diagnostic. Only for labels the
# model could test; a label with too few loci has no callable window, and a
# page per such label would be a run of empty placeholders. They are named on
# the key page instead.
.label_pages <- function(windows_df, opts, plot_species, labels) {
  unlist(lapply(labels, function(lbl) {
    windows <- dplyr::filter(windows_df, .data$label == lbl)
    list(plot_manhattan(windows, opts$pvalue_threshold, plot_species, lbl),
         plot_qq(windows, plot_species, lbl))
  }), recursive = FALSE)
}

.write_plots <- function(args, inputs, result, group_col) {
  log_section("Phase 6: plots")
  use_retroseek_style()
  opts <- inputs$opts
  plot_species <- display_species(inputs$species, inputs$config$species)
  callable <- result$windows %>%
    dplyr::group_by(.data$label) %>%
    dplyr::summarise(tested = any(!is.na(.data$qval_nb)), .groups = "drop")
  untested <- sort(callable$label[!callable$tested])
  pages <- c(
    .label_pages(result$windows, opts, plot_species, callable$label[callable$tested]),
    list(plot_karyotype(inputs$seqlengths, result$hotspots, plot_species),
         plot_summary_panel(result$hotspots, inputs$seqlengths, plot_species),
         plot_hotspot_composition(result$hotspots, plot_species,
                                  sprintf("The %s tier, grouped by %s.", opts$input,
                                          group_col)))
  )
  key <- key_page(
    sprintf("Integration hotspots in %s", plot_species),
    paste(
      sprintf(paste("Windows of the genome holding more %s loci than a negative binomial",
                    "model expects, merged into hotspots (q below %s). %d hotspots were",
                    "called. The composition page shows what each is made of."),
              opts$input, format(opts$pvalue_threshold), length(result$hotspots)),
      if (length(untested)) {
        sprintf("Too few loci to test: %s.", paste(display_label(untested), collapse = ", "))
      }),
    colours = stats::setNames(unname(.STRUCTURE_COLOUR),
                              display_label(names(.STRUCTURE_COLOUR))),
    pages = page_titles(pages))
  save_stage_pdf(c(list(key), pages),
                 file.path(args$pdf_output_dir, paste0(inputs$species, ".hotspots.pdf")))
}

# -----------------------------------------------------------------------------
# main(): the six phases. Run through run_main() so warnings are logged and
# every ending is recorded the same way (ADR-021).
# -----------------------------------------------------------------------------
main <- function(args) {
  inputs <- .load_inputs(args)
  layout <- .build_windows(inputs)

  log_section("Phase 3: per-group NB GLM")
  group_col <- .group_column(inputs$hits, inputs$opts$group_by)
  groups <- .split_hits(inputs$hits, group_col)
  log_info("%d group(s): %s", length(groups), paste(names(groups), collapse = ", "))
  scans <- lapply(stats::setNames(nm = names(groups)), function(label) {
    .scan_group(label, groups[[label]], layout, inputs$opts)
  })

  result <- .combine_groups(scans, inputs, group_col)
  .write_outputs(args, inputs, layout, result, scans)
  .write_plots(args, inputs, result, group_col)
  log_ok("%s hotspots in %s windows, from %s events",
         format(length(result$hotspots), big.mark = ","),
         format(nrow(result$windows), big.mark = ","),
         format(length(inputs$hits), big.mark = ","))
}

run_main(function() main(args))
