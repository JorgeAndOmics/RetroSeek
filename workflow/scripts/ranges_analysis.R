# =============================================================================
# ranges_analysis.R — orchestrator
# =============================================================================
# Phase-separated pipeline that integrates tBLASTn results with LTRdigest
# annotations. The heavy lifting is in workflow/scripts/range_analysis/*.R;
# this file only argument-parses, sources the modules, and wires the data
# flow phase by phase. Each section is intentionally short so the high-level
# data shape is visible at a glance.
#
# Phases:
#   1. Load inputs (config, BLAST parquet, LTRdigest GFF3, probes, FASTA).
#   2. Build the BLAST GRanges and apply quality / length filters.
#   3. Reduce overlapping ranges per (probe x virus|label) -> composite metrics.
#   4. Globally reduce per probe across virus|label groupings.
#   5. Process LTRdigest into retrotransposons, domains-with-probes, and
#      flanking LTRs.
#   6. Identify candidate hits (LTR-overlapping) and valid hits (probe matches
#      a Pfam-domain probe in the same retrotransposon).
#   7. Build the per-row plot dataframe + attach probe_category to all tracks.
#   8. Export GFF3 / BED6 / parquet / overlap matrix / manifest.

suppressMessages({
  library(argparse)
  library(GenomicRanges)
  library(plyranges)
  library(dplyr)
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source modules
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "range_aggregation_strategies.R"))
source(file.path(.script_dir, "range_analysis", "io.R"))
source(file.path(.script_dir, "range_analysis", "granges_build.R"))
source(file.path(.script_dir, "range_analysis", "filtering.R"))
source(file.path(.script_dir, "range_analysis", "reductions.R"))
source(file.path(.script_dir, "range_analysis", "validation.R"))
source(file.path(.script_dir, "range_analysis", "plot_dataframe.R"))
source(file.path(.script_dir, "range_analysis", "stage_dataframe.R"))
source(file.path(.script_dir, "range_analysis", "exporters.R"))


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
parser <- ArgumentParser(description = "Process tBLASTn and LTRdigest integration overlaps")
parser$add_argument("--fasta",                    required = TRUE)
parser$add_argument("--blast",                    required = TRUE)
parser$add_argument("--ltrdigest",                required = TRUE)
parser$add_argument("--probes",                   required = TRUE)
parser$add_argument("--probe_dict",               required = TRUE,
                    help = paste("Post-fetch probe_dict parquet (with",
                                 "genbank_seq populated). Used to compute",
                                 "per-hit query_coverage."))
parser$add_argument("--config",                   required = TRUE)
parser$add_argument("--original_ranges",          required = TRUE)
parser$add_argument("--original_ranges_reduced",  required = TRUE)
parser$add_argument("--candidate_ranges",         required = TRUE)
parser$add_argument("--candidate_ranges_reduced", required = TRUE)
parser$add_argument("--valid_ranges",             required = TRUE)
parser$add_argument("--valid_ranges_reduced",     required = TRUE)
parser$add_argument("--flanking_ltr_ranges",      required = TRUE)
parser$add_argument("--overlap_matrix",           required = TRUE)
parser$add_argument("--plot_dataframe",           required = TRUE)
# Middle-stage dataframes (homology + LTR integration) consumed by
# stage_plot_generator.R. Optional so non-Snakemake callers / tests can omit
# them — emitted only when the path is supplied (mirrors --manifest).
parser$add_argument("--stage_hits_dataframe",     required = FALSE)
parser$add_argument("--stage_ltr_dataframe",      required = FALSE)
parser$add_argument("--stage_reduced_dataframe",  required = FALSE)
parser$add_argument("--manifest",                 required = FALSE)
args <- parser$parse_args()


# ----------------------------------------------------------------------------
# Pipeline instrumentation
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}
.counts <- list()
record_count <- function(key, value) { .counts[[key]] <<- value }

message(paste0("Processing ranges for ",
               tools::file_path_sans_ext(basename(args$blast)), "..."))


# ----------------------------------------------------------------------------
# Phase 1. Load inputs
# ----------------------------------------------------------------------------
log_section("Phase 1: loading config, FASTA, BLAST parquet, LTRdigest, probes")
config        <- read_config(args$config)
opts          <- read_pipeline_options(config)
chrom_lengths <- load_chrom_lengths(args$fasta)
blast_df      <- load_blast_parquet(args$blast)
ltr_data      <- load_ltrdigest_gff3(args$ltrdigest)
probes        <- load_probes(args$probes)
probe_lengths <- load_probe_lengths(args$probe_dict)
record_count("raw_blast_hits",       nrow(blast_df))
record_count("ltrdigest_features",   length(ltr_data))


# ----------------------------------------------------------------------------
# Phase 2. Build the BLAST GRanges + filter
# ----------------------------------------------------------------------------
log_section("Phase 2: building BLAST GRanges and applying filters")
gr <- build_blast_gr(blast_df, probe_lengths = probe_lengths)
.pre_n <- length(gr)
gr <- filter_blast_gr(gr, opts$probe_min_length, opts$bitscore_threshold, opts$identity_threshold)
gr <- attach_min_gapwidth(gr, opts$probe_min_length)
record_count("filtered_blast_hits",  length(gr))
message(sprintf("   %d of %d hits kept", length(gr), .pre_n))


# ----------------------------------------------------------------------------
# Phase 3. First reduction (per probe x virus | label)
# ----------------------------------------------------------------------------
log_section("Phase 3: first reduction (per probe x virus|label)")
gr_virus <- reduce_first(gr, opts$merge_option, opts)
gr_virus <- attach_probe_id(gr_virus)
record_count("first_reduced_ranges", length(gr_virus))


# ----------------------------------------------------------------------------
# Phase 4. Global reduction (per probe across groupings)
# ----------------------------------------------------------------------------
log_section("Phase 4: global reduction (per probe across groupings)")
gr_global <- reduce_global(gr_virus, opts)
gr_global <- attach_probe_id(gr_global)
record_count("global_reduced_ranges", length(gr_global))


# ----------------------------------------------------------------------------
# Phase 5. LTRdigest processing (retros + domains + flanking LTRs)
# ----------------------------------------------------------------------------
log_section("Phase 5: extracting retrotransposons, domains, flanking LTRs")
domain_map        <- build_domain_map(opts$domains)
retrotransposons  <- extract_retrotransposons(ltr_data, resize_bp = opts$ltr_resize)
domains_w_probes  <- extract_domains_with_probes(ltr_data, domain_map)
flanking_ltrs     <- extract_flanking_ltrs(ltr_data)
record_count("retrotransposons",      length(retrotransposons))
record_count("domains_with_probes",   length(domains_w_probes))
record_count("flanking_ltrs",         length(flanking_ltrs))


# ----------------------------------------------------------------------------
# Phase 6. Candidate + valid hits
# ----------------------------------------------------------------------------
log_section("Phase 6: identifying candidate + valid hits")
candidate_hits           <- find_candidate_hits(gr_virus,  retrotransposons)
candidate_hits_reduced   <- find_candidate_hits(gr_global, retrotransposons)
valid_hits               <- find_valid_hits(candidate_hits,         retrotransposons,
                                            domains_w_probes, opts$agg_concat_separator)
valid_hits_reduced       <- find_valid_hits(candidate_hits_reduced, retrotransposons,
                                            domains_w_probes, opts$agg_concat_separator)
record_count("candidate_ranges",           length(candidate_hits))
record_count("candidate_ranges_reduced",   length(candidate_hits_reduced))
record_count("valid_ranges",               length(valid_hits))
record_count("valid_ranges_reduced",       length(valid_hits_reduced))


# ----------------------------------------------------------------------------
# Phase 7. Plot dataframe + probe-category tagging
# ----------------------------------------------------------------------------
log_section("Phase 7: plot dataframe + probe_category tagging")
# Plot dataframe is built from the valid-reduced tier (the final, LTR-integrated
# output) rather than gr_virus (the homology-only original tier). The middle
# stage is visualised separately by stage_plot_generator.R.
plot_df <- build_plot_dataframe(
  valid_hits_reduced, probes$df_sum, opts$main_probes,
  opts$agg_virus, opts$agg_concat_separator
)
gr_virus               <- attach_probe_category(gr_virus,             opts$main_probes, opts$agg_concat_separator)
gr_global              <- attach_probe_category(gr_global,            opts$main_probes, opts$agg_concat_separator)
candidate_hits         <- attach_probe_category(candidate_hits,       opts$main_probes, opts$agg_concat_separator)
candidate_hits_reduced <- attach_probe_category(candidate_hits_reduced, opts$main_probes, opts$agg_concat_separator)
valid_hits             <- attach_probe_category(valid_hits,           opts$main_probes, opts$agg_concat_separator)
valid_hits_reduced     <- attach_probe_category(valid_hits_reduced,   opts$main_probes, opts$agg_concat_separator)


# ----------------------------------------------------------------------------
# Phase 8. Export tracks + tables + manifest
# ----------------------------------------------------------------------------
log_section("Phase 8: exporting tracks, BED6, parquet, manifest")
gen_ver <- resolve_generator_version()

track_exporter(gr_virus,               args$original_ranges,          gen_ver)
track_exporter(gr_global,              args$original_ranges_reduced,  gen_ver)
bed_exporter(  gr_global,              sub("\\.gff3$", ".bed", args$original_ranges_reduced))

track_exporter(candidate_hits,         args$candidate_ranges,         gen_ver)
track_exporter(candidate_hits_reduced, args$candidate_ranges_reduced, gen_ver)
bed_exporter(  candidate_hits_reduced, sub("\\.gff3$", ".bed", args$candidate_ranges_reduced))

track_exporter(valid_hits,             args$valid_ranges,             gen_ver)
track_exporter(valid_hits_reduced,     args$valid_ranges_reduced,     gen_ver)
bed_exporter(  valid_hits_reduced,     sub("\\.gff3$", ".bed", args$valid_ranges_reduced))

track_exporter(flanking_ltrs,          args$flanking_ltr_ranges,      gen_ver)

overlap_matrix_exporter(gr_virus, candidate_hits, valid_hits, args$overlap_matrix)
arrow::write_parquet(plot_df, args$plot_dataframe)

# Middle-stage dataframes — emitted only when paths are supplied.
if (!is.null(args$stage_hits_dataframe)) {
  arrow::write_parquet(
    build_stage_hits_df(gr_virus, retrotransposons, candidate_hits, valid_hits),
    args$stage_hits_dataframe
  )
}
if (!is.null(args$stage_ltr_dataframe)) {
  arrow::write_parquet(
    build_stage_ltr_df(retrotransposons, flanking_ltrs, domains_w_probes,
                       ltr_data, gr_virus),
    args$stage_ltr_dataframe
  )
}
if (!is.null(args$stage_reduced_dataframe)) {
  arrow::write_parquet(
    build_stage_reduced_df(gr_global),
    args$stage_reduced_dataframe
  )
}

if (!is.null(args$manifest)) {
  emit_manifest(args, .counts, gen_ver, opts, args$manifest)
}

log_section("Done")
