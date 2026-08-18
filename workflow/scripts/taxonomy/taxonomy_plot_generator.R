# =============================================================================
# taxonomy_plot_generator.R
# =============================================================================
# Builds RetroSeek's taxonomy panel: per-locus ERV taxon-call composition and
# resolution plots, derived from the per-genome `<genome>.loci.parquet` tables
# that taxonomy_classify_loci.py writes to data/tables/taxonomy_classification/.
# Because every plot is computed from those same loci tables, the panel is
# concordant with the tables by construction.
#
# Plots:
#   1. taxon_composition      - per-species stacked counts of confident (axis-resolved) taxon calls.
#   2. rank_resolution        - per-species resolved-rank distribution (genus /
#                               subfamily / family / unclassified).
#   3. method_mix             - how taxon calls were made (placement / lca / presence).
#   4. erv_class_composition  - per-species Class I/II/III composition (the
#                               literature anchor: bats+human Class I, mouse Class II).
#   5. mosaic sub-panel       - recombination views over mosaic loci:
#                               mosaic_alluvial (gene -> taxon flows), mosaic_burden
#                               (per-species mosaic fraction), mosaic_taxon_pairs
#                               (recombination-partner heatmap), mosaic_gene_discordance
#                               (odd-one-out gene), mosaic_composition_by_species.
# Plus the confidence / evidence / domain-tier / structure-class panels (18 PNGs total).
#
# Shared infrastructure (empty_plot, add_titles, save_plot) is reused from
# plot2sort/*.R - not duplicated - so the panel matches the existing plots.
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard keeps
# the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
  library(ggsci)        # Scientific colour palettes
  library(ggalluvial)   # Alluvial flows for the mosaic plot
  library(GenomicRanges)  # catalog reconciliation (ltr-flanked-precedence overlap)
  library(patchwork)      # tree | bars composition for the tree-attached panels
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source shared modules (same idiom as the other
# plot generators, so the panel reuses one set of helpers).
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- tryCatch(sys.frame(i), error = function(e) NULL)
    if (is.null(fr)) next
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))  # empty_plot, add_titles
source(file.path(.script_dir, "..", "plot2sort", "io.R"))       # save_plot


# ----------------------------------------------------------------------------
# Pipeline instrumentation - same idiom as the other plot generators (save_plot
# calls log_section, so it must exist in the global env at call time).
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}


# ----------------------------------------------------------------------------
# Load every per-genome loci table in `input_dir`, tagging each row with its
# species (the filename stem). Returns one tidy data frame (empty if none).
# ----------------------------------------------------------------------------
load_loci <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.loci\\.parquet$", full.names = TRUE)
  if (length(files) == 0L) return(tibble())
  frames <- lapply(files, function(f) {
    df <- as_tibble(arrow::read_parquet(f))
    if (nrow(df) == 0L) return(NULL)
    df$species <- sub("\\.loci$", "", tools::file_path_sans_ext(basename(f)))
    df
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) return(tibble())
  bind_rows(frames)
}


# Load every per-genome orphan table (the recovered non-LTR tier). Same schema
# as the loci tables (source == "orphan"); the on-disk file keeps its historical
# `.orphans.parquet` name. Empty if none.
load_orphans <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.orphans\\.parquet$", full.names = TRUE)
  if (length(files) == 0L) return(tibble())
  frames <- lapply(files, function(f) {
    df <- as_tibble(arrow::read_parquet(f))
    if (nrow(df) == 0L) return(NULL)
    df$species <- sub("\\.orphans$", "", tools::file_path_sans_ext(basename(f)))
    df
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) return(tibble())
  bind_rows(frames)
}


# Reconcile the unified catalog to a fully non-overlapping record set with
# LTR-FLANKED PRECEDENCE: within each (species, seqname), drop any orphan locus
# whose span overlaps an LTR-flanked locus. LTR-flanked loci are LTR-confirmed;
# orphan loci are proximity-inferred, and a proximity cluster can bridge OVER an
# LTR-flanked provirus (its member hits flank the element). Where they collide the
# LTR-confirmed call wins. Orphan-orphan and LTR-flanked/LTR-flanked are already
# non-overlapping (clustering + one-per-element), so this only resolves the
# cross-tier edge. Dropped orphans remain in the per-genome .orphans table.
reconcile_catalog <- function(combined) {
  if (nrow(combined) == 0L) return(combined)
  src <- as.character(combined$source)
  anch_i <- which(src == "ltr-flanked")
  orph_i <- which(src != "ltr-flanked")
  if (length(anch_i) == 0L || length(orph_i) == 0L) return(combined)
  gr <- GenomicRanges::GRanges(
    seqnames = paste(combined$species, combined$seqname, sep = "|"),
    ranges   = IRanges::IRanges(suppressWarnings(as.integer(combined$start)),
                                suppressWarnings(as.integer(combined$end)))
  )
  ov <- GenomicRanges::findOverlaps(gr[orph_i], gr[anch_i], ignore.strand = TRUE)
  drop <- orph_i[unique(S4Vectors::queryHits(ov))]
  combined[setdiff(seq_len(nrow(combined)), drop), , drop = FALSE]
}


# Build the tidy classification report from the combined (LTR-flanked + orphan)
# loci frame: counts by genus, by confidence, by method, plus mosaic and
# integration totals - split by `source` so the two tiers stay distinguishable.
# Returns a long tibble (source, dimension, level, count); empty-safe.
build_report <- function(combined) {
  if (nrow(combined) == 0L) {
    return(tibble(
      source = character(), dimension = character(),
      level = character(), count = integer()
    ))
  }
  cols <- c("source", "dimension", "level", "count")
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  per_source <- function(df, src) {
    taxon <- df %>% filter(.data$resolved == "True") %>%
      count(level = .data$taxon_call, name = "count") %>%
      mutate(dimension = "taxon")
    conf <- df %>%
      count(level = .data$confidence_tag, name = "count") %>%
      mutate(dimension = "confidence")
    method <- df %>% filter(.data$resolved == "True") %>%
      count(level = .data$method, name = "count") %>%
      mutate(dimension = "method")
    summary <- tibble(
      dimension = "summary",
      level     = c("mosaic", "integrations"),
      count     = c(sum(df$is_mosaic == "True"), nrow(df))
    )
    bind_rows(taxon, conf, method, summary) %>%
      mutate(source = src) %>%
      select(all_of(cols))
  }
  combined %>%
    group_split(.data$source) %>%
    lapply(function(df) per_source(df, df$source[1])) %>%
    bind_rows()
}


# ----------------------------------------------------------------------------
# Shared categorical palettes (ADR-009). Levels are declared in full and used
# with drop = FALSE so a tier/class that is absent in one genome still keeps its
# colour and legend slot - and the per-hit `positional` mode's extra non_domain
# level never breaks a plot built on membership-mode data.
# ----------------------------------------------------------------------------
.SOURCE_FILL       <- c(`ltr-flanked` = "#1F78B4", orphan = "#33A02C")
.DOMAIN_TIER_LEVELS <- c("domain_selected", "domain_unlisted", "non_domain")
.DOMAIN_TIER_FILL   <- c(domain_selected = "#1B9E77",
                         domain_unlisted = "#D95F02",
                         non_domain      = "#999999")
.STRUCTURE_LEVELS  <- c("full", "partial", "gene")
.STRUCTURE_FILL    <- c(full = "#1B9E77", partial = "#D95F02", gene = "#7570B3")


# ----------------------------------------------------------------------------
# Plot builders. Each takes the combined loci frame and returns a ggplot (or
# empty_plot() when there is nothing to show), so the orchestrator stays flat.
# ----------------------------------------------------------------------------

# Confident (axis-resolved) calls, stacked per species and coloured by taxon.
# ----------------------------------------------------------------------------
# Panel registry - ONE declaration, TWO consumers
# ----------------------------------------------------------------------------
# main() below and taxonomy_segments.R both iterate this, so a plot added here
# is offered to both by construction. Previously each kept its own hand-written
# emit list, and they drifted: the per-segment subset stayed at 3 plots while
# this panel grew to 24.
#
# Fields:
#   file    output PNG name
#   build   function(d, ctx) -> ggplot. ctx carries `tree_dir` and
#           `confidence_min` for the handful of builders that take them, so
#           every entry has the same shape.
#   data    "loci"     = LTR-flanked only (the taxon-founded assembly)
#           "combined" = both tiers
#           Losing this distinction would quietly mix orphans into the
#           composition and mosaic panels, which deliberately exclude them.
#   n_x     "species" scales canvas + tick labels to the species count;
#           NULL keeps the fixed config canvas (histograms, heatmaps, alluvials)
#   axis    "y" routes through auto_dims instead - the tree panels, whose height
#           grows with tip count and whose y axis is supplied by the tree
#   segment  FALSE where the plot is DEGENERATE for a single segment:
#           erv_class is constant within a genus (measured: exactly 1.00
#           distinct values), and the taxonomy cladogram collapses to one tip.
panel_registry <- function() {
  list(
    list(file = "taxon_composition.png", build = function(d, ctx) taxon_composition_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    list(file = "rank_resolution.png", build = function(d, ctx) rank_resolution_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    list(file = "method_mix.png", build = function(d, ctx) method_mix_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    # erv_class is a function of genus, so within one segment it is constant.
    list(file = "erv_class_composition.png", build = function(d, ctx) erv_class_composition_plot(d), data = "loci", n_x = "species", axis = "x", segment = FALSE),
    list(file = "mosaic_alluvial.png", build = function(d, ctx) mosaic_alluvial_plot(d), data = "loci", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "mosaic_burden.png", build = function(d, ctx) mosaic_burden_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    list(file = "mosaic_taxon_pairs.png", build = function(d, ctx) mosaic_taxon_pairs_plot(d), data = "loci", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "mosaic_gene_discordance.png", build = function(d, ctx) mosaic_gene_discordance_plot(d), data = "loci", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "mosaic_composition_by_species.png", build = function(d, ctx) mosaic_composition_by_species_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    list(file = "confidence.png", build = function(d, ctx) confidence_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    list(file = "confidence_count.png", build = function(d, ctx) confidence_count_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    list(file = "confidence_gradient.png", build = function(d, ctx) confidence_gradient_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    list(file = "evidence_depth.png", build = function(d, ctx) evidence_depth_plot(d), data = "combined", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "confidence_density.png", build = function(d, ctx) confidence_density_plot(d, ctx$confidence_min), data = "combined", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "confidence_vs_evidence.png", build = function(d, ctx) confidence_vs_evidence_plot(d), data = "combined", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "structure_by_tier.png", build = function(d, ctx) structure_by_tier_plot(d), data = "combined", n_x = NULL, axis = "x", segment = TRUE),
    list(file = "source_yield.png", build = function(d, ctx) source_yield_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    list(file = "taxon_by_source.png", build = function(d, ctx) taxon_by_source_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    list(file = "domain_tier_composition.png", build = function(d, ctx) domain_tier_composition_plot(d), data = "loci", n_x = "species", axis = "x", segment = TRUE),
    list(file = "structure_class_composition.png", build = function(d, ctx) structure_class_composition_plot(d), data = "combined", n_x = "species", axis = "x", segment = TRUE),
    # Tree panels: height grows with tip count, and the TREE supplies the tip
    # labels, so these route through auto_dims rather than the categorical axis.
    # The two taxonomy-cladogram trees are one tip for a single segment.
    list(file = "taxon_confidence_tree.png", build = function(d, ctx) taxon_confidence_tree_plot(d, ctx$tree_dir), data = "combined", n_x = "taxa", axis = "y", segment = FALSE),
    list(file = "taxon_tier_tree.png", build = function(d, ctx) taxon_tier_tree_plot(d, ctx$tree_dir), data = "combined", n_x = "taxa", axis = "y", segment = FALSE),
    list(file = "species_confidence_tree.png", build = function(d, ctx) species_confidence_tree_plot(d, ctx$tree_dir), data = "combined", n_x = "species", axis = "y", segment = TRUE),
    list(file = "species_composition_tree.png", build = function(d, ctx) species_composition_tree_plot(d, ctx$tree_dir), data = "combined", n_x = "species", axis = "y", segment = TRUE)
  )
}


# Entries a per-segment panel should render, per `plots.segment_panel`:
#   full     every entry that is meaningful for one segment (the default)
#   curated  the small legacy subset - who carries it, how sure, what shape
#   none     tables only
# Unknown values fall back to `full` rather than silently rendering nothing.
.CURATED_SEGMENT_PLOTS <- c("taxon_composition.png", "confidence_gradient.png",
                            "structure_class_composition.png")

segment_panel <- function(registry, mode = "full") {
  if (identical(mode, "none")) return(list())
  meaningful <- Filter(function(e) isTRUE(e$segment), registry)
  if (identical(mode, "curated")) {
    return(Filter(function(e) e$file %in% .CURATED_SEGMENT_PLOTS, meaningful))
  }
  meaningful
}


# Numeric companions for the evidence/confidence builders. The loci tables and
# the catalog store EVERY column as a string, so any builder needing a number
# derives it from here rather than coercing on the spot. Guarded so an all-empty
# input stays well-formed.
#
# Shared rather than inlined in main(): taxonomy_segments.R reuses these builders
# without running main(), and when it did not derive these columns
# `confidence_gradient_plot` hit its `confidence_num` guard and returned the
# empty placeholder for every segment (see test-taxonomy_segments.R).
add_numeric_companions <- function(df) {
  if (nrow(df) == 0L) return(df)
  df %>%
    mutate(
      confidence_num   = suppressWarnings(as.numeric(.data$confidence)),
      n_hits           = suppressWarnings(as.integer(.data$n_blastx_hits)),
      completeness_num = suppressWarnings(as.numeric(.data$completeness))
    )
}


taxon_composition_plot <- function(loci) {
  d <- loci %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("no confident taxon calls"))
  counts <- d %>% count(.data$species, .data$taxon_call, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon_call)) +
    geom_col() +
    scale_fill_igv() +
    labs(x = NULL, y = "loci (axis-resolved)", fill = "taxon") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "ERV taxon composition", "Confident (axis-resolved) calls per species")
}

# Resolved-rank distribution per species (how deep the call goes).
rank_resolution_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  order_lvl <- c("genus", "subfamily", "family", "none")
  d <- loci %>%
    mutate(rank = factor(ifelse(.data$rank %in% order_lvl, .data$rank, "none"),
                         levels = order_lvl)) %>%
    count(.data$species, .data$rank, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$rank)) +
    geom_col(position = "fill") +
    scale_fill_npg() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "resolved rank") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Taxonomic resolution", "Resolved rank per species")
}

# How each taxon call was made (placement vs weighted-LCA vs presence).
method_mix_plot <- function(loci) {
  d <- loci %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("no confident taxon calls"))
  counts <- d %>% count(.data$species, .data$method, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$method)) +
    geom_col() +
    scale_fill_aaas() +
    labs(x = NULL, y = "taxon calls", fill = "method") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Taxon-call method mix", "Placement vs weighted-LCA vs presence")
}

# ERV class (I/II/III) composition per species - the literature anchor.
erv_class_composition_plot <- function(loci) {
  d <- loci %>% filter(nzchar(.data$erv_class))
  if (nrow(d) == 0L) return(empty_plot("no ERV-class assignments"))
  counts <- d %>% count(.data$species, .data$erv_class, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$erv_class)) +
    geom_col(position = "fill") +
    scale_fill_jco() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of classified loci", fill = "ERV class") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "ERV class composition",
             "Class I / II / III per species (Jern/Blomberg)")
}

# Gene -> taxon flows across mosaic loci (loci whose member genes disagree).
# mosaic_composition is packed as "GENE:Taxon;GENE:Taxon"; unpack to flows.
mosaic_alluvial_plot <- function(loci) {
  d <- loci %>% filter(.data$is_mosaic == "True", nzchar(.data$mosaic_composition))
  if (nrow(d) == 0L) return(empty_plot("no mosaic loci"))
  flows <- d %>%
    mutate(.locus = row_number()) %>%
    separate_rows("mosaic_composition", sep = ";") %>%
    separate("mosaic_composition", into = c("gene", "taxon"),
             sep = ":", fill = "right", extra = "merge") %>%
    filter(nzchar(.data$gene), nzchar(.data$taxon))
  if (nrow(flows) == 0L) return(empty_plot("no mosaic loci"))
  counts <- flows %>% count(.data$gene, .data$taxon, name = "n")
  p <- ggplot(counts,
              aes(axis1 = .data$gene, axis2 = .data$taxon, y = .data$n)) +
    geom_alluvium(aes(fill = .data$taxon)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("gene", "taxon"), expand = c(0.1, 0.1)) +
    scale_fill_igv() +
    labs(y = "mosaic-locus gene calls", fill = "taxon") +
    theme_bw()
  add_titles(p, "Mosaic composition", "Per-gene taxon calls within mosaic loci")
}

# Unpack mosaic_composition ("GENE:Taxon;GENE:Taxon") to a long gene->taxon frame
# with a per-locus id (.locus). Shared by the mosaic sub-panel builders below.
.unpack_mosaic <- function(loci) {
  loci %>%
    filter(.data$is_mosaic == "True", nzchar(.data$mosaic_composition)) %>%
    mutate(.locus = dplyr::row_number()) %>%
    separate_rows("mosaic_composition", sep = ";") %>%
    separate("mosaic_composition", into = c("gene", "taxon"),
             sep = ":", fill = "right", extra = "merge") %>%
    filter(nzchar(.data$gene), nzchar(.data$taxon))
}

# Mosaic burden: fraction of LTR-flanked loci per species whose member genes
# disagree on their axis taxon (a recombination-load proxy). Orphans are
# single-gene and never mosaic, so this reads the LTR-flanked loci frame.
mosaic_burden_plot <- function(loci) {
  if (nrow(loci) == 0L || !"is_mosaic" %in% names(loci)) return(empty_plot("no loci"))
  counts <- loci %>%
    mutate(kind = ifelse(.data$is_mosaic == "True", "mosaic", "single-lineage")) %>%
    count(.data$species, .data$kind, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$kind)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c(mosaic = "#D95F02", `single-lineage` = "#7570B3")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of ltr-flanked loci", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Mosaic burden",
             "Fraction of loci with discordant gene taxa, per species")
}

# Recombination-partner heatmap: within mosaic loci, how often each unordered
# pair of axis taxa co-occurs (which lineages recombine, e.g. BetaxGamma). Upper
# triangle via the taxon.x < taxon.y filter on a per-locus self-join.
mosaic_taxon_pairs_plot <- function(loci) {
  flows <- .unpack_mosaic(loci) %>% distinct(.data$.locus, .data$taxon)
  if (nrow(flows) == 0L) return(empty_plot("no mosaic loci"))
  pairs <- flows %>%
    dplyr::inner_join(flows, by = ".locus", relationship = "many-to-many") %>%
    filter(.data$taxon.x < .data$taxon.y) %>%
    count(.data$taxon.x, .data$taxon.y, name = "n")
  if (nrow(pairs) == 0L) return(empty_plot("no co-occurring taxon pairs"))
  p <- ggplot(pairs, aes(x = .data$taxon.x, y = .data$taxon.y, fill = .data$n)) +
    geom_tile(colour = "grey80") +
    geom_text(aes(label = .data$n), size = 3, fontface = "bold") +
    scale_fill_viridis_c(trans = "log10") +
    labs(x = NULL, y = NULL, fill = "mosaic loci") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Recombination partners",
             "Co-occurring axis-taxon pairs within mosaic loci")
}

# Per-gene discordance: within each mosaic locus, the majority (consensus) taxon
# is the backbone; genes calling a different taxon are the recombinant signal.
# Fraction discordant per gene surfaces env-capture (ENV usually highest).
mosaic_gene_discordance_plot <- function(loci) {
  flows <- .unpack_mosaic(loci)
  if (nrow(flows) == 0L) return(empty_plot("no mosaic loci"))
  consensus <- flows %>%
    count(.data$.locus, .data$taxon, name = "n") %>%
    group_by(.data$.locus) %>%
    dplyr::slice_max(.data$n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(".locus", consensus = "taxon")
  disc <- flows %>%
    dplyr::left_join(consensus, by = ".locus") %>%
    group_by(.data$gene) %>%
    summarise(frac = mean(.data$taxon != .data$consensus), n = dplyr::n(),
              .groups = "drop")
  p <- ggplot(disc, aes(x = stats::reorder(.data$gene, -.data$frac), y = .data$frac)) +
    geom_col(fill = "#D95F02") +
    geom_text(aes(label = .data$n), vjust = -0.3, size = 3) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
    labs(x = "gene", y = "fraction discordant vs locus consensus") +
    theme_bw()
  add_titles(p, "Per-gene discordance",
             "How often each gene breaks from its locus's consensus taxon (n = gene occurrences)")
}

# Per-species mosaic composition: within mosaic loci only, the stacked axis-taxon
# mix per species - which genera drive the chimeras in each host.
mosaic_composition_by_species_plot <- function(loci) {
  flows <- .unpack_mosaic(loci)
  if (nrow(flows) == 0L) return(empty_plot("no mosaic loci"))
  flows <- collapse_long_tail(flows, "taxon", top_n = 20)
  counts <- flows %>% count(.data$species, .data$taxon, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon)) +
    geom_col(position = "fill") +
    scale_fill_igv() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of mosaic gene calls", fill = "taxon") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Mosaic composition by species",
             "Axis-taxon mix within mosaic loci, per host")
}

# Confidence-tag composition per species (HC/LC), faceted by tier so LTR-flanked
# loci and recovered orphans are both visible. Reads the combined frame.
confidence_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no classified loci"))
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  counts <- combined %>%
    count(.data$species, .data$source, .data$confidence_tag, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$confidence_tag)) +
    geom_col(position = "fill") +
    facet_wrap(~ .data$source) +
    scale_fill_manual(values = c(HC = "#1B9E77", LC = "#D95F02")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "confidence") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Call confidence", "High vs low confidence (< confidence_min) per species")
}

# Stacked COUNT bar of loci by confidence tag (HC/LC) per species, faceted by tier.
# Unlike confidence_plot (proportional), this shows absolute counts.
confidence_count_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no classified loci"))
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  counts <- combined %>%
    count(.data$species, .data$source, .data$confidence_tag, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$confidence_tag)) +
    geom_col(position = "stack") +
    facet_wrap(~ .data$source) +
    scale_fill_manual(values = c(HC = "#1B9E77", LC = "#D95F02")) +
    labs(x = NULL, y = "loci", fill = "confidence") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Confidence composition (counts)",
             "HC vs LC locus counts per species, by tier")
}

# Stacked COUNT bar with a continuous confidence GRADIENT: confidence binned
# (0.05 steps), each species' bar stacked low->high and coloured by the bin
# midpoint on a sequential viridis scale. Reveals the intra-HC spread that the
# binary HC/LC view collapses (most calls are HC).
confidence_gradient_plot <- function(combined) {
  if (nrow(combined) == 0L || !"confidence_num" %in% names(combined)) {
    return(empty_plot("no confidence values"))
  }
  d <- combined %>% filter(!is.na(.data$confidence_num))
  if (nrow(d) == 0L) return(empty_plot("no confidence values"))
  if (!"source" %in% names(d)) d$source <- "ltr-flanked"
  brks <- seq(0, 1, by = 0.05)
  d <- d %>% mutate(
    bin = cut(.data$confidence_num, breaks = brks, include.lowest = TRUE, right = FALSE),
    mid = brks[as.integer(.data$bin)] + 0.025
  )
  counts <- d %>% count(.data$species, .data$source, .data$bin, .data$mid, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n,
                          group = .data$bin, fill = .data$mid)) +
    geom_col(position = position_stack(reverse = TRUE), colour = NA) +
    facet_wrap(~ .data$source) +
    scale_fill_viridis_c(name = "confidence", limits = c(0, 1)) +
    labs(x = NULL, y = "loci") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Confidence distribution (gradient)",
             "Per-species locus counts stacked by confidence bin")
}


# ----------------------------------------------------------------------------
# Tree-attached confidence panels (ADR-011)
# ----------------------------------------------------------------------------
# A tree beside the bars turns an alphabetical list into a phylogenetic
# statement: you can see whether ERV burden tracks host relatedness, or which
# viral lineages dominate where. The tree is drawn from coordinates precomputed
# by tree_layout.py (Bio.Phylo), so no R tree library is needed; patchwork
# aligns the two panels on a shared y scale.
#
# Softer sequential palette than viridis - these panels are meant to be read
# next to a tree, where saturated colour fights the topology for attention.
.SOFT_CONFIDENCE <- c("#efe7d6", "#dbe3cf", "#bfd4c6", "#9cc0bf", "#7aa7ad", "#5f88a1")

# Read a tree coordinate CSV written by tree_layout.py. Missing/empty file =>
# NULL, which the builders below turn into an explanatory placeholder rather
# than an error: a study with no species tree configured is a normal state.
read_tree_part <- function(dir, name, part) {
  f <- file.path(dir, sprintf("%s.tree_%s.csv", name, part))
  if (!file.exists(f)) return(NULL)
  df <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE))
  if (nrow(df) == 0L) NULL else df
}


# Horizontal confidence-gradient bars whose y order is fixed by a tree, with the
# tree drawn alongside. `key` is the catalog column the tips correspond to
# (`species` or `taxon_call`) - level-agnostic: the taxon tree's tips are
# whatever ranks the calls resolved to.
tree_confidence_plot <- function(combined, tree_dir, tree_name, key,
                                 title, subtitle) {
  tips <- read_tree_part(tree_dir, tree_name, "tips")
  segs <- read_tree_part(tree_dir, tree_name, "segments")
  if (is.null(tips)) {
    return(empty_plot(sprintf("no %s tree available", tree_name)))
  }
  if (nrow(combined) == 0L || !"confidence_num" %in% names(combined) ||
      !key %in% names(combined)) {
    return(empty_plot("no confidence values"))
  }
  d <- combined %>% filter(!is.na(.data$confidence_num))
  d <- d[as.character(d[[key]]) %in% tips$tip, , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("no loci matching the tree tips"))
  if (!"source" %in% names(d)) d$source <- "ltr-flanked"

  brks <- seq(0, 1, by = 0.05)
  d <- d %>% mutate(
    bin = cut(.data$confidence_num, breaks = brks, include.lowest = TRUE,
              right = FALSE),
    mid = brks[as.integer(.data$bin)] + 0.025,
    .y  = tips$y[match(as.character(.data[[key]]), tips$tip)]
  )
  counts <- d %>% count(.data$.y, .data$source, .data$bin, .data$mid, name = "n")
  ylim <- c(0.4, nrow(tips) + 0.6)

  bars <- ggplot(counts, aes(x = .data$n, y = .data$.y,
                             group = .data$bin, fill = .data$mid)) +
    geom_col(position = position_stack(reverse = TRUE), colour = NA,
             orientation = "y") +
    facet_wrap(~ .data$source) +
    scale_fill_gradientn(colours = .SOFT_CONFIDENCE, limits = c(0, 1),
                         name = "confidence") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = "loci (stacked low -> high confidence)", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  tree <- ggplot() +
    { if (!is.null(segs)) {
        geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend),
                     colour = "grey45", linewidth = 0.4, lineend = "round")
      } } +
    geom_text(data = tips, aes(x = .data$x, y = .data$y, label = .data$tip),
              hjust = -0.05, size = 3, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.9))) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_void()

  p <- patchwork::wrap_plots(tree, bars, widths = c(1.2, 3))
  # add_titles() styles a ggplot; a patchwork needs its annotation instead.
  p + patchwork::plot_annotation(
    title = title, subtitle = subtitle,
    theme = theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle   = element_text(hjust = 0.5, size = 11),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )
}


# Viral lineages on the y axis, ordered by the reference taxonomy cladogram.
taxon_confidence_tree_plot <- function(combined, tree_dir) {
  tree_confidence_plot(
    combined, tree_dir, "taxon", "taxon_call",
    "ERV confidence by viral lineage",
    "Reference taxonomy cladogram - stacked by confidence, split by tier"
  )
}

# Counts/composition sibling of tree_confidence_plot (ADR-014). Same tree
# machinery, but the bars show what each tip is MADE OF rather than how
# confident it is: `fill_col` is any categorical catalog column. Keeping the two
# separate rather than adding a mode flag keeps each one readable.
tree_composition_plot <- function(combined, tree_dir, tree_name, key, fill_col,
                                  title, subtitle) {
  tips <- read_tree_part(tree_dir, tree_name, "tips")
  segs <- read_tree_part(tree_dir, tree_name, "segments")
  if (is.null(tips)) {
    return(empty_plot(sprintf("no %s tree available", tree_name)))
  }
  if (nrow(combined) == 0L || !key %in% names(combined) ||
      !fill_col %in% names(combined)) {
    return(empty_plot(sprintf("no %s values", fill_col)))
  }
  # The ADR-011 tip-label trap: a locus keyed by a name the tree does not carry
  # must be dropped, not silently drawn at whatever y position match() returns.
  d <- combined[as.character(combined[[key]]) %in% tips$tip, , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("no loci matching the tree tips"))

  d$.y <- tips$y[match(as.character(d[[key]]), tips$tip)]
  counts <- d %>%
    count(.data$.y, .fill = as.character(.data[[fill_col]]), name = "n")
  ylim <- c(0.4, nrow(tips) + 0.6)

  bars <- ggplot(counts, aes(x = .data$n, y = .data$.y, fill = .data$.fill)) +
    geom_col(position = position_stack(reverse = TRUE), colour = NA,
             orientation = "y") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    scale_fill_manual(values = igv_unlimited_palette(
      length(unique(counts$.fill))), name = fill_col) +
    labs(x = "loci", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  tree <- ggplot() +
    { if (!is.null(segs)) {
        geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend),
                     colour = "grey45", linewidth = 0.4, lineend = "round")
      } } +
    geom_text(data = tips, aes(x = .data$x, y = .data$y, label = .data$tip),
              hjust = -0.05, size = 3, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.9))) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_void()

  patchwork::wrap_plots(tree, bars, widths = c(1.2, 3)) +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.subtitle   = element_text(hjust = 0.5, size = 11),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )
}


# Which viral lineages each host carries, ordered by host relatedness. Reads
# directly as "did related hosts keep related viruses?" - the catalog-side
# counterpart to the placement co-phylogeny (ADR-014).
species_composition_tree_plot <- function(combined, tree_dir) {
  tree_composition_plot(
    combined, tree_dir, "species", "species", "segment",
    "Viral lineage composition by host species",
    "Host phylogeny (input.species_tree) - loci per lineage"
  )
}

# Per-lineage evidence tiers. A lineage that is almost entirely orphan is one
# whose structural evidence has eroded away - a decay signal that the pooled
# counts hide.
taxon_tier_tree_plot <- function(combined, tree_dir) {
  tree_composition_plot(
    combined, tree_dir, "taxon", "taxon_call", "source",
    "Evidence tier by viral lineage",
    "Reference taxonomy cladogram - LTR-flanked (LTR-confirmed) vs orphan"
  )
}


# Host species on the y axis, ordered by the user-supplied species phylogeny.
species_confidence_tree_plot <- function(combined, tree_dir) {
  tree_confidence_plot(
    combined, tree_dir, "species", "species",
    "ERV confidence by host species",
    "Host phylogeny (input.species_tree) - stacked by confidence, split by tier"
  )
}

# Bucket a per-locus blastx hit count into ordered evidence bands. Pure helper
# (unit-tested): 0 / 1 / 2-5 / 6+. Robust to numeric (non-integer) input.
bucket_evidence <- function(n) {
  n <- as.integer(n)
  out <- dplyr::case_when(
    is.na(n) ~ NA_character_,
    n <= 0L  ~ "0",
    n == 1L  ~ "1",
    n <= 5L  ~ "2-5",
    TRUE     ~ "6+"
  )
  factor(out, levels = c("0", "1", "2-5", "6+"))
}

# Per-locus blastx evidence depth; the zero bin is the candidate-novel-retrovirus
# pile. Pseudo-log x so the long tail of well-supported loci stays readable.
evidence_depth_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no classified loci"))
  d <- combined %>% mutate(is_novel = .data$n_hits == 0L)
  n_novel <- sum(d$is_novel, na.rm = TRUE)
  p <- ggplot(d, aes(x = .data$n_hits, fill = .data$is_novel)) +
    geom_histogram(binwidth = 1, colour = "grey30", linewidth = 0.2) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_fill_manual(values = c("FALSE" = "#4575B4", "TRUE" = "#D73027"),
                      labels = c("FALSE" = "has homology", "TRUE" = "novel (0 hits)")) +
    facet_wrap(~ .data$source, scales = "free_y") +
    labs(x = "blastx hits per locus (pseudo-log)", y = "loci", fill = NULL) +
    theme_bw()
  add_titles(p, "Blastx evidence depth",
             sprintf("Per-locus homology - %d novel candidates (0 hits)", n_novel))
}

# Raw confidence distribution split by method, with the HC/LC threshold line.
# A histogram (not a KDE): confidence values pile up at discrete points - most
# at exactly 1.000 - so a density estimate smears mass past the [0,1] domain and
# drops single-value groups (e.g. presence, always 1.0). Faceting by method keeps
# the very different scales legible; the dashed line marks confidence_min.
confidence_density_plot <- function(combined, confidence_min = 0.5) {
  d <- combined %>% filter(!is.na(.data$confidence_num))
  if (nrow(d) == 0L) return(empty_plot("no confidence values"))
  p <- ggplot(d, aes(x = .data$confidence_num, fill = .data$method)) +
    geom_histogram(binwidth = 0.05, boundary = 0, colour = "grey30", linewidth = 0.15) +
    geom_vline(xintercept = confidence_min, linetype = "dashed", colour = "grey20") +
    facet_wrap(~ .data$method, scales = "free_y", ncol = 1) +
    scale_fill_aaas(guide = "none") +
    scale_x_continuous(limits = c(-0.02, 1.02)) +
    labs(x = "call confidence", y = "loci") +
    theme_bw()
  add_titles(p, "Confidence calibration",
             sprintf("Confidence by method; dashed = confidence_min (%.2f)", confidence_min))
}

# Call confidence across blastx evidence-depth buckets - does more homology mean
# a more confident call?
confidence_vs_evidence_plot <- function(combined) {
  d <- combined %>%
    filter(!is.na(.data$confidence_num), !is.na(.data$n_hits)) %>%
    mutate(bucket = bucket_evidence(.data$n_hits))
  if (nrow(d) == 0L) return(empty_plot("no data"))
  p <- ggplot(d, aes(x = .data$bucket, y = .data$confidence_num)) +
    geom_boxplot(fill = "#74ADD1", outlier.size = 0.6) +
    labs(x = "blastx hits per locus", y = "call confidence") +
    theme_bw()
  add_titles(p, "Confidence vs evidence", "Call confidence across blastx hit-count buckets")
}

# Structural completeness by tier - how many main genes each locus carries,
# LTR-flanked proviruses vs recovered orphans. (Replaces a novel-vs-classified
# view: 0-hit "novel" loci are essentially absent here - they are domain-validated
# so they have homology - so that comparison was empty. This populated view is
# the useful one: it shows orphans are structurally simpler, mostly single
# markers, while LTR-flanked loci carry more of the gag/pol/env complement.)
# Counts are shown as a fraction within each tier so the two tiers' very
# different sizes don't swamp the comparison.
structure_by_tier_plot <- function(combined) {
  if (nrow(combined) == 0L || !"completeness_num" %in% names(combined)) {
    return(empty_plot("no loci"))
  }
  d <- combined %>% filter(!is.na(.data$completeness_num))
  if (nrow(d) == 0L) return(empty_plot("no completeness data"))
  p <- ggplot(d, aes(x = .data$completeness_num, y = after_stat(.data$density),
                     fill = .data$source)) +
    geom_histogram(binwidth = 0.1, boundary = 0, position = "identity",
                   alpha = 0.5, colour = "grey40", linewidth = 0.15) +
    scale_fill_manual(values = .SOURCE_FILL) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "completeness (fraction of main genes present)",
         y = "within-tier density", fill = "tier") +
    theme_bw()
  add_titles(p, "Structural completeness by tier",
             "LTR-flanked proviruses carry more genes; orphans are mostly single markers")
}

# Yield boost from the orphans tier - loci recovered per tier, per species.
source_yield_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no loci"))
  counts <- combined %>% count(.data$species, .data$source, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$source)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = .SOURCE_FILL) +
    labs(x = NULL, y = "loci", fill = "tier") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "LTR-flanked vs orphan yield", "Loci recovered per tier, per species")
}

# Taxon composition split by tier - surfaces taxa present only in the orphan
# tier (the novel-lineage check). Long taxon tail folded via collapse_long_tail.
taxon_by_source_plot <- function(combined) {
  d <- combined %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("no confident taxon calls"))
  d <- collapse_long_tail(d, "taxon_call", top_n = 20)
  counts <- d %>% count(.data$species, .data$source, .data$taxon_call, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon_call)) +
    geom_col() +
    facet_wrap(~ .data$source, scales = "free_x") +
    scale_fill_igv() +
    labs(x = NULL, y = "axis-resolved loci", fill = "taxon") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Taxon composition by tier",
             "Taxa per species, ltr-flanked vs orphan (new-lineage check)")
}


# Per-provirus domain-tier composition on LTR-flanked loci: the recall the new
# labelling preserves - domain_selected (config-matched Pfam) vs domain_unlisted
# (a Pfam domain, just not in the curated set) vs non_domain (LTR-flanked purely by
# position). Fraction within species so genome size doesn't swamp the mix.
domain_tier_composition_plot <- function(loci) {
  if (nrow(loci) == 0L || !"domain_tier" %in% names(loci)) {
    return(empty_plot("no ltr-flanked loci"))
  }
  d <- loci %>%
    mutate(domain_tier = factor(.data$domain_tier, levels = .DOMAIN_TIER_LEVELS))
  counts <- d %>% count(.data$species, .data$domain_tier, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$domain_tier)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = .DOMAIN_TIER_FILL, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of ltr-flanked loci", fill = "domain tier") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Domain-tier composition",
             "Per-provirus domain support of ltr-flanked loci (recall preserved)")
}


# Discrete structural class (full / partial / gene) per species, faceted by tier.
# The first catalogued full-vs-partial-vs-single-gene view (ADR-009).
structure_class_composition_plot <- function(combined) {
  if (nrow(combined) == 0L || !"structure_class" %in% names(combined)) {
    return(empty_plot("no loci"))
  }
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  d <- combined %>%
    mutate(structure_class = factor(.data$structure_class, levels = .STRUCTURE_LEVELS))
  counts <- d %>% count(.data$species, .data$source, .data$structure_class, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$structure_class)) +
    geom_col(position = "fill") +
    facet_wrap(~ .data$source) +
    scale_fill_manual(values = .STRUCTURE_FILL, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "structure") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Structural class composition",
             "Full / partial / gene per species, by tier")
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Generate RetroSeek taxonomy panel (per-locus taxon-call plots)"
  )
  parser$add_argument("--input", required = TRUE,
                      help = "Directory with per-genome <genome>.loci.parquet tables.")
  parser$add_argument("--output", required = TRUE,
                      help = "Directory to save output plots.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  parser$add_argument("--report_csv", required = TRUE,
                      help = "Output path for the tidy classification report CSV.")
  parser$add_argument("--tree_dir", required = FALSE, default = "",
                      help = paste("Directory of tree coordinate CSVs from",
                                   "tree_layout.py. Absent/empty renders the",
                                   "tree panels as placeholders."))
  parser$add_argument("--catalog_csv", required = TRUE,
                      help = paste("Output path for the unified authoritative ERV",
                                   "catalog CSV (ltr-flanked proviruses + clustered",
                                   "orphan loci, one non-overlapping record each)."))
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi       <- cfg$plots$dpi    %||% 300
  plot_height    <- cfg$plots$height %||% 12
  plot_width     <- cfg$plots$width  %||% 15
  confidence_min <- cfg$classification$confidence_min %||% 0.5
  # Canvas growth per extra category past the base, and the hard ceiling in
  # inches (at 300 dpi, 60in is ~18,000 px - the practical PNG limit).
  per_stratum    <- cfg$plots$per_stratum %||% 0.18
  max_dim        <- cfg$plots$max_dim     %||% 60

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
  log_section(sprintf("RetroSeek taxonomy plot generation (output: %s)", args$output))

  loci <- load_loci(args$input)
  orphans <- load_orphans(args$input)
  # Relabel genome stems to the config display names (species: map) for all axes.
  if (nrow(loci) > 0L) loci$species <- relabel_species(loci$species, cfg$species)
  if (nrow(orphans) > 0L) {
    orphans$species <- relabel_species(orphans$species, cfg$species)
  }
  combined <- bind_rows(loci, orphans)
  combined <- add_numeric_companions(combined)
  n_species <- length(unique(combined$species))
  log_section(sprintf("Loaded %d ltr-flanked loci + %d recovered orphans across %d species",
                      nrow(loci), nrow(orphans), n_species))

  # `n_x` = number of categories on the x axis (species, for the per-species
  # panels). Supplying it scales BOTH the canvas and the tick-label text/angle,
  # so a 102-genome study stays readable; omitting it keeps the fixed config
  # canvas for plots whose x axis is not per-species (histograms, heatmaps,
  # alluvials). Applied here rather than in each builder so the whole panel's
  # scaling decisions are visible in one place.
  emit <- function(name, plot, n_x = NULL) {
    if (!is.null(n_x)) plot <- scale_categorical_axis(plot, n_x, axis = "x",
                                                      base_w = plot_width,
                                                      base_h = plot_height,
                                                      per_stratum = per_stratum,
                                                      cap = max_dim)
    save_plot(name, plot, args$output, dims = attr(plot, "intended_dims"),
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  # Every panel entry is declared once in panel_registry(); this loop is the
  # only place the full panel is emitted. `data` picks the tier scope, `axis`
  # picks the sizing rule - the tree panels grow on y via auto_dims, because the
  # TREE supplies the tip labels and scale_categorical_axis() would re-enable
  # axis.text.y and print row indices next to it.
  tree_dir <- args$tree_dir %||% ""
  n_taxa <- length(unique(combined$taxon_call))
  ctx <- list(tree_dir = tree_dir, confidence_min = confidence_min)
  for (e in panel_registry()) {
    d <- if (identical(e$data, "loci")) loci else combined
    p <- e$build(d, ctx)
    if (identical(e$axis, "y")) {
      n_tips <- if (identical(e$n_x, "taxa")) n_taxa else n_species
      dims <- auto_dims(n_tips, axis = "y", base_w = plot_width,
                        base_h = plot_height, per_stratum = per_stratum,
                        cap = max_dim)
      save_plot(e$file, p, args$output, dims = dims,
                base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
    } else {
      emit(e$file, p, if (identical(e$n_x, "species")) n_species else NULL)
    }
  }

  # Tidy report: counts by taxon / confidence / method + mosaic + integrations,
  # split by tier. Concordant with the plots (same combined frame).
  dir.create(dirname(args$report_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(build_report(combined), args$report_csv)

  # Unified authoritative catalog: every locus as ONE non-overlapping record -
  # LTR-flanked proviruses (LTR-confirmed) + clustered orphan loci (proximity-
  # inferred), `source` keeping the confidence gradient explicit. The single
  # "this is what we found at this location, and here's everything about it" table.
  catalog_cols <- c(
    "species", "source", "seqname", "start", "end", "strand",
    "taxon_call", "rank", "segment", "segment_rank",
    "resolved", "confidence", "confidence_tag", "erv_class",
    "structure_class", "domain_tier", "oversized", "canonical_order",
    "completeness", "n_main_genes", "genes_present", "is_mosaic",
    "mosaic_composition", "n_blastx_hits", "method", "id"
  )
  catalog <- reconcile_catalog(combined) %>% dplyr::select(dplyr::any_of(catalog_cols))
  if (nrow(catalog) > 0L && all(c("species", "seqname", "start") %in% names(catalog))) {
    catalog <- catalog %>%
      dplyr::arrange(.data$species, .data$seqname,
                     suppressWarnings(as.integer(.data$start)))
  }
  dir.create(dirname(args$catalog_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(catalog, args$catalog_csv)

  log_section(sprintf("Done - wrote 24 PNGs to %s + report %s + catalog %s",
                      args$output, args$report_csv, args$catalog_csv))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript taxonomy_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
