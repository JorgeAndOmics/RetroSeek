# =============================================================================
# taxonomy_plot_generator.R
# =============================================================================
# Builds RetroSeek's taxonomy panel: per-locus ERV taxon-call composition and
# resolution plots, derived from the per-genome `<genome>.loci.parquet` tables
# that taxonomy_classify_loci.py writes to data/tables/taxonomy_classification/.
# Because every plot is computed from those same loci tables, the panel is
# concordant with the tables by construction.
#
# Output is ONE PDF (results/plots/classification/taxonomy/taxonomy.pdf): a key
# page, then one page per panel_registry() entry, in registry order: lineage
# composition per host and per tier, call depth and method, ERV class, tier
# yield, structure and domain support, confidence in four views, the two
# lineage trees, and the mosaic (recombination) pages.
#
# Pages follow the house style (plot2sort/style.R, docs/visual_style.md): hosts
# on rows in the host-tree order with the tree beside them, fixed genus colours,
# readable species names.
#
# Shared infrastructure (style.R, helpers.R, tree_axis.R) is reused from
# plot2sort/. `testthat` sources this file; the `if (sys.nframe() == 0L) main()`
# guard keeps the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
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
source(file.path(.script_dir, "..", "utils", "log.R"))  # line contract, run_main (ADR-021)
source(file.path(.script_dir, "..", "plot2sort", "style.R"))  # palette, theme, labels, stage PDFs
source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))  # empty_plot, add_titles
source(file.path(.script_dir, "..", "plot2sort", "tree_axis.R"))  # species rows, host tree


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
    ranges   = IRanges::IRanges(as.integer(combined$start), as.integer(combined$end))
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
# Category levels (ADR-009). Declared in full and used with drop = FALSE, so a
# level absent from one genome still keeps its colour and legend slot. The
# colours themselves come from plot2sort/style.R.
# ----------------------------------------------------------------------------
.DOMAIN_TIER_LEVELS <- c("domain_selected", "domain_unlisted", "non_domain")
.STRUCTURE_LEVELS   <- c("full", "partial", "gene")
.RANK_LEVELS        <- c("genus", "subfamily", "family", "none")

# Facet strips name tiers and methods in words ("ltr-flanked" -> "LTR-flanked").
.word_strips <- ggplot2::as_labeller(display_label)


# ----------------------------------------------------------------------------
# Panel registry: ONE declaration, TWO consumers
# ----------------------------------------------------------------------------
# main() below and taxonomy_segments.R both render this list, so a page added
# here reaches both by construction. Previously each kept its own hand-written
# list, and they drifted: the per-segment subset stayed at 3 plots while this
# panel grew to 24.
#
# Fields:
#   name     page identifier (tests, the curated segment subset)
#   build    function(d, ctx) -> page. ctx is panel_ctx() (plot2sort/tree_axis.R):
#            the host tree and config order for species rows, the taxon tree
#            directory and confidence_min, so every entry has the same shape.
#   data     "loci"     = LTR-flanked only (the taxon-founded assembly)
#            "combined" = both tiers
#            Losing this distinction would quietly mix orphans into the
#            composition and mosaic pages, which deliberately exclude them.
#   segment  FALSE where the page is DEGENERATE for a single segment:
#            erv_class is constant within a genus (measured: exactly 1.00
#            distinct values), and the taxonomy cladogram collapses to one tip.
#
# Registry order is page order in taxonomy.pdf.
panel_registry <- function() {
  list(
    list(name = "taxon_composition", build = function(d, ctx) taxon_composition_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "lineage_composition", build = function(d, ctx) lineage_composition_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "taxon_by_source", build = function(d, ctx) taxon_by_source_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "rank_resolution", build = function(d, ctx) rank_resolution_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "method_mix", build = function(d, ctx) method_mix_plot(d, ctx), data = "loci", segment = TRUE),
    # erv_class is a function of genus, so within one segment it is constant.
    list(name = "erv_class_composition", build = function(d, ctx) erv_class_composition_plot(d, ctx), data = "loci", segment = FALSE),
    list(name = "source_yield", build = function(d, ctx) source_yield_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "structure_class_composition", build = function(d, ctx) structure_class_composition_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "structure_by_tier", build = function(d, ctx) structure_by_tier_plot(d), data = "combined", segment = TRUE),
    list(name = "domain_tier_composition", build = function(d, ctx) domain_tier_composition_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "confidence", build = function(d, ctx) confidence_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "confidence_count", build = function(d, ctx) confidence_count_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "confidence_gradient", build = function(d, ctx) confidence_gradient_plot(d, ctx), data = "combined", segment = TRUE),
    list(name = "confidence_density", build = function(d, ctx) confidence_density_plot(d, ctx$confidence_min), data = "combined", segment = TRUE),
    list(name = "evidence_depth", build = function(d, ctx) evidence_depth_plot(d), data = "combined", segment = TRUE),
    list(name = "confidence_vs_evidence", build = function(d, ctx) confidence_vs_evidence_plot(d), data = "combined", segment = TRUE),
    # The taxonomy cladogram is a single tip within one segment.
    list(name = "taxon_confidence_tree", build = function(d, ctx) taxon_confidence_tree_plot(d, ctx$tree_dir), data = "combined", segment = FALSE),
    list(name = "taxon_tier_tree", build = function(d, ctx) taxon_tier_tree_plot(d, ctx$tree_dir), data = "combined", segment = FALSE),
    list(name = "mosaic_burden", build = function(d, ctx) mosaic_burden_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "mosaic_composition_by_species", build = function(d, ctx) mosaic_composition_by_species_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "mosaic_alluvial", build = function(d, ctx) mosaic_alluvial_plot(d), data = "loci", segment = TRUE),
    list(name = "mosaic_taxon_pairs", build = function(d, ctx) mosaic_taxon_pairs_plot(d), data = "loci", segment = TRUE),
    list(name = "mosaic_gene_discordance", build = function(d, ctx) mosaic_gene_discordance_plot(d), data = "loci", segment = TRUE)
  )
}


# Entries a per-segment panel should render, per `plots.segment_panel`:
#   full     every entry that is meaningful for one segment (the default)
#   curated  the small legacy subset: who carries it, how sure, what shape
#   none     tables only
# Unknown values fall back to `full` rather than silently rendering nothing.
.CURATED_SEGMENT_PLOTS <- c("taxon_composition", "confidence_gradient",
                            "structure_class_composition")

segment_panel <- function(registry, mode = "full") {
  if (identical(mode, "none")) return(list())
  meaningful <- Filter(function(e) isTRUE(e$segment), registry)
  if (identical(mode, "curated")) {
    return(Filter(function(e) e$name %in% .CURATED_SEGMENT_PLOTS, meaningful))
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
      # Blank for loci without a call; NA is the intended reading.
      confidence_num   = suppressWarnings(as.numeric(.data$confidence)),  # blank -> NA
      n_hits           = suppressWarnings(as.integer(.data$n_blastx_hits)),  # blank -> NA
      completeness_num = suppressWarnings(as.numeric(.data$completeness))  # blank -> NA
    )
}


# ----------------------------------------------------------------------------
# Plot builders. Each takes the loci frame (and the panel ctx where species go on
# rows) and returns a page, or empty_plot() when there is nothing to show.
# Species-on-x builders hand their plot to on_rows(), which flips it so species
# are rows in the canonical order, beside the host tree when one is configured.
# ----------------------------------------------------------------------------

# Confident calls per host, coloured by lineage.
taxon_composition_plot <- function(loci, ctx = NULL) {
  d <- loci %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("No confident taxon calls"))
  counts <- d %>% count(.data$species, .data$taxon_call, name = "n")
  counts$taxon_call <- taxon_factor(counts$taxon_call, counts$n)
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon_call)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_taxon(counts$taxon_call, counts$n) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "LTR-flanked loci", fill = NULL)
  p <- add_titles(p, "Viral lineages per host",
                  "LTR-flanked loci with a confident call, by the lineage they resolve to.")
  on_rows(p, counts$species, ctx)
}

# Every locus of both tiers by the lineage it rolls up to at the segment rank,
# so unresolved calls are one honest "unassigned" bar rather than scattered
# higher ranks. Reads as "did related hosts keep related viruses?" (ADR-014).
lineage_composition_plot <- function(combined, ctx = NULL) {
  if (nrow(combined) == 0L || !"segment" %in% names(combined)) {
    return(empty_plot("No segmented loci"))
  }
  counts <- combined %>% count(.data$species, segment = as.character(.data$segment),
                               name = "n")
  counts$segment <- taxon_factor(counts$segment, counts$n)
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$segment)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_taxon(counts$segment, counts$n) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Loci, both tiers", fill = NULL)
  p <- add_titles(p, "Viral lineages per host, both tiers",
                  "Every locus by the lineage it rolls up to at the segment rank.")
  on_rows(p, counts$species, ctx)
}

# Lineages split by tier: a lineage seen only among orphans is either new or has
# lost its structural evidence (the new-lineage check).
taxon_by_source_plot <- function(combined, ctx = NULL) {
  d <- combined %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("No confident taxon calls"))
  d <- collapse_long_tail(d, "taxon_call", top_n = 20)
  counts <- d %>% count(.data$species, .data$source, .data$taxon_call, name = "n")
  counts$taxon_call <- taxon_factor(counts$taxon_call, counts$n)
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon_call)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    facet_wrap(~ .data$source, labeller = .word_strips) +
    scale_fill_taxon(counts$taxon_call, counts$n) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Loci with a confident call", fill = NULL)
  p <- add_titles(p, "Viral lineages per tier",
                  paste("Confident calls by lineage, per host and tier. A lineage seen",
                        "only among orphans is either new or has lost its LTRs."))
  on_rows(p, counts$species, ctx)
}

# How deep each call goes, from genus down to no resolution at all.
rank_resolution_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>%
    mutate(rank = factor(ifelse(.data$rank %in% .RANK_LEVELS, .data$rank, "none"),
                         levels = .RANK_LEVELS)) %>%
    count(.data$species, .data$rank, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$rank)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = .RANK_COLOUR, labels = display_label, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of LTR-flanked loci", fill = "Resolved to")
  p <- add_titles(p, "How deep the calls go",
                  "The rank each LTR-flanked locus resolves to, per host.")
  on_rows(p, d$species, ctx)
}

# How each confident call was made.
method_mix_plot <- function(loci, ctx = NULL) {
  d <- loci %>% filter(.data$resolved == "True")
  if (nrow(d) == 0L) return(empty_plot("No confident taxon calls"))
  counts <- d %>% count(.data$species, .data$method, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$method)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = category_colours(sort(unique(counts$method))),
                      labels = display_label) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Confident calls", fill = NULL)
  p <- add_titles(p, "How each call was made",
                  "Phylogenetic placement or weighted LCA of blastx hits, per host.")
  on_rows(p, counts$species, ctx)
}

# ERV class composition per host: the literature anchor (bats and human mostly
# Class I, mouse Class II).
erv_class_composition_plot <- function(loci, ctx = NULL) {
  d <- loci %>% filter(nzchar(.data$erv_class))
  if (nrow(d) == 0L) return(empty_plot("No ERV class assignments"))
  counts <- d %>% count(.data$species, .data$erv_class, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$erv_class)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = .ERV_CLASS_COLOUR) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of classified loci", fill = NULL)
  p <- add_titles(p, "ERV classes per host",
                  paste("Class I is gamma-like, Class II beta-like, Class III",
                        "spuma-like (Jern and Blomberg). Each wears its genus colour."))
  on_rows(p, counts$species, ctx)
}

# Loci recovered per tier, per host.
source_yield_plot <- function(combined, ctx = NULL) {
  if (nrow(combined) == 0L) return(empty_plot("No loci"))
  counts <- combined %>% count(.data$species, .data$source, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$source)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = .TIER_COLOUR, labels = display_label) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Loci", fill = NULL)
  p <- add_titles(p, "Loci per tier",
                  "LTR-flanked elements and recovered orphans, per host.")
  on_rows(p, counts$species, ctx)
}

# Full / partial / single-gene loci per host, by tier (ADR-009).
structure_class_composition_plot <- function(combined, ctx = NULL) {
  if (nrow(combined) == 0L || !"structure_class" %in% names(combined)) {
    return(empty_plot("No loci"))
  }
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  d <- combined %>%
    mutate(structure_class = factor(.data$structure_class, levels = .STRUCTURE_LEVELS))
  counts <- d %>% count(.data$species, .data$source, .data$structure_class, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$structure_class)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    facet_wrap(~ .data$source, labeller = .word_strips) +
    scale_fill_manual(values = .STRUCTURE_COLOUR, labels = display_label, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of loci", fill = NULL)
  p <- add_titles(p, "Structural class per host",
                  "Full, partial and single-gene loci, by tier.")
  on_rows(p, counts$species, ctx)
}

# Structural completeness by tier: how many main genes each locus carries. Shown
# as a density within each tier so the tiers' very different sizes do not swamp
# the comparison. (It replaced a novel-versus-classified view that was empty:
# domain-validated loci always have homology.)
structure_by_tier_plot <- function(combined) {
  if (nrow(combined) == 0L || !"completeness_num" %in% names(combined)) {
    return(empty_plot("No loci"))
  }
  d <- combined %>% filter(!is.na(.data$completeness_num))
  if (nrow(d) == 0L) return(empty_plot("No completeness data"))
  # Completeness is a fraction of a handful of main genes, so it takes only a
  # few values: side-by-side bars per value, not overlapping histograms.
  shares <- d %>%
    count(.data$source, completeness = round(.data$completeness_num, 2), name = "n") %>%
    group_by(.data$source) %>%
    mutate(share = .data$n / sum(.data$n)) %>%
    ungroup()
  steps <- sort(unique(shares$completeness))
  shares$completeness <- factor(scales::percent(shares$completeness, accuracy = 1),
                                levels = scales::percent(steps, accuracy = 1))
  p <- ggplot(shares, aes(x = .data$completeness, y = .data$share, fill = .data$source)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = .TIER_COLOUR, labels = display_label) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Main genes present", y = "Share of the tier's loci", fill = NULL) +
    theme(panel.grid.major.x = element_blank())
  add_titles(p, "Structural completeness by tier",
             "LTR-flanked elements carry more of their genes; orphans are mostly single markers.")
}

# Domain support of LTR-flanked loci: the recall the labelling preserves.
domain_tier_composition_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L || !"domain_tier" %in% names(loci)) {
    return(empty_plot("No LTR-flanked loci"))
  }
  d <- loci %>%
    mutate(domain_tier = factor(.data$domain_tier, levels = .DOMAIN_TIER_LEVELS))
  counts <- d %>% count(.data$species, .data$domain_tier, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$domain_tier)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = .DOMAIN_TIER_COLOUR, labels = display_label, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of LTR-flanked loci", fill = NULL)
  p <- add_titles(p, "Domain support of LTR-flanked loci",
                  paste("A retroviral Pfam domain, another Pfam domain, or none: an",
                        "element placed by its position between LTRs alone."))
  on_rows(p, counts$species, ctx)
}

# High versus low confidence, as shares (confidence_plot) or counts
# (confidence_count_plot), faceted by tier.
.confidence_bars <- function(combined, position, y_label, title, subtitle, ctx) {
  if (nrow(combined) == 0L) return(empty_plot("No classified loci"))
  if (!"source" %in% names(combined)) combined$source <- "ltr-flanked"
  counts <- combined %>%
    count(.data$species, .data$source, .data$confidence_tag, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$confidence_tag)) +
    geom_col(position = position, width = 0.7) +
    facet_wrap(~ .data$source, labeller = .word_strips) +
    scale_fill_manual(values = .CONFIDENCE_COLOUR, labels = display_label) +
    labs(x = NULL, y = y_label, fill = NULL)
  p <- p + scale_y_continuous(labels = if (inherits(position, "PositionFill")) {
    scales::percent
  } else {
    scales::comma
  })
  on_rows(add_titles(p, title, subtitle), counts$species, ctx)
}

confidence_plot <- function(combined, ctx = NULL) {
  .confidence_bars(combined, position_fill(reverse = TRUE), "Share of loci",
                   "Call confidence per host",
                   "Share of loci whose call is high or low confidence, by tier.", ctx)
}

confidence_count_plot <- function(combined, ctx = NULL) {
  .confidence_bars(combined, position_stack(reverse = TRUE), "Loci",
                   "Call confidence per host, as counts",
                   "Loci whose call is high or low confidence, by tier.", ctx)
}

# Loci stacked by confidence in 0.05 steps: the spread inside "high confidence"
# that the two-way split hides (most calls are high confidence).
confidence_gradient_plot <- function(combined, ctx = NULL) {
  if (nrow(combined) == 0L || !"confidence_num" %in% names(combined)) {
    return(empty_plot("No confidence values"))
  }
  d <- combined %>% filter(!is.na(.data$confidence_num))
  if (nrow(d) == 0L) return(empty_plot("No confidence values"))
  if (!"source" %in% names(d)) d$source <- "ltr-flanked"
  brks <- seq(0, 1, by = 0.05)
  d <- d %>% mutate(
    bin = cut(.data$confidence_num, breaks = brks, include.lowest = TRUE, right = FALSE),
    mid = brks[as.integer(.data$bin)] + 0.025
  )
  counts <- d %>% count(.data$species, .data$source, .data$bin, .data$mid, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n,
                          group = .data$bin, fill = .data$mid)) +
    geom_col(position = position_stack(reverse = TRUE), colour = NA, width = 0.7) +
    facet_wrap(~ .data$source, labeller = .word_strips) +
    scale_fill_ramp(limits = c(0, 1), name = "Confidence") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Loci")
  p <- add_titles(p, "Confidence distribution per host",
                  "Loci stacked from low to high confidence, in steps of 0.05, by tier.")
  on_rows(p, counts$species, ctx)
}

# Raw confidence distribution split by method, with the high/low threshold. A
# histogram, not a density estimate: confidence values pile up at discrete
# points (most at exactly 1.000), so a density smears mass past [0, 1] and drops
# single-value groups. One facet per method keeps their very different scales
# legible.
confidence_density_plot <- function(combined, confidence_min = 0.5) {
  d <- combined %>% filter(!is.na(.data$confidence_num))
  if (nrow(d) == 0L) return(empty_plot("No confidence values"))
  p <- ggplot(d, aes(x = .data$confidence_num)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = .DATA_COLOUR,
                   colour = .PAPER, linewidth = 0.2) +
    geom_vline(xintercept = confidence_min, linetype = "dashed", colour = .INK_SOFT) +
    facet_wrap(~ .data$method, scales = "free_y", ncol = 1, labeller = .word_strips) +
    # Zoom rather than limit the scale: scale limits would drop the edge bins.
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Call confidence", y = "Loci")
  add_titles(p, "Confidence calibration",
             sprintf(paste("Call confidence by method. Dashed line: the %.2f threshold",
                           "between high and low confidence."), confidence_min))
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

# Per-locus blastx evidence depth; the zero bin is the candidate novel retrovirus
# pile. Pseudo-log x so the long tail of well-supported loci stays readable.
evidence_depth_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("No classified loci"))
  d <- combined %>% mutate(is_novel = .data$n_hits == 0L)
  n_novel <- sum(d$is_novel, na.rm = TRUE)
  p <- ggplot(d, aes(x = .data$n_hits, fill = .data$is_novel)) +
    # `bins`, not `binwidth`: the bins are cut on the transformed (log) scale.
    geom_histogram(bins = 40, colour = NA) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 10, 100, 1000, 10000),
                       labels = scales::comma) +
    scale_fill_manual(values = c("FALSE" = .GREY_MID, "TRUE" = .DATA_COLOUR),
                      labels = c("FALSE" = "Has homology", "TRUE" = "No blastx hit")) +
    facet_wrap(~ .data$source, scales = "free_y", labeller = .word_strips) +
    labs(x = "Blastx hits per locus (pseudo-log scale)", y = "Loci", fill = NULL)
  add_titles(p, "Blastx evidence per locus",
             sprintf(paste("Hits per locus, by tier. %s loci have no hit at all:",
                           "candidate novel retroviruses."), scales::comma(n_novel)))
}

# Call confidence across blastx evidence-depth buckets: does more homology mean
# a more confident call? Buckets are ordered, so they take the ordinal ramp.
confidence_vs_evidence_plot <- function(combined) {
  d <- combined %>%
    filter(!is.na(.data$confidence_num), !is.na(.data$n_hits)) %>%
    mutate(bucket = bucket_evidence(.data$n_hits))
  if (nrow(d) == 0L) return(empty_plot("No data"))
  p <- ggplot(d, aes(x = .data$bucket, y = .data$confidence_num, fill = .data$bucket)) +
    geom_boxplot(outlier.size = 0.6, outlier.colour = .INK_SOFT, colour = .INK_SOFT,
                 show.legend = FALSE) +
    scale_fill_manual(values = stats::setNames(seq_colours(5)[2:5], levels(d$bucket))) +
    labs(x = "Blastx hits per locus", y = "Call confidence")
  add_titles(p, "Confidence against evidence",
             "Call confidence by the number of blastx hits supporting the locus.")
}


# ----------------------------------------------------------------------------
# Lineage trees (ADR-011, ADR-014)
# ----------------------------------------------------------------------------
# The reference taxonomy cladogram beside the bars turns a list of lineages into
# a phylogenetic statement. The tree comes from coordinates precomputed by
# tree_layout.py (Bio.Phylo), so no R tree library is needed; tree_column() and
# compose_with_tree() in plot2sort/tree_axis.R draw and align it.

# Horizontal confidence bars whose rows are fixed by a tree. `key` is the column
# the tips correspond to (`taxon_call`); the function is level-agnostic, since
# the taxon tree's tips are whatever ranks the calls resolved to.
tree_confidence_plot <- function(combined, tree_dir, tree_name, key,
                                 title, subtitle) {
  tips <- read_tree_part(tree_dir, tree_name, "tips")
  segs <- read_tree_part(tree_dir, tree_name, "segments")
  if (is.null(tips)) {
    return(empty_plot(sprintf("No %s tree available", tree_name)))
  }
  if (nrow(combined) == 0L || !"confidence_num" %in% names(combined) ||
      !key %in% names(combined)) {
    return(empty_plot("No confidence values"))
  }
  d <- combined %>% filter(!is.na(.data$confidence_num))
  d <- d[as.character(d[[key]]) %in% tips$tip, , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("No loci matching the tree tips"))
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
             orientation = "y", width = 0.7) +
    facet_wrap(~ .data$source, labeller = .word_strips) +
    scale_fill_ramp(limits = c(0, 1), name = "Confidence") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = "Loci, stacked from low to high confidence", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  compose_with_tree(tree_column(tips, segs, ylim),
                    bars + labs(title = title, subtitle = subtitle))
}

# Counts sibling of tree_confidence_plot: the bars show what each tip is MADE OF
# rather than how confident it is. `fill_col` is any categorical column and
# `colours` its named colours (category colours when NULL). Kept separate rather
# than adding a mode flag, so each stays readable.
tree_composition_plot <- function(combined, tree_dir, tree_name, key, fill_col,
                                  title, subtitle, colours = NULL) {
  tips <- read_tree_part(tree_dir, tree_name, "tips")
  segs <- read_tree_part(tree_dir, tree_name, "segments")
  if (is.null(tips)) {
    return(empty_plot(sprintf("No %s tree available", tree_name)))
  }
  if (nrow(combined) == 0L || !key %in% names(combined) ||
      !fill_col %in% names(combined)) {
    return(empty_plot(sprintf("No %s values", fill_col)))
  }
  # The ADR-011 tip-label trap: a locus keyed by a name the tree does not carry
  # must be dropped, not silently drawn at whatever y position match() returns.
  d <- combined[as.character(combined[[key]]) %in% tips$tip, , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("No loci matching the tree tips"))

  d$.y <- tips$y[match(as.character(d[[key]]), tips$tip)]
  counts <- d %>%
    count(.data$.y, .fill = as.character(.data[[fill_col]]), name = "n")
  ylim <- c(0.4, nrow(tips) + 0.6)
  if (is.null(colours)) colours <- category_colours(sort(unique(counts$.fill)))

  bars <- ggplot(counts, aes(x = .data$n, y = .data$.y, fill = .data$.fill)) +
    geom_col(position = position_stack(reverse = TRUE), colour = NA,
             orientation = "y", width = 0.7) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    scale_x_continuous(labels = scales::comma) +
    scale_fill_manual(values = colours, labels = display_label, name = NULL) +
    labs(x = "Loci", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  compose_with_tree(tree_column(tips, segs, ylim),
                    bars + labs(title = title, subtitle = subtitle))
}

taxon_confidence_tree_plot <- function(combined, tree_dir) {
  tree_confidence_plot(
    combined, tree_dir, "taxon", "taxon_call",
    "Confidence by viral lineage",
    paste("Lineages in the order of the reference taxonomy. Loci stacked from low",
          "to high confidence, by tier.")
  )
}

# Per-lineage evidence tiers. A lineage that is almost entirely orphan is one
# whose structural evidence has eroded away, a decay signal the pooled counts hide.
taxon_tier_tree_plot <- function(combined, tree_dir) {
  tree_composition_plot(
    combined, tree_dir, "taxon", "taxon_call", "source",
    "Tier by viral lineage",
    paste("Lineages in the order of the reference taxonomy. A lineage found almost",
          "only as orphans has lost its structural evidence."),
    colours = .TIER_COLOUR
  )
}


# ----------------------------------------------------------------------------
# Mosaic pages: loci whose genes resolve to different lineages
# ----------------------------------------------------------------------------

# Unpack mosaic_composition ("GENE:Taxon;GENE:Taxon") to a long gene/taxon frame
# with a per-locus id (.locus). Shared by the mosaic builders below.
.unpack_mosaic <- function(loci) {
  loci %>%
    filter(.data$is_mosaic == "True", nzchar(.data$mosaic_composition)) %>%
    mutate(.locus = dplyr::row_number()) %>%
    separate_rows("mosaic_composition", sep = ";") %>%
    separate("mosaic_composition", into = c("gene", "taxon"),
             sep = ":", fill = "right", extra = "merge") %>%
    filter(nzchar(.data$gene), nzchar(.data$taxon))
}

# Mosaic burden: the share of LTR-flanked loci per host whose genes disagree on
# their lineage (a recombination-load proxy). Orphans are single-gene and never
# mosaic, so this reads the LTR-flanked frame.
mosaic_burden_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L || !"is_mosaic" %in% names(loci)) return(empty_plot("No loci"))
  counts <- loci %>%
    mutate(kind = ifelse(.data$is_mosaic == "True", "mosaic", "single_lineage")) %>%
    count(.data$species, .data$kind, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$kind)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(mosaic = .DATA_COLOUR, single_lineage = .GREY_OTHER),
                      labels = display_label) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of LTR-flanked loci", fill = NULL)
  p <- add_titles(p, "Mosaic loci per host",
                  paste("Share of LTR-flanked loci whose genes resolve to different",
                        "lineages, a proxy for recombination."))
  on_rows(p, counts$species, ctx)
}

# Which lineages drive the chimeras in each host.
mosaic_composition_by_species_plot <- function(loci, ctx = NULL) {
  flows <- .unpack_mosaic(loci)
  if (nrow(flows) == 0L) return(empty_plot("No mosaic loci"))
  flows <- collapse_long_tail(flows, "taxon", top_n = 20)
  counts <- flows %>% count(.data$species, .data$taxon, name = "n")
  counts$taxon <- taxon_factor(counts$taxon, counts$n)
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$taxon)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_taxon(counts$taxon, counts$n) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of gene calls in mosaic loci", fill = NULL)
  p <- add_titles(p, "Lineages inside mosaic loci, per host",
                  "The lineage each gene of a mosaic locus resolves to.")
  on_rows(p, counts$species, ctx)
}

# Gene-to-lineage flows across mosaic loci.
mosaic_alluvial_plot <- function(loci) {
  flows <- .unpack_mosaic(loci)
  if (nrow(flows) == 0L) return(empty_plot("No mosaic loci"))
  counts <- flows %>% count(.data$gene, .data$taxon, name = "n")
  p <- ggplot(counts,
              aes(axis1 = .data$gene, axis2 = .data$taxon, y = .data$n)) +
    geom_alluvium(aes(fill = .data$taxon), alpha = 0.75, width = 0.25) +
    geom_stratum(fill = .PAPER, colour = .GREY_MID, width = 0.25) +
    # Lineage and gene names both take italics. Strata under 1.5% of the calls
    # are too thin to hold a label without overprinting their neighbours.
    geom_text(stat = "stratum",
              aes(label = after_stat(ifelse(prop > 0.015, as.character(stratum), ""))),
              size = 3.2, family = .FONT, fontface = "italic", colour = .INK) +
    scale_x_discrete(limits = c("Gene", "Lineage"), expand = c(0.1, 0.1)) +
    scale_fill_taxon(counts$taxon, counts$n) +
    labs(x = NULL, y = "Gene calls in mosaic loci", fill = NULL) +
    theme(panel.grid.major.x = element_blank())
  add_titles(p, "Genes and the lineages they call in mosaic loci",
             "Each flow is one gene of a mosaic locus, running to the lineage it resolves to.")
}

# Recombination partners: within mosaic loci, how often each unordered pair of
# lineages co-occurs (upper triangle via taxon.x < taxon.y on a self-join).
mosaic_taxon_pairs_plot <- function(loci) {
  flows <- .unpack_mosaic(loci) %>% distinct(.data$.locus, .data$taxon)
  if (nrow(flows) == 0L) return(empty_plot("No mosaic loci"))
  pairs <- flows %>%
    dplyr::inner_join(flows, by = ".locus", relationship = "many-to-many") %>%
    filter(.data$taxon.x < .data$taxon.y) %>%
    count(.data$taxon.x, .data$taxon.y, name = "n")
  if (nrow(pairs) == 0L) return(empty_plot("No co-occurring lineage pairs"))
  pairs$ink <- ink_on_ramp(pairs$n, trans = log10)
  p <- ggplot(pairs, aes(x = .data$taxon.x, y = .data$taxon.y, fill = .data$n)) +
    geom_tile(colour = .PAPER, linewidth = 0.6) +
    geom_text(aes(label = scales::comma(.data$n), colour = .data$ink), size = 3.2,
              family = .FONT) +
    scale_colour_identity() +
    scale_fill_ramp(trans = "log10", labels = scales::comma) +
    scale_x_discrete(labels = taxon_labels) +
    scale_y_discrete(labels = taxon_labels) +
    labs(x = NULL, y = NULL, fill = "Mosaic loci") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  add_titles(p, "Which lineages recombine",
             "Pairs of lineages found together in one mosaic locus. Numbers are loci.")
}

# Per-gene discordance: within each mosaic locus the majority lineage is the
# backbone, and a gene calling a different lineage is the recombinant signal.
mosaic_gene_discordance_plot <- function(loci) {
  flows <- .unpack_mosaic(loci)
  if (nrow(flows) == 0L) return(empty_plot("No mosaic loci"))
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
  p <- ggplot(disc, aes(x = stats::reorder(.data$gene, .data$frac), y = .data$frac)) +
    geom_col(fill = .DATA_COLOUR, width = 0.65) +
    geom_text(aes(label = scales::comma(.data$n)), hjust = -0.2, size = 3.2,
              family = .FONT) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Share of calls differing from the locus majority") +
    theme(panel.grid.major.y = element_blank())
  add_titles(p, "Which genes break from their locus",
             paste("Share of each gene's calls that differ from its locus's majority",
                   "lineage. Numbers are gene calls."))
}


# ----------------------------------------------------------------------------
# The key page: what taxonomy.pdf shows and what its colours mean.
# ----------------------------------------------------------------------------
# Tiers first, then every lineage the study resolved, in its fixed colour.
taxonomy_key_colours <- function(combined) {
  tiers <- stats::setNames(unname(.TIER_COLOUR[c("ltr-flanked", "orphan")]),
                           display_label(c("ltr-flanked", "orphan")))
  if (!"segment" %in% names(combined) || nrow(combined) == 0L) return(tiers)
  segments <- taxon_levels(combined$segment)
  c(tiers, stats::setNames(unname(taxon_colours(segments)[segments]),
                           display_label(segments)))
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
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "The stage PDF: a key page, then one page per panel entry.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  parser$add_argument("--report_csv", required = TRUE,
                      help = "Output path for the tidy classification report CSV.")
  parser$add_argument("--species_tree_dir", required = FALSE, default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--tree_dir", required = FALSE, default = "",
                      help = paste("Directory of tree coordinate CSVs from",
                                   "tree_layout.py. Absent/empty renders the",
                                   "tree panels as placeholders."))
  parser$add_argument("--catalog_csv", required = TRUE,
                      help = paste("Output path for the unified authoritative ERV",
                                   "catalog CSV (ltr-flanked proviruses + clustered",
                                   "orphan loci, one non-overlapping record each)."))
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "taxonomy_plot_generator")

  use_retroseek_style()
  cfg <- yaml::read_yaml(args$config)
  log_section(sprintf("RetroSeek taxonomy plot generation (output: %s)", args$out_pdf))

  loci <- load_loci(args$input)
  orphans <- load_orphans(args$input)
  # Readable species names (the config `species:` values) for every page; the
  # catalog below carries them too, as it always has.
  if (nrow(loci) > 0L) loci$species <- display_species(loci$species, cfg$species)
  if (nrow(orphans) > 0L) {
    orphans$species <- display_species(orphans$species, cfg$species)
  }
  combined <- bind_rows(loci, orphans)
  combined <- add_numeric_companions(combined)
  log_section(sprintf("Loaded %d ltr-flanked loci + %d recovered orphans across %d species",
                      nrow(loci), nrow(orphans), length(unique(combined$species))))

  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "", tree_dir = args$tree_dir %||% "")
  pages <- render_panel(panel_registry(), loci, combined, ctx)
  key <- key_page(
    "ERV taxonomy",
    paste("What lineage each ERV locus belongs to, how it was called and how sure the",
          "call is, for LTR-flanked elements and recovered orphans. Hosts are rows in",
          "the order of the host tree; each viral lineage keeps one colour across every",
          "RetroSeek figure, and lineages not resolved to a genus are grey."),
    colours = taxonomy_key_colours(combined), pages = page_titles(pages))
  # One page height for the whole PDF: tall enough for the most rows any page has.
  n_rows <- max(length(ctx$species_order), length(unique(combined$species)),
                length(unique(combined$taxon_call)))
  save_stage_pdf(c(list(key), pages), args$out_pdf,
                 height = page_height_for(n_rows, per_species = cfg$plots$per_stratum %||% 0.18))

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
      dplyr::arrange(.data$species, .data$seqname, as.integer(.data$start))
  }
  dir.create(dirname(args$catalog_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(catalog, args$catalog_csv)

  log_info("wrote %s (%d pages), report %s, catalog %s",
           args$out_pdf, length(pages) + 1L, args$report_csv, args$catalog_csv)
  log_ok("catalog of %s loci (%s LTR-flanked, %s orphans), %s pages",
         format(nrow(catalog), big.mark = ","), format(nrow(loci), big.mark = ","),
         format(nrow(orphans), big.mark = ","), format(length(pages) + 1L, big.mark = ","))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript taxonomy_plot_generator.R`.
# run_main() logs how the job ended (ADR-021).
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) run_main(main)
