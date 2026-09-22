# =============================================================================
# plot2sort/style.R - the one visual system for every RetroSeek figure
# =============================================================================
# Every colour, font, text rule and output format decision lives here, so a
# figure made anywhere in the pipeline looks like it belongs to the same study.
# The principles, and the reasons for them, are written up for people in
# docs/visual_style.md; this file is their executable form.
#
# Four rules drive everything below:
#
#   1. One colour, one meaning. A recurring concept (a tier, a fate, a genus) has
#      exactly one colour in every plot, and no colour stands for two concepts.
#      Plotting code never writes a hex value: it asks this file.
#   2. Colour-blind safe. The categorical palette is Paul Tol's "muted" set,
#      chosen because it keeps every pair distinguishable under deuteranopia,
#      protanopia and tritanopia (minimum CIELAB distance 15.8; the Futurama
#      palette this replaces fell to 4.7). test-style.R re-checks it.
#   3. Ordered things get ordered colours. Ordinal categories (structure class,
#      confidence) use a light-to-dark ramp, so "more" reads as darker before the
#      legend is consulted.
#   4. Readable text. IBM Plex Sans (shipped inside the conda env, so every machine
#      renders the same), left-aligned sentence-case titles, readable species
#      names, plain-word labels, and no dash used as punctuation.
#
# Species-axis behaviour (rows, canonical order, the host tree) lives beside this
# in plot2sort/tree_axis.R.

# ---------------------------------------------------------------------------
# Colour
# ---------------------------------------------------------------------------

# Paul Tol's "muted" qualitative palette, named for readability at call sites.
.TOL_MUTED <- c(
  rose   = "#CC6677",
  indigo = "#332288",
  sand   = "#DDCC77",
  green  = "#117733",
  cyan   = "#88CCEE",
  wine   = "#882255",
  teal   = "#44AA99",
  olive  = "#999933",
  purple = "#AA4499"
)

# Neutrals. "Other" and "not resolved" are always grey, so grey never means a
# real category.
.GREY_OTHER  <- "#DDDDDD"
.GREY_MID    <- "#BBBBBB"
.INK         <- "#222222"
.INK_SOFT    <- "#5A5A5A"
.GRID        <- "#EBEBEB"
.PAPER       <- "#FFFFFF"
.WARNING_INK <- "#A12A2A"   # caveat captions only

# The three ways a retroviral integration can be seen in this pipeline. The
# pipeline's tiers and the solo stage's fates are the same three situations under
# two names, so each pair shares one colour and a reader learns them once.
.TIER_COLOUR <- c(
  `ltr-flanked` = .TOL_MUTED[["indigo"]],   # an intact element
  orphan        = .TOL_MUTED[["teal"]],     # coding sequence without an LTR pair
  `solo-ltr`    = .TOL_MUTED[["rose"]]      # a lone LTR, provirus excised
)
.FATE_COLOUR <- c(
  intact_flank       = .TIER_COLOUR[["ltr-flanked"]],
  mono_ltr_at_orphan = .TIER_COLOUR[["orphan"]],
  solo               = .TIER_COLOUR[["solo-ltr"]]
)

# One sequential scale (light to indigo) and one diverging scale (rose through
# white to indigo), both built from palette colours so they sit with the rest.
seq_colours <- function(n) {
  grDevices::colorRampPalette(c("#F3F1F9", .TOL_MUTED[["indigo"]]))(max(n, 1L))
}
div_colours <- function(n) {
  grDevices::colorRampPalette(
    c(.TOL_MUTED[["rose"]], "#F7F7F7", .TOL_MUTED[["indigo"]]))(max(n, 1L))
}

# Ordinal categories on the sequential ramp, most complete darkest.
.STRUCTURE_COLOUR <- stats::setNames(seq_colours(4)[4:2], c("full", "partial", "gene"))

# Viral genera. The six orthoretrovirus genera each get one fixed colour, chosen
# from the palette entries not already carrying a tier meaning. The spumavirus
# genera are rare here (133 of ~42,000 loci) and form their own subfamily, so they
# share one warm stone hue in shades rather than using up distinct colours.
# Anything not resolved to a genus is grey.
.GENUS_COLOUR <- c(
  Betaretrovirus    = .TOL_MUTED[["green"]],
  Gammaretrovirus   = .TOL_MUTED[["wine"]],
  Epsilonretrovirus = .TOL_MUTED[["sand"]],
  Lentivirus        = .TOL_MUTED[["purple"]],
  Alpharetrovirus   = .TOL_MUTED[["olive"]],
  Deltaretrovirus   = .TOL_MUTED[["cyan"]]
)
.SPUMA_SHADES <- c("#7A6A58", "#9C8B77", "#BCAE9C", "#D8CEC1")

# Fixed colours for any set of taxon names, in any plot. Orthoretrovirus genera
# get their own colour; spumavirus genera (names ending "spumavirus") get the
# stone shades in alphabetical order, so the assignment never depends on which
# plot is drawing them; everything else (higher ranks, unresolved, Other) is grey.
taxon_colours <- function(taxa) {
  taxa <- unique(as.character(taxa))
  out <- stats::setNames(rep(.GREY_OTHER, length(taxa)), taxa)
  hit <- taxa %in% names(.GENUS_COLOUR)
  out[hit] <- .GENUS_COLOUR[taxa[hit]]
  spuma <- sort(taxa[grepl("spumavirus$", taxa, ignore.case = TRUE)])
  if (length(spuma)) {
    shades <- grDevices::colorRampPalette(.SPUMA_SHADES)(max(length(spuma), 4L))
    out[spuma] <- shades[seq_along(spuma)]
  }
  higher <- taxa %in% c("Orthoretrovirinae", "Spumaretrovirinae", "Retroviridae")
  out[higher] <- .GREY_MID
  out
}

# Legend order for taxa: most abundant first, and anything that is not a real
# taxon (unassigned, unclassified, Other) last, so the legend reads from what
# matters to what is left over. `weights` are counts aligned with `values`.
taxon_levels <- function(values, weights = NULL) {
  values <- as.character(values)
  if (is.null(weights)) weights <- rep(1, length(values))
  totals <- tapply(weights, values, sum)
  ranked <- names(sort(totals, decreasing = TRUE))
  leftover <- grepl(.NOT_A_TAXON, ranked, ignore.case = TRUE)
  c(ranked[!leftover], ranked[leftover])
}


# Colours for an open-ended categorical variable with no fixed meaning (probes,
# viruses). Palette order, then interpolation if a plot has more levels than the
# palette has colours, which a long-tail collapse should normally prevent.
category_colours <- function(levels) {
  levels <- unique(as.character(levels))
  n <- length(levels)
  cols <- if (n <= length(.TOL_MUTED)) unname(.TOL_MUTED[seq_len(n)])
          else grDevices::colorRampPalette(unname(.TOL_MUTED))(n)
  cols[grepl("^Other", levels)] <- .GREY_OTHER
  stats::setNames(cols, levels)
}

# ---------------------------------------------------------------------------
# Words
# ---------------------------------------------------------------------------

# Readable species names. The config `species:` map holds stem -> display name;
# a genome missing from it falls back to its stem with underscores as spaces, so
# a plot never shows a file name. File names on disk keep their stems.
display_species <- function(values, species_map = NULL) {
  values <- as.character(values)
  vapply(values, function(v) {
    nm <- if (!is.null(species_map)) species_map[[v]] else NULL
    if (is.null(nm) || !nzchar(as.character(nm))) gsub("_", " ", v) else as.character(nm)
  }, character(1), USE.NAMES = FALSE)
}

# Identifiers as plain words. Anything not listed falls back to underscores as
# spaces with the first letter capitalised, so an unlisted value is still
# readable rather than raw.
.LABELS <- c(
  `ltr-flanked`        = "LTR-flanked",
  orphan               = "Orphan",
  `solo-ltr`           = "Solo LTR",
  intact_flank         = "Flank of an intact element",
  mono_ltr_at_orphan   = "MonoLTR at an orphan",
  solo                 = "Solo LTR",
  unassigned_at_genus  = "Unassigned at genus",
  unclassified         = "Unclassified",
  domain_selected      = "Retroviral domain",
  domain_unlisted      = "Other domain",
  non_domain           = "No domain",
  full                 = "Full",
  partial              = "Partial",
  gene                 = "Single gene",
  HC                   = "High confidence",
  LC                   = "Low confidence",
  placement            = "Phylogenetic placement",
  weighted_lca         = "Weighted LCA",
  presence             = "Gene presence",
  no_intact            = "No intact member",
  with_intact          = "Has an intact member",
  no_solo              = "No solos"
)
display_label <- function(values) {
  values <- as.character(values)
  vapply(values, function(v) {
    if (!is.na(v) && v %in% names(.LABELS)) return(.LABELS[[v]])
    if (is.na(v) || !nzchar(v)) return(v)
    words <- gsub("_", " ", v)
    paste0(toupper(substr(words, 1, 1)), substr(words, 2, nchar(words)))
  }, character(1), USE.NAMES = FALSE)
}

# Legend or axis labels with taxa in italics (ICTV and binomial convention) and
# everything that is not a taxon name upright. Returns plotmath expressions,
# which ggplot accepts anywhere it accepts labels.
.NOT_A_TAXON <- "^(Other|Unassigned|Unclassified|unassigned|unclassified|No |Not )"
italic_labels <- function(values) {
  lapply(as.character(values), function(v) {
    if (grepl(.NOT_A_TAXON, v)) bquote(.(v)) else bquote(italic(.(v)))
  })
}

# ---------------------------------------------------------------------------
# Type and theme
# ---------------------------------------------------------------------------

.FONT <- "IBM Plex Sans"

# The single theme. Light and quiet: text does the work, grid lines only guide.
theme_retroseek <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = .FONT) +
    ggplot2::theme(
      text               = ggplot2::element_text(colour = .INK),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 4,
                                            hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = .INK_SOFT,
                                            hjust = 0, lineheight = 1.15,
                                            margin = ggplot2::margin(b = 10)),
      plot.caption  = ggplot2::element_text(size = base_size - 2.5, colour = .INK_SOFT,
                                            hjust = 0, margin = ggplot2::margin(t = 8)),
      axis.title    = ggplot2::element_text(size = base_size - 1, colour = .INK_SOFT),
      axis.text     = ggplot2::element_text(size = base_size - 1.5, colour = .INK),
      panel.grid.major = ggplot2::element_line(colour = .GRID, linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text    = ggplot2::element_text(size = base_size - 1, hjust = 0,
                                            face = "bold", colour = .INK),
      legend.position = "bottom",
      legend.title  = ggplot2::element_text(size = base_size - 1.5, colour = .INK_SOFT),
      legend.text   = ggplot2::element_text(size = base_size - 1.5),
      plot.background = ggplot2::element_rect(fill = .PAPER, colour = NA),
      plot.margin   = ggplot2::margin(14, 30, 12, 14)
    )
}

# The same theme on a blank canvas, for pages with no data axes: trees, sankeys.
# Titles, legend and paper stay the house ones.
theme_retroseek_blank <- function(base_size = 11) {
  theme_retroseek(base_size) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   # The child, not `panel.grid`: theme_retroseek sets the child,
                   # and a set child outranks a blank parent.
                   panel.grid.major = ggplot2::element_blank())
}

# Make every plot built in this R session use the theme and font unless a builder
# overrides a detail. Called once by each generator after sourcing this file.
use_retroseek_style <- function() {
  ggplot2::theme_set(theme_retroseek())
  ggplot2::update_geom_defaults("text", list(family = .FONT, colour = .INK))
  ggplot2::update_geom_defaults("label", list(family = .FONT))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Output: one multi-page PDF per stage
# ---------------------------------------------------------------------------

# A4 landscape. Pages are one size per document (a PDF device cannot vary it),
# so a stage with many species passes a taller height for its whole document.
.PAGE_WIDTH  <- 11.69
.PAGE_HEIGHT <- 8.27

page_height_for <- function(n_species, base = .PAGE_HEIGHT, per_species = 0.18,
                            threshold = 20, cap = 40) {
  min(cap, base + max(0, n_species - threshold) * per_species)
}

# Write plots as pages of one vector PDF. cairo_pdf rather than pdf() because it
# embeds the actual font (pdf() would substitute Helvetica).
save_stage_pdf <- function(pages, path, width = .PAGE_WIDTH, height = .PAGE_HEIGHT) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE,
                       family = .FONT)
  on.exit(grDevices::dev.off())
  for (p in pages) print(p)
  invisible(path)
}

# The first page of every stage PDF: what the stage shows, what its colours mean,
# and what the following pages are. Makes each PDF readable on its own.
#   title        stage name
#   description  a sentence or two, plain words
#   colours      named character vector: label -> colour (may be empty)
#   pages        character vector of page titles, in order
key_page <- function(title, description, colours = character(0), pages = character(0)) {
  # The description sits in the subtitle slot, so the theme lays it out above
  # the panel however long it is; the two blocks below start at the panel top.
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(title = title,
                  subtitle = paste(strwrap(description, 100), collapse = "\n")) +
    ggplot2::theme_void(base_family = .FONT) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0,
                                         margin = ggplot2::margin(b = 12)),
      plot.subtitle = ggplot2::element_text(size = 12, colour = .INK, hjust = 0,
                                            lineheight = 1.2,
                                            margin = ggplot2::margin(b = 28)),
      plot.background = ggplot2::element_rect(fill = .PAPER, colour = NA),
      plot.margin = ggplot2::margin(28, 28, 24, 28)
    )
  heading <- function(x, label) {
    ggplot2::annotate("text", x = x, y = 1, label = label, hjust = 0, vjust = 1,
                      size = 4.4, fontface = "bold", family = .FONT, colour = .INK)
  }
  if (length(colours)) {
    ys <- 1 - 0.07 * seq_along(colours) - 0.02
    p <- p + heading(0, "Colours") +
      ggplot2::annotate("rect", xmin = 0, xmax = 0.025, ymin = ys - 0.022,
                        ymax = ys + 0.022, fill = unname(colours)) +
      ggplot2::annotate("text", x = 0.035, y = ys, label = names(colours), hjust = 0,
                        size = 4, family = .FONT, colour = .INK)
  }
  if (length(pages)) {
    x <- if (length(colours)) 0.5 else 0
    listing <- paste(sprintf("%2d. %s", seq_along(pages) + 1L, pages), collapse = "\n")
    p <- p + heading(x, "Pages") +
      ggplot2::annotate("text", x = x, y = 0.93, label = listing, hjust = 0, vjust = 1,
                        size = 3.8, family = .FONT, colour = .INK_SOFT, lineheight = 1.3)
  }
  p
}
