# tree_axis.R - species on rows, in one canonical order, with the host tree.
#
# Every RetroSeek plot that compares genomes puts them on ROWS (the house style,
# docs/visual_style.md), except the slide-shaped highlights.pdf, which puts them
# on columns: species_columns() and on_columns() at the end of this file. Rows
# read naturally for long italic binomials, never need rotated labels, and let
# the host phylogeny sit directly beside them. The same
# species is therefore in the same place on every page of every stage.
#
# The canonical order:
#   with a host tree configured (input.species_tree)  the tree's tip order;
#   without one                                        the config `species:` order.
#
# The tree itself is ADR-011's coordinate bridge: species_tree_layout.py lays it out and
# writes `species.tree_tips.csv` (tip, x, y) and `species.tree_segments.csv`; this
# file only draws those coordinates. No R tree package is involved.

# No directory at all (an unset option) means no tree, same as an empty file.
.no_dir <- function(dir) is.null(dir) || length(dir) == 0L || !nzchar(dir)

# Read a tree coordinate CSV written by tree_layout.py. A missing or empty file
# means no tree is configured, which is a normal state, so this returns NULL for
# the caller to handle rather than failing.
read_tree_part <- function(dir, name, part) {
  if (.no_dir(dir)) return(NULL)
  f <- file.path(dir, sprintf("%s.tree_%s.csv", name, part))
  if (!file.exists(f)) return(NULL)
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (nrow(df) == 0L) NULL else df
}

# A tree as list(tips, segments), or NULL when it has no tips (none configured).
read_tree <- function(dir, name) {
  tips <- read_tree_part(dir, name, "tips")
  if (is.null(tips)) return(NULL)
  list(tips = tips, segments = read_tree_part(dir, name, "segments"))
}

# The host tree as list(tips, segments), or NULL when none is configured.
read_species_tree <- function(dir) read_tree(dir, "species")

# Every species, bottom row first (the order a flipped discrete axis draws them),
# in the canonical order. ALL canonical species are returned, not only those with
# data: a species with nothing to show keeps its row, empty, so each species sits
# on the same row on every page. Species in the data but outside the canonical
# list (not on the tree, or not in the config) are appended at the top.
species_order <- function(present, tree = NULL, fallback_order = NULL) {
  present <- unique(as.character(present))
  canonical <- if (!is.null(tree)) {
    tree$tips$tip[order(tree$tips$y)]
  } else if (!is.null(fallback_order)) {
    rev(fallback_order)
  } else {
    rev(sort(present))
  }
  c(canonical, setdiff(present, canonical))
}


# The tree column: grey topology, italic species names as the row labels.
# `ylim` is shared with the data panel so the two line up row for row.
tree_column <- function(tips, segs, ylim) {
  ggplot() + {
    if (!is.null(segs)) {
      geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                    xend = .data$xend, yend = .data$yend),
                   colour = .GREY_MID, linewidth = 0.45, lineend = "round")
    }
  } +
    geom_text(data = tips, aes(x = .data$x, y = .data$y, label = .data$tip),
              hjust = -0.06, size = 3.4, colour = .INK, fontface = "italic",
              family = .FONT) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 1.1))) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_void(base_family = .FONT)
}


# Compose a tree column with a data panel under one title. The panel's own title
# and subtitle are lifted onto the composition, so the two stack in the right
# order (a patchwork otherwise draws its annotation above the panel's title).
compose_with_tree <- function(tree, panel, widths = c(1.1, 3)) {
  title <- panel$labels$title
  subtitle <- panel$labels$subtitle
  panel <- panel + labs(title = NULL, subtitle = NULL)
  patchwork::wrap_plots(tree, panel, widths = widths) +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme_retroseek()
    )
}


# Put a plot's species on rows, in canonical order, beside the host tree when one
# is configured.
#
#   panel           a finished ggplot whose discrete `axis` ("x" or "y") is species
#   species         the species (display names) the plot shows
#   tree            read_species_tree(dir), or NULL
#   fallback_order  display names top to bottom (the config order), used when
#                   there is no tree
#   axis            which aesthetic carries species in the panel. "x" panels are
#                   flipped; "y" panels (a heatmap already on rows) are not.
#
# Without a tree the species labels stay on the panel, in italics. With a tree,
# the tree prints them and the panel's own labels are switched off.
species_rows <- function(panel, species, tree = NULL, fallback_order = NULL,
                         axis = "x") {
  levels <- species_order(species, tree, fallback_order)
  species_text <- theme(axis.text.y = element_text(face = "italic", hjust = 1),
                        axis.title.y = element_blank())
  panel <- if (axis == "x") {
    panel + scale_x_discrete(limits = levels) + coord_flip() + species_text +
      # After the flip the value axis is horizontal: undo any rotation the builder
      # applied for species-on-x, and keep grid lines only along the values.
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
            panel.grid.major.y = element_blank())
  } else {
    # Already on rows: only the species axis is touched, the x axis is the
    # builder's own (a heatmap's categories, for instance).
    panel + scale_y_discrete(limits = levels) + species_text
  }
  if (is.null(tree)) return(panel)

  not_on_tree <- setdiff(species, tree$tips$tip)
  if (length(not_on_tree)) {
    # The tree cannot place these, so keep plain labelled rows and say why.
    return(panel + labs(caption = paste("No host tree drawn: not on the tree:",
                                        paste(sort(not_on_tree), collapse = ", "))))
  }
  # Every tip keeps its tree row, so segments line up without re-indexing.
  tips <- tree$tips[order(tree$tips$y), , drop = FALSE]
  ylim <- c(0.4, nrow(tips) + 0.6)
  panel <- panel + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  compose_with_tree(tree_column(tips, tree$segments, ylim), panel)
}


# species_rows() for a panel-registry builder. Every registry passes its builders
# one `ctx` list; the host tree and the config order travel in it as
# ctx$species_tree (read_species_tree()) and ctx$species_order (display names,
# config order). With no ctx, as in a unit test, rows fall back to sorted order.
on_rows <- function(panel, species, ctx = NULL, axis = "x") {
  species_rows(panel, species, tree = ctx$species_tree,
               fallback_order = ctx$species_order, axis = axis)
}


# The `ctx` list a panel registry hands to its builders: the host tree and the
# config order for species rows, plus the two settings some builders need.
#   cfg               the parsed config
#   species_tree_dir  species_tree_layout.py output; "" or missing means no tree
#   tree_dir          the taxon tree coordinates (tree_layout.py), for the
#                     lineage trees
panel_ctx <- function(cfg, species_tree_dir = "", tree_dir = "") {
  species <- cfg$species
  confidence_min <- cfg$classification$confidence_min
  list(
    species_tree   = read_species_tree(species_tree_dir),
    species_order  =
      if (length(species)) display_species(names(species), species) else NULL,
    tree_dir       = tree_dir,
    confidence_min = if (is.null(confidence_min)) 0.5 else confidence_min
  )
}


# Small multiples, one row per species, the first species of the canonical order
# on top: the species-on-rows rule for plots whose x axis is a measurement rather
# than species (histograms). `panel`'s data must carry a `species` column; it is
# turned into a factor in canonical order, which is what facet_grid lays out.
# A species with no data has no row here, since a facet cannot be empty.
species_facets <- function(panel, ctx = NULL, scales = "free_y") {
  top_first <- rev(species_order(panel$data$species, ctx$species_tree,
                                 ctx$species_order))
  panel$data$species <- factor(panel$data$species, levels = top_first)
  panel +
    facet_grid(rows = vars(.data$species), switch = "y", scales = scales) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "italic"),
          strip.placement = "outside")
}


# ---------------------------------------------------------------------------
# Species as COLUMNS under the host tree: the slide-shaped highlights pages
# (docs/visual_style.md), the one exception to species on rows. A wide page holds
# more hosts across than down, and a column keeps its host from top to bottom
# through every stacked panel. The results deck kit uses the same functions.
# ---------------------------------------------------------------------------

# The column order and the tree to draw above it.
#   levels  species, left to right: the tree's tips top row first (else the
#           config order, else alphabetical), then any species outside that
#           order, alphabetically, at the right end
#   tree    the tree, or NULL
#   note    names the species the tree could not place, or NULL
species_columns <- function(species, tree = NULL, fallback_order = NULL) {
  species <- unique(as.character(species))
  # Tips are stored bottom row first (y = 1); the top row becomes the left column.
  canonical <- if (!is.null(tree)) rev(tree$tips$tip[order(tree$tips$y)]) else
    if (!is.null(fallback_order)) fallback_order else sort(species)
  extra <- sort(setdiff(species, canonical))
  note <- if (!is.null(tree) && length(extra)) {
    paste("Not on the host tree, shown at the right:", paste(extra, collapse = ", "))
  }
  list(levels = c(canonical, extra), tree = tree, note = note)
}

# The tree turned on its side: root at the top, each tip directly above its
# column. A tip on row y sits in column n_tips + 1 - y, the same flip that puts
# the top row on the left.
.tree_on_top <- function(tree, n_columns) {
  n_tips <- nrow(tree$tips)
  ggplot(tree$segments) +
    geom_segment(aes(x = n_tips + 1 - .data$y, xend = n_tips + 1 - .data$yend,
                     y = -.data$x, yend = -.data$xend),
                 colour = .GREY_MID, linewidth = 0.45, lineend = "round") +
    scale_x_continuous(limits = c(0.5, n_columns + 0.5), expand = c(0, 0)) +
    theme_void(base_family = .FONT)
}

# Species names under the bottom panel: as large as the column width allows,
# within a readable range, and turned upright past a dozen species.
.column_label_style <- function(n, width_in) {
  size <- max(4.5, min(10, width_in * 72 / n * 0.62))
  if (n <= 12) return(element_text(face = "italic", size = size, colour = .INK))
  element_text(face = "italic", size = size, colour = .INK, angle = 90, hjust = 1,
               vjust = 0.5)
}

# Stack the tree (when there is one) above panels whose x axis is species. Only
# the bottom panel prints the names; one legend is collected below. The tree
# takes the same share of every page, so it does not jump between pages.
#   panels           ggplots with species on a discrete x axis
#   axis             species_columns()
#   title, subtitle  the page's, set over the whole composition
#   width_in         the page width, which sizes the species labels
on_columns <- function(panels, axis, title = NULL, subtitle = NULL,
                       width_in = .SLIDE_WIDTH) {
  tree_share <- 0.16
  n <- length(axis$levels)
  last <- length(panels)
  panels <- lapply(seq_along(panels), function(i) {
    p <- panels[[i]] +
      scale_x_discrete(limits = axis$levels, drop = FALSE,
                       expand = expansion(add = 0.5)) +
      theme(axis.title.x = element_blank())
    if (i < last) return(p + theme(axis.text.x = element_blank()))
    p + theme(axis.text.x = .column_label_style(n, width_in))
  })
  heights <- rep((1 - tree_share) / last, last)
  if (!is.null(axis$tree)) {
    panels <- c(list(.tree_on_top(axis$tree, n)), panels)
    heights <- c(tree_share, heights)
  }
  page <- patchwork::wrap_plots(panels, ncol = 1, heights = heights,
                                guides = "collect") &
    theme(legend.position = "bottom")
  # One annotation for the whole page: its title (which page_titles() reads for
  # the key page), its subtitle, and the note on species the tree lacks.
  page + patchwork::plot_annotation(title = title, subtitle = subtitle,
                                    caption = axis$note, theme = theme_retroseek())
}

# A column panel of loci per host, stacked by `fill`, from pre-counted rows
# (species, the fill column, n): counting once keeps a large study's pages
# small. A total sits above each bar while the hosts are few enough for the
# numbers not to collide.
column_count_panel <- function(counts, fill, colours, labels, n_species, title = NULL) {
  totals <- counts %>% dplyr::count(.data$species, wt = .data$n, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data[[fill]])) +
    # show.legend keeps a key for a class with no loci (a study without orphans).
    geom_col(width = 0.72, position = position_stack(reverse = TRUE),
             show.legend = TRUE) +
    scale_fill_manual(values = colours, breaks = names(colours), labels = labels,
                      limits = names(colours), drop = FALSE, name = NULL) +
    scale_y_count(expand = expansion(mult = c(0, 0.22))) +
    labs(y = "Loci", title = title) +
    theme(panel.grid.major.x = element_blank())
  if (n_species > 30L || !nrow(totals)) return(p)
  p + geom_text(data = totals, aes(x = .data$species, y = .data$n,
                                   label = scales::comma(.data$n)),
                inherit.aes = FALSE, vjust = -0.4, size = 3, colour = .INK_SOFT)
}

# The column panel of a lineage with no loci: it keeps its place and says so.
no_loci_column_panel <- function(n_columns, title = NULL) {
  ggplot() +
    annotate("text", x = (n_columns + 1) / 2, y = 0.5, label = "No loci",
             colour = .INK_SOFT, size = 3.5, family = .FONT) +
    labs(title = title, y = "Loci") +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
}
