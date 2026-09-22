# tree_axis.R - species on rows, in one canonical order, with the host tree.
#
# Every RetroSeek plot that compares genomes puts them on ROWS (the house style,
# docs/visual_style.md). Rows read naturally for long italic binomials, never need
# rotated labels, and let the host phylogeny sit directly beside them. The same
# species is therefore in the same place on every page of every stage.
#
# The canonical order:
#   with a host tree configured (input.species_tree)  the tree's tip order;
#   without one                                        the config `species:` order.
#
# The tree itself is ADR-011's coordinate bridge: tree_layout.py lays it out and
# writes `species.tree_tips.csv` (tip, x, y) and `species.tree_segments.csv`; this
# file only draws those coordinates. No R tree package is involved.

# Read a tree coordinate CSV written by tree_layout.py. A missing or empty file
# means no tree is configured, which is a normal state, so this returns NULL for
# the caller to handle rather than failing.
read_tree_part <- function(dir, name, part) {
  # No directory at all (an unset option) means no tree, same as an empty file.
  if (is.null(dir) || length(dir) == 0L || !nzchar(dir)) return(NULL)
  f <- file.path(dir, sprintf("%s.tree_%s.csv", name, part))
  if (!file.exists(f)) return(NULL)
  df <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE))
  if (nrow(df) == 0L) NULL else df
}

# The host tree as list(tips, segments), or NULL when none is configured.
read_species_tree <- function(dir) {
  if (is.null(dir) || !nzchar(dir)) return(NULL)
  tips <- read_tree_part(dir, "species", "tips")
  if (is.null(tips)) return(NULL)
  list(tips = tips, segments = read_tree_part(dir, "species", "segments"))
}

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
  ggplot() +
    { if (!is.null(segs)) {
        geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend),
                     colour = .GREY_MID, linewidth = 0.45, lineend = "round")
      } } +
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
      theme = theme_retroseek() + theme(plot.margin = margin(14, 18, 12, 14))
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
species_rows <- function(panel, species, tree = NULL, fallback_order = NULL, axis = "x") {
  levels <- species_order(species, tree, fallback_order)
  panel <- panel +
    { if (axis == "x") scale_x_discrete(limits = levels) else scale_y_discrete(limits = levels) } +
    { if (axis == "x") coord_flip() } +
    theme(
      # After the flip the value axis is horizontal: undo any label rotation a
      # builder applied for species-on-x, and keep grid lines only along values.
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.text.y = element_text(face = "italic", hjust = 1),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank()
    )
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
