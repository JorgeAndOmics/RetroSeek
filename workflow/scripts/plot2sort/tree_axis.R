# tree_axis.R - attach a phylogeny to any panel whose categorical axis is species.
#
# ADR-011 built the hard part already: tree_layout.py parses the Newick, lays it
# out, and writes plain coordinate CSVs (`tip, x, y` and `x, y, xend, yend`). No R
# tree package is needed or wanted - the ADR chose this coordinate bridge over
# ape/ggtree deliberately, to keep the environment lean.
#
# What lives here is the other half: drawing that coordinate table as a tree
# column, and composing it with a data panel so the rows line up. Two panels in
# taxonomy_plot_generator.R each carried their own copy of the drawing code; they
# now share this one, and any other panel can gain a tree by asking for it.
#
# Why the tree variants are ADDITIONAL files rather than a change to the existing
# panels: a tree axis puts tips on ROWS, but almost every species panel in this
# project puts species on the x axis. Rewriting them all to be horizontal would
# change every figure a reader is used to, to no benefit where the tree adds
# nothing. So `<panel>_tree.png` sits beside `<panel>.png` and the original is
# untouched.

# Read a tree coordinate CSV written by tree_layout.py. A missing or empty file
# means no tree is configured, which is a normal state (input.species_tree
# defaults to empty), so it returns NULL for the caller to handle rather than
# failing.
read_tree_part <- function(dir, name, part) {
  f <- file.path(dir, sprintf("%s.tree_%s.csv", name, part))
  if (!file.exists(f)) return(NULL)
  df <- suppressWarnings(readr::read_csv(f, show_col_types = FALSE))
  if (nrow(df) == 0L) NULL else df
}


# The tree column itself: segments for the topology, tip labels on the right.
#
# `ylim` is shared with the data panel so the two line up row for row; it is the
# caller's job to pass the same value to both. theme_void() is deliberate - the
# tree supplies the tip labels, so the data panel must NOT also print them.
tree_column <- function(tips, segs, ylim) {
  ggplot() +
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
}


# Compose a tree column with a data panel, titled as one figure.
#
# add_titles() styles a plain ggplot; a patchwork needs plot_annotation instead,
# which is why this does not simply call the former.
compose_with_tree <- function(tree, panel, title, subtitle,
                              widths = c(1.2, 3)) {
  patchwork::wrap_plots(tree, panel, widths = widths) +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.subtitle   = element_text(hjust = 0.5, size = 11),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )
}


# Turn a finished species-on-x panel into a tree-axis variant.
#
# The panel is not rebuilt. Three things are done to it:
#   1. its discrete x scale is re-limited to the tree's tip order, which both
#      orders the categories and drops any species the tree does not carry;
#   2. coord_flip() turns that axis vertical, so categories become rows;
#   3. its own category labels are switched off, because the tree prints them.
#
# Returns NULL when there is no tree or no overlap, so the caller can skip the
# variant rather than emit a misleading empty figure.
#
# `dropped` species are named in the subtitle rather than silently vanishing:
# tree_layout.py already warns about them in the log, but a reader looking at the
# figure would otherwise never know a genome was missing from it.
attach_tree_axis <- function(panel, tips, segs, present_species,
                             title, subtitle) {
  if (is.null(tips) || is.null(panel)) return(NULL)
  on_tree <- intersect(as.character(present_species), tips$tip)
  if (!length(on_tree)) return(NULL)

  ordered_tips <- tips[order(tips$y), , drop = FALSE]
  dropped <- setdiff(as.character(present_species), tips$tip)
  if (length(dropped)) {
    subtitle <- sprintf("%s\nNot on the tree, so not shown: %s",
                        subtitle, paste(sort(dropped), collapse = ", "))
  }

  ylim <- c(0.4, nrow(ordered_tips) + 0.6)
  flipped <- panel +
    scale_x_discrete(limits = ordered_tips$tip) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  compose_with_tree(tree_column(ordered_tips, segs, ylim), flipped,
                    title, subtitle)
}
