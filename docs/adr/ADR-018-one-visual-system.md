# ADR-018: One visual system for every figure, one PDF per stage

- **Status**: Accepted
- **Date**: 2026-09-22
- **Deciders**: Jorge González García
- **Builds on**: [ADR-011](ADR-011-phylogeny-aware-plots-and-rank-segmentation.md) (the tree coordinate bridge)

## Context

Nothing governed RetroSeek's figures, and they looked like it. A survey of every
plotting script found:

- **About thirty colour schemes**: seven ggsci families, five hand-rolled
  gradients and fifteen hard-coded sets. One colour stood for many concepts (one
  orange meant eight things), and one concept wore many colours (the orphan tier
  was green, grey or gold depending on the page). The Futurama palette in use
  failed colour-blind simulation (minimum tritan distance 4.7).
- **Three stock themes, no set font**, bold body text, centred titles, and four
  competing species orderings.
- **Dashes as punctuation** in about twenty titles and seven "->" arrows, and file
  stems leaking into titles and legends.
- **Around a hundred PNGs per run** spread over directories, with `_tree.png`
  duplicates for the phylogeny-attached variants.

## Decision

Every figure goes through one style module, `workflow/scripts/plot2sort/style.R`,
and one species-axis module, `workflow/scripts/plot2sort/tree_axis.R`:

- **Colour**: Paul Tol's muted palette (colour-blind safe, minimum distance 15.8)
  with fixed semantic maps: the three evidence tiers and the three solo-LTR fates
  share three colours; each orthoretrovirus genus has one colour and the
  spumaviruses fixed stone shades; ERV classes wear their genus colour; viruses
  are shades of their genus; ordinal categories and continuous values use one
  light-to-indigo ramp. Plotting code never writes a colour.
- **Type**: IBM Plex Sans, shipped in the conda environment
  (`font-ttf-ibm-plex-sans`) and embedded in every PDF.
- **Words**: readable species names from the config, one dictionary for
  identifiers, left-aligned sentence-case titles, no dash as punctuation.
- **Species on rows**, in one canonical order (host tree, else config), every
  configured species keeping its row, the host tree beside them. A new early rule,
  `species_tree_layout`, lays the tree out from the config alone so the first
  plotting stage can use it.
- **Output**: one multi-page A4-landscape PDF per stage, or per genome where the
  pages are per genome, opening with a key page. The README demo figures stay PNG.

`tests/unit/test_visual_style.py` and `workflow/tests/testthat/test-style.R` guard
the rules. The system is described for people in [visual_style.md](../visual_style.md).

## Consequences

- Positive:
  - A reader who has learnt one stage can read any other: a genus, a tier or a
    species looks the same and sits in the same place everywhere.
  - Colour-blind readers can use every figure.
  - One PDF per stage replaces around a hundred PNGs; the key page makes each
    readable on its own, and vector output scales for publication.
  - New figures inherit the style by construction, and the guards catch drift.
- Negative:
  - A new dependency, the IBM Plex Sans font package.
  - The per-plot PNGs are gone, so a single figure for a slide is now a page to
    export rather than a file to copy.
  - All pages of a stage share one page height, the tallest any page needs.
  - Sets with no fixed meaning (13 probes, LTR families) borrow the palette's
    hues, which can coincide with a genus colour on another page; their legend
    disambiguates. Nine colour-blind-safe colours cannot cover every set.
- Neutral:
  - The config keys that sized PNG canvases (`plots.dpi`, `width`, `height`,
    `max_dim`) are removed; `plots.per_stratum` now grows the page height past
    20 genomes.

## Alternatives considered

- **Keep PNGs and only unify the palette**: fixes colour but not the scattered
  output, the missing font or the rotated species labels.
- **A PDF per plot**: keeps the file sprawl without the key page that explains
  the colours once.
- **Okabe-Ito or viridis as the categorical palette**: Okabe-Ito has eight
  colours, one short of what the tiers and genera need together; viridis is
  sequential, not categorical.
- **An R tree package (ggtree) for species rows**: rejected by ADR-011's choice
  of a coordinate bridge; the existing Python layout already feeds every tree.

## Revisit trigger

- A study whose genera or categories outgrow the nine palette colours in one figure.
- A journal requiring a figure format other than vector PDF.
- The circle-plot stage being repaired: it is broken and not yet in the house style.

## References

- [visual_style.md](../visual_style.md)
- Paul Tol, "Colour schemes and templates", SRON technical note.
- [ADR-011](ADR-011-phylogeny-aware-plots-and-rank-segmentation.md): the tree coordinate bridge.
