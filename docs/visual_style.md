# Visual style

Every RetroSeek figure follows one visual system, so a reader who has learnt one
stage's pages can read any other. This page explains the system and the reasons
for it. Its executable form is `workflow/scripts/plot2sort/style.R` (colour,
words, type, output) and `workflow/scripts/plot2sort/tree_axis.R` (species on
rows); plotting code asks those two files instead of deciding for itself.

## Four rules

1. **One colour, one meaning.** A recurring concept (an evidence tier, a solo-LTR
   fate, a viral genus) has exactly one colour in every figure, and no two fixed
   meanings share a colour. Plotting code never writes a colour value. Sets
   with no fixed meaning (probes, methods, LTR families) borrow palette hues,
   always with a legend on the page: nine colour-blind-safe colours cannot give
   every open-ended set its own.
2. **Colour-blind safe.** The categorical palette is Paul Tol's "muted" set. Every
   pair stays distinguishable under deuteranopia, protanopia and tritanopia
   (minimum CIELAB distance 15.8; the Futurama palette it replaced fell to 4.7).
3. **Ordered things get ordered colours.** Ordinal categories take a light-to-dark
   ramp, so "more" reads as darker before the legend is consulted.
4. **Readable text.** One typeface, left-aligned sentence-case titles, readable
   species names, plain words for identifiers, and no dash used as punctuation.

## Colour

### The palette

| Name | Hex | Fixed meaning |
|---|---|---|
| indigo | `#332288` | LTR-flanked tier, flank of an intact element; the dark end of every ramp |
| teal | `#44AA99` | orphan tier, monoLTR at an orphan |
| rose | `#CC6677` | solo LTR |
| green | `#117733` | *Betaretrovirus* (and ERV Class II) |
| wine | `#882255` | *Gammaretrovirus* (and ERV Class I) |
| sand | `#DDCC77` | *Epsilonretrovirus* |
| purple | `#AA4499` | *Lentivirus* |
| olive | `#999933` | *Alpharetrovirus* |
| cyan | `#88CCEE` | *Deltaretrovirus* |

Neutrals: light grey `#DDDDDD` for "Other" and anything not resolved, mid grey
`#BBBBBB` for higher ranks and for the unselected half of a two-way split, ink
`#222222` and soft ink `#5A5A5A` for text.

### Semantic maps

- **Tiers and fates are the same three situations**, so they share colours: an
  LTR-flanked element and the flank of an intact element are indigo; an orphan
  and a monoLTR at an orphan are teal; a solo LTR is rose.
- **Viral genera** keep one colour each everywhere. The spumavirus genera are rare
  in mammals and form their own subfamily, so they share a stone hue in fixed
  shades (a genus keeps its shade whichever others share the page). Higher ranks
  (*Orthoretrovirinae*, *Retroviridae*) are mid grey; unassigned and unclassified
  loci are light grey.
- **ERV classes** wear the colour of the genus they resemble: Class I
  gamma-like, Class II beta-like, Class III spuma-like.
- **Viruses** (the probeset's virus names) are shades of their genus colour, the
  most abundant virus of a genus taking the genus colour itself, so a virus reads
  as a member of its genus. The shades are fixed once per stage.
- **Ordinal categories** take the sequential ramp, strongest evidence darkest:
  structural class (full, partial, single gene), call confidence, resolved rank,
  domain support, and the call method (placement dark, weighted LCA light). A
  level meaning "no evidence" is grey, not the lightest step. A page drawn with a
  fixed class list keeps a legend key for every class, present or not.
- **Continuous values** (counts, shares, confidence) take the same sequential
  ramp, from `#F3F1F9` to indigo; numbers printed on a tile switch to white on
  the dark end.
- **A single data series** with no category meaning is a mid step of the ramp.
  In a two-way split, the subset of interest takes that colour and the rest grey.
- **Open-ended categories** with no fixed meaning (probes, LTR families) take
  the palette in an order that leaves the tier colours for last.
  Past nine levels the nine keep their colours and extra levels take lighter
  tints. Within a stage, a probe keeps one colour on every page. These borrowed
  hues can coincide with a genus colour on another page, so they always come
  with a legend.
- **Flows** in alluvial plots are coloured by the axis that carries a meaning
  (probe or lineage), never by host: hosts have no colour of their own.

## Type

- **IBM Plex Sans** everywhere, installed into the conda environment
  (`font-ttf-ibm-plex-sans`), so every machine renders the same figure and every
  PDF embeds the font.
- Titles bold, left-aligned, sentence case. Subtitles regular and grey: what the
  page shows and how to read it, in one or two sentences.
- **Italics** for species binomials, viral taxa and gene symbols (ICTV and
  nomenclature convention). "Other", "Unassigned" and similar stay upright, and
  so does a virus name ("Murine leukemia virus") in a mixed list.
- **Counts are whole numbers**: every count axis uses `scale_y_count()`, with
  thousands separators and no fractional ticks.

## Words

- **Readable species names.** A plot shows the config `species:` value for a
  genome, falling back to the file stem with underscores as spaces. File names on
  disk keep their stems.
- **Plain words for identifiers.** One dictionary in `style.R` maps identifiers to
  words (`ltr-flanked` to "LTR-flanked", `mono_ltr_at_orphan` to "MonoLTR at an
  orphan", `lca` to "Weighted LCA"); anything unlisted falls back to underscores
  as spaces, capitalised.
- **What a page is about leads its subtitle** ("*Homo sapiens*. Every LTR began
  as ..."), so titles stay identical across genomes and subsets. A species leads
  in italics (`add_titles(species = )`); a probe set or segment leads upright
  (`subset_label`). plotmath, which sets the italics, cannot break lines, so a
  subtitle with a line break keeps an upright lead.
- **No dash as punctuation.** Never a spaced hyphen, an arrow, an em dash or an
  en dash between words: use a colon, a comma, "to", or a new sentence. Number
  ranges read "2 to 5". Hyphens inside compound words ("LTR-flanked") are words,
  not punctuation.

## Species on rows

Every figure that compares genomes puts them on **rows**, never on the x axis.
Rows read naturally for long italic binomials, need no rotated labels, and let
the host phylogeny sit directly beside them.

- **One canonical order**: the host tree's tip order when a tree is configured
  (`input.species_tree`), otherwise the config `species:` order. The early rule
  `species_tree_layout` lays the tree out from the config alone, so every stage,
  including the first plotting stage, reads the same order.
- **Every configured species keeps its row**, empty if it has no data, so a
  species sits on the same row on every page.
- **The tree is drawn beside the rows** when configured, and prints the names.
  If any species is missing from the tree, the rows keep their labels and the
  page says which species the tree could not place.
- Plots whose x axis is a measurement (histograms, funnels) become **small
  multiples, one per species**, in the same order.
- **The one exception is `highlights.pdf`**: a few 16:9 pages of headline
  findings, with hosts as **columns** under the host tree
  (`species_columns()` and `on_columns()` in `tree_axis.R`), since a wide page
  holds more hosts across than down. Anything beyond the headlines belongs in
  its stage PDF.

## Output

- **One multi-page vector PDF per stage**, or per genome where the pages are per
  genome. `results/plots/` is one level deep: each stage's PDF sits directly in
  it, and only the stages that write one file per genome, per segment or per
  tree get a folder.

  | Stage | File (under `results/plots/`) |
  |---|---|
  | Homology | `homology.pdf` |
  | Integration | `integration.pdf` |
  | Taxonomy | `taxonomy.pdf` |
  | Structure | `structure.pdf` |
  | Loss | `loss.pdf` |
  | Segments | `segments.pdf` (all segments), `segments/<segment>.pdf` |
  | Solo LTRs | `solo_ltr.pdf` (all genomes), `solo_ltr/<genome>.solo_ltr.pdf` |
  | Hotspots | `hotspots/<genome>.hotspots.pdf` |
  | Highlights | `highlights.pdf` (16:9, hosts as columns) |
  | Placement | `placement/<genome>.<tier>.<gene>.tree.svg` |

- **A4 landscape.** A page grows taller past 20 species, by `plots.per_stratum`
  inches per species, so rows keep their room in a large study.
- **Every PDF opens with a key page**: what the stage shows, what its colours
  mean (taxa in italics), and the list of pages, read from the pages themselves,
  in columns of at most 14 so any number fits. A title naming a taxon or species
  sets it in bold italics (`key_title()`).
- **No empty pages.** A page whose tier has no loci is left out, and the key page
  says so; a segment with fewer than 5 loci gets one page listing them instead
  of the full panel.
- The README's demo figures are the one exception: PNG, because GitHub renders
  images inline.
- The placement heat-trees are drawn by gappa, given house colours: branches
  with no placement mass in mid grey, so the tree keeps its shape, darkening to
  indigo for the heaviest.

## Adding a plot

1. Build it with `theme_retroseek()` (set for the session by
   `use_retroseek_style()`); take colours from `style.R` (`.TIER_COLOUR`,
   `scale_fill_taxon()`, `.STRUCTURE_COLOUR`, `scale_fill_ramp()`,
   `category_colours()`), never a hex value; a count axis through
   `scale_y_count()`.
2. Title and subtitle through `add_titles()` (a genome's page passes
   `species =`); labels through `display_label()`, species through
   `display_species()`, taxa through `taxon_labels()`.
3. If genomes are an axis, map them to x and pass the plot to `on_rows()` (or
   `species_facets()` for a histogram) with the panel's `ctx`.
4. Add the page to its stage's page list or panel registry; the key page picks
   it up.
5. `make check` runs the guards: `tests/unit/test_visual_style.py` (no stray
   colour, dash, hyphen range or foreign palette) and `workflow/tests/testthat/test-style.R`
   (colour-blind distances, fixed genus colours, ordinal ramps).
