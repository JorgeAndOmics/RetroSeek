# -----------------------------------------------------------------------------
# range_analysis / erv_assembly.R
# -----------------------------------------------------------------------------
# Assemble "ERV-like" candidates from the UNREDUCED valid tier (`valid_hits`).
#
# Each valid hit is one main-or-accessory probe locus. A real ERV is a chain of
# co-located, co-oriented main genes (e.g. GAG -> POL -> ENV). This module
# chains adjacent MAIN-probe loci of the same group (virus / label / none) on
# the same strand into composite candidates: full ERVs when every main probe is
# present, or the longest recoverable fragment otherwise (e.g. GAG+POL with ENV
# absent). The result is purely additive — isolated single-gene loci are never
# emitted here; they remain in the valid output unchanged.
#
# Why manual gap-chaining rather than plyranges::reduce_ranges_directed: the
# reduce idiom is built for *metric aggregation across overlapping ranges* and
# discards per-contributor identity. We need the opposite — to keep each
# constituent locus as a first-class child feature (for GFF3 Parent linkage),
# to know each probe's genomic order (for the canonical-order check), and to
# enforce a ">=2 distinct main probes" predicate during clustering. A
# single-linkage walk over a per-group, start-sorted set is shorter, exact, and
# deterministic for those goals.
#
# Public:
#   assemble_erv_like(valid_hits, main_probes, group_by, max_join_distance,
#                     require_canonical_order, completeness_threshold, agg)
#       -> list(parents = <GRanges>, children = <GRanges>,
#               dropped_noncanonical = <integer>)
#   build_erv_like_df(parents) -> per-genome tibble (feeds the erv-like plots)

suppressMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

# Null-coalescing helper (io.R defines its own copy locally inside
# read_pipeline_options; redefine here so this module is testable in isolation).
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x


# Normalise any aggregation-strategy mcol into a list-of-tokens (one character
# vector per range). Handles CharacterList (`list` strategy), separator-joined
# strings (`concatenate`), and plain single values (`best`/`first`/...). Empty
# and NA tokens are dropped.
.col_as_list <- function(col, sep) {
  if (is.null(col)) return(list())
  clean <- function(x) {
    x <- trimws(as.character(x))
    x[!is.na(x) & nzchar(x)]
  }
  if (methods::is(col, "List") || is.list(col)) {
    lapply(col, clean)
  } else {
    lapply(strsplit(as.character(col), sep, fixed = TRUE), clean)
  }
}


# Collapse a token list back to a single separator-joined string per range
# (sorted-unique). Used for the parent's virus/label summary.
.tokens_to_string <- function(token_list, sep) {
  vapply(token_list, function(x) {
    if (length(x) == 0L) return(NA_character_)
    paste(sort(unique(x)), collapse = sep)
  }, character(1))
}


# Single-linkage gap chaining over start-sorted intervals. `starts`/`ends` must
# already be sorted by `starts` ascending. Two consecutive intervals join when
# the gap to the running cluster end is <= `max_gap` (inclusive; overlaps give a
# negative gap and always join). Returns an integer cluster id per interval.
.chain_by_gap <- function(starts, ends, max_gap) {
  n <- length(starts)
  if (n == 0L) return(integer(0))
  cid <- integer(n)
  cid[1] <- 1L
  running_end <- ends[1]
  for (i in seq_len(n)[-1]) {
    gap <- starts[i] - running_end - 1L
    cid[i] <- if (gap <= max_gap) cid[i - 1L] else cid[i - 1L] + 1L
    running_end <- max(running_end, ends[i])
  }
  cid
}


# Main probes in 5'->3' genomic order, de-duplicated (first occurrence kept).
# On the '-' strand the gene's 5'->3' direction runs against genomic
# coordinates, so we order by descending start.
.observed_order <- function(probes, starts, strand) {
  if (length(probes) == 0L) return(character(0))
  ord <- if (identical(strand, "-")) order(-starts) else order(starts)
  p <- probes[ord]
  p[!duplicated(p)]
}


# Is the observed 5'->3' probe order consistent with the canonical order given
# by `main_probes`? Map each observed probe to its main_probes index; canonical
# means strictly increasing indices. `.observed_order` already made the sequence
# 5'->3' for +/-, so both check "increasing". For '*' (strand unknown) accept
# either direction. Fewer than two probes is trivially ordered.
.check_canonical <- function(observed, main_probes, strand) {
  idx <- match(observed, main_probes)
  if (anyNA(idx)) return(NA)
  if (length(idx) < 2L) return(TRUE)
  increasing <- all(diff(idx) > 0L)
  if (identical(strand, "*")) increasing || all(diff(idx) < 0L) else increasing
}


# Typed-empty results so the orchestrator's exporters / Snakemake outputs are
# always satisfied even with no candidates.
.empty_erv_parents <- function() {
  gr <- GenomicRanges::GRanges()
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    ID = character(0), type = character(0),
    probes_present = character(0), n_main_present = integer(0),
    n_main_expected = integer(0), completeness_fraction = numeric(0),
    is_full = logical(0), observed_order = character(0),
    is_canonical = logical(0), virus = character(0), label = character(0),
    species = character(0), n_loci = integer(0),
    max_bitscore = numeric(0), max_identity = numeric(0),
    min_evalue = numeric(0), query_coverage = numeric(0)
  )
  gr
}

.empty_erv_children <- function() {
  gr <- GenomicRanges::GRanges()
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    ID = character(0), type = character(0), Parent = character(0),
    source_id = character(0), probe = character(0),
    virus = character(0), label = character(0),
    max_bitscore = numeric(0), max_identity = numeric(0),
    min_evalue = numeric(0), query_coverage = numeric(0)
  )
  gr
}

.empty_erv_result <- function() {
  list(parents = .empty_erv_parents(), children = .empty_erv_children(),
       dropped_noncanonical = 0L)
}


# Assemble ERV-like candidates. See module header for the contract.
assemble_erv_like <- function(valid_hits, main_probes, group_by,
                              max_join_distance, require_canonical_order,
                              completeness_threshold, agg) {
  sep <- agg$agg_concat_separator %||% "; "
  n_expected <- length(main_probes)
  if (length(valid_hits) < 2L || n_expected == 0L) return(.empty_erv_result())

  mc <- S4Vectors::mcols(valid_hits)

  # 1. Keep only loci carrying at least one MAIN probe; record the matched set.
  probe_tokens <- .col_as_list(mc$probe, sep)
  matched_main <- lapply(probe_tokens, function(p) intersect(p, main_probes))
  keep <- lengths(matched_main) > 0L
  if (sum(keep) < 2L) return(.empty_erv_result())

  hits         <- valid_hits[keep]
  matched_main <- matched_main[keep]
  mc           <- S4Vectors::mcols(hits)

  # Per-locus scalars (valid_hits carries the first-reduction composites).
  num <- function(x) if (is.null(x)) rep(NA_real_, length(hits)) else as.numeric(x)
  loc <- list(
    seqnames       = as.character(GenomicRanges::seqnames(hits)),
    start          = BiocGenerics::start(hits),
    end            = BiocGenerics::end(hits),
    strand         = as.character(BiocGenerics::strand(hits)),
    source_id      = if (!is.null(mc$ID)) as.character(mc$ID) else paste0("hit_", seq_along(hits)),
    probe_str      = .tokens_to_string(probe_tokens[keep], sep),
    virus_str      = .tokens_to_string(.col_as_list(mc$virus, sep), sep),
    label_str      = .tokens_to_string(.col_as_list(mc$label, sep), sep),
    species        = if (!is.null(mc$species)) as.character(mc$species) else NA_character_,
    max_bitscore   = num(mc$max_bitscore),
    max_identity   = num(mc$max_identity),
    min_evalue     = num(mc$min_evalue),
    query_coverage = num(mc$query_coverage)
  )

  # 2. Explode each locus across the group-by tokens so co-located but
  #    different-group probes seed SEPARATE candidate streams ("retain both").
  group_tokens <- switch(group_by,
    virus = .col_as_list(mc$virus, sep),
    label = .col_as_list(mc$label, sep),
    as.list(rep("", length(hits)))            # "none": one empty group
  )
  # A locus with no group token (e.g. NA virus under group_by=virus) can't be
  # placed in a group, so it cannot chain — drop it.
  group_tokens <- lapply(group_tokens, function(g) if (length(g) == 0L) NA_character_ else g)
  times <- lengths(group_tokens)
  ridx  <- rep(seq_along(hits), times)           # instance -> locus row
  inst <- data.frame(
    locus  = ridx,
    group  = unlist(group_tokens, use.names = FALSE),
    seqkey = paste(loc$seqnames[ridx], loc$strand[ridx],
                   unlist(group_tokens, use.names = FALSE), sep = "\r"),
    start  = loc$start[ridx],
    end    = loc$end[ridx],
    strand = loc$strand[ridx],
    stringsAsFactors = FALSE
  )
  inst <- inst[!is.na(inst$group), , drop = FALSE]
  if (nrow(inst) < 2L) return(.empty_erv_result())

  # 3. Chain within each (seqnames, strand, group) key, start-sorted. Offset the
  #    per-key cluster ids into a single global namespace.
  inst <- inst[order(inst$seqkey, inst$start, inst$end), , drop = FALSE]
  inst$cid <- NA_integer_
  offset <- 0L
  for (k in unique(inst$seqkey)) {
    rows <- which(inst$seqkey == k)
    local_cid <- .chain_by_gap(inst$start[rows], inst$end[rows], max_join_distance)
    inst$cid[rows] <- local_cid + offset
    offset <- offset + max(local_cid)
  }

  # 4 + 5. One candidate per cluster with >=2 DISTINCT main probes; compute
  #        descriptors (always — even when later filtered by canonical order).
  clusters <- lapply(split(seq_len(nrow(inst)), inst$cid), function(rows) {
    lrows  <- inst$locus[rows]                     # locus indices in `hits`
    strand <- inst$strand[rows][1]
    present <- unique(unlist(matched_main[lrows]))
    if (length(present) < 2L) return(NULL)
    present <- present[order(match(present, main_probes))]   # canonical order

    # Probe-per-position for observed order: expand each locus by its main probes.
    pm     <- matched_main[lrows]
    probes <- unlist(pm)
    starts <- rep(inst$start[rows], lengths(pm))
    observed <- .observed_order(probes, starts, strand)

    frac <- length(present) / n_expected
    # The grouped attribute carries the group value (this candidate *is* the HIV
    # candidate); the ungrouped attribute carries the union of constituents.
    union_tok <- function(strs) paste(sort(unique(unlist(
      strsplit(strs, sep, fixed = TRUE)))), collapse = sep)
    gval <- inst$group[rows][1]
    list(
      seqnames = loc$seqnames[lrows][1],
      start    = min(inst$start[rows]),
      end      = max(inst$end[rows]),
      strand   = strand,
      inst_rows = rows,
      probes_present = paste(present, collapse = sep),
      n_main_present = length(present),
      completeness_fraction = frac,
      is_full = frac >= completeness_threshold,
      observed_order = paste(observed, collapse = sep),
      is_canonical = .check_canonical(observed, main_probes, strand),
      virus = if (identical(group_by, "virus")) gval else union_tok(loc$virus_str[lrows]),
      label = if (identical(group_by, "label")) gval else union_tok(loc$label_str[lrows]),
      species = loc$species[lrows][1],
      n_loci = length(rows),
      max_bitscore   = max(loc$max_bitscore[lrows],   na.rm = TRUE),
      max_identity   = max(loc$max_identity[lrows],   na.rm = TRUE),
      min_evalue     = min(loc$min_evalue[lrows],     na.rm = TRUE),
      query_coverage = min(1, sum(loc$query_coverage[lrows], na.rm = TRUE))
    )
  })
  clusters <- Filter(Negate(is.null), clusters)
  if (length(clusters) == 0L) return(.empty_erv_result())

  # Canonical-order filter: count drops for the canonical-vs-rearranged plot,
  # then optionally discard non-canonical candidates from the output.
  is_canon <- vapply(clusters, function(c) isTRUE(c$is_canonical), logical(1))
  dropped_noncanonical <- if (require_canonical_order) sum(!is_canon) else 0L
  if (require_canonical_order) clusters <- clusters[is_canon]
  if (length(clusters) == 0L) {
    res <- .empty_erv_result(); res$dropped_noncanonical <- dropped_noncanonical
    return(res)
  }

  # Deterministic ID assignment: sort candidates by genomic position.
  ord <- order(
    vapply(clusters, `[[`, character(1), "seqnames"),
    vapply(clusters, `[[`, numeric(1),   "start"),
    vapply(clusters, `[[`, numeric(1),   "end")
  )
  clusters <- unname(clusters[ord])   # drop split()'s cid names so mcols stay unnamed
  ids <- sprintf("ERV_like_%d", seq_along(clusters))

  # 5. Parent GRanges.
  parents <- GenomicRanges::GRanges(
    seqnames = vapply(clusters, `[[`, character(1), "seqnames"),
    ranges   = IRanges::IRanges(
      start = vapply(clusters, `[[`, numeric(1), "start"),
      end   = vapply(clusters, `[[`, numeric(1), "end")
    ),
    strand   = vapply(clusters, `[[`, character(1), "strand")
  )
  S4Vectors::mcols(parents) <- S4Vectors::DataFrame(
    ID = ids, type = "erv_like",
    probes_present  = vapply(clusters, `[[`, character(1), "probes_present"),
    n_main_present  = as.integer(vapply(clusters, `[[`, numeric(1), "n_main_present")),
    n_main_expected = rep(as.integer(n_expected), length(clusters)),
    completeness_fraction = vapply(clusters, `[[`, numeric(1), "completeness_fraction"),
    is_full         = vapply(clusters, `[[`, logical(1), "is_full"),
    observed_order  = vapply(clusters, `[[`, character(1), "observed_order"),
    is_canonical    = vapply(clusters, function(c) isTRUE(c$is_canonical), logical(1)),
    virus = vapply(clusters, `[[`, character(1), "virus"),
    label = vapply(clusters, `[[`, character(1), "label"),
    species = vapply(clusters, `[[`, character(1), "species"),
    n_loci = as.integer(vapply(clusters, `[[`, numeric(1), "n_loci")),
    max_bitscore   = vapply(clusters, `[[`, numeric(1), "max_bitscore"),
    max_identity   = vapply(clusters, `[[`, numeric(1), "max_identity"),
    min_evalue     = vapply(clusters, `[[`, numeric(1), "min_evalue"),
    query_coverage = vapply(clusters, `[[`, numeric(1), "query_coverage")
  )

  # 6. Child member GRanges — the actual retained loci, parented to the
  #    candidate. Built in one shot (rather than per-cluster c()) so combining
  #    single-chromosome blocks never triggers a Seqinfo-merge warning.
  member_rows <- lapply(clusters, `[[`, "inst_rows")
  n_members   <- lengths(member_rows)
  flat_rows   <- unlist(member_rows, use.names = FALSE)
  flat_parent <- rep(ids, n_members)
  flat_mindex <- unlist(lapply(n_members, seq_len), use.names = FALSE)
  lrows       <- inst$locus[flat_rows]
  children <- GenomicRanges::GRanges(
    seqnames = loc$seqnames[lrows],
    ranges   = IRanges::IRanges(start = inst$start[flat_rows], end = inst$end[flat_rows]),
    strand   = inst$strand[flat_rows]
  )
  S4Vectors::mcols(children) <- S4Vectors::DataFrame(
    ID = sprintf("%s.m%d", flat_parent, flat_mindex),
    type = "erv_like_member", Parent = flat_parent,
    source_id = loc$source_id[lrows], probe = loc$probe_str[lrows],
    virus = loc$virus_str[lrows], label = loc$label_str[lrows],
    max_bitscore = loc$max_bitscore[lrows], max_identity = loc$max_identity[lrows],
    min_evalue = loc$min_evalue[lrows], query_coverage = loc$query_coverage[lrows]
  )

  list(parents = parents, children = children,
       dropped_noncanonical = as.integer(dropped_noncanonical))
}


# Per-genome erv-like table (one row per candidate). Feeds the erv-like plots.
build_erv_like_df <- function(parents) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0), ID = character(0),
    virus = character(0), label = character(0), species = character(0),
    probes_present = character(0), n_main_present = integer(0),
    n_main_expected = integer(0), completeness_fraction = numeric(0),
    is_full = logical(0), observed_order = character(0),
    is_canonical = logical(0), n_loci = integer(0),
    max_bitscore = numeric(0), max_identity = numeric(0),
    query_coverage = numeric(0)
  )
  if (length(parents) == 0L) return(empty)

  m <- S4Vectors::mcols(parents)
  tibble::tibble(
    seqnames = as.character(GenomicRanges::seqnames(parents)),
    start    = as.integer(BiocGenerics::start(parents)),
    end      = as.integer(BiocGenerics::end(parents)),
    width    = as.integer(BiocGenerics::width(parents)),
    strand   = as.character(BiocGenerics::strand(parents)),
    ID       = as.character(m$ID),
    virus    = as.character(m$virus),
    label    = as.character(m$label),
    species  = as.character(m$species),
    probes_present        = as.character(m$probes_present),
    n_main_present        = as.integer(m$n_main_present),
    n_main_expected       = as.integer(m$n_main_expected),
    completeness_fraction = as.numeric(m$completeness_fraction),
    is_full               = as.logical(m$is_full),
    observed_order        = as.character(m$observed_order),
    is_canonical          = as.logical(m$is_canonical),
    n_loci                = as.integer(m$n_loci),
    max_bitscore          = as.numeric(m$max_bitscore),
    max_identity          = as.numeric(m$max_identity),
    query_coverage        = as.numeric(m$query_coverage)
  )
}


# Per-member table (one row per constituent locus of each candidate). `gap_to_prev`
# is the bp gap to the running-max end of the preceding members within the same
# candidate (the exact chaining gap, so it never exceeds max_join_distance); the
# first member of each candidate has NA. Feeds the inter-probe gap (tuning) plot.
build_erv_like_members_df <- function(children) {
  empty <- tibble::tibble(
    Parent = character(0), ID = character(0), probe = character(0),
    seqnames = character(0), start = integer(0), end = integer(0),
    gap_to_prev = integer(0), virus = character(0), label = character(0)
  )
  if (length(children) == 0L) return(empty)

  m <- S4Vectors::mcols(children)
  df <- data.frame(
    Parent   = as.character(m$Parent),
    ID       = as.character(m$ID),
    probe    = as.character(m$probe),
    seqnames = as.character(GenomicRanges::seqnames(children)),
    start    = as.integer(BiocGenerics::start(children)),
    end      = as.integer(BiocGenerics::end(children)),
    virus    = as.character(m$virus),
    label    = as.character(m$label),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$Parent, df$start), , drop = FALSE]
  gap <- rep(NA_integer_, nrow(df))
  for (p in unique(df$Parent)) {
    idx <- which(df$Parent == p)
    if (length(idx) < 2L) next
    st <- df$start[idx]; en <- df$end[idx]
    running_end <- en[1]
    for (j in 2:length(idx)) {
      gap[idx[j]] <- as.integer(st[j] - running_end - 1L)
      running_end <- max(running_end, en[j])
    }
  }
  df$gap_to_prev <- gap
  tibble::as_tibble(df)
}
