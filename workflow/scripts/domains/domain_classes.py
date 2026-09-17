# =============================================================================
# domain_classes.py
# =============================================================================
# Interpretation of Pfam domain hits: the curated class table and the two derived
# columns that read it.
#
# Deliberately free of heavy imports and of any dependency on the scanner, because
# both the scanner (workflow/scripts/domains/scan_domains.py) and the classifier
# (workflow/scripts/taxonomy/taxonomy_classify_loci.py) import it. The scanner
# already imports the classifier for locus grouping, so putting these here is what
# keeps that from becoming a cycle.
#
# The class of a family is a biological judgement recorded in
# data/config/pfam_domain_classes.tsv, not something derivable from the sequence.
# It replaces the substring regexes in the retired `config.domains` block, which
# matched `ase` against `Transposase_22` (an L1 ORF1p domain, 29,081 instances)
# while missing `rve`, `RVP`, `IN_DBD_C` and `GP41` entirely.
# =============================================================================

from __future__ import annotations

from pathlib import Path
from typing import Any

# Ordered most to least informative about retroviral identity.
CLASS_RANK = {
    "retroviral_diagnostic": 0,
    "retroelement_shared": 1,
    "non_ltr": 2,
    "dna_transposon": 3,
    "other": 4,
}

# Families whose presence makes a locus `domain_selected`. Includes
# retroelement_shared: RT alone is genuinely ambiguous (42% of RT-bearing elements
# sit in L1 company against 32.6% retroviral) and `domain_evidence` records that
# ambiguity, but it is still protein-coding retroelement evidence.
SELECTED_CLASSES = frozenset({"retroviral_diagnostic", "retroelement_shared"})

# A family the table does not list is interpreted, never dropped.
DEFAULT_CLASS = "other"


def load_classes(classes_tsv: Path) -> dict[str, str]:
    """Map unversioned Pfam accession to curated class."""
    classes: dict[str, str] = {}
    with classes_tsv.open(encoding="utf-8") as fh:
        next(fh, None)  # header
        for raw in fh:
            row = raw.strip()
            if not row or row.startswith("#"):
                continue
            acc, _name, cls = row.split("\t")[:3]
            classes[acc.split(".")[0]] = cls
    return classes


def evidence_for(classes: set[str]) -> str:
    """The strongest class present, or `none` when nothing was found."""
    if not classes:
        return "none"
    return min(classes, key=lambda c: CLASS_RANK.get(c, CLASS_RANK[DEFAULT_CLASS]))


def tier_for(classes: set[str]) -> str:
    """ADR-009 `domain_tier`, recomputed from curated classes instead of regexes.

    Keeps the original three values and their ordering; only the mechanism moved.
    """
    if classes & SELECTED_CLASSES:
        return "domain_selected"
    if classes:
        return "domain_unlisted"
    return "non_domain"


def summarise(
    hits: list[dict[str, Any]], classes: dict[str, str]
) -> dict[str, dict[str, str]]:
    """Collapse hits to one row per locus.

    Keyed on the SET of distinct families, never on hit counts: `gt ltrdigest`
    chains fragments of one model into a single feature and `hmmsearch` does not,
    so a count-based statistic would differ between them for reasons unrelated to
    biology.

    A family absent from the curated table is classed `other` and still listed in
    `domain_names`. Dropping it would reintroduce the silent loss this feature
    exists to remove.
    """
    names: dict[str, set[str]] = {}
    found: dict[str, set[str]] = {}
    for hit in hits:
        locus = hit["locus_id"]
        names.setdefault(locus, set()).add(hit["pfam_name"])
        found.setdefault(locus, set()).add(classes.get(hit["pfam_acc"], DEFAULT_CLASS))
    return {
        locus: {
            "domain_evidence": evidence_for(found[locus]),
            "domain_names": ",".join(sorted(names[locus])),
            "domain_tier": tier_for(found[locus]),
        }
        for locus in names
    }
