"""Unit tests for workflow/scripts/domains/domain_classes.py.

These pin the meaning of the curated table: which class wins when a locus carries
several, how `domain_tier` is rebuilt from classes rather than name regexes, and
the rule that a family missing from the table is still recorded rather than
silently dropped.
"""

from pathlib import Path

import domain_classes as scan_domains


def test_evidence_picks_the_strongest_class_present() -> None:
    assert (
        scan_domains.evidence_for({"retroelement_shared", "retroviral_diagnostic"})
        == "retroviral_diagnostic"
    )
    assert scan_domains.evidence_for({"non_ltr", "other"}) == "non_ltr"
    assert scan_domains.evidence_for(set()) == "none"


def test_tier_selected_requires_a_retroviral_or_shared_family() -> None:
    assert scan_domains.tier_for({"retroviral_diagnostic"}) == "domain_selected"
    assert scan_domains.tier_for({"retroelement_shared"}) == "domain_selected"


def test_tier_unlisted_when_families_are_present_but_uninformative() -> None:
    assert scan_domains.tier_for({"non_ltr"}) == "domain_unlisted"
    assert scan_domains.tier_for({"other", "dna_transposon"}) == "domain_unlisted"


def test_tier_non_domain_only_when_nothing_was_found() -> None:
    assert scan_domains.tier_for(set()) == "non_domain"


def test_summarise_keys_on_distinct_families_not_hit_counts(tmp_path: Path) -> None:
    """Two hits of the same family summarise identically to one hit of it.

    Keying on distinct families keeps the result independent of fragment
    chaining.
    """
    classes = {"PF00665": "retroviral_diagnostic"}
    one = [{"locus_id": "L", "pfam_acc": "PF00665", "pfam_name": "rve"}]
    two = [*one, {"locus_id": "L", "pfam_acc": "PF00665", "pfam_name": "rve"}]
    assert scan_domains.summarise(one, classes) == scan_domains.summarise(two, classes)


def test_summarise_reports_sorted_distinct_names(tmp_path: Path) -> None:
    classes = {"PF00665": "retroviral_diagnostic", "PF00078": "retroelement_shared"}
    hits = [
        {"locus_id": "L", "pfam_acc": "PF00078", "pfam_name": "RVT_1"},
        {"locus_id": "L", "pfam_acc": "PF00665", "pfam_name": "rve"},
    ]
    got = scan_domains.summarise(hits, classes)["L"]
    assert got["domain_names"] == "RVT_1,rve"
    assert got["domain_evidence"] == "retroviral_diagnostic"
    assert got["domain_tier"] == "domain_selected"


def test_family_absent_from_the_table_is_recorded_but_classed_other() -> None:
    """A hit we cannot interpret must still appear in domain_names.

    Dropping it would reintroduce the silent-loss bug this feature exists to
    remove.
    """
    got = scan_domains.summarise(
        [{"locus_id": "L", "pfam_acc": "PF99999", "pfam_name": "Mystery"}], {}
    )["L"]
    assert got["domain_names"] == "Mystery"
    assert got["domain_evidence"] == "other"
    assert got["domain_tier"] == "domain_unlisted"


def test_load_classes_maps_accession_to_class(tmp_path: Path) -> None:
    p = tmp_path / "c.tsv"
    p.write_text(
        "pfam_acc\tpfam_name\tclass\n"
        "PF00665\trve\tretroviral_diagnostic\n"
        "PF00225\tKinesin\tother\n"
    )
    assert scan_domains.load_classes(p) == {
        "PF00665": "retroviral_diagnostic",
        "PF00225": "other",
    }
