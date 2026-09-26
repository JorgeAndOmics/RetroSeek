"""Unit tests for workflow/scripts/domains/subset_pfam.py.

The subsetter is the only thing standing between a curated accession list and a
2.2 GB Pfam flat file, so the tests pin the two behaviours that matter: it must
select records by ACCESSION (ignoring the version suffix Pfam appends), and it
must fail loudly naming anything it could not find rather than silently emitting
a short library.
"""

from pathlib import Path

import pytest
import subset_pfam


def _record(name: str, acc: str) -> str:
    """A minimal but structurally faithful HMMER3 record."""
    return (
        "HMMER3/f [3.4 | Aug 2023]\n"
        f"NAME  {name}\n"
        f"ACC   {acc}\n"
        f"DESC  a test model called {name}\n"
        "LENG  10\n"
        "GA    25.00 25.00;\n"
        "HMM   A C D E\n"
        "//\n"
    )


@pytest.fixture
def hmm_file(tmp_path: Path) -> Path:
    p = tmp_path / "Pfam-A.hmm"
    p.write_text(
        _record("rve", "PF00665.33")
        + _record("RVT_1", "PF00078.31")
        + _record("Kinesin", "PF00225.28")
    )
    return p


@pytest.fixture
def classes_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "classes.tsv"
    p.write_text(
        "pfam_acc\tpfam_name\tclass\n"
        "PF00665\trve\tretroviral_diagnostic\n"
        "PF00078\tRVT_1\tretroelement_shared\n"
    )
    return p


def test_wanted_accessions_strips_header_and_keeps_order(classes_tsv: Path) -> None:
    assert subset_pfam.wanted_accessions(classes_tsv) == {"PF00665", "PF00078"}


def test_subset_selects_only_requested_records(
    hmm_file: Path, classes_tsv: Path, tmp_path: Path
) -> None:
    out = tmp_path / "subset.hmm"
    n = subset_pfam.subset_pfam(hmm_file, classes_tsv, out)
    assert n == 2
    text = out.read_text()
    assert "NAME  rve" in text
    assert "NAME  RVT_1" in text
    assert "Kinesin" not in text  # the host model must not survive


def test_subset_matches_versioned_accessions(
    hmm_file: Path, classes_tsv: Path, tmp_path: Path
) -> None:
    """Pfam writes `ACC PF00665.33`; the table carries the unversioned key."""
    out = tmp_path / "subset.hmm"
    subset_pfam.subset_pfam(hmm_file, classes_tsv, out)
    assert "ACC   PF00665.33" in out.read_text()


def test_subset_emits_complete_records(
    hmm_file: Path, classes_tsv: Path, tmp_path: Path
) -> None:
    """Every emitted record must run from the HMMER3 banner to its `//`."""
    out = tmp_path / "subset.hmm"
    subset_pfam.subset_pfam(hmm_file, classes_tsv, out)
    text = out.read_text()
    assert text.count("HMMER3/f") == 2
    assert text.count("//") == 2
    assert text.startswith("HMMER3/f")
    assert text.rstrip().endswith("//")


def test_missing_accession_fails_loudly(hmm_file: Path, tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    bad.write_text("pfam_acc\tpfam_name\tclass\nPF99999\tnope\tother\n")
    out = tmp_path / "subset.hmm"
    with pytest.raises(subset_pfam.PipelineError, match="PF99999"):
        subset_pfam.subset_pfam(hmm_file, bad, out)


def test_missing_accession_message_says_why_and_how_to_fix(
    hmm_file: Path, tmp_path: Path
) -> None:
    """The bare list of accessions once left a user deleting files by hand.

    The message must name the likely cause (an older Pfam release) and the command
    that fixes it.
    """
    bad = tmp_path / "bad.tsv"
    bad.write_text("pfam_acc\tpfam_name\tclass\nPF99999\tnope\tother\n")
    with pytest.raises(subset_pfam.PipelineError) as caught:
        subset_pfam.subset_pfam(hmm_file, bad, tmp_path / "subset.hmm")
    message = str(caught.value)
    assert "older Pfam release" in message
    assert "input.pfam_release" in message
    assert "--download-hmm" in message


def test_missing_accessions_lists_what_the_library_lacks(hmm_file: Path) -> None:
    wanted = {"PF00665", "PF99999", "PF88888"}
    assert subset_pfam.missing_accessions(hmm_file, wanted) == ["PF88888", "PF99999"]


def test_missing_accessions_is_empty_when_all_are_present(hmm_file: Path) -> None:
    assert subset_pfam.missing_accessions(hmm_file, {"PF00665", "PF00078"}) == []


def test_empty_table_is_an_error_not_an_empty_library(
    hmm_file: Path, tmp_path: Path
) -> None:
    empty = tmp_path / "empty.tsv"
    empty.write_text("pfam_acc\tpfam_name\tclass\n")
    out = tmp_path / "subset.hmm"
    with pytest.raises(subset_pfam.PipelineError, match="no accessions"):
        subset_pfam.subset_pfam(hmm_file, empty, out)


def test_name_map_covers_every_model(hmm_file: Path, tmp_path: Path) -> None:
    """The name map covers every model in the full Pfam library.

    It must resolve LTRdigest names for families that are not in the curated
    subset.
    """
    out = tmp_path / "names.tsv"
    n = subset_pfam.write_name_map(hmm_file, out)
    assert n == 3
    rows = dict(line.split("\t") for line in out.read_text().splitlines())
    assert rows["rve"] == "PF00665"
    assert rows["Kinesin"] == "PF00225"  # present even though never curated
