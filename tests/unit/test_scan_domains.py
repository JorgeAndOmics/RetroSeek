"""Unit tests for workflow/scripts/domains/scan_domains.py.

Three things are easy to get wrong here and all three have bitten this project
or a neighbouring one before:

* `--domtblout` field order. Field 1 is the TARGET sequence and field 4 is the
  QUERY model; swapping them silently produces plausible nonsense.
* Six-frame translation. Losing the reverse frames systematically drops domains
  on the minus strand.
* Deriving a tier from a SET of families rather than a hit count. Keying on
  counts would make results depend on whether the search chained fragments,
  which `hmmsearch` does not do and `gt ltrdigest` does.
"""

from pathlib import Path

import pytest
import scan_domains


# ------------------------------------------------------------ translation
def test_six_frames_returns_all_six() -> None:
    frames = scan_domains.six_frame_translations("ATGGCCATTGTAATGGGCCGC")
    assert sorted(frames) == [-3, -2, -1, 1, 2, 3]


def test_frame_one_is_a_plain_translation() -> None:
    # ATG GCC ATT -> M A I
    assert scan_domains.six_frame_translations("ATGGCCATT")[1] == "MAI"


def test_negative_frames_translate_the_reverse_complement() -> None:
    # revcomp(ATGGCCATT) = AATGGCCAT -> AAT GGC CAT -> N G H
    assert scan_domains.six_frame_translations("ATGGCCATT")[-1] == "NGH"


def test_short_sequence_yields_empty_strings_not_an_error() -> None:
    frames = scan_domains.six_frame_translations("AT")
    assert all(p == "" for p in frames.values())


def test_query_fasta_roundtrips_locus_ids_containing_pipes(tmp_path: Path) -> None:
    """Locus ids carry `|` already (`{id}|{gene}` elsewhere), so the frame suffix
    must be split off from the right."""
    fasta = tmp_path / "q.faa"
    n = scan_domains.write_query_fasta(
        {"chr1|LTR_retrotransposon2": "ATGGCCATT"}, fasta
    )
    assert n == 6
    names = [ln[1:] for ln in fasta.read_text().splitlines() if ln.startswith(">")]
    assert "chr1|LTR_retrotransposon2|f1" in names
    assert scan_domains.locus_of("chr1|LTR_retrotransposon2|f-2") == (
        "chr1|LTR_retrotransposon2"
    )


# --------------------------------------------------------------- parsing
_DOMTBL = """\
# target name        accession   tlen query name  accession   qlen   E-value  score
locusA|f1            -            120 rve         PF00665.33   102   1.2e-30  105.4   0.1   1 1 1e-33 2e-30 104.9 0.1 1 100 5 106
locusA|f-2           -            120 RVT_1       PF00078.31   250   3.0e-10   40.2   0.0   1 1 2e-13 4e-10  39.8 0.0 1 200 3 210
locusB|f3            -             90 Kinesin     PF00225.28   340   1.0e-05   25.0   0.0   1 1 9e-08 2e-05  24.6 0.0 1 300 2 305
#
"""


@pytest.fixture
def domtbl(tmp_path: Path) -> Path:
    p = tmp_path / "hits.domtbl"
    p.write_text(_DOMTBL)
    return p


def test_parse_reads_target_from_field1_and_model_from_field4(domtbl: Path) -> None:
    hits = scan_domains.parse_domtblout(domtbl)
    assert len(hits) == 3
    first = hits[0]
    assert first["locus_id"] == "locusA"
    assert first["pfam_name"] == "rve"
    assert first["pfam_acc"] == "PF00665"  # version stripped
    assert first["frame"] == 1


def test_parse_skips_comment_lines(domtbl: Path) -> None:
    assert all(
        not h["pfam_name"].startswith("#") for h in scan_domains.parse_domtblout(domtbl)
    )


def test_parse_of_empty_file_returns_no_hits(tmp_path: Path) -> None:
    p = tmp_path / "none.domtbl"
    p.write_text("# all comments\n#\n")
    assert scan_domains.parse_domtblout(p) == []
