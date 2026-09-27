"""Tests for tabular.py: the shared row, GFF3 feature and attribute readers."""

from __future__ import annotations

from pathlib import Path

import pytest

from log import PipelineError
from tabular import gff3_attributes, gff3_features, tab_rows


def test_splits_on_tabs_without_the_newline() -> None:
    assert list(tab_rows(["a\tb\tc\n"], 3)) == [["a", "b", "c"]]


def test_comment_lines_are_skipped() -> None:
    assert list(tab_rows(["##gff-version 3\n", "# a\tb\tc\n", "x\ty\tz\n"], 3)) == [
        ["x", "y", "z"]
    ]


def test_short_blank_and_tabless_rows_are_skipped() -> None:
    lines = ["a\tb\n", "\n", "   \n", "no tabs\n", "a\tb\tc\n"]
    assert list(tab_rows(lines, 3)) == [["a", "b", "c"]]


def test_longer_rows_pass_and_spaces_are_kept() -> None:
    assert list(tab_rows([" a \tb\tc\td"], 3)) == [[" a ", "b", "c", "d"]]


def test_empty_input_yields_nothing() -> None:
    assert list(tab_rows([], 1)) == []


# ---- gff3_attributes: the one reader of GFF3 column 9 ----


def test_attributes_become_a_dict() -> None:
    assert gff3_attributes("ID=h1;probe=pol;Parent=LTR_retrotransposon7") == {
        "ID": "h1",
        "probe": "pol",
        "Parent": "LTR_retrotransposon7",
    }


def test_a_key_is_never_found_inside_a_longer_key() -> None:
    # A pattern search for "probe=" would read `subprobe` here.
    assert gff3_attributes("subprobe=ENV;probe=POL")["probe"] == "POL"
    assert "probe" not in gff3_attributes("subprobe=ENV")


def test_a_value_keeps_its_equals_signs() -> None:
    assert gff3_attributes("label=a=b") == {"label": "a=b"}


def test_empty_entries_and_values_are_ignored_and_spaces_trimmed() -> None:
    assert gff3_attributes(" ID = h1 ;;flag;empty=; Parent=p\n") == {
        "ID": "h1",
        "Parent": "p",
    }


def test_a_repeated_key_keeps_its_first_value() -> None:
    # GFF3 allows a key once per feature; a repeat is malformed, so pick one rule.
    assert gff3_attributes("ID=first;ID=second") == {"ID": "first"}


def test_an_empty_column_has_no_attributes() -> None:
    assert gff3_attributes("") == {}
    assert gff3_attributes(".") == {}


# ---- gff3_features: the strict reader of the pipeline's own GFF3 tracks ----

_ROW = "chr1\tsrc\tLTR_retrotransposon\t100\t200\t.\t+\t.\tID=e1"


def _track(tmp_path: Path, *lines: str) -> Path:
    path = tmp_path / "track.gff3"
    path.write_text("".join(f"{line}\n" for line in lines))
    return path


def test_feature_rows_come_with_integer_coordinates(tmp_path: Path) -> None:
    rows = list(gff3_features(_track(tmp_path, "##gff-version 3", "", _ROW)))
    assert [(f[2], start, end) for f, start, end in rows] == [
        ("LTR_retrotransposon", 100, 200)
    ]


def test_reading_stops_at_the_fasta_section(tmp_path: Path) -> None:
    track = _track(tmp_path, _ROW, "##FASTA", ">chr1", "ACGT")
    assert len(list(gff3_features(track))) == 1


@pytest.mark.parametrize(
    "bad", ["chr1\tsrc\tx\t1\t2", "chr1\tsrc\tx\tone\t2\t.\t+\t.\tID=a", "no tabs"]
)
def test_a_row_that_is_not_a_feature_names_its_line(tmp_path: Path, bad: str) -> None:
    with pytest.raises(PipelineError, match=r"track\.gff3, line 3"):
        list(gff3_features(_track(tmp_path, "# note", _ROW, bad)))
