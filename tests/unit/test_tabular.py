"""Tests for tabular.tab_rows, the shared first step of every GFF3/TSV reader."""

from __future__ import annotations

from tabular import tab_rows


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
