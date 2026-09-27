"""Rows of tab-separated text: GFF3 tracks and the pipeline's own tables.

Every reader of these files needs the same first step: split each line on tabs,
skip ``#`` comment lines, and skip rows too short to hold the columns the caller
reads. Doing that in one place keeps the readers about their own columns. The
GFF3 tracks the pipeline writes are read stricter, through ``gff3_features``, and
all GFF3 readers share the parser of column 9, the ``key=value`` attributes.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path

from log import PipelineError


def tab_rows(lines: Iterable[str], min_fields: int) -> Iterator[list[str]]:
    """Yield each line split on tabs, without its newline.

    Lines starting with ``#`` are skipped, and so are rows with fewer than
    ``min_fields`` fields, which covers blank and malformed lines. Values are not
    stripped: a field keeps any spaces it was written with.
    """
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= min_fields:
            yield fields


def gff3_attributes(column: str) -> dict[str, str]:
    """The ``key=value`` attributes of a GFF3 column 9, as a dict.

    Keys and values are stripped of spaces, and an entry without a value (or the
    GFF3 placeholder ``.``) adds nothing. A value keeps any ``=`` after the first.
    GFF3 allows a key once per feature; if one repeats anyway, its first value
    is kept. Values are returned as written: percent-escapes are the caller's.
    """
    attributes: dict[str, str] = {}
    for entry in column.split(";"):
        key, _, value = entry.partition("=")
        if value.strip():
            attributes.setdefault(key.strip(), value.strip())
    return attributes


def _coordinates(fields: list[str]) -> tuple[int, int] | None:
    """Start and end of a GFF3 row, or None when it is not a feature row."""
    if len(fields) < 9:
        return None
    try:
        return int(fields[3]), int(fields[4])
    except ValueError:
        return None


def gff3_features(gff3: Path) -> Iterator[tuple[list[str], int, int]]:
    """Every feature row of a GFF3 track the pipeline wrote, with its start and end.

    Comment and blank lines are skipped, and reading stops at a ``##FASTA``
    directive: the GFF3 spec puts sequences, not features, after it. Any other
    row must be a feature row, nine tab-separated columns with whole-number
    coordinates. The pipeline writes these tracks, so a row that is not means a
    damaged file, and skipping it would quietly lose an element, an arm or a
    locus.

    Raises:
        PipelineError: Naming the file and line of the first row that is not a
            feature row.
    """
    with gff3.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            coordinates = _coordinates(fields)
            if coordinates is None:
                raise PipelineError(
                    f"{gff3}, line {number}: not a GFF3 feature row (nine "
                    "tab-separated columns, whole-number start and end)",
                    hint="regenerate the file with the rule that writes it",
                )
            yield fields, *coordinates
