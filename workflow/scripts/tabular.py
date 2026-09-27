"""Rows of tab-separated text: GFF3 tracks and the pipeline's own tables.

Every reader of these files needs the same first step: split each line on tabs,
skip ``#`` comment lines, and skip rows too short to hold the columns the caller
reads. Doing that in one place keeps the readers about their own columns. GFF3
readers also share the parser of column 9, the ``key=value`` attributes.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator


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
