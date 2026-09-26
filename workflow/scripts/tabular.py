"""Rows of tab-separated text: GFF3 tracks and the pipeline's own tables.

Every reader of these files needs the same first step: split each line on tabs,
skip ``#`` comment lines, and skip rows too short to hold the columns the caller
reads. Doing that in one place keeps the readers about their own columns.
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
