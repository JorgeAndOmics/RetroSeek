"""Unit tests for ``workflow/scripts/ltr_retriever_prefilter.py``.

Phase 0.1 establishes regression coverage for ``_parse_des``: the .des
files produced by ``gt suffixerator`` end with a binary 0xFF-padded
trailer that crashes naive text-mode reading. The reader must tolerate
the trailer and stop cleanly at it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ltr_retriever_prefilter import _parse_des


def test_parse_des_basic_text_file(tmp_path: Path) -> None:
    """Plain text .des file: one chromosome name per line, returned in order."""
    des = tmp_path / "basic.des"
    des.write_bytes(b"chr1\nchr2\nchr3\n")
    assert _parse_des(des) == ["chr1", "chr2", "chr3"]


def test_parse_des_handles_binary_trailer(tmp_path: Path) -> None:
    """Real-world .des files end with a 0xFF-padded trailer.

    Reproduces the failure mode that motivated the binary-mode reader:
    a UTF-8 text reader chokes on 0xFF mid-stream. The parser must
    treat the first 0xFF byte as end-of-data and return only the
    chromosome names that preceded it.
    """
    des = tmp_path / "with_trailer.des"
    des.write_bytes(b"chr1\nchr2\nchr3\n" + b"\xff" * 16)
    assert _parse_des(des) == ["chr1", "chr2", "chr3"]


def test_parse_des_strips_to_first_token(tmp_path: Path) -> None:
    """LTRharvest -seqids uses only the first whitespace-delimited token.

    A full FASTA header like ``chr1 dna:chromosome chromosome:GRCh38:1``
    must be reduced to ``chr1`` so the chromosome lookup matches the
    GFF3 seqid column.
    """
    des = tmp_path / "with_descriptions.des"
    des.write_bytes(b"chr1 dna:chromosome chromosome:GRCh38:1\nchr2 some other text\n")
    assert _parse_des(des) == ["chr1", "chr2"]


def test_parse_des_missing_file_raises(tmp_path: Path) -> None:
    """An absent .des file must raise FileNotFoundError, not silently empty out."""
    missing = tmp_path / "does_not_exist.des"
    with pytest.raises(FileNotFoundError, match="descriptor file not found"):
        _parse_des(missing)
