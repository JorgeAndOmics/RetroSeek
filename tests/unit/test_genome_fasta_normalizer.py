"""Unit tests for ``workflow/scripts/genome_fasta_normalizer.py``.

Phase 1 scaffold only — the normalizer script lands in Phase 4. This
file exists ahead of time so the entry-point name is locked in.
"""

from __future__ import annotations

import importlib.util

import pytest


@pytest.mark.xfail(
    reason="genome_fasta_normalizer.py not implemented yet (Phase 4)",
    strict=True,
)
def test_genome_fasta_normalizer_module_imports() -> None:
    """The script must be importable as ``genome_fasta_normalizer``."""
    spec = importlib.util.find_spec("genome_fasta_normalizer")
    assert spec is not None
