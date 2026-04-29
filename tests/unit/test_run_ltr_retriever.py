"""Unit tests for ``workflow/scripts/run_ltr_retriever.py``.

Phase 1 scaffold only — the runner script lands in Phase 3. This file
exists ahead of time so the entry-point name (``run_ltr_retriever``) is
locked in, and so Phase 3's TDD work has somewhere obvious to add tests.
"""

from __future__ import annotations

import importlib.util

import pytest


@pytest.mark.xfail(
    reason="run_ltr_retriever.py not implemented yet (Phase 3)",
    strict=True,
)
def test_run_ltr_retriever_module_imports() -> None:
    """The script must be importable as ``run_ltr_retriever``."""
    spec = importlib.util.find_spec("run_ltr_retriever")
    assert spec is not None
