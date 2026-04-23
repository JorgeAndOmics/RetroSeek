"""Shared pytest fixtures for RetroSeek tests.

Fixtures defined here are available to every test file under ``tests/``.
Add helpers here when two or more tests would otherwise duplicate setup.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pytest


@pytest.fixture
def project_root() -> Path:
    """Absolute path to the repository root.

    Useful for tests that need to resolve paths relative to the project
    (e.g., ``data/config/config.yaml``) without relying on CWD.
    """
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def fixtures_dir() -> Path:
    """Absolute path to ``tests/fixtures/``."""
    return Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def isolated_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Run the test in a clean, temporary working directory.

    Yields the temporary path; restores the original CWD on teardown.
    """
    monkeypatch.chdir(tmp_path)
    yield tmp_path
