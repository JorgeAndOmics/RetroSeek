"""Smoke test: the project layout is importable and key files exist.

This is scaffolding so ``make test-py`` has something to run before real
tests are written. Remove or replace once real unit tests arrive.
"""

from __future__ import annotations

from pathlib import Path


def test_project_root_has_key_files(project_root: Path) -> None:
    """Confirm the repo has the files this test suite expects.

    Guards against accidentally running tests from the wrong working
    directory or after a refactor that moves core config files.
    """
    assert (project_root / "workflow" / "Snakefile").is_file()
    assert (project_root / "data" / "config" / "config.yaml").is_file()
    assert (project_root / "data" / "config" / "schema.yaml").is_file()
    assert (project_root / "data" / "config" / "environment.yml").is_file()
    assert (project_root / "RetroSeek").is_file()


def test_python_scripts_directory_is_populated(project_root: Path) -> None:
    """Every Python script under workflow/scripts should be importable by path."""
    scripts = list((project_root / "workflow" / "scripts").glob("*.py"))
    assert scripts, "workflow/scripts/ should contain Python modules"
