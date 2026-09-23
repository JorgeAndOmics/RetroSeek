"""Tests for the root `./RetroSeek` shim (ADR-020).

The shim once threw away the pipeline's exit code, so every run looked like a
success. These tests run it as a real process: help works without a config, and
a failing pipeline makes the command fail.
"""

import importlib.util
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

SHIM = Path(__file__).resolve().parents[2] / "RetroSeek"


def _run(*argv: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SHIM), *argv],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def test_help_needs_no_config_and_lists_the_phases() -> None:
    result = _run("-h")
    assert result.returncode == 0
    for phase in ("Setup:", "Indexing:", "Discovery:", "Analysis:", "Figures:"):
        assert phase in result.stdout


def test_no_stage_prints_help_and_succeeds() -> None:
    result = _run("-skp")
    assert result.returncode == 0
    assert "Analysis:" in result.stdout


def test_a_failing_pipeline_fails_the_command(tmp_path: Path) -> None:
    """The pipeline cannot even read a missing config; the shim must say so."""
    result = _run("--classify", "-skp", "--configfile", str(tmp_path / "absent.yaml"))
    assert result.returncode != 0


@pytest.fixture
def shim():
    """The shim as a module (it has no .py suffix, hence the explicit loader)."""
    loader = SourceFileLoader("retroseek_shim", str(SHIM))
    spec = importlib.util.spec_from_loader("retroseek_shim", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_configfile_becomes_absolute(shim, tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    argv, config = shim.absolute_configfile(
        ["--classify", "--configfile", "c.yaml", "-n"]
    )
    assert config == str(tmp_path / "c.yaml")
    assert argv == ["--classify", "--configfile", str(tmp_path / "c.yaml"), "-n"]


def test_configfile_with_equals_becomes_absolute(
    shim, tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.chdir(tmp_path)
    argv, config = shim.absolute_configfile(["--configfile=c.yaml"])
    assert argv == [f"--configfile={tmp_path / 'c.yaml'}"]
    assert config == str(tmp_path / "c.yaml")


def test_no_configfile_means_no_override(shim) -> None:
    assert shim.absolute_configfile(["--classify"]) == (["--classify"], None)
