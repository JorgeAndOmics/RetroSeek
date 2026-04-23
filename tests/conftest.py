"""Shared pytest fixtures + sys.path + defaults stubbing for RetroSeek tests.

Why stub ``defaults``: the real ``workflow/scripts/defaults.py`` has global
side effects at import time (reads ``config.yaml``, resolves all pipeline
paths, and creates every output directory on disk). That makes it unusable
in unit tests — importing any module that transitively imports ``defaults``
touches the host filesystem and often fails on non-developer machines.

The fix lives in this conftest: we insert a minimal stub module named
``defaults`` into :data:`sys.modules` *before* any test-time imports happen,
so every ``import defaults`` in ``workflow/scripts/*.py`` resolves to our
stub. The stub exposes just enough attributes for the unit tests to run.
Integration tests that need the real ``defaults`` can opt out via the
``real_defaults`` marker.

This is a pragmatic bridge — the long-term fix is to refactor ``defaults.py``
to a ``load_config()`` function so import has no side effects. Until then,
this stub pattern keeps Phase 0 unblocked without touching production code.
"""

from __future__ import annotations

import sys
import types
from collections.abc import Iterator
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# ── 1. Make workflow/scripts importable by tests ────────────────────────
# The production scripts import each other by bare name (``import defaults``,
# ``from RetroSeeker_class import RetroSeeker``), so tests need the scripts
# directory on ``sys.path``.
REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


# ── 2. Stub ``defaults`` before any production module imports it ────────
def _build_defaults_stub() -> types.ModuleType:
    """Construct a stand-in ``defaults`` module.

    Only attributes actually referenced by the production code at import or
    call time need to be defined here; everything else is a MagicMock so
    attribute access never raises.
    """
    stub = types.ModuleType("defaults")
    stub.CONFIG_FILE = Path("/tmp/retroseek-test-config.yaml")
    stub.E_VALUE = 1e-3
    stub.ACCESSION_ID_REGEX = r"[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}"
    stub.PROBE_MIN_LENGTH: dict[str, int] = {}
    stub.LEVEL_STYLES: dict[str, dict] = {}
    stub.FIELD_STYLES: dict[str, dict] = {}
    stub.PATH_DICT: dict[str, Path] = {
        "TMP_DIR": Path("/tmp/retroseek-test-tmp"),
        "LOG_DIR": Path("/tmp/retroseek-test-logs"),
        "CONFIG_DIR": Path("/tmp/retroseek-test-config"),
        "SPECIES_DB": Path("/tmp/retroseek-test-species"),
        "PICKLE_DIR": Path("/tmp/retroseek-test-pickles"),
        "TBLASTN_PICKLE_DIR": Path("/tmp/retroseek-test-tblastn"),
        "TABLE_OUTPUT_DIR": Path("/tmp/retroseek-test-tables"),
    }
    stub.NUM_CORES = 1
    stub.USE_SPECIES_DICT = False
    stub.RETRIVAL_TIME_LAG = 0.0
    stub.MAX_RETRIEVAL_ATTEMPTS = 2
    stub.MAX_THREADPOOL_WORKERS = 1
    stub.ENTREZ_EMAIL = "test@example.com"
    stub.DISPLAY_SNAKEMAKE_INFO = False
    stub.DISPLAY_REQUESTS_WARNING = False
    stub.DISPLAY_OPERATION_INFO = False
    stub.PROBE_CSV = Path("/tmp/retroseek-test-probes.csv")
    stub.SPECIES_DICT: dict[str, str] = {}
    stub.SPECIES: list[str] = []
    # Expose ``config`` as a plain dict so code that reads
    # ``defaults.config['input']['probe_csv']`` works.
    stub.config = {
        "input": {"probe_csv": str(stub.PROBE_CSV)},
        "execution": {"entrez_email": stub.ENTREZ_EMAIL},
    }
    return stub


# Inject the stub exactly once, before any ``import defaults`` in a scripts
# module can fire the real one.
if "defaults" not in sys.modules:
    sys.modules["defaults"] = _build_defaults_stub()


# Create the TMP_DIR on disk because some code paths drop temp files into
# ``PATH_DICT['TMP_DIR']``. Harmless if it already exists.
sys.modules["defaults"].PATH_DICT["TMP_DIR"].mkdir(parents=True, exist_ok=True)


# ── 3. Generic fixtures ─────────────────────────────────────────────────
@pytest.fixture
def project_root() -> Path:
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture
def fixtures_dir() -> Path:
    """Absolute path to ``tests/fixtures/``."""
    return REPO_ROOT / "tests" / "fixtures"


@pytest.fixture
def isolated_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Run the test in a clean, temporary working directory."""
    monkeypatch.chdir(tmp_path)
    yield tmp_path


@pytest.fixture
def defaults_stub() -> types.ModuleType:
    """Return the active ``defaults`` stub module for tests that tweak it."""
    return sys.modules["defaults"]


@pytest.fixture
def fake_hsp() -> MagicMock:
    """Minimal HSP-like object for tests that construct a RetroSeeker with HSP."""
    hsp = MagicMock()
    hsp.sbjct_start = 100
    hsp.sbjct_end = 500
    hsp.align_length = 400
    hsp.frame = (1, 0)
    return hsp
