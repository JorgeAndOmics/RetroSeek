"""Guards for the console convention (ADR-021, docs/console_style.md).

Every script writes the same plain line through log.py or log.R, and only the
launcher draws the screen. These checks catch the old habits coming back: a
private logging setup, a bare print, a copied log_section, warnings switched
off, or raw stderr writes that the launcher cannot place.
"""

import re
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
PY_FILES = sorted(SCRIPTS.rglob("*.py"))
R_FILES = sorted(SCRIPTS.rglob("*.R"))


def _hits(files: list[Path], pattern: str, allowed: set[str]) -> list[str]:
    regex = re.compile(pattern)
    found = []
    for path in files:
        if path.name in allowed:
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if regex.search(line) and not line.lstrip().startswith("#"):
                found.append(f"{path.relative_to(SCRIPTS)}:{number}: {line.strip()}")
    return found


@pytest.mark.parametrize(
    ("pattern", "why"),
    [
        (r"logging\.basicConfig|coloredlogs|colored_logging", "use log.job_logging"),
        (r"(?<![\w.])print\(", "use a logger; only console.py draws the screen"),
        (r"sys\.stderr\.write", "log it, or use external.run_tool for tool output"),
    ],
)
def test_python_scripts_follow_the_convention(pattern: str, why: str) -> None:
    hits = _hits(PY_FILES, pattern, allowed={"console.py", "log.py"})
    assert not hits, f"{why}:\n" + "\n".join(hits)


@pytest.mark.parametrize(
    ("pattern", "why", "allowed"),
    [
        (
            r"log_section\s*<-\s*function",
            "log_section lives in utils/log.R only",
            set(),
        ),
        (
            r"options\(\s*warn\s*=",
            "run_main routes warnings; do not switch them off",
            set(),
        ),
        # style.R prints ggplot pages onto the PDF device: that draws, not logs.
        (r"(?<![\w.])print\(", "use log_info; print() writes to stdout", {"style.R"}),
    ],
)
def test_r_scripts_follow_the_convention(
    pattern: str, why: str, allowed: set[str]
) -> None:
    hits = _hits(R_FILES, pattern, allowed={"log.R", *allowed})
    assert not hits, f"{why}:\n" + "\n".join(hits)


def test_every_suppressed_warning_says_why() -> None:
    """A suppressed warning hides a signal; the same line must say why that is safe."""
    hits = [
        hit
        for hit in _hits(R_FILES, r"suppressWarnings\(", allowed=set())
        if "#" not in hit.split("suppressWarnings(", 1)[1]
    ]
    assert not hits, "add a same-line reason:\n" + "\n".join(hits)
