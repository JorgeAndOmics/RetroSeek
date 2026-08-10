"""Guard: every tracked text file must be ASCII-only.

This is not style policing. `yaml::read_yaml` -> `readLines` HALTS on a byte that
is invalid in the read encoding and returns the file TRUNCATED at that point, with
no error. An em-dash in a `config.yaml` comment - above `probe_min_length` - once
truncated the parse there and produced "0 of N hits kept" for every one of 102 bat
genomes. Nothing failed; the run just silently produced nothing.

ASCII-only makes that class of failure impossible, in configs and everywhere else.

The only exemptions are proper nouns: people, labs and cited authors. A name is
data, not punctuation, and none sits on a parsing path. They are matched as
SUBSTRINGS, so a stray accented letter outside a name is still caught.

NOTE: the allowlist below is written with \\u escapes so that THIS file stays
pure ASCII and satisfies its own invariant.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

# Real people / labs / cited authors - preserved verbatim.
ALLOWED_PROPER_NOUNS = (
    "Gonz\u00e1lez",  # Gonzalez  - maintainer (ADR Deciders)
    "Garc\u00eda",  # Garcia    - maintainer (ADR Deciders)
    "N\u00ed Leathlobhair",  # Ni        - lab (README)
    "M\u00f6lder",  # Molder    - Snakemake citation (README)
    "Gon\u00e7alves",  # Goncalves - citation (04_FINDINGS)
)


def _tracked_text_files() -> list[Path]:
    """Every git-tracked file that decodes as text."""
    out = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    files = []
    for rel in out:
        p = REPO_ROOT / rel
        if not p.is_file():
            continue
        try:
            p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue  # binary / unreadable - not our concern
        files.append(p)
    return files


def _strip_allowed(text: str) -> str:
    for noun in ALLOWED_PROPER_NOUNS:
        text = text.replace(noun, "")
    return text


def test_tracked_files_are_ascii_only() -> None:
    """No tracked text file may contain a non-ASCII character.

    Fails with file:line:char so the offender is fixable at a glance. If you hit
    this: replace the character with an ASCII equivalent that PRESERVES MEANING
    (`>=` not `>`, `->` not `-`), never delete it.
    """
    offenders: list[str] = []
    for path in _tracked_text_files():
        rel = path.relative_to(REPO_ROOT)
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            bad = {c for c in _strip_allowed(line) if ord(c) > 0x7F}
            if bad:
                shown = ", ".join(f"{c!r} (U+{ord(c):04X})" for c in sorted(bad))
                offenders.append(f"{rel}:{lineno}: {shown}")

    assert not offenders, (
        f"{len(offenders)} line(s) contain non-ASCII characters. A non-ASCII byte "
        "can truncate a config parse mid-file and silently zero a whole run "
        "(see gotchas #29). Replace with a meaning-preserving ASCII equivalent:\n  "
        + "\n  ".join(offenders[:40])
    )


def test_config_files_are_strictly_ascii_with_no_exemptions() -> None:
    """Configs get no proper-noun exemption at all - they are the parsing path."""
    offenders = []
    for path in _tracked_text_files():
        rel = str(path.relative_to(REPO_ROOT))
        if not rel.startswith("data/config/"):
            continue
        text = path.read_text(encoding="utf-8")
        if any(ord(c) > 0x7F for c in text):
            offenders.append(rel)
    assert not offenders, f"config files must be pure ASCII: {offenders}"


@pytest.mark.parametrize("noun", ALLOWED_PROPER_NOUNS)
def test_allowlist_entries_are_still_used(noun: str) -> None:
    """A stale allowlist entry silently widens the exemption - drop it when the
    name it protects is gone."""
    used = any(noun in p.read_text(encoding="utf-8") for p in _tracked_text_files())
    assert used, f"allowlist entry {noun!r} no longer appears in the repo; remove it"
