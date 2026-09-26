"""Guards for the house style of every RetroSeek figure (docs/visual_style.md).

Three rules are cheap to break by accident and tedious to spot by eye, so they
are checked mechanically:

1. No em dash or en dash anywhere in the project, tracked or not (local notes,
   notebooks and CLAUDE.md included). The owner's rule, without exceptions.
2. No colour written as a hex literal outside plot2sort/style.R, which is the
   single source of every colour. One documented exception: placement_figures.py
   mirrors the few house colours gappa needs (a separate test pins the mirror to
   style.R).
3. No dash or arrow used as punctuation inside the string literals of the
   plotting scripts: a page says "Removed: flanks", never "- flanks" or "a -> b".
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"

_SKIP_DIRS = {
    ".git",
    ".snakemake",
    "__pycache__",
    ".ruff_cache",
    ".mypy_cache",
    ".pytest_cache",
    ".Rproj.user",
}
_DASHES = ("\u2014", "\u2013")  # em dash, en dash

_HEX_ALLOWED = {
    SCRIPTS / "plot2sort" / "style.R",
    SCRIPTS / "taxonomy" / "placement_figures.py",
}
_HEX = re.compile(r"""["']#[0-9A-Fa-f]{6}(?:[0-9A-Fa-f]{2})?["']""")

# Scripts whose string literals end up on a page (titles, labels, key pages) or
# in the log beside them.
_PLOTTING_SCRIPTS = [
    *sorted((SCRIPTS / "plot2sort").glob("*.R")),
    *sorted((SCRIPTS / "stage_plot_generator").glob("*.R")),
    *sorted((SCRIPTS / "taxonomy").glob("*.R")),
    *sorted((SCRIPTS / "solo_ltr").glob("*.R")),
    SCRIPTS / "hotspot" / "plots.R",
    SCRIPTS / "plot2sort.R",
    SCRIPTS / "stage_plot_generator.R",
    SCRIPTS / "demo_figures.R",
    SCRIPTS / "hotspot_detector.R",
]


def _text_files(root: Path) -> list[Path]:
    """Every UTF-8 text file under `root`, skipping VCS and tool caches."""
    files = []
    for path in root.rglob("*"):
        if any(part in _SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\x00" in data[:4096]:
            continue  # binary
        try:
            data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        files.append(path)
    return files


def find_dashes(root: Path) -> list[str]:
    """`path:line` for every line under `root` holding an em or en dash."""
    found = []
    for path in _text_files(root):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if any(dash in line for dash in _DASHES):
                found.append(f"{path.relative_to(root)}:{number}")
    return found


def r_string_literals(line: str) -> list[str]:
    """The string literals on one line of R, ignoring anything after a comment.

    A small scanner rather than a regex: a `#` inside a string (a colour, a
    format) is not a comment, and a quote inside a comment is not a string.
    """
    literals: list[str] = []
    quote: str | None = None
    current: list[str] = []
    i = 0
    while i < len(line):
        char = line[i]
        if quote:
            if char == "\\":
                current.append(line[i : i + 2])
                i += 2
                continue
            if char == quote:
                literals.append("".join(current))
                quote, current = None, []
            else:
                current.append(char)
        elif char == "#":
            break
        elif char in "\"'":
            quote = char
        i += 1
    return literals


def test_find_dashes_catches_a_planted_dash(tmp_path: Path) -> None:
    """The scanner itself must fail on a dash, or the repo-wide test proves nothing."""
    (tmp_path / "clean.md").write_text("a plain hyphen-word and a colon: fine\n")
    (tmp_path / "dirty.md").write_text("one line\nan aside \u2014 like this\n")
    (tmp_path / "range.md").write_text("pages 3\u20135\n")
    assert sorted(find_dashes(tmp_path)) == ["dirty.md:2", "range.md:1"]


def test_no_em_or_en_dash_anywhere_in_the_project() -> None:
    found = find_dashes(ROOT)
    assert not found, "em/en dashes found (rewrite with a colon, comma or 'to'):\n" + (
        "\n".join(found[:50])
    )


def test_no_hex_colour_outside_style_r() -> None:
    offenders = []
    for path in [*SCRIPTS.rglob("*.R"), *SCRIPTS.rglob("*.py")]:
        if path in _HEX_ALLOWED:
            continue
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if _HEX.search(line):
                offenders.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()}")
    assert not offenders, "colours belong in plot2sort/style.R:\n" + "\n".join(
        offenders
    )


def test_no_foreign_palette_or_theme_in_plotting_scripts() -> None:
    """ggsci palettes, brewer/viridis scales and stock themes bypass the house style."""
    banned = re.compile(
        r"library\(ggsci\)|ggsci::|pal_futurama|scale_(?:fill|colou?r)_"
        r"(?:brewer|viridis_[cd]|distiller|npg|aaas|jco|igv|futurama)|"
        r"theme_(?:bw|minimal|classic|light|gr[ae]y)\(|\"gr[ae]y\d+\""
    )
    offenders = []
    for path in _PLOTTING_SCRIPTS:
        if path == SCRIPTS / "plot2sort" / "style.R":
            continue
        for number, line in enumerate(path.read_text().splitlines(), 1):
            code = line.split("#", 1)[0] if not line.lstrip().startswith("#") else ""
            if banned.search(code):
                offenders.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()}")
    assert not offenders, "use style.R's scales and theme:\n" + "\n".join(offenders)


@pytest.mark.parametrize("path", _PLOTTING_SCRIPTS, ids=lambda p: p.name)
def test_no_dash_or_arrow_as_punctuation_in_plot_text(path: Path) -> None:
    offenders = [
        f"{number}: {literal!r}"
        for number, line in enumerate(path.read_text().splitlines(), 1)
        for literal in r_string_literals(line)
        if " - " in literal or "->" in literal or literal.startswith("- ")
    ]
    assert not offenders, f"{path.name}: dash or arrow as punctuation:\n" + "\n".join(
        offenders
    )


def test_r_string_scanner_ignores_comments_and_keeps_hashes_in_strings() -> None:
    line = 'x <- c("#332288", "a - b") # a comment - with "quotes"'
    assert r_string_literals(line) == ["#332288", "a - b"]
