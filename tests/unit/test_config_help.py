"""Unit tests for ``workflow/scripts/config_help.py``.

The help engine parses ``docs/configuration.md`` (the canonical prose source)
and reads live values from ``data/config/config.yaml``. These tests run
against the real repository files so they double as a drift guard: if a key
the CLI advertises stops resolving, a test here fails.

``conftest.py`` puts ``workflow/scripts`` on ``sys.path`` and stubs
``defaults``; ``config_help`` itself imports neither, so a plain
``import config_help`` is enough.
"""

from __future__ import annotations

import pytest

import config_help

# Representative keys spanning bare section-children, nested-map (dotted)
# keys, and the previously-undocumented `seed`.
REPRESENTATIVE_KEYS = [
    "e_value",
    "merge_option",
    "aggregation.virus",
    "erv_like.group_by",
    "source_scn",
    "seed",
]

# The 12 top-level config sections, all of which the overview must list.
TOP_LEVEL_SECTIONS = [
    "blast",
    "genome_tools",
    "parameters",
    "ltr_retriever",
    "logging",
    "plots",
    "execution",
    "input",
    "display",
    "root",
    "domains",
    "species",
]


@pytest.mark.parametrize("key", REPRESENTATIVE_KEYS)
def test_representative_keys_resolve(key: str) -> None:
    """Each representative key renders a multi-line block with a doc link."""
    out = config_help.render(key)
    assert key in out
    assert "docs/configuration.md#" in out
    # A real field block has a header plus at least the footer link.
    assert len(out.splitlines()) >= 2
    assert not out.startswith("Unknown")
    assert not out.startswith("Ambiguous")


def test_enum_type_is_rendered() -> None:
    """A strict-enum field shows its allowed values and default."""
    out = config_help.render("merge_option")
    assert "virus | label" in out
    assert "default: virus" in out


def test_seed_drift_is_closed() -> None:
    """`seed` (the field missing from the old docs) now resolves with prose."""
    out = config_help.render("seed")
    assert "default: 67" in out
    assert "manifest" in out.lower()


def test_list_mode_covers_all_sections() -> None:
    """`render(None)` lists every top-level section."""
    out = config_help.render(None)
    for section in TOP_LEVEL_SECTIONS:
        assert section in out, f"section {section!r} missing from overview"
    assert "config-help KEY" in out


def test_bare_key_resolves_via_suffix() -> None:
    """A bare nested key resolves to its dotted doc entry."""
    out = config_help.render("group_by")
    assert "erv_like.group_by" in out


def test_dotted_and_bare_forms_agree() -> None:
    """The dotted form and its unique bare suffix resolve to the same field."""
    dotted = config_help.render("erv_like.group_by")
    bare = config_help.render("group_by")
    assert dotted == bare


def test_unknown_key_suggests_close_match() -> None:
    """A typo'd key is rejected with a close suggestion."""
    out = config_help.render("merge_optino")
    assert out.startswith("Unknown")
    assert "merge_option" in out


def test_render_has_no_filesystem_side_effects() -> None:
    """Rendering must not create directories (unlike importing ``defaults``)."""
    before = set(config_help.DEFAULT_CONFIG_PATH.parent.iterdir())
    config_help.render(None)
    config_help.render("e_value")
    after = set(config_help.DEFAULT_CONFIG_PATH.parent.iterdir())
    assert before == after
