"""Render config-field help for the ``RetroSeek --config-help`` CLI.

The single source of prose is [`docs/configuration.md`](../../docs/configuration.md):
this module parses its per-field Markdown tables (the ``Key | Type | Default |
Meaning`` blocks) so the terminal help can never drift from the documentation.
Type/enum/range shown to the user therefore mirror the doc, and the *current*
value is read live from ``data/config/config.yaml``.

The module is deliberately self-contained — it must **not** ``import defaults``,
whose import has filesystem side effects (it resolves and creates every
pipeline output directory). A ``--config-help`` lookup must stay read-only.

Table conventions parsed (doc-specific, by column count):

* 4 columns ``Key | Type | Default | Meaning`` — the standard field table.
* 3 columns ``Key | Type | Default`` — e.g. the ``display`` section.
* 2 columns ``Key | Meaning`` — e.g. the ``root`` section.

Only tables whose header row's first cell is exactly ``Key`` are treated as
field tables; everything else (e.g. the aggregation-strategy *vocabulary*
table) is ignored. Nested-map keys are addressed dotted in the doc
(``aggregation.virus``); section-level keys are bare (``e_value``).
"""

from __future__ import annotations

import difflib
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DOCS_PATH = _REPO_ROOT / "docs" / "configuration.md"
DEFAULT_CONFIG_PATH = _REPO_ROOT / "data" / "config" / "config.yaml"

_WRAP_WIDTH = 76
_LIST_KEY_WIDTH = 30
_LIST_MEANING_WIDTH = 60


@dataclass(frozen=True)
class FieldDoc:
    """One documented config field, extracted from a Markdown table row."""

    key: str
    """Key as written in the doc (bare like ``e_value`` or dotted like
    ``aggregation.virus``)."""
    type: str
    default: str
    meaning: str
    section: str
    """Top-level (``##``) section the field lives under, e.g. ``parameters``."""
    anchor: str
    """Nearest-heading anchor used for the ``docs/configuration.md#…`` link."""


@dataclass(frozen=True)
class DocIndex:
    """Parsed view of ``docs/configuration.md``."""

    fields: dict[str, FieldDoc]
    sections: list[str]
    section_anchor: dict[str, str]


# ── Markdown helpers ────────────────────────────────────────────────────


def _strip_md(text: str) -> str:
    """Reduce inline Markdown to plain text suitable for a terminal."""
    text = text.replace("\\|", "|")
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = text.replace("**", "")
    text = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", text)
    return text.strip()


def _anchor(heading: str) -> str:
    """Approximate GitHub's heading-to-anchor slug."""
    slug = heading.strip().lower().replace("`", "")
    slug = re.sub(r"[^\w\s-]", "", slug)
    return re.sub(r"\s+", "-", slug)


def _split_row(line: str) -> list[str]:
    """Split a Markdown table row into stripped cells (on *unescaped* pipes).

    Cells may contain ``\\|`` to mean a literal pipe (e.g. enum types like
    ``virus \\| label``); those must not be treated as column delimiters.
    """
    inner = line.strip()
    inner = inner.removeprefix("|").removesuffix("|")
    return [cell.strip() for cell in re.split(r"(?<!\\)\|", inner)]


def _is_separator(cells: list[str]) -> bool:
    """True for a ``|---|---|`` style separator row."""
    return all(cell and set(cell) <= {"-", ":", " "} for cell in cells)


# ── Parsing ─────────────────────────────────────────────────────────────


def parse_docs(path: Path = DEFAULT_DOCS_PATH) -> DocIndex:
    """Parse the configuration reference into a :class:`DocIndex`."""
    fields: dict[str, FieldDoc] = {}
    sections: list[str] = []
    section_anchor: dict[str, str] = {}
    current_section = ""
    current_anchor = ""
    in_field_table = False

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()

        heading = re.match(r"(#+)\s+(.*)", line)
        if heading:
            level, text = len(heading.group(1)), heading.group(2).strip()
            current_anchor = _anchor(text)
            if level == 2:
                current_section = _strip_md(text)
                if current_section not in section_anchor:
                    sections.append(current_section)
                    section_anchor[current_section] = current_anchor
            in_field_table = False
            continue

        if not line.strip():
            in_field_table = False
            continue

        if not line.lstrip().startswith("|"):
            continue

        cells = _split_row(line)
        if _is_separator(cells):
            continue
        if cells[0] == "Key":  # header row of a field table
            in_field_table = True
            continue
        if not in_field_table or not cells[0].startswith("`"):
            continue

        key = cells[0].strip("`").strip()
        if len(cells) >= 4:
            type_, default, meaning = cells[1], cells[2], "|".join(cells[3:])
        elif len(cells) == 3:
            type_, default, meaning = cells[1], cells[2], ""
        else:  # 2-column: Key | Meaning
            type_, default, meaning = "", "", cells[1]

        fields[key] = FieldDoc(
            key=key,
            type=_strip_md(type_),
            default=_strip_md(default),
            meaning=_strip_md(meaning),
            section=current_section,
            anchor=current_anchor,
        )

    return DocIndex(fields=fields, sections=sections, section_anchor=section_anchor)


def flatten_config(data: Any, prefix: str = "") -> dict[str, Any]:
    """Flatten nested config dicts to dotted keys → leaf value."""
    out: dict[str, Any] = {}
    if isinstance(data, dict):
        for key, value in data.items():
            dotted = f"{prefix}.{key}" if prefix else str(key)
            if isinstance(value, dict):
                out.update(flatten_config(value, dotted))
            else:
                out[dotted] = value
    return out


def load_config_values(path: Path = DEFAULT_CONFIG_PATH) -> dict[str, Any]:
    """Load and flatten ``config.yaml``; empty dict if it cannot be read."""
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return {}
    return flatten_config(data)


# ── Value resolution + formatting ───────────────────────────────────────


def _lookup_value(doc_key: str, flat: dict[str, Any]) -> Any | None:
    """Resolve a doc key's live value by dotted-suffix match; None if ambiguous."""
    if doc_key in flat:
        return flat[doc_key]
    suffix = f".{doc_key}"
    matches = [value for path, value in flat.items() if path.endswith(suffix)]
    return matches[0] if len(matches) == 1 else None


def _fmt(value: Any) -> str:
    """Format a config value the way it reads in YAML."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return ", ".join(str(item) for item in value)
    return str(value)


# ── Rendering ───────────────────────────────────────────────────────────


def _render_field(fd: FieldDoc, flat: dict[str, Any]) -> str:
    """Format a single field's help block."""
    signature_parts = [part for part in (fd.type,) if part]
    if fd.default:
        signature_parts.append(f"default: {fd.default}")
    header = fd.key + (f"  ({', '.join(signature_parts)})" if signature_parts else "")

    out = [header]

    value = _lookup_value(fd.key, flat)
    if (
        value is not None
        and not isinstance(value, list)
        and fd.default
        and _fmt(value) != fd.default
    ):
        out.append(f"  set to: {_fmt(value)}  (differs from default)")

    if fd.meaning:
        out.extend("  " + line for line in textwrap.wrap(fd.meaning, _WRAP_WIDTH))

    out.append(f"  Full reference: docs/configuration.md#{fd.anchor}")
    return "\n".join(out)


def _render_list(index: DocIndex) -> str:
    """Format the all-fields overview grouped by section."""
    by_section: dict[str, list[FieldDoc]] = {name: [] for name in index.sections}
    for fd in index.fields.values():
        by_section.setdefault(fd.section, []).append(fd)

    out = ["RetroSeek configuration fields  (source: docs/configuration.md)", ""]
    for name in index.sections:
        out.append(name)
        section_fields = by_section.get(name, [])
        if not section_fields:
            out.append("    (structured section — see docs/configuration.md)")
        for fd in section_fields:
            summary = fd.meaning.split(". ")[0]
            summary = textwrap.shorten(summary, _LIST_MEANING_WIDTH, placeholder="…")
            out.append(f"    {fd.key:<{_LIST_KEY_WIDTH}} {summary}".rstrip())
        out.append("")
    out.append("Run 'RetroSeek --config-help KEY' for a single field.")
    return "\n".join(out)


def _resolve(key: str, index: DocIndex) -> FieldDoc | None:
    """Find a field by exact, case-insensitive, or dotted-suffix match."""
    if key in index.fields:
        return index.fields[key]
    lowered = {name.lower(): name for name in index.fields}
    if key.lower() in lowered:
        return index.fields[lowered[key.lower()]]
    suffix = f".{key}"
    candidates = [name for name in index.fields if name.endswith(suffix)]
    if len(candidates) == 1:
        return index.fields[candidates[0]]
    return None


def _render_unknown(key: str, index: DocIndex) -> str:
    """Format the not-found message, with close suggestions."""
    suffix = f".{key}"
    candidates = sorted(name for name in index.fields if name.endswith(suffix))
    if candidates:
        joined = ", ".join(candidates)
        return (
            f"Ambiguous config key: {key!r}. Did you mean one of: {joined}?\n"
            "Run 'RetroSeek --config-help' to list all fields."
        )

    near = difflib.get_close_matches(key, list(index.fields), n=5, cutoff=0.4)
    message = f"Unknown config key: {key!r}."
    if near:
        message += " Did you mean: " + ", ".join(near) + "?"
    message += "\nRun 'RetroSeek --config-help' to list all fields."
    return message


def render(
    key: str | None,
    docs_path: Path = DEFAULT_DOCS_PATH,
    config_path: Path = DEFAULT_CONFIG_PATH,
) -> str:
    """Render help text for ``key`` (or the full overview when ``key`` is None)."""
    index = parse_docs(docs_path)
    if key is None:
        return _render_list(index)

    fd = _resolve(key, index)
    if fd is not None:
        return _render_field(fd, load_config_values(config_path))

    if key in index.section_anchor:
        return (
            f"{key}\n"
            "  Structured section — no per-field entries.\n"
            f"  Full reference: docs/configuration.md#{index.section_anchor[key]}"
        )

    return _render_unknown(key, index)
