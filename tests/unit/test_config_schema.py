"""The configs the repo ships validate against ``data/config/schema.yaml``.

The launcher's preflight rejects a config the schema does not accept, so a
committed example that fails it would fail for anyone who copies it. The
test fixture once kept a ``domains:`` block for months after the schema
dropped it (the curated Pfam class table replaced it).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yamale

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA = REPO_ROOT / "data" / "config" / "schema.yaml"


@pytest.mark.parametrize(
    "config",
    ["data/config/config.yaml", "tests/fixtures/example_test_config.yaml"],
)
def test_shipped_config_matches_the_schema(config: str) -> None:
    schema = yamale.make_schema(str(SCHEMA))
    data = yamale.make_data(str(REPO_ROOT / config))
    yamale.validate(schema, data)  # raises YamaleError naming each bad field


@pytest.mark.parametrize(("release", "valid"), [("38.2", True), ("38x2", False)])
def test_pfam_release_must_be_major_dot_minor(
    release: str, valid: bool, tmp_path: Path
) -> None:
    config = tmp_path / "config.yaml"
    text = (REPO_ROOT / "data/config/config.yaml").read_text(encoding="utf-8")
    assert "pfam_release: '38.2'" in text  # the line this test rewrites
    config.write_text(
        text.replace("pfam_release: '38.2'", f"pfam_release: '{release}'"),
        encoding="utf-8",
    )
    schema = yamale.make_schema(str(SCHEMA))
    data = yamale.make_data(str(config))
    if valid:
        yamale.validate(schema, data)
    else:
        with pytest.raises(yamale.YamaleError, match="pfam_release"):
            yamale.validate(schema, data)
