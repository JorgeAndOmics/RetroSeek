# ruff: noqa: PLC0415
# Per-test imports of `validator` are deliberate: the conftest defaults
# stub must take effect before the module under test is imported, and
# every test then triggers a fresh import to avoid sharing global state
# between cases.

"""Unit tests for validator.py.

Covers two prior defects:

1. ``green_light`` handled invalid user input by recursing but not
   returning — the top-level caller then saw the implicit ``None`` and
   treated it as "abort", even though the user may have simply mistyped.
2. The module used to self-import (``import validator``) at the top,
   which is pointless and confusing. The import should be gone.

Also verifies the positive / negative user-input paths and the
``all_valid=False`` short-circuit, so future refactors don't regress.
"""

from __future__ import annotations

import importlib
from unittest.mock import patch


class TestGreenLight:
    """Interactive confirmation behaviour."""

    def test_returns_false_when_validation_failed(self) -> None:
        """If ``all_valid`` is False, green_light returns False without prompting."""
        import validator as v

        with patch("builtins.input") as mocked_input:
            result = v.green_light(all_valid=False)
        assert result is False
        mocked_input.assert_not_called()

    def test_returns_true_on_yes(self) -> None:
        """Typing 'Y' proceeds."""
        import validator as v

        with patch("builtins.input", return_value="Y"):
            assert v.green_light(all_valid=True) is True

    def test_returns_true_on_default_empty(self) -> None:
        """Hitting enter (empty string) is treated as 'Y'."""
        import validator as v

        with patch("builtins.input", return_value=""):
            assert v.green_light(all_valid=True) is True

    def test_lowercase_y_also_proceeds(self) -> None:
        """'y' is accepted (case-insensitive)."""
        import validator as v

        with patch("builtins.input", return_value="y"):
            assert v.green_light(all_valid=True) is True

    def test_returns_false_on_no(self) -> None:
        """'N' aborts."""
        import validator as v

        with patch("builtins.input", return_value="N"):
            assert v.green_light(all_valid=True) is False

    def test_invalid_then_yes_returns_true(self) -> None:
        """Invalid input re-prompts; valid response on the retry is honoured.

        Regression: previously the recursive call's return was discarded,
        so the function fell through to an implicit ``None`` regardless of
        the second prompt. Now it must forward the recursive result.
        """
        import validator as v

        with patch("builtins.input", side_effect=["maybe", "Y"]):
            result = v.green_light(all_valid=True)
        assert result is True

    def test_invalid_then_no_returns_false(self) -> None:
        """Invalid input → retry → 'N' → returns False."""
        import validator as v

        with patch("builtins.input", side_effect=["blarg", "N"]):
            assert v.green_light(all_valid=True) is False


class TestNoSelfImport:
    """The module must not import itself at the top level."""

    def test_module_has_no_self_import(self) -> None:
        """Regression guard: ``import validator`` shouldn't live inside validator.py."""

        source = importlib.util.find_spec("validator")
        assert source is not None, "validator module should be findable on sys.path"
        from pathlib import Path

        text = Path(source.origin).read_text(encoding="utf-8")
        lines = [line.strip() for line in text.splitlines()]
        assert "import validator" not in lines, (
            "validator.py should not import itself — remove the dead import"
        )
