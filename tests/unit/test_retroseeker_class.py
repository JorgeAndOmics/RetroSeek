"""Unit tests for RetroSeeker_class.

These tests pin down the invariants the class is supposed to uphold:

- Identifier-based equality: two RetroSeekers with the same ``identifier``
  compare equal regardless of other fields. The current implementation
  checks ``isinstance(other, Object)`` against an undefined ``Object`` — so
  ``a == b`` raises :class:`NameError` the moment either side is a real
  :class:`RetroSeeker`. The fix is to compare against ``RetroSeeker``.
- Hash consistency: objects that compare equal must hash equal (Python
  contract required for dict / set correctness).
- Non-RetroSeeker comparands return ``False``, never raise.
"""

from __future__ import annotations

import pytest

# ``conftest.py`` puts ``workflow/scripts`` on sys.path and stubs ``defaults``
# before this import runs, so the real RetroSeeker_class can be imported
# without triggering defaults.py's filesystem side effects.
from RetroSeeker_class import RetroSeeker


def _make(identifier: str, label: str = "toy_label") -> RetroSeeker:
    """Build a minimally-populated RetroSeeker for comparison tests."""
    return RetroSeeker(
        label=label,
        virus="toy_virus",
        abbreviation="TOY",
        species="Toyus_simplex",
        probe="GAG",
        accession="TOY_001",
        identifier=identifier,
    )


class TestRetroSeekerEquality:
    """Identifier-based equality invariants."""

    def test_same_identifier_compare_equal(self) -> None:
        """Two RetroSeekers with the same identifier are equal."""
        a = _make(identifier="abc123")
        b = _make(identifier="abc123", label="different_label")
        assert a == b

    def test_different_identifier_compare_unequal(self) -> None:
        """Two RetroSeekers with different identifiers are not equal."""
        a = _make(identifier="abc123")
        b = _make(identifier="xyz789")
        assert a != b

    def test_equality_is_symmetric(self) -> None:
        """``a == b`` implies ``b == a``."""
        a = _make(identifier="sym01")
        b = _make(identifier="sym01")
        assert (a == b) == (b == a)

    def test_non_retroseeker_comparand_returns_false(self) -> None:
        """Comparing against a non-RetroSeeker returns False, never raises."""
        rs = _make(identifier="solo")
        assert (rs == "not a retroseeker") is False
        assert (rs == 42) is False
        assert (rs == None) is False  # noqa: E711 — explicit identity-check intent


class TestRetroSeekerHashContract:
    """Hash / equality contract required for dict and set usage."""

    def test_equal_objects_have_equal_hashes(self) -> None:
        """Python's contract: ``a == b`` implies ``hash(a) == hash(b)``."""
        a = _make(identifier="hash01")
        b = _make(identifier="hash01", label="different")
        assert hash(a) == hash(b)

    def test_retroseeker_usable_as_dict_key(self) -> None:
        """RetroSeekers can be used as dict keys; same identifier collapses."""
        a = _make(identifier="key01")
        b = _make(identifier="key01", label="different")
        d: dict[RetroSeeker, int] = {a: 1}
        d[b] = 2
        assert len(d) == 1
        assert d[a] == 2

    def test_retroseeker_usable_in_set(self) -> None:
        """RetroSeekers deduplicate in a set by identifier."""
        a = _make(identifier="set01")
        b = _make(identifier="set01", label="different")
        c = _make(identifier="set02")
        assert len({a, b, c}) == 2


class TestRetroSeekerRepr:
    """The string form is stable and contains the identifying fields."""

    def test_repr_is_comma_separated_of_core_fields(self) -> None:
        """``__repr__`` returns a comma-joined summary of core attributes."""
        rs = _make(identifier="rep01", label="repr_lbl")
        text = repr(rs)
        for token in ["repr_lbl", "toy_virus", "TOY", "Toyus_simplex", "GAG", "TOY_001", "rep01"]:
            assert token in text, f"expected {token!r} in repr, got {text!r}"


class TestRetroSeekerIsComplete:
    """``is_complete`` should require all three of alignment, HSP, genbank."""

    def test_is_complete_false_when_all_none(self) -> None:
        """Default-constructed RetroSeeker is not complete."""
        rs = _make(identifier="inc01")
        assert rs.is_complete() is False

    def test_is_complete_false_when_missing_any(self) -> None:
        """Missing any one of alignment/HSP/genbank keeps the record incomplete."""
        rs = _make(identifier="inc02")
        rs.alignment = object()
        rs.HSP = object()
        # ``genbank`` still None.
        assert rs.is_complete() is False

    def test_is_complete_true_when_all_three_set(self) -> None:
        """All three of alignment/HSP/genbank populated ⇒ complete."""
        rs = _make(identifier="cmp01")
        rs.alignment = object()
        rs.HSP = object()
        rs.genbank = object()
        assert rs.is_complete() is True


# Sanity check the test plumbing itself: the stub ``defaults`` was applied
# before ``RetroSeeker_class`` imported, so RetroSeeker is the real one.
def test_imported_class_is_the_real_retroseeker() -> None:
    """Guard against an accidental mock leak from conftest."""
    rs = _make(identifier="id01")
    assert rs.__class__.__name__ == "RetroSeeker"
    assert rs.__class__.__module__ == "RetroSeeker_class"
