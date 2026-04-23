"""Characterization tests for utils.py.

Pin down the current observable behavior of ``utils`` before refactoring
``os.path`` usage to :class:`pathlib.Path`. Every function must continue
to accept ``str`` arguments (the Snakefile passes string paths), with
:class:`Path` also accepted as an ergonomic upgrade.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import dill
import pytest

import utils


class TestPicklerUnpicklerRoundtrip:
    """``pickler`` writes what ``unpickler`` reads — for both str and Path."""

    def test_string_path_roundtrip(self, tmp_path: Path) -> None:
        """Canonical caller: Snakefile-style string paths."""
        payload = {"a": 1, "b": [1, 2, 3]}
        utils.pickler(payload, str(tmp_path), "payload.pkl")
        loaded = utils.unpickler(str(tmp_path), "payload.pkl")
        assert loaded == payload

    def test_path_object_roundtrip(self, tmp_path: Path) -> None:
        """Ergonomic caller: pass Path objects directly."""
        payload = [1, 2, 3]
        utils.pickler(payload, tmp_path, "list.pkl")
        loaded = utils.unpickler(tmp_path, "list.pkl")
        assert loaded == payload

    def test_uses_dill_not_stdlib_pickle(self, tmp_path: Path) -> None:
        """Pickles must be dill-compatible (production code relies on dill)."""
        payload = lambda x: x + 1  # noqa: E731 — lambdas require dill, not pickle
        utils.pickler(payload, tmp_path, "lam.pkl")
        with (tmp_path / "lam.pkl").open("rb") as handle:
            reloaded = dill.load(handle)
        assert reloaded(41) == 42

    def test_unpickler_missing_file_raises(self, tmp_path: Path) -> None:
        """Missing file → Exception with path in message (current contract)."""
        with pytest.raises(Exception, match="Failed to unpickle"):
            utils.unpickler(tmp_path, "does-not-exist.pkl")


class TestDirectoryFileRetriever:
    """Lists files (not subdirectories) in a directory."""

    def test_returns_file_basenames(self, tmp_path: Path) -> None:
        """Returns file names (not subdirectory names)."""
        (tmp_path / "a.txt").write_text("hi")
        (tmp_path / "b.txt").write_text("hi")
        (tmp_path / "subdir").mkdir()
        (tmp_path / "subdir" / "nested.txt").write_text("hi")

        result = utils.directory_file_retriever(tmp_path)
        assert sorted(result) == ["a.txt", "b.txt"]

    def test_accepts_string_path(self, tmp_path: Path) -> None:
        """Caller can pass a ``str`` path."""
        (tmp_path / "only.txt").write_text("hi")
        result = utils.directory_file_retriever(str(tmp_path))
        assert result == ["only.txt"]

    def test_empty_directory_returns_empty_list(self, tmp_path: Path) -> None:
        """Empty dir → empty list, not None."""
        assert utils.directory_file_retriever(tmp_path) == []


class TestDirectoryContentEraser:
    """Deletes files in the given directory; leaves subdirs alone."""

    def test_removes_files(self, tmp_path: Path) -> None:
        """Plain files under the dir are removed."""
        (tmp_path / "doomed.txt").write_text("bye")
        (tmp_path / "alsodoomed.txt").write_text("bye")
        utils.directory_content_eraser(tmp_path)
        assert list(tmp_path.iterdir()) == []

    def test_leaves_subdirectories_intact(self, tmp_path: Path) -> None:
        """Subdirectories are not removed (current contract)."""
        (tmp_path / "keep").mkdir()
        (tmp_path / "doomed.txt").write_text("bye")
        utils.directory_content_eraser(tmp_path)
        remaining = [p.name for p in tmp_path.iterdir()]
        assert "keep" in remaining
        assert "doomed.txt" not in remaining

    def test_accepts_string_path(self, tmp_path: Path) -> None:
        """String path input works."""
        (tmp_path / "x.txt").write_text("x")
        utils.directory_content_eraser(str(tmp_path))
        assert list(tmp_path.iterdir()) == []


class TestIncompleteDictCleaner:
    """Filters an object dict by ``.is_complete()``."""

    def test_keeps_complete_drops_incomplete(self) -> None:
        """Only records where ``is_complete()`` is truthy survive."""
        complete = MagicMock(); complete.is_complete.return_value = True
        incomplete = MagicMock(); incomplete.is_complete.return_value = False
        result = utils.incomplete_dict_cleaner(
            {"keep": complete, "drop": incomplete}
        )
        assert "keep" in result
        assert "drop" not in result

    def test_empty_in_empty_out(self) -> None:
        """No inputs → no outputs."""
        assert utils.incomplete_dict_cleaner({}) == {}


class TestRandomStringGenerator:
    """6-char identifiers; upper-case ASCII + digits only."""

    def test_length_matches_requested(self) -> None:
        """Output has the requested character count."""
        assert len(utils.random_string_generator(6)) == 6
        assert len(utils.random_string_generator(12)) == 12

    def test_only_upper_ascii_and_digits(self) -> None:
        """Output uses only ``A-Z`` and ``0-9``."""
        allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
        for _ in range(100):
            assert set(utils.random_string_generator(8)) <= allowed
