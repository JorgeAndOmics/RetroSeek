# ruff: noqa: PLC0415
# Per-test imports of `seq_utils` are deliberate (conftest defaults
# stub must take effect first).

"""Unit tests for seq_utils.

Focus: the Entrez kwargs that ``gb_fetcher`` builds are free of the
deprecated ``expand_by`` variable. The old code tried to pad the hit range
via ``instance.HSP.sbjct_start + expand_by`` where ``expand_by`` was never
defined — any call with a populated HSP raised :class:`NameError`, which
was then swallowed by the retry/``except Exception`` block and counted as
"fetch failed" instead of surfacing the real error. Range expansion belongs
in ``ranges_analysis.R`` via ``GenomicRanges::resize()`` — not here.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from RetroSeeker_class import RetroSeeker


def _make_instance_with_hsp(*, start: int = 100, end: int = 500) -> RetroSeeker:
    """Return a RetroSeeker whose ``HSP.sbjct_*`` attributes are populated."""
    rs = RetroSeeker(
        label="toy_label",
        virus="toy_virus",
        abbreviation="TOY",
        species="Toyus_simplex",
        probe="POL",
        accession="TOY_POL_001",
        identifier="seq001",
    )
    hsp = MagicMock()
    hsp.sbjct_start = start
    hsp.sbjct_end = end
    hsp.frame = (1, 0)  # so __set_HSP doesn't fail on strand derivation
    rs.set_HSP(hsp)
    return rs


class TestGbFetcherKwargs:
    """What kwargs does gb_fetcher actually hand to Entrez.efetch?"""

    def test_kwargs_include_seq_start_and_seq_stop_unchanged(self) -> None:
        """With an HSP, gb_fetcher passes the raw sbjct coordinates.

        Not ``sbjct_start + expand_by``, not ``sbjct_end + expand_by`` — the
        expansion feature is deprecated and its implementation was broken
        (``expand_by`` was an undefined name).
        """
        from seq_utils import gb_fetcher

        instance = _make_instance_with_hsp(start=100, end=500)

        # Patch Entrez.efetch so the test doesn't hit the network. We record
        # what kwargs the real gb_fetcher handed in.
        fake_handle = MagicMock()
        fake_handle.__enter__ = MagicMock(return_value=fake_handle)
        fake_handle.__exit__ = MagicMock(return_value=False)
        fake_handle.read = MagicMock(return_value="LOCUS dummy")

        with (
            patch("seq_utils.Entrez.efetch", return_value=fake_handle) as efetch,
            patch.object(instance, "set_genbank"),
        ):
            gb_fetcher(instance=instance, online_database="nuccore")

        assert efetch.called, "gb_fetcher should have called Entrez.efetch"
        kwargs = efetch.call_args.kwargs
        assert kwargs["seq_start"] == 100, (
            f"expected raw sbjct_start=100, got {kwargs.get('seq_start')!r} — "
            "are we still applying a deprecated expand_by offset?"
        )
        assert kwargs["seq_stop"] == 500
        assert kwargs["id"] == "TOY_POL_001"
        assert kwargs["rettype"] == "gb"
        assert kwargs["retmode"] == "text"
        assert kwargs["db"] == "nuccore"

    def test_kwargs_omit_seq_range_when_no_hsp(self) -> None:
        """Without an HSP, gb_fetcher should not pass seq_start / seq_stop."""
        from seq_utils import gb_fetcher

        rs = RetroSeeker(
            label="l",
            virus="v",
            abbreviation="A",
            species="S",
            probe="P",
            accession="NOHSP_001",
            identifier="noh01",
        )
        assert rs.HSP is None

        fake_handle = MagicMock()
        fake_handle.__enter__ = MagicMock(return_value=fake_handle)
        fake_handle.__exit__ = MagicMock(return_value=False)
        fake_handle.read = MagicMock(return_value="LOCUS dummy")

        with (
            patch("seq_utils.Entrez.efetch", return_value=fake_handle) as efetch,
            patch.object(rs, "set_genbank"),
        ):
            gb_fetcher(instance=rs, online_database="nuccore")

        kwargs = efetch.call_args.kwargs
        assert "seq_start" not in kwargs
        assert "seq_stop" not in kwargs

    def test_sbjct_start_below_one_is_clamped(self) -> None:
        """``max(1, sbjct_start)`` guard on the lower bound is preserved."""
        from seq_utils import gb_fetcher

        instance = _make_instance_with_hsp(start=-10, end=500)

        fake_handle = MagicMock()
        fake_handle.__enter__ = MagicMock(return_value=fake_handle)
        fake_handle.__exit__ = MagicMock(return_value=False)
        fake_handle.read = MagicMock(return_value="LOCUS dummy")

        with (
            patch("seq_utils.Entrez.efetch", return_value=fake_handle) as efetch,
            patch.object(instance, "set_genbank"),
        ):
            gb_fetcher(instance=instance, online_database="nuccore")

        assert efetch.call_args.kwargs["seq_start"] == 1


class TestGbFetcherNoUndefinedNameError:
    """Regression guard: gb_fetcher must not reference undefined variables.

    A NameError was previously raised inside the ``if instance.HSP:`` block,
    caught by the bare ``except Exception``, retried, and then silently
    returned the instance without a GenBank record. That masked the real
    error, which the following test ensures stays dead.
    """

    def test_no_nameerror_raised_on_hsp_path(self) -> None:
        """The HSP-populated path completes without NameError.

        We patch Entrez to succeed; the only way NameError could surface is
        from the function body itself.
        """
        from seq_utils import gb_fetcher

        instance = _make_instance_with_hsp(start=50, end=150)

        fake_handle = MagicMock()
        fake_handle.__enter__ = MagicMock(return_value=fake_handle)
        fake_handle.__exit__ = MagicMock(return_value=False)
        fake_handle.read = MagicMock(return_value="LOCUS dummy")

        with (
            patch("seq_utils.Entrez.efetch", return_value=fake_handle),
            patch.object(instance, "set_genbank"),
        ):
            result = gb_fetcher(instance=instance, online_database="nuccore")

        # Successful path returns the (mutated) instance.
        assert result is instance
