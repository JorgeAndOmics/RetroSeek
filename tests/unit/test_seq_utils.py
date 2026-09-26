# ruff: noqa: PLC0415
# Per-test imports of `seq_utils` are deliberate (conftest defaults
# stub must take effect first).

"""Unit tests for seq_utils.

Focus: the Entrez kwargs that ``gb_fetcher`` builds are free of the
deprecated ``expand_by`` variable. The old code tried to pad the hit range
via ``instance.HSP.sbjct_start + expand_by`` where ``expand_by`` was never
defined - any call with a populated HSP raised :class:`NameError`, which
was then swallowed by the retry/``except Exception`` block and counted as
"fetch failed" instead of surfacing the real error. Range expansion belongs
in ``ranges_analysis.R`` via ``GenomicRanges::resize()`` - not here.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
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

        Not ``sbjct_start + expand_by``, not ``sbjct_end + expand_by`` - the
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
            f"expected raw sbjct_start=100, got {kwargs.get('seq_start')!r} - "
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


class TestBlastFailuresStopTheJob:
    """A failed BLAST used to become a silently missing probe (ADR-021).

    ``blaster`` returned None on failure, ``blaster_parser`` swallowed every
    exception, and ``_blast_task`` caught the rest, so a genome's pickle could lack
    a probe's hits with nothing but a log line to show for it.
    """

    @staticmethod
    def _fake_blast(tmp_path, body: str):
        import stat

        tool = tmp_path / "fake_tblastn"
        tool.write_text(f"#!/bin/sh\n{body}\n")
        tool.chmod(tool.stat().st_mode | stat.S_IEXEC)
        return str(tool)

    @staticmethod
    def _instance(tmp_path):
        query = tmp_path / "q.fa"
        query.write_text(">q\nMKV\n")
        instance = MagicMock()
        instance.get_fasta.return_value = str(query)
        instance.probe = "POL"
        return instance

    def test_a_failing_blast_raises(self, tmp_path) -> None:
        import pytest

        import seq_utils
        from log import PipelineError

        tool = self._fake_blast(tmp_path, "echo 'BLAST Database error' >&2; exit 3")
        with pytest.raises(PipelineError, match="exit code 3"):
            seq_utils.blaster(self._instance(tmp_path), tool, tmp_path, "Toyus", 1)

    def test_an_empty_blast_output_raises(self, tmp_path) -> None:
        """outfmt 11 (ASN.1) is never empty, even with no hits: empty means broken."""
        import pytest

        import seq_utils
        from log import PipelineError

        tool = self._fake_blast(tmp_path, "exit 0")
        with pytest.raises(PipelineError, match="no output"):
            seq_utils.blaster(self._instance(tmp_path), tool, tmp_path, "Toyus", 1)

    def test_an_unreadable_archive_raises(self) -> None:
        import pytest

        import seq_utils
        from log import PipelineError

        with pytest.raises(PipelineError, match="blast_formatter"):
            seq_utils.blaster_parser("not an ASN.1 archive", MagicMock(), "Toyus")

    def test_blast_task_does_not_swallow_failures(self) -> None:
        import pytest

        import seq_utils
        from log import PipelineError

        with (
            patch.object(seq_utils, "blaster", side_effect=PipelineError("broken")),
            pytest.raises(PipelineError),
        ):
            seq_utils._blast_task(MagicMock(), "tblastn", "Toyus", "/db", 1)


class TestBlasterParserReadsEveryHsp:
    """blaster_parser turns every HSP of every hit into one RetroSeeker.

    blast_formatter is replaced by a stand-in that returns a small fixed XML
    report, so the test needs no BLAST install and no real archive.
    """

    XML = Path(__file__).resolve().parents[1] / "fixtures" / "blast" / "two_hits.xml"

    def test_one_object_per_hsp_with_query_metadata(self) -> None:
        import seq_utils

        query = RetroSeeker(
            label="ALV",
            virus="Avian leukosis virus",
            abbreviation="ALV",
            species=None,
            probe=" POL ",
            accession="Q1",
            identifier="x",
        )
        report = SimpleNamespace(stdout=self.XML.read_text())
        with patch.object(seq_utils, "run_tool", return_value=report):
            hits = seq_utils.blaster_parser("archive", query, "Toyus_toyus")

        by_key = sorted(hits.items())
        assert [key.split("-")[0] for key, _ in by_key] == ["CM1.1", "CM1.1", "CM2.1"]
        assert all(len(key.split("-")[1]) == 6 for key, _ in by_key)
        rows = sorted(
            (o.accession, o.HSP.sbjct_start, o.strand, o.species, o.probe, o.label)
            for _, o in by_key
        )
        assert rows == [
            ("CM1.1", 100, "+", "Toyus_toyus", "POL", "ALV"),
            ("CM1.1", 900, "-", "Toyus_toyus", "POL", "ALV"),
            ("CM2.1", 10, "+", "Toyus_toyus", "POL", "ALV"),
        ]
        assert all(o.alignment.hit_def.startswith(o.accession) for o in hits.values())

    def test_a_repeated_random_identifier_does_not_drop_a_hit(self) -> None:
        """Keys are {accession}-{6 random characters}; a repeated draw on one
        accession used to replace the earlier hit without a word."""
        import seq_utils
        import utils

        query = RetroSeeker(
            label="ALV",
            virus="Avian leukosis virus",
            abbreviation="ALV",
            species=None,
            probe="POL",
            accession="Q1",
            identifier="x",
        )
        report = SimpleNamespace(stdout=self.XML.read_text())
        draws = iter(["AAAAAA", "AAAAAA", "BBBBBB", "CCCCCC"])
        with (
            patch.object(seq_utils, "run_tool", return_value=report),
            patch.object(utils, "random_string_generator", lambda n: next(draws)),
        ):
            hits = seq_utils.blaster_parser("archive", query, "Toyus_toyus")
        assert sorted(hits) == ["CM1.1-AAAAAA", "CM1.1-BBBBBB", "CM2.1-CCCCCC"]
        assert all(key.endswith(o.identifier) for key, o in hits.items())

    def test_identifiers_stay_unique_across_the_probes_of_a_genome(self) -> None:
        """Every probe of a genome shares one key space: two probes hitting the
        same chromosome must not reuse an identifier either."""
        import seq_utils
        import utils

        def probe(name: str) -> RetroSeeker:
            return RetroSeeker(
                label=name,
                virus=name,
                abbreviation=name,
                species=None,
                probe="POL",
                accession="Q",
                identifier=name,
            )

        report = SimpleNamespace(stdout=self.XML.read_text())
        draws = iter(["A", "B", "C", "A", "B", "D", "E", "F"])
        with (
            patch.object(seq_utils, "blaster", return_value="archive"),
            patch.object(seq_utils, "run_tool", return_value=report),
            patch.object(utils, "random_string_generator", lambda n: next(draws)),
        ):
            hits = seq_utils.blast_executor(
                {"p1": probe("P1"), "p2": probe("P2")}, "tblastn", "db", 1, "Toyus"
            )
        identifiers = [o.identifier for o in hits.values()]
        assert len(hits) == 6
        assert sorted(identifiers) == ["A", "B", "C", "D", "E", "F"]
