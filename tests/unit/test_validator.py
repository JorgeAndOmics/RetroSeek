# ruff: noqa: PLC0415
# Per-test imports of `validator` are deliberate: the conftest defaults
# stub must take effect before the module under test is imported, and
# every test then triggers a fresh import to avoid sharing global state
# between cases.

"""Unit tests for validator.py.

Covers two prior defects:

1. ``green_light`` handled invalid user input by recursing but not
   returning - the top-level caller then saw the implicit ``None`` and
   treated it as "abort", even though the user may have simply mistyped.
2. The module used to self-import (``import validator``) at the top,
   which is pointless and confusing. The import should be gone.

Also verifies the positive / negative user-input paths and the
``all_valid=False`` short-circuit, so future refactors don't regress.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path
from unittest.mock import patch

import pytest


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
        """Invalid input -> retry -> 'N' -> returns False."""
        import validator as v

        with patch("builtins.input", side_effect=["blarg", "N"]):
            assert v.green_light(all_valid=True) is False


class TestUnattendedRuns:
    """Prompts must not abort a run that has nobody to answer them.

    Regression: ``validate_ncbi_key`` and ``green_light`` called bare
    ``input()``. With no terminal attached (CI, an agent, ``nohup``) that
    raises ``EOFError`` and killed the whole pipeline before Snakemake was
    ever reached.
    """

    def test_ask_returns_the_typed_answer(self) -> None:
        """With a terminal present, the answer is used, not the default."""
        import validator as v

        with patch("builtins.input", return_value="N"):
            assert v.ask("Proceed [Y/n]: ", default="Y") == "N"

    def test_ask_falls_back_to_default_without_stdin(self) -> None:
        """No stdin means the default answer, not an exception."""
        import validator as v

        with patch("builtins.input", side_effect=EOFError):
            assert v.ask("Proceed [Y/n]: ", default="Y") == "Y"

    def test_green_light_proceeds_without_stdin(self) -> None:
        """An unattended run continues once validation has passed."""
        import validator as v

        with patch("builtins.input", side_effect=EOFError):
            assert v.green_light(all_valid=True) is True

    def test_green_light_still_aborts_when_validation_failed(self) -> None:
        """The fallback must not turn a failed validation into a go-ahead."""
        import validator as v

        with patch("builtins.input", side_effect=EOFError):
            assert v.green_light(all_valid=False) is False

    def test_ncbi_key_prompt_is_skipped_without_stdin(self) -> None:
        """A missing API key downgrades to a warning instead of crashing."""
        import validator as v

        with (
            patch.dict("os.environ", {}, clear=True),
            patch("builtins.input", side_effect=EOFError),
        ):
            v.validate_ncbi_key()
            assert "NCBI_API_KEY" not in os.environ


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
            "validator.py should not import itself - remove the dead import"
        )


class TestPreflight:
    """The fast checks that always run, whatever -skp says (ADR-020).

    They catch, in seconds, what used to fail hours into a run: a tool missing
    from the environment, or a Pfam library older than the curated table.
    """

    @staticmethod
    def _hmm(path: Path, *accessions: str) -> Path:
        path.write_text(
            "".join(f"HMMER3/f\nNAME  m{a}\nACC   {a}.1\n//\n" for a in accessions)
        )
        return path

    @staticmethod
    def _classes(path: Path, *accessions: str) -> Path:
        rows = "".join(f"{a}\tname\tother\n" for a in accessions)
        path.write_text("pfam_acc\tpfam_name\tclass\n" + rows)
        return path

    def test_missing_tools_names_only_the_absent_ones(self) -> None:
        import validator as v

        assert v.missing_tools(["sh", "retroseek-no-such-tool"]) == [
            "retroseek-no-such-tool"
        ]

    def test_pfam_problem_is_none_when_the_library_has_every_family(
        self, tmp_path: Path
    ) -> None:
        import validator as v

        hmm = self._hmm(tmp_path / "Pfam-A.hmm", "PF00001", "PF00002")
        classes = self._classes(tmp_path / "c.tsv", "PF00001")
        assert v.pfam_problem(hmm, classes) is None

    def test_pfam_problem_explains_an_old_library(self, tmp_path: Path) -> None:
        import validator as v

        hmm = self._hmm(tmp_path / "Pfam-A.hmm", "PF00001")
        classes = self._classes(tmp_path / "c.tsv", "PF00001", "PF29843")
        problem = v.pfam_problem(hmm, classes)
        assert problem is not None
        assert "PF29843" in problem
        assert "--download-hmm" in problem

    def test_pfam_problem_is_none_before_the_first_download(
        self, tmp_path: Path
    ) -> None:
        """No library yet: the downloader will fetch the pinned release."""
        import validator as v

        classes = self._classes(tmp_path / "c.tsv", "PF00001")
        assert v.pfam_problem(tmp_path / "absent.hmm", classes) is None

    def test_preflight_fails_on_a_missing_tool(self) -> None:
        import stages
        import validator as v

        fake = stages.Stage(
            "--x", "Analysis", ("r",), "x", tools=("retroseek-no-such-tool",)
        )
        with patch.object(v, "yaml_validator", return_value=True):
            assert v.preflight([fake]) is False

    def test_preflight_passes_when_nothing_is_missing(self) -> None:
        import stages
        import validator as v

        fake = stages.Stage("--x", "Analysis", ("r",), "x", tools=("sh",))
        with patch.object(v, "yaml_validator", return_value=True):
            assert v.preflight([fake]) is True


class TestRetiredKeys:
    """A config written for an older RetroSeek must say what replaced each key.

    The schema alone would only say "unexpected key", leaving the user to guess.
    """

    def test_retired_display_switches_name_their_replacement(self) -> None:
        import validator as v

        config = {"display": {"display_snakemake_info": True, "verbosity": "normal"}}
        messages = v.retired_key_messages(config)
        assert len(messages) == 1
        assert "display.display_snakemake_info" in messages[0]
        assert "display.verbosity" in messages[0]

    def test_retired_logging_block_is_named(self) -> None:
        import validator as v

        messages = v.retired_key_messages({"logging": {"level_styles": {}}})
        assert len(messages) == 1
        assert "logging" in messages[0]

    def test_retired_circle_plot_threshold_is_named(self) -> None:
        """The circle-plot stage was removed; its display cutoff went with it."""
        import validator as v

        config = {"plots": {"circle_plot_bitscore_threshold": 0}}
        messages = v.retired_key_messages(config)
        assert len(messages) == 1
        assert "plots.circle_plot_bitscore_threshold" in messages[0]

    def test_a_current_config_has_no_retired_keys(self) -> None:
        import validator as v

        assert v.retired_key_messages({"display": {"verbosity": "quiet"}}) == []

    def test_a_retired_key_set_to_false_still_counts(self) -> None:
        """The template shipped `display_snakemake_info: false`."""
        import validator as v

        assert (
            len(v.retired_key_messages({"display": {"display_snakemake_info": False}}))
            == 1
        )


class TestGenomeFasta:
    """The file validation reads for a genome is the one the normalizer rule links.

    The launcher used to rename `.fasta`/`.fna`/`.fas` files to `.fa` in place
    before every run (replacing an existing `.fa` of the same name). The
    `genome_fasta_normalizer_setup` rule does that job without touching the
    source files, so validation now reads what that rule will link to.
    """

    @staticmethod
    def _fasta(path: Path) -> None:
        path.write_text(">chr1\nACGT\n")

    def test_a_fna_only_genome_is_read_from_the_fna(self, tmp_path: Path) -> None:
        import validator as v

        self._fasta(tmp_path / "Toyus.fna")
        assert v.genome_fasta(tmp_path, "Toyus") == str(tmp_path / "Toyus.fna")

    def test_an_existing_fa_wins(self, tmp_path: Path) -> None:
        import validator as v

        self._fasta(tmp_path / "Toyus.fa")
        self._fasta(tmp_path / "Toyus.fna")
        assert v.genome_fasta(tmp_path, "Toyus") == str(tmp_path / "Toyus.fa")

    def test_a_missing_genome_names_the_fa_path(self, tmp_path: Path) -> None:
        import validator as v

        assert v.genome_fasta(tmp_path, "Toyus") == str(tmp_path / "Toyus.fa")

    def test_ambiguous_genomes_names_only_the_ambiguous_ones(
        self, tmp_path: Path
    ) -> None:
        import validator as v

        self._fasta(tmp_path / "Toyus.fna")
        self._fasta(tmp_path / "Toyus.fasta")
        self._fasta(tmp_path / "Clean.fna")
        problems = v.ambiguous_genomes(tmp_path, ["Toyus", "Clean", "Absent"])
        assert len(problems) == 1
        assert "Toyus" in problems[0]

    def test_a_link_contradicting_its_neighbour_stops_the_preflight(
        self, tmp_path: Path
    ) -> None:
        """The launcher is the only place that can stop a stale link being used.

        Snakemake never reruns the normalizer while `{genome}.fa` exists.
        """
        import validator as v

        self._fasta(tmp_path / "Toyus.fasta")
        self._fasta(tmp_path / "Toyus.fna")
        (tmp_path / "Toyus.fa").symlink_to(tmp_path / "Toyus.fasta")
        with (
            patch.object(v, "yaml_validator", return_value=True),
            patch.object(v.defaults, "SPECIES", ["Toyus"]),
            patch.dict(v.defaults.PATH_DICT, {"SPECIES_DB": tmp_path}),
        ):
            assert v.preflight([]) is False

    def test_the_source_files_are_never_renamed(self, tmp_path: Path) -> None:
        import validator as v

        self._fasta(tmp_path / "Toyus.fasta")
        v.genome_fasta(tmp_path, "Toyus")
        assert sorted(p.name for p in tmp_path.iterdir()) == ["Toyus.fasta"]


def v_probe_names(path: Path) -> set[str] | None:
    import validator as v

    return v.probe_names(path)


class TestGeneLists:
    """The three gene lists are checked before a run (ADR-022).

    ``main_probes`` (which genes are main), ``classification.gene_priority``
    (which gene's call wins) and ``parameters.gene_order`` (5' to 3'). R and
    Python read them differently when a name repeats or is not upper case, so
    those stop the run; a missing list only warns, because it falls back to
    ``main_probes``.
    """

    @staticmethod
    def _config(main: list[str], order: object = None, priority: object = None) -> dict:
        parameters: dict[str, object] = {"main_probes": main}
        classification: dict[str, object] = {}
        if order is not None:
            parameters["gene_order"] = order
        if priority is not None:
            classification["gene_priority"] = priority
        return {"parameters": parameters, "classification": classification}

    def test_a_clean_config_has_no_errors_and_no_notes(self) -> None:
        import validator as v

        config = self._config(["POL", "GAG"], ["GAG", "POL"], ["POL", "GAG"])
        assert v.gene_list_errors(config) == []
        assert v.gene_list_notes(config, {"POL", "GAG", "ENV"}) == []

    def test_a_repeated_name_is_an_error(self) -> None:
        import validator as v

        (error,) = v.gene_list_errors(self._config(["POL", "GAG", "POL"], [], []))
        assert "parameters.main_probes" in error
        assert "POL" in error

    def test_a_name_that_is_not_upper_case_is_an_error(self) -> None:
        import validator as v

        (error,) = v.gene_list_errors(self._config(["POL"], ["Gag", "POL"], ["POL"]))
        assert "parameters.gene_order" in error
        assert "Gag" in error

    def test_a_missing_list_is_noted_with_its_fallback(self) -> None:
        import validator as v

        notes = v.gene_list_notes(self._config(["POL", "GAG"]), {"POL", "GAG"})
        assert len(notes) == 2
        assert any("parameters.gene_order" in note for note in notes)
        assert any("classification.gene_priority" in note for note in notes)
        assert all("POL, GAG" in note for note in notes)

    def test_a_listed_name_that_is_no_probe_is_noted(self) -> None:
        import validator as v

        config = self._config(["POL"], ["GAG", "PRO", "POL"], ["POL"])
        (note,) = v.gene_list_notes(config, {"POL", "GAG"})
        assert "PRO" in note
        assert "parameters.gene_order" in note

    def test_no_probe_table_means_no_name_check(self) -> None:
        import validator as v

        config = self._config(["POL"], ["GAG", "POL"], ["POL"])
        assert v.gene_list_notes(config, None) == []

    def test_probe_names_are_read_upper_cased(self, tmp_path: Path) -> None:
        import validator as v

        csv_file = tmp_path / "probes.csv"
        csv_file.write_text(
            "Label,Name,Abbreviation,Probe,Accession\nA,V,Va,Pr160,P1\n"
        )
        assert v.probe_names(csv_file) == {"PR160"}
        assert v.probe_names(tmp_path / "absent.csv") is None

    def test_probe_names_are_not_trimmed(self, tmp_path: Path) -> None:
        # The pipeline upper-cases a probe name and keeps its spaces, so "POL "
        # is a different probe from "POL" and the check must say so.
        csv_file = tmp_path / "probes.csv"
        csv_file.write_text("Label,Name,Abbreviation,Probe,Accession\nA,V,Va,POL ,P1\n")
        assert v_probe_names(csv_file) == {"POL "}

    def test_a_damaged_probe_table_gives_no_names_instead_of_a_crash(
        self, tmp_path: Path
    ) -> None:
        ragged = tmp_path / "ragged.csv"
        ragged.write_text("Label,Name,Abbreviation,Probe,Accession\nA,V\n")
        latin1 = tmp_path / "latin1.csv"
        latin1.write_bytes("Label,Probe\nGonz\xe1lez,POL\n".encode("latin-1"))
        no_column = tmp_path / "no_column.csv"
        no_column.write_text("Label,Name\nA,V\n")
        assert v_probe_names(ragged) == set()
        assert v_probe_names(latin1) is None
        assert v_probe_names(no_column) is None

    def test_an_error_stops_the_preflight_and_a_note_does_not(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        import validator as v

        def run(config: dict) -> bool:
            with (
                patch.object(v, "yaml_validator", return_value=True),
                patch.object(v.defaults, "config", config),
                patch.object(v, "probe_names", return_value=None),
            ):
                return v.preflight([])

        with caplog.at_level("WARNING"):
            assert run(self._config(["POL", "GAG"])) is True
        assert "parameters.gene_order" in caplog.text
        assert run(self._config(["POL", "POL"], ["POL"], ["POL"])) is False
