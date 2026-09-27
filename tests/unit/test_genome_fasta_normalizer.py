"""Unit tests for ``workflow/scripts/genome_fasta_normalizer.py``.

The normalizer canonicalises whatever FASTA extension upstream provided
(.fa / .fna / .fasta / .ffn) into ``{genome}.fa`` via symlink. Tests
cover the extension-preference policy, the ambiguity-refusal rule, and
the idempotency guarantees.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from genome_fasta_normalizer import (
    EXT_PREFERENCE,
    normalize,
    pick_canonical_source,
)


# ---------------------------------------------------------------------
# pick_canonical_source - extension preference
# ---------------------------------------------------------------------
def _write_fasta(path: Path) -> None:
    """Write a minimal 2-line FASTA file."""
    path.write_text(">chr1\nACGTACGT\n")


def test_pick_canonical_source_prefers_dot_fa(tmp_path: Path) -> None:
    """If both ``Toyus.fa`` and ``Toyus.fna`` exist, ``.fa`` wins."""
    _write_fasta(tmp_path / "Toyus.fa")
    _write_fasta(tmp_path / "Toyus.fna")
    assert pick_canonical_source(tmp_path, "Toyus") == tmp_path / "Toyus.fa"


def test_pick_canonical_source_falls_back_to_fna(tmp_path: Path) -> None:
    """``.fna`` only - NCBI Datasets common case."""
    _write_fasta(tmp_path / "Toyus.fna")
    assert pick_canonical_source(tmp_path, "Toyus") == tmp_path / "Toyus.fna"


def test_pick_canonical_source_falls_back_to_fasta(tmp_path: Path) -> None:
    """``.fasta`` only - manual-download convention."""
    _write_fasta(tmp_path / "Toyus.fasta")
    assert pick_canonical_source(tmp_path, "Toyus") == tmp_path / "Toyus.fasta"


def test_pick_canonical_source_falls_back_to_ffn(tmp_path: Path) -> None:
    """``.ffn`` only - older GenBank exports."""
    _write_fasta(tmp_path / "Toyus.ffn")
    assert pick_canonical_source(tmp_path, "Toyus") == tmp_path / "Toyus.ffn"


def test_pick_canonical_source_raises_when_none_present(tmp_path: Path) -> None:
    """No matching FASTA -> clear error listing the four extensions tried."""
    with pytest.raises(FileNotFoundError) as excinfo:
        pick_canonical_source(tmp_path, "Toyus")
    msg = str(excinfo.value)
    for ext in EXT_PREFERENCE:
        assert ext in msg


def test_pick_canonical_source_refuses_on_ambiguous_pair(tmp_path: Path) -> None:
    """``.fna`` + ``.fasta`` (no ``.fa``) -> refuse, demand user disambiguation.

    A `SPECIES_DB` directory containing two non-`.fa` variants of the same
    genome could plausibly hold a working `.fna` (NCBI Datasets) plus an
    older `.fasta` from manual download. Silent preference would mask user
    error; explicit refusal forces an informed decision.
    """
    _write_fasta(tmp_path / "Toyus.fna")
    _write_fasta(tmp_path / "Toyus.fasta")
    with pytest.raises(RuntimeError, match=r"ambiguous|multiple"):
        pick_canonical_source(tmp_path, "Toyus")


# ---------------------------------------------------------------------
# normalize - symlink creation, idempotency, FASTA-shape validation
# ---------------------------------------------------------------------
def test_normalize_creates_canonical_fa_via_symlink(tmp_path: Path) -> None:
    """``.fna`` only -> symlink ``Toyus.fa -> Toyus.fna``."""
    src = tmp_path / "Toyus.fna"
    _write_fasta(src)
    output = tmp_path / "Toyus.fa"

    result = normalize(tmp_path, "Toyus", output)

    assert result == output
    assert output.is_symlink()
    assert output.resolve() == src.resolve()


def test_normalize_idempotent_when_canonical_correct(tmp_path: Path) -> None:
    """Re-running on a directory with a correct symlink is a no-op."""
    src = tmp_path / "Toyus.fna"
    _write_fasta(src)
    output = tmp_path / "Toyus.fa"
    normalize(tmp_path, "Toyus", output)

    # Re-run; must not raise, must still point at the same target.
    normalize(tmp_path, "Toyus", output)

    assert output.is_symlink()
    assert output.resolve() == src.resolve()


def test_normalize_replaces_stale_symlink(tmp_path: Path) -> None:
    """A pre-existing symlink to a deleted target gets replaced cleanly."""
    src = tmp_path / "Toyus.fna"
    _write_fasta(src)
    output = tmp_path / "Toyus.fa"
    output.symlink_to(tmp_path / "does_not_exist")

    normalize(tmp_path, "Toyus", output)

    assert output.resolve() == src.resolve()


def test_normalize_when_fa_is_already_a_real_file(tmp_path: Path) -> None:
    """A real ``Toyus.fa`` (not a symlink) is left untouched.

    Protects existing pipeline state where the genome arrived as a true
    `.fa` file (e.g., from `genome_downloader_setup`'s prior outputs).
    """
    real = tmp_path / "Toyus.fa"
    _write_fasta(real)
    original_bytes = real.read_bytes()

    normalize(tmp_path, "Toyus", real)

    assert real.is_file()
    assert not real.is_symlink()
    assert real.read_bytes() == original_bytes


def test_normalize_validates_first_byte_is_gt(tmp_path: Path) -> None:
    """A misnamed non-FASTA file (e.g. tarball) is rejected with a clear error."""
    bogus = tmp_path / "Toyus.fna"
    bogus.write_bytes(b"\x1f\x8b\x08\x00binary garbage")  # gzip magic
    output = tmp_path / "Toyus.fa"

    with pytest.raises(RuntimeError, match="FASTA"):
        normalize(tmp_path, "Toyus", output)


# ---------------------------------------------------------------------
# A working `.fa` link is only as good as what it points at
# ---------------------------------------------------------------------
def test_a_link_to_the_wrong_file_is_refused_not_kept(tmp_path: Path) -> None:
    """``Toyus.fa -> Toyus.fasta`` while ``Toyus.fna`` also exists: ambiguous.

    The link used to win the preference order and resolve to itself, so an old
    link survived next to the file that should back it. Which of the two is the
    genome is the user's call, so the normalizer refuses and changes nothing.
    """
    old = tmp_path / "Toyus.fasta"
    _write_fasta(old)
    new = tmp_path / "Toyus.fna"
    _write_fasta(new)
    output = tmp_path / "Toyus.fa"
    output.symlink_to(old.resolve())

    with pytest.raises(RuntimeError, match=r"Toyus\.fa links to .*Toyus\.fna"):
        normalize(tmp_path, "Toyus", output)
    assert output.resolve() == old.resolve()


def test_a_link_elsewhere_with_no_variant_beside_it_is_kept(tmp_path: Path) -> None:
    """A hand-made ``Toyus.fa`` link to a genome stored elsewhere is honoured."""
    elsewhere = tmp_path / "store" / "toyus_assembly.fa"
    elsewhere.parent.mkdir()
    _write_fasta(elsewhere)
    species = tmp_path / "species"
    species.mkdir()
    output = species / "Toyus.fa"
    output.symlink_to(elsewhere)

    assert pick_canonical_source(species, "Toyus") == output
    normalize(species, "Toyus", output)
    assert output.resolve() == elsewhere.resolve()


def test_a_link_to_its_only_variant_names_that_variant(tmp_path: Path) -> None:
    src = tmp_path / "Toyus.fna"
    _write_fasta(src)
    (tmp_path / "Toyus.fa").symlink_to(src.resolve())
    assert pick_canonical_source(tmp_path, "Toyus") == src
