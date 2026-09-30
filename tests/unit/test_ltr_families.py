"""Tests for the LTR family builder (solo_ltr/ltr_families.py, ADR-023).

A family is a group of bait arms at least ``identity`` identical over the whole
shorter arm, found by cd-hit-est. What can go wrong quietly: identifiers that
change between runs of the same input, a representative read from the wrong
line of the cluster file, an element counted twice because it has two arms, and
a genus majority that depends on dictionary order.
"""

from __future__ import annotations

from pathlib import Path

import ltr_families as lf
import pytest

from log import PipelineError

# cd-hit-est's cluster file, as it writes it: the representative carries "*".
_CLSTR = """\
>Cluster 0
0\t999nt, >chr2|LTR_retrotransposon5|L... *
1\t998nt, >chr2|LTR_retrotransposon5|R... at +/99.80%
2\t549nt, >chr1|LTR_retrotransposon1|L... at -/96.54%
>Cluster 1
0\t600nt, >chr1|LTR_retrotransposon2|L... *
>Cluster 2
0\t610nt, >chr1|LTR_retrotransposon3|R... *
1\t500nt, >chr1|LTR_retrotransposon2|R... at +/88.10%
"""

# The bait BED (0-based starts), one line per arm.
_BED = """\
chr1\t99\t649\tchr1|LTR_retrotransposon1|L\t.\t+
chr1\t999\t1600\tchr1|LTR_retrotransposon2|L\t.\t+
chr1\t4999\t5500\tchr1|LTR_retrotransposon2|R\t.\t+
chr1\t8999\t9610\tchr1|LTR_retrotransposon3|R\t.\t+
chr2\t199\t1198\tchr2|LTR_retrotransposon5|L\t.\t+
chr2\t8999\t9997\tchr2|LTR_retrotransposon5|R\t.\t+
"""


@pytest.fixture
def clstr(tmp_path: Path) -> Path:
    path = tmp_path / "bait.clstr"
    path.write_text(_CLSTR)
    return path


@pytest.fixture
def arms(tmp_path: Path) -> dict[str, lf.Arm]:
    path = tmp_path / "bait.bed"
    path.write_text(_BED)
    return lf.read_bait_bed(path)


# ---- names ----


@pytest.mark.parametrize(
    ("genome", "code"),
    [("Mus_musculus", "Mmus"), ("Molossus_molossus", "Mmol"), ("Homo_sapiens", "Hsap")],
)
def test_the_species_code_is_one_genus_letter_and_three_species_letters(
    genome: str, code: str
) -> None:
    assert lf.species_code(genome) == code


def test_a_genome_name_without_a_species_part_is_an_error() -> None:
    with pytest.raises(PipelineError, match="Genus_species"):
        lf.species_code("Mus")


def test_the_element_is_read_from_the_right_of_the_arm_name() -> None:
    # Pooled arms carry a genome prefix; the element is still the second field
    # from the right.
    pooled = lf.Arm("Mus_musculus|chr1|LTR_retrotransposon7|R", "chr1", 1, 2)
    assert pooled.element == "LTR_retrotransposon7"


def test_the_bait_bed_is_read_one_based(arms: dict[str, lf.Arm]) -> None:
    arm = arms["chr1|LTR_retrotransposon1|L"]
    assert (arm.seqname, arm.start, arm.end, arm.element) == (
        "chr1",
        100,
        649,
        "LTR_retrotransposon1",
    )


# ---- cd-hit-est ----


@pytest.mark.parametrize(
    ("identity", "word"), [(0.80, 5), (0.85, 6), (0.88, 7), (0.90, 8), (0.95, 10)]
)
def test_the_word_size_follows_the_identity(identity: float, word: int) -> None:
    assert lf.word_size(identity) == word


def test_an_identity_below_what_cd_hit_accepts_is_an_error() -> None:
    with pytest.raises(PipelineError, match=r"0\.8"):
        lf.word_size(0.75)


def test_the_command_counts_identity_over_the_whole_shorter_arm(
    tmp_path: Path,
) -> None:
    cmd = lf.cdhit_command(tmp_path / "a.fna", tmp_path / "out", 0.80, 4)
    flags = dict(zip(cmd[1::2], cmd[2::2], strict=False))
    assert cmd[0] == "cd-hit-est"
    assert flags["-c"] == "0.8"
    assert flags["-n"] == "5"
    assert flags["-G"] == "1"  # global identity: the whole shorter arm
    assert flags["-r"] == "1"  # both strands: bait arms are not oriented
    assert flags["-d"] == "0"  # keep full names
    assert flags["-T"] == "4"


def test_the_cluster_file_gives_representatives_and_members(clstr: Path) -> None:
    clusters = lf.parse_clstr(clstr)
    assert [rep for rep, _ in clusters] == [
        "chr2|LTR_retrotransposon5|L",
        "chr1|LTR_retrotransposon2|L",
        "chr1|LTR_retrotransposon3|R",
    ]
    assert len(clusters[0][1]) == 3


# ---- families ----


def test_families_are_named_by_size_largest_first(
    clstr: Path, arms: dict[str, lf.Arm]
) -> None:
    families = lf.name_families(lf.parse_clstr(clstr), arms, "Toyu")
    assert [(f.name, f.representative, len(f.members)) for f in families] == [
        ("Toyu_F001", "chr2|LTR_retrotransposon5|L", 3),
        ("Toyu_F002", "chr1|LTR_retrotransposon3|R", 2),
        ("Toyu_F003", "chr1|LTR_retrotransposon2|L", 1),
    ]


def test_families_of_one_size_are_ordered_by_position(
    tmp_path: Path, arms: dict[str, lf.Arm]
) -> None:
    clstr = tmp_path / "tie.clstr"
    clstr.write_text(
        ">Cluster 0\n0\t600nt, >chr1|LTR_retrotransposon2|L... *\n"
        ">Cluster 1\n0\t550nt, >chr1|LTR_retrotransposon1|L... *\n"
    )
    families = lf.name_families(lf.parse_clstr(clstr), arms, "Toyu")
    assert [f.representative for f in families] == [
        "chr1|LTR_retrotransposon1|L",  # starts at 100, before 1000
        "chr1|LTR_retrotransposon2|L",
    ]


def test_a_summary_counts_elements_once_and_flags_split_pairs(
    clstr: Path, arms: dict[str, lf.Arm]
) -> None:
    families = lf.name_families(lf.parse_clstr(clstr), arms, "Toyu")
    genus = {
        ("chr2", "LTR_retrotransposon5"): "Betaretrovirus",
        ("chr1", "LTR_retrotransposon1"): "Gammaretrovirus",
        ("chr1", "LTR_retrotransposon2"): "Gammaretrovirus",
        ("chr1", "LTR_retrotransposon3"): "Gammaretrovirus",
    }
    similarity = {
        ("chr2", "LTR_retrotransposon5"): 99.0,
        ("chr1", "LTR_retrotransposon1"): 91.0,
    }
    rows = {
        r["ltr_family"]: r for r in lf.summary_rows(families, arms, genus, similarity)
    }
    first = rows["Toyu_F001"]
    assert (first["n_arms"], first["n_elements"]) == (3, 2)
    # Two elements, one of each genus: the tie goes to the name first in byte order.
    assert first["majority_genus"] == "Betaretrovirus"
    assert first["genus_purity"] == pytest.approx(0.5)
    assert first["median_arm_similarity"] == pytest.approx(95.0)
    assert first["split_elements"] == 0
    # Element 2 has one arm in F002 and the other in F003.
    assert rows["Toyu_F002"]["split_elements"] == 1
    assert rows["Toyu_F003"]["split_elements"] == 1
    assert rows["Toyu_F003"]["median_arm_similarity"] == ""


def test_the_genus_table_counts_elements_per_family(
    clstr: Path, arms: dict[str, lf.Arm]
) -> None:
    families = lf.name_families(lf.parse_clstr(clstr), arms, "Toyu")
    genus = {("chr2", "LTR_retrotransposon5"): "Betaretrovirus"}
    rows = lf.genus_rows(families, arms, genus)
    assert {
        "ltr_family": "Toyu_F001",
        "genus": "Betaretrovirus",
        "n_elements": 1,
    } in rows
    assert {"ltr_family": "Toyu_F001", "genus": "", "n_elements": 1} in rows


def test_an_arm_missing_from_the_bait_bed_is_an_error(tmp_path: Path) -> None:
    clstr = tmp_path / "x.clstr"
    clstr.write_text(">Cluster 0\n0\t500nt, >chrX|LTR_retrotransposon9|L... *\n")
    with pytest.raises(PipelineError, match="chrX"):
        lf.name_families(lf.parse_clstr(clstr), {}, "Toyu")


@pytest.mark.integration
def test_a_real_run_puts_two_near_identical_arms_in_one_family(tmp_path: Path) -> None:
    seq = "ACGTTGCAGGCTTACCGATCGGATCCTAGGCTAACGTTAGCCTAGGATCCGATCGGTAAGCCTGCAACGT" * 5
    fna = tmp_path / "bait.fna"
    fna.write_text(f">c|e1|L\n{seq}\n>c|e1|R\n{seq[:-1]}A\n")
    clusters = lf.run_cdhit(fna, tmp_path / "work", 0.80, 1)
    assert len(clusters) == 1
