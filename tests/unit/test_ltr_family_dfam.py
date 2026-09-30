"""Tests for the optional Dfam labelling of LTR families (solo_ltr/ltr_family_dfam.py).

Each family's representative arm is searched with every Dfam curated model
(nhmmer, gathering thresholds). The label is evidence, never a filter. What can go
wrong quietly: the best hit taken by E-value order in the file rather than by
score, a minus-strand alignment giving a negative coverage, and a family with no
hit silently dropped from the table.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import ltr_family_dfam as dfam
import pytest

_HEADERS = """\
HMMER3/f [3.4 | Aug 2023]
NAME  MER4A
ACC   DF000000161.4
DESC  MER4A LTR of an ERV1 family
LENG  454
GA    25.00;
TC    27.00;
NC    20.00;
CT    Type; LTR;
CT    SubType; ERV1;
MS    TaxId:9606 TaxName:Homo sapiens
HMM          A        C        G        T
  COMPO   1.38629  1.38629  1.38629  1.38629
//
HMMER3/f [3.4 | Aug 2023]
NAME  IAPLTR1_Mm
ACC   DF000000470.4
LENG  337
GA    30.00;
CT    Type; LTR;
CT    SubType; ERVK;
HMM          A        C        G        T
//
"""

# nhmmer --tblout: the target is our representative, the query the Dfam model.
_TBLOUT = """\
# target name  accession  query name  accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
Mmus_F001  -  MER4A  DF000000161.4  1  300  11  310  1  320  400  +  1e-30  90.0  0.1  -
Mmus_F001  -  IAPLTR1_Mm  DF000000470.4  1  337  1  390  1  390  400  +  1e-80  250.0  0.2  -
Mmus_F002  -  MER4A  DF000000161.4  5  200  300  101  300  100  300  -  1e-12  40.0  0.0  -
"""


def test_headers_give_each_models_accession_and_class(tmp_path: Path) -> None:
    hmm = tmp_path / "dfam.hmm"
    hmm.write_text(_HEADERS)
    models = dfam.read_model_headers(hmm)
    assert models["MER4A"] == ("DF000000161.4", "LTR/ERV1")
    assert models["IAPLTR1_Mm"] == ("DF000000470.4", "LTR/ERVK")


def test_each_representatives_best_hit_is_the_highest_score(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    best = dfam.best_hits(tbl)
    assert best["Mmus_F001"].model == "IAPLTR1_Mm"
    assert best["Mmus_F001"].coverage == pytest.approx(390 / 400)


def test_a_minus_strand_hit_has_a_positive_coverage(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    assert dfam.best_hits(tbl)["Mmus_F002"].coverage == pytest.approx(200 / 300)


def test_every_family_gets_a_row_hit_or_not(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    hmm = tmp_path / "dfam.hmm"
    hmm.write_text(_HEADERS)
    rows = dfam.label_rows(
        ["Mmus_F001", "Mmus_F002", "Mmus_F003"],
        dfam.best_hits(tbl),
        dfam.read_model_headers(hmm),
    )
    assert [r["dfam_name"] for r in rows] == ["IAPLTR1_Mm", "MER4A", ""]
    assert rows[0]["dfam_class"] == "LTR/ERVK"
    assert rows[2]["dfam_evalue"] == ""


def test_representatives_are_written_under_their_family_names(tmp_path: Path) -> None:
    families = tmp_path / "Toyus_toyus.ltr_families.csv"
    families.write_text(
        "arm,seqname,start,end,element,ltr_family,representative\n"
        "c|e1|L,c,1,8,e1,Ttoy_F001,True\n"
        "c|e1|R,c,20,27,e1,Ttoy_F001,False\n"
    )
    bait = tmp_path / "Toyus_toyus.bait.fna"
    bait.write_text(">c|e1|L\nACGTACGT\n>c|e1|R\nTTTTTTTT\n")
    out = tmp_path / "reps.fna"
    names = dfam.write_representatives([(families, bait)], out)
    assert names == ["Ttoy_F001"]
    assert out.read_text() == ">Ttoy_F001\nACGTACGT\n"


@pytest.mark.integration
def test_nhmmer_finds_a_representative_with_a_toy_model(tmp_path: Path) -> None:
    seq = "ACGTTGCAGGCTTACCGATCGGATCCTAGGCTAACGTTAGCCTAGGATCCGATCGGTAAGCCTGCAACGT"
    stockholm = tmp_path / "toy.sto"
    stockholm.write_text(f"# STOCKHOLM 1.0\ntoy1 {seq}\ntoy2 {seq}\n//\n")
    hmm = tmp_path / "toy.hmm"
    subprocess.run(
        ["hmmbuild", "--dna", "-n", "TOY", str(hmm), str(stockholm)],
        check=True,
        capture_output=True,
    )
    # Dfam models carry gathering thresholds and class lines; add them.
    text = hmm.read_text().replace("\nNSEQ", "\nGA    10.00;\nCT    Type; LTR;\nNSEQ")
    hmm.write_text(text)
    reps = tmp_path / "reps.fna"
    reps.write_text(f">Ttoy_F001\n{seq}\n")
    best = dfam.best_hits(dfam.run_nhmmer(hmm, reps, tmp_path / "work", 1))
    assert best["Ttoy_F001"].model == "TOY"
