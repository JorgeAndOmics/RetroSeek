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

from log import PipelineError

# Headers as Dfam 4.0 writes them: the RepeatMasker class is in the CC lines, and
# CT holds the long classification path. Not every model has a GA line.
_HEADERS = """\
HMMER3/f [3.4 | Aug 2023]
NAME  MER4A
ACC   DF000000161.4
LENG  454
TH    TaxId:9606; TaxName:Homo_sapiens; GA:15.00; TC:27.00; NC:14.00; fdr:0.002;
CT    Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1;
CC         Type: LTR
CC         SubType: ERV1
HMM          A        C        G        T
  COMPO   1.38629  1.38629  1.38629  1.38629
//
HMMER3/f [3.4 | Aug 2023]
NAME  IAPLTR1_Mm
ACC   DF000004162.1
GA    30.94;
CT    Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;
CC         Type: LTR
CC         SubType: ERVK
HMM          A        C        G        T
//
HMMER3/f [3.4 | Aug 2023]
NAME  L1_Mus1
ACC   DF000000001.1
CC         Type: LINE
HMM          A        C        G        T
//
"""

# nhmmer --tblout: the target is our representative, the query the Dfam model.
_TBLOUT = """\
# target name  accession  query name  accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
Mus_musculus|Mmus_F001  -  MER4A  DF000000161.4  1  300  11  310  1  320  400  +  1e-30  90.0  0.1  -
Mus_musculus|Mmus_F001  -  IAPLTR1_Mm  DF000004162.1  1  337  1  390  1  390  400  +  1e-80  250.0  0.2  -
Mus_musculus|Mmus_F002  -  MER4A  DF000000161.4  5  200  300  101  300  100  300  -  1e-12  40.0  0.0  -
"""


def test_headers_give_each_models_accession_and_class(tmp_path: Path) -> None:
    hmm = tmp_path / "dfam.hmm"
    hmm.write_text(_HEADERS)
    models = dfam.read_model_headers(hmm, {"MER4A", "IAPLTR1_Mm", "L1_Mus1"})
    assert models["MER4A"] == ("DF000000161.4", "LTR/ERV1")
    assert models["IAPLTR1_Mm"] == ("DF000004162.1", "LTR/ERVK")
    assert models["L1_Mus1"] == ("DF000000001.1", "LINE")


def test_only_the_wanted_models_are_read(tmp_path: Path) -> None:
    hmm = tmp_path / "dfam.hmm"
    hmm.write_text(_HEADERS)
    assert set(dfam.read_model_headers(hmm, {"IAPLTR1_Mm"})) == {"IAPLTR1_Mm"}


def test_each_representatives_best_hit_is_the_highest_score(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    best = dfam.best_hits(tbl)
    assert best["Mus_musculus|Mmus_F001"].model == "IAPLTR1_Mm"
    assert best["Mus_musculus|Mmus_F001"].coverage == pytest.approx(390 / 400)


def test_a_minus_strand_hit_has_a_positive_coverage(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    assert dfam.best_hits(tbl)["Mus_musculus|Mmus_F002"].coverage == pytest.approx(
        200 / 300
    )


def test_every_family_gets_a_row_hit_or_not(tmp_path: Path) -> None:
    tbl = tmp_path / "hits.tbl"
    tbl.write_text(_TBLOUT)
    hmm = tmp_path / "dfam.hmm"
    hmm.write_text(_HEADERS)
    rows = dfam.label_rows(
        ["Mus_musculus|Mmus_F001", "Mus_musculus|Mmus_F002", "Mus_musculus|Mmus_F003"],
        dfam.best_hits(tbl),
        dfam.read_model_headers(hmm, {"MER4A", "IAPLTR1_Mm"}),
    )
    assert [r["dfam_name"] for r in rows] == ["IAPLTR1_Mm", "MER4A", ""]
    assert (rows[0]["genome"], rows[0]["ltr_family"]) == ("Mus_musculus", "Mmus_F001")
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
    names = dfam.write_representatives([("Toyus_toyus", families, bait)], out)
    # Named by genome too: two genomes can share a family code.
    assert names == ["Toyus_toyus|Ttoy_F001"]
    assert out.read_text() == ">Toyus_toyus|Ttoy_F001\nACGTACGT\n"


def test_a_representative_missing_from_the_bait_stops_the_job(tmp_path: Path) -> None:
    families = tmp_path / "families.csv"
    families.write_text(
        "arm,seqname,start,end,element,ltr_family,representative\n"
        "c|e9|L,c,1,8,e9,Ttoy_F001,True\n"
    )
    bait = tmp_path / "bait.fna"
    bait.write_text(">c|e1|L\nACGTACGT\n")
    with pytest.raises(PipelineError, match="Ttoy_F001"):
        dfam.write_representatives(
            [("Toyus_toyus", families, bait)], tmp_path / "r.fna"
        )


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
    reps = tmp_path / "reps.fna"
    reps.write_text(f">Toyus_toyus|Ttoy_F001\n{seq}\n")
    best = dfam.best_hits(dfam.run_nhmmer(hmm, reps, tmp_path / "work", 1))
    assert best["Toyus_toyus|Ttoy_F001"].model == "TOY"
