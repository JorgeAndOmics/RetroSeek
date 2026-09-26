"""Tests for taxonomy_build_tree: the reference-alignment quality line."""

from __future__ import annotations

from pathlib import Path

import taxonomy_build_tree as tbt


class TestAlignmentQuality:
    """mean_pident samples the first 40 informative pairs in alignment order."""

    @staticmethod
    def _write(path: Path, seqs: list[str]) -> Path:
        path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
        return path

    def test_value_and_flag(self, tmp_path: Path) -> None:
        afa = self._write(tmp_path / "a.afa", ["AAAA", "AAAT", "--AA"])
        # pairs: (0,1) 3/4, (0,2) 2/2, (1,2) 1/2 -> mean 75%
        assert tbt.alignment_quality(afa) == (
            "n=3 cols=4 gap=17% mean_pident=75.0% -> OK"
        )

    def test_only_the_first_40_informative_pairs_count(self, tmp_path: Path) -> None:
        # 10 sequences give 45 pairs; the 40 sampled are the first 40 in (i, j)
        # order. With the odd one last, 6 of them include it ((0,9) to (5,9)):
        # 34/40 = 85%. With it first, 9 do ((0,1) to (0,9)): 31/40 = 77.5%.
        same, odd = "AAAAAAAAAA", "TTTTTTTTTT"
        last = self._write(tmp_path / "last.afa", [same] * 9 + [odd])
        first = self._write(tmp_path / "first.afa", [odd] + [same] * 9)
        assert "mean_pident=85.0%" in tbt.alignment_quality(last)
        assert "mean_pident=77.5%" in tbt.alignment_quality(first)
