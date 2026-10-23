"""The launcher's note that Dfam is off (validator.dfam_filter_notes, ADR-026).

With Dfam off, LTR elements whose arms are LINE, SINE or DNA stay in the LTR set and
their arms seed the solo-LTR search, which then counts every L1 or Alu copy as a solo
LTR. The run still works, so this is a warning, never a stop. What can go wrong
quietly: warning on runs that never reach ranges analysis, staying silent when the
switch is simply absent from the config, or repeating the note once per stage.
"""

from __future__ import annotations

import stages
import validator as v


def _stage(flag: str, *targets: str) -> stages.Stage:
    return stages.Stage(flag, "Analysis", tuple(targets), "x")


def _config(dfam: object) -> dict[str, object]:
    families = {} if dfam is None else {"dfam": dfam}
    return {"solo_ltr": {"families": families}}


def test_warns_when_ranges_analysis_runs_with_dfam_off():
    chosen = [_stage("--ranges-analysis", "ranges_analysis")]
    (note,) = v.dfam_filter_notes(chosen, _config(False))
    assert "solo_ltr.families.dfam" in note
    assert "LINE" in note


def test_silent_when_dfam_is_on():
    chosen = [_stage("--ranges-analysis", "ranges_analysis")]
    assert v.dfam_filter_notes(chosen, _config(True)) == []


def test_silent_when_ranges_analysis_is_not_part_of_the_run():
    chosen = [_stage("--classify", "taxonomy_classify")]
    assert v.dfam_filter_notes(chosen, _config(False)) == []


def test_an_unset_switch_counts_as_off():
    chosen = [_stage("--ranges-analysis", "ranges_analysis")]
    assert len(v.dfam_filter_notes(chosen, _config(None))) == 1
    assert len(v.dfam_filter_notes(chosen, {})) == 1


def test_one_note_however_many_stages_are_chosen():
    chosen = [
        _stage("--ranges-analysis", "ranges_analysis"),
        _stage("--downstream", "ranges_analysis", "taxonomy_classify"),
    ]
    assert len(v.dfam_filter_notes(chosen, _config(False))) == 1
