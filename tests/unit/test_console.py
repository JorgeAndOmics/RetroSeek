"""Unit tests for workflow/scripts/console.py, the launcher's renderer (ADR-021).

The launcher reads one merged stream from Snakemake and every job. It recognises
three kinds of line and nothing else: our line contract, Snakemake's progress line
and Snakemake's error block. These tests pin that parsing on recorded text, the
verbosity filter, and the summary built from what streamed past.
"""

import io
import logging
import os
import signal
import sys
import threading

import pytest

import console

OURS = "14:02:40 WARN solo_finder Mus_musculus | 3 baits under 300 bp"


def test_parse_event_reads_our_line() -> None:
    event = console.parse_event(OURS)
    assert event == console.Event(
        "14:02:40", "WARN", "solo_finder", "Mus_musculus", "3 baits under 300 bp"
    )


@pytest.mark.parametrize(
    "line",
    [
        "Select jobs to execute...",
        "[Wed Sep 23 14:48:28 2026]",
        "rule ranges_analysis_setup:",
        "14:02:40 NOTE step g | not a level we use",
        "",
    ],
)
def test_parse_event_ignores_everything_else(line: str) -> None:
    assert console.parse_event(line) is None


def test_parse_progress() -> None:
    assert console.parse_progress("212 of 340 steps (62%) done") == (212, 340)
    assert console.parse_progress("Finished jobid: 3 (Rule: x)") is None


@pytest.mark.parametrize(
    ("verbosity", "level", "shown"),
    [
        ("quiet", "ERROR", True),
        ("quiet", "WARN", False),
        ("quiet", "OK", False),
        ("normal", "OK", True),
        ("normal", "WARN", True),
        ("normal", "INFO", False),
        ("verbose", "DEBUG", True),
    ],
)
def test_shown_follows_verbosity(verbosity: str, level: str, shown: bool) -> None:
    event = console.Event("t", level, "s", "g", "m")
    assert console.shown(event, verbosity) is shown


SNAKEMAKE_ERROR = """\
[Wed Sep 23 15:10:02 2026]
Error in rule solo_finder_setup:
    message: None
    jobid: 12
    input: /data/study/x.tsv
    output: /data/study/y.csv
    log: /data/study/logs/solo_finder/Mus_musculus.log (check log file(s) for error details)
    shell:
        python scripts/solo_ltr/solo_finder.py
        (command exited with non-zero exit status 1)
    wildcards: genome=Mus_musculus

"""


def _feed(tally: console.Tally, text: str) -> None:
    for line in text.splitlines():
        tally.feed(line)


def test_tally_reads_a_snakemake_error_block() -> None:
    tally = console.Tally()
    _feed(tally, SNAKEMAKE_ERROR)
    assert tally.failures == [
        console.Failure(
            "solo_finder_setup",
            "Mus_musculus",
            "/data/study/logs/solo_finder/Mus_musculus.log",
        )
    ]


def test_tally_keeps_the_last_error_line_of_the_failed_job() -> None:
    tally = console.Tally()
    tally.feed("15:10:01 ERROR solo_finder_setup Mus_musculus | bait file is empty")
    tally.feed("15:10:01 HINT solo_finder_setup Mus_musculus | not a level")
    _feed(tally, SNAKEMAKE_ERROR)
    assert tally.last_error("solo_finder_setup", "Mus_musculus") == "bait file is empty"


def test_tally_groups_warnings_that_differ_only_in_numbers() -> None:
    tally = console.Tally()
    tally.feed("10:00:00 WARN solo_finder Mus_musculus | 3 baits under 300 bp")
    tally.feed("10:00:01 WARN solo_finder Homo_sapiens | 12 baits under 300 bp")
    tally.feed("10:00:02 WARN hotspot Homo_sapiens | no windows")
    assert tally.warning_groups() == [
        ("solo_finder", "N baits under N bp", 2),
        ("hotspot", "no windows", 1),
    ]


def test_tally_tracks_progress() -> None:
    tally = console.Tally()
    tally.feed("3 of 10 steps (30%) done")
    tally.feed("4 of 10 steps (40%) done")
    assert tally.progress == (4, 10)


def test_summary_text_names_status_warnings_and_failures() -> None:
    tally = console.Tally()
    tally.feed("15:10:01 WARN hotspot Homo_sapiens | no windows")
    tally.feed("15:10:01 ERROR solo_finder_setup Mus_musculus | bait file is empty")
    _feed(tally, SNAKEMAKE_ERROR)
    text = console.summary_lines(
        tally, status="failed", elapsed_s=3725, run_log="/r.log"
    )
    joined = "\n".join(text)
    assert "failed" in joined
    assert "1h02m" in joined
    assert "hotspot: no windows (1)" in joined
    assert "solo_finder_setup Mus_musculus: bait file is empty" in joined
    assert "/data/study/logs/solo_finder/Mus_musculus.log" in joined
    assert "rerun the same command" in joined


def test_render_writes_plain_text_without_a_terminal() -> None:
    """Piped or redirected output must carry no colour codes."""
    out = io.StringIO()
    screen = console.Screen(file=out, verbosity="normal")
    screen.show(console.parse_event(OURS))
    assert "\x1b[" not in out.getvalue()
    # On screen the genome reads as a name; files keep the stem.
    assert "Mus musculus" in out.getvalue()
    assert "3 baits under 300 bp" in out.getvalue()


def test_no_color_turns_colour_off(monkeypatch) -> None:
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setenv("FORCE_COLOR", "1")
    out = io.StringIO()
    console.Screen(file=out, verbosity="normal").show(console.parse_event(OURS))
    assert "\x1b[" not in out.getvalue()


def test_elapsed_formats() -> None:
    assert console.elapsed(59) == "59s"
    assert console.elapsed(125) == "2m05s"
    assert console.elapsed(3725) == "1h02m"


def test_stream_routes_lines_and_returns_the_exit_code(tmp_path) -> None:
    """A stand-in for Snakemake: our line, a progress line, a tool line, exit 2."""
    script = (
        "import sys\n"
        "print('10:00:00 OK ranges Homo_sapiens | 26,499 elements')\n"
        "print('1 of 2 steps (50%) done')\n"
        "print('mafft chatter', file=sys.stderr)\n"
        "sys.exit(2)\n"
    )
    out = io.StringIO()
    screen = console.Screen(file=out, verbosity="normal")
    tally = console.Tally()
    run_log = tmp_path / "runs" / "r.log"

    code = console.stream([sys.executable, "-c", script], screen, tally, run_log)

    assert code == 2
    assert tally.progress == (1, 2)
    shown = out.getvalue()
    assert "26,499 elements" in shown
    assert "mafft chatter" not in shown  # not ours: run log only at normal
    written = run_log.read_text()
    assert "mafft chatter" in written
    assert "26,499 elements" in written


def test_launcher_messages_reach_the_screen_and_the_run_log(tmp_path) -> None:
    out = io.StringIO()
    screen = console.Screen(file=out, verbosity="normal")
    run_log = tmp_path / "r.log"
    handler = console.ScreenHandler(screen, run_log)
    record = logging.LogRecord(
        "x", logging.ERROR, __file__, 1, "Preflight failed", None, None
    )
    handler.emit(record)
    info = logging.LogRecord(
        "x", logging.INFO, __file__, 1, "config is valid", None, None
    )
    handler.emit(info)

    assert "Preflight failed" in out.getvalue()
    assert "config is valid" not in out.getvalue()  # INFO: screen at verbose only
    written = run_log.read_text()
    assert " ERROR launcher all | Preflight failed" in written
    assert " INFO launcher all | config is valid" in written


SNAKEMAKE9_ERROR = """\
Error in rule make:
    message: None
    jobid: 2
    output: out/Bad_genome.txt
    log: logs/make/Bad_genome.log (check log file(s) for error details)
    shell:
        python -c "..."
        (command exited with non-zero exit code)
"""


def test_tally_reads_the_genome_from_the_log_path() -> None:
    """Snakemake 9 prints no wildcards line; the log layout names the genome."""
    tally = console.Tally()
    _feed(tally, SNAKEMAKE9_ERROR)
    tally.finish()
    assert tally.failures == [
        console.Failure("make", "Bad_genome", "logs/make/Bad_genome.log")
    ]


def test_tally_counts_a_failure_once_although_snakemake_repeats_it() -> None:
    """The block is printed when the job fails and again in the closing summary."""
    tally = console.Tally()
    _feed(tally, SNAKEMAKE9_ERROR + "Select jobs to execute...\n" + SNAKEMAKE9_ERROR)
    tally.finish()
    assert len(tally.failures) == 1


def test_an_interrupt_keeps_reading_until_snakemake_has_stopped(tmp_path) -> None:
    """Ctrl-C reaches Snakemake too; while it stops its jobs, its last lines must
    still be read (a full pipe would stall it) and the run reported as 130."""
    script = (
        "import sys, time\n"
        "print('started', flush=True)\n"
        "time.sleep(1.5)\n"
        "print('10:00:02 ERROR make g | stopped by the interrupt', flush=True)\n"
    )
    timer = threading.Timer(0.7, lambda: os.kill(os.getpid(), signal.SIGINT))
    timer.start()
    tally = console.Tally()
    run_log = tmp_path / "r.log"
    code = console.stream(
        [sys.executable, "-c", script],
        console.Screen(file=io.StringIO(), verbosity="normal"),
        tally,
        run_log,
    )
    timer.join()
    assert code == 130
    assert "stopped by the interrupt" in run_log.read_text()


def test_an_interrupted_run_is_told_to_continue_not_to_fix() -> None:
    lines = console.summary_lines(console.Tally(), "interrupted", 5.0, "r.log")
    assert lines[-1] == "Rerun the same command to continue: finished work is kept."
