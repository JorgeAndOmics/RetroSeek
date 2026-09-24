"""Unit tests for workflow/scripts/log.py, the console line contract (ADR-021).

Every script in the pipeline writes the same plain line to stderr:

    HH:MM:SS LEVEL step genome | message

The launcher parses these lines while many jobs write at once, so the tests pin
the parts that matter to it: the exact shape, one prefix per physical line, the
level names, and which levels reach the terminal at each verbosity.
"""

import logging
import re
from pathlib import Path

import pytest

import log

LINE = re.compile(r"^\d\d:\d\d:\d\d (DEBUG|INFO|OK|WARN|ERROR) \S+ \S+ \| .*$")


@pytest.fixture(autouse=True)
def _clean_root():
    """setup_logging configures the root logger; leave it as we found it."""
    root = logging.getLogger()
    saved = root.handlers[:], root.level
    yield
    for handler in root.handlers:
        handler.close()
    root.handlers, root.level = saved[0], saved[1]


def _record(level: int, message: str) -> logging.LogRecord:
    return logging.LogRecord("x", level, __file__, 1, message, None, None)


def test_line_shape() -> None:
    text = log.LineFormatter("solo_finder", "Mus_musculus").format(
        _record(logging.WARNING, "3 baits under 300 bp")
    )
    assert LINE.match(text)
    assert text.endswith(" WARN solo_finder Mus_musculus | 3 baits under 300 bp")


def test_no_genome_reads_all() -> None:
    text = log.LineFormatter("hotspot", None).format(_record(logging.INFO, "done"))
    assert " INFO hotspot all | done" in text


@pytest.mark.parametrize(
    ("level", "name"),
    [
        (logging.DEBUG, "DEBUG"),
        (logging.INFO, "INFO"),
        (log.OK, "OK"),
        (logging.WARNING, "WARN"),
        (logging.ERROR, "ERROR"),
        (logging.CRITICAL, "ERROR"),
    ],
)
def test_level_names(level: int, name: str) -> None:
    text = log.LineFormatter("s", "g").format(_record(level, "m"))
    assert f" {name} s g | m" in text


def test_every_physical_line_carries_the_prefix() -> None:
    """A multi-line message must not leave bare lines the launcher cannot place."""
    text = log.LineFormatter("s", "g").format(_record(logging.INFO, "one\ntwo\nthree"))
    lines = text.split("\n")
    assert len(lines) == 3
    assert all(LINE.match(line) for line in lines)


@pytest.mark.parametrize(
    ("verbosity", "expected"),
    [
        ("quiet", log.OK),
        ("normal", log.OK),
        ("verbose", logging.DEBUG),
        (None, logging.INFO),
    ],
)
def test_console_threshold_follows_verbosity(verbosity, expected: int) -> None:
    assert log.console_threshold(verbosity) == expected


def test_unknown_verbosity_is_an_error() -> None:
    with pytest.raises(ValueError, match="loud"):
        log.console_threshold("loud")


def test_setup_writes_info_to_the_file_and_only_ok_up_to_stderr(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    job_log = tmp_path / "logs" / "rule" / "g.log"
    log.setup_logging("rule", "g", job_log)
    logger = logging.getLogger("t")
    logger.info("detail")
    logger.log(log.OK, "headline")
    logger.warning("careful")

    err = capsys.readouterr().err
    assert "detail" not in err
    assert "OK rule g | headline" in err
    assert "WARN rule g | careful" in err

    written = job_log.read_text()
    assert "INFO rule g | detail" in written
    assert "OK rule g | headline" in written


def test_setup_appends_to_an_existing_job_log(tmp_path: Path, monkeypatch) -> None:
    """A rule's shell may write tool output to the log before the script starts."""
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    job_log = tmp_path / "g.log"
    job_log.write_text("mafft said hello\n")
    log.setup_logging("rule", "g", job_log)
    logging.getLogger("t").info("after")
    assert job_log.read_text().startswith("mafft said hello\n")


def test_job_logging_reads_the_job_from_its_log_path(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    """LOG_DIR/<step>/<genome>.log names the job: scripts shared by two rules
    (taxonomy_classify_loci.py, solo_plots.R) still report the right step."""
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    log.job_logging(tmp_path / "taxonomy_orphans" / "Mus_musculus.log", "fallback")
    logging.getLogger("t").warning("w")
    assert "WARN taxonomy_orphans Mus_musculus | w" in capsys.readouterr().err


def test_job_logging_all_means_every_genome(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    log.job_logging(tmp_path / "hotspot_detector" / "all.log", "fallback")
    logging.getLogger("t").warning("w")
    assert "WARN hotspot_detector all | w" in capsys.readouterr().err


def test_job_logging_without_a_path_uses_the_fallback_step(capsys, monkeypatch) -> None:
    monkeypatch.delenv("RETROSEEK_VERBOSITY", raising=False)
    log.job_logging(None, "subset_pfam")
    logging.getLogger("t").info("by hand")
    assert "INFO subset_pfam all | by hand" in capsys.readouterr().err


def test_run_main_known_error_is_one_line_with_its_fix(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    job_log = tmp_path / "solo_finder" / "g.log"
    log.job_logging(job_log, "x")

    def main() -> None:
        raise log.PipelineError(
            "bait file is empty", hint="run --solo-ltr-detector again"
        )

    with pytest.raises(SystemExit) as caught:
        log.run_main(main)
    assert caught.value.code == 1
    err = [line for line in capsys.readouterr().err.splitlines() if " ERROR " in line]
    assert err == [err[0]]
    assert err[0].endswith("| bait file is empty. Fix: run --solo-ltr-detector again")
    assert "finished: failed" in job_log.read_text()


def test_run_main_unexpected_error_keeps_the_traceback_in_the_job_log(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    job_log = tmp_path / "solo_finder" / "g.log"
    log.job_logging(job_log, "x")

    def main() -> None:
        {}["missing"]

    with pytest.raises(SystemExit) as caught:
        log.run_main(main)
    assert caught.value.code == 1
    err = capsys.readouterr().err
    assert "Traceback" not in err  # the screen gets one line
    assert "KeyError" in err
    written = job_log.read_text()
    assert "Traceback" in written
    assert "finished: failed" in written


def test_run_main_success_writes_the_footer(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("RETROSEEK_VERBOSITY", "normal")
    job_log = tmp_path / "step" / "g.log"
    log.job_logging(job_log, "x")
    log.run_main(lambda: None)
    written = job_log.read_text()
    assert "started:" in written
    assert "finished: done" in written
