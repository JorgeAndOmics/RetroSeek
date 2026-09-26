# =============================================================================
# log.py
# =============================================================================
# The console line contract for every Python script (ADR-021).
#
# Each message is one plain line on stderr, and in the job's log file:
#
#   14:02:40 WARN solo_finder Mus_musculus | 3 baits under 300 bp
#   time     level step       genome         message
#
# Scripts never colour anything. The launcher reads these lines from every job
# at once, draws them (colour, symbols, the progress bar) and writes the run log;
# see console.py. Run a script by hand and the same lines are plain and readable.
#
# Levels: DEBUG, INFO, OK (the one headline a job prints when it finishes), WARN
# and ERROR. Which of them reach stderr depends on RETROSEEK_VERBOSITY, set by
# the launcher from `display.verbosity`:
#
#   quiet / normal  OK, WARN, ERROR   (INFO and DEBUG stay in the job log)
#   verbose         everything
#   unset           INFO and up, for scripts run by hand
#
# Why lines are safe with many writers: Linux writes up to 4 KB to a pipe in one
# piece, and a handler writes each record with one call, so lines from parallel
# jobs interleave but never tear.
#
# A script's entry point is usually just:
#
#   job_logging(args.log, "subset_pfam")   # the job, from its Snakemake log path
#   run_main(main)                          # errors become one line + exit 1
#
# and a failure the script can explain is raised as
# `PipelineError("what went wrong", hint="what to do")`.
# =============================================================================

"""The console line contract and job logging for every Python script (ADR-021)."""

from __future__ import annotations

import logging
import os
import sys
import time
from collections.abc import Callable
from pathlib import Path

OK = 25
logging.addLevelName(OK, "OK")

VERBOSITIES = ("quiet", "normal", "verbose")

# Python's names for the two levels we spell differently.
_LEVEL_NAMES = {logging.WARNING: "WARN", logging.CRITICAL: "ERROR"}


class LineFormatter(logging.Formatter):
    """Formats a record as `HH:MM:SS LEVEL step genome | message`, per line.

    :param step: the rule or script name, one word.
    :param genome: the genome the job works on, or None for all genomes.
    """

    def __init__(self, step: str, genome: str | None) -> None:
        super().__init__()
        self.step = step
        self.genome = genome or "all"

    def format(self, record: logging.LogRecord) -> str:
        """Return the record as contract lines, one prefixed line per message line.

        A traceback, when the record carries one, is appended to the message.
        """
        level = _LEVEL_NAMES.get(record.levelno, record.levelname)
        stamp = time.strftime("%H:%M:%S", time.localtime(record.created))
        prefix = f"{stamp} {level} {self.step} {self.genome} | "
        text = record.getMessage()
        if record.exc_info:
            text = f"{text}\n{self.formatException(record.exc_info)}"
        # Every physical line gets the prefix, so the launcher can place each one.
        return "\n".join(prefix + line for line in text.splitlines() or [""])


def console_threshold(verbosity: str | None) -> int:
    """The lowest level written to stderr at `verbosity` (None: run by hand)."""
    if verbosity is None:
        return logging.INFO
    if verbosity not in VERBOSITIES:
        raise ValueError(
            f"verbosity {verbosity!r} is not one of {', '.join(VERBOSITIES)}"
        )
    return logging.DEBUG if verbosity == "verbose" else OK


def setup_logging(
    step: str, genome: str | None = None, log_file: Path | None = None
) -> None:
    """Send every logger in this process through the line contract.

    Configures the root logger, so library modules that only call
    `logging.getLogger(__name__)` follow the same format.

    :param step: the rule or script name shown on each line.
    :param genome: the genome this job works on, or None for all.
    :param log_file: the job log (Snakemake's `log:`), appended to; its folder is
        created. Without it, messages go to stderr only.
    """
    verbosity = os.environ.get("RETROSEEK_VERBOSITY")
    formatter = LineFormatter(step, genome)

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(console_threshold(verbosity))
    console.setFormatter(formatter)
    handlers: list[logging.Handler] = [console]

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        to_file = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        to_file.setLevel(logging.DEBUG if verbosity == "verbose" else logging.INFO)
        to_file.setFormatter(formatter)
        handlers.append(to_file)

    root = logging.getLogger()
    for old in root.handlers[:]:
        root.removeHandler(old)
        old.close()
    root.setLevel(logging.DEBUG)
    for handler in handlers:
        root.addHandler(handler)


def job_logging(log_file: Path | None, step: str) -> None:
    """Set up logging for a pipeline job from its Snakemake log path.

    Log paths are LOG_DIR/<step>/<genome>.log (`all.log` for jobs over every
    genome), so the path alone names the job: a script shared by two rules
    reports the right one. Run by hand without a path, lines go to stderr under
    `step`. Writes a "started" line with the command into the job log.
    """
    if log_file is None:
        setup_logging(step)
    else:
        genome = None if log_file.stem == "all" else log_file.stem
        setup_logging(log_file.parent.name, genome, log_file)
    logging.getLogger(__name__).info("started: %s", " ".join(sys.argv))


class PipelineError(Exception):
    """A failure the script understands: what went wrong, and what fixes it.

    :param message: what went wrong, in the user's terms.
    :param hint: the action that usually fixes it, or "" if there is none.
    """

    def __init__(self, message: str, hint: str = "") -> None:
        super().__init__(message)
        self.hint = hint

    def __str__(self) -> str:
        message = super().__str__().rstrip(".")
        return f"{message}. Fix: {self.hint}" if self.hint else f"{message}."


def run_main(main: Callable[[], None]) -> None:
    """Run a script's `main` so every ending is recorded the same way.

    A PipelineError becomes one ERROR line (message and fix) and exit 1. Any
    other exception becomes one ERROR line naming it; its full traceback goes to
    the job log (and to the screen only at verbose). The job log ends with
    "finished: done" or "finished: failed" and the elapsed time.
    """
    logger = logging.getLogger(__name__)
    started = time.monotonic()

    def footer(status: str) -> None:
        logger.info("finished: %s in %.1f s", status, time.monotonic() - started)

    try:
        main()
    except PipelineError as error:
        logger.error("%s", error)
        footer("failed")
        sys.exit(1)
    except (SystemExit, KeyboardInterrupt):
        footer("stopped")
        raise
    except Exception as error:
        # INFO so the traceback lands in the job log without flooding the screen.
        logger.info("traceback of the error below:", exc_info=True)
        logger.error(
            "unexpected %s: %s (traceback in the job log)", type(error).__name__, error
        )
        footer("failed")
        sys.exit(1)
    footer("done")
