# =============================================================================
# console.py
# =============================================================================
# What the launcher shows while a run goes (ADR-021).
#
# Snakemake and every job write into one pipe; the launcher reads it line by line
# (many writers, one reader, no threads) and does three things with each line:
#
#   1. writes it, uncoloured, to the run log (LOG_DIR/runs/<time>.log);
#   2. counts it: warnings, errors, failed jobs, progress (the Tally);
#   3. shows it on screen if the verbosity asks for it (the Screen).
#
# It recognises three kinds of line and nothing else: our own line contract
# (log.py), Snakemake's "N of M steps (P%) done", and Snakemake's "Error in rule"
# block. Anything else is Snakemake's or a tool's own output: always in the run
# log, on screen only at `verbose`.
#
# Colour and symbols come from rich, used here and nowhere else. rich turns them
# off by itself when the output is not a terminal; NO_COLOR turns them off too.
# =============================================================================

"""What the launcher shows and records while a run goes (ADR-021)."""

from __future__ import annotations

import logging
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from rich.console import Console
from rich.progress import BarColumn, Progress, TaskID, TextColumn, TimeElapsedColumn

from log import LineFormatter

_EVENT = re.compile(
    r"^(\d\d:\d\d:\d\d) (DEBUG|INFO|OK|WARN|ERROR) (\S+) (\S+) \| (.*)$"
)
_PROGRESS = re.compile(r"^(\d+) of (\d+) steps \(\d+%\) done")
_ERROR_IN_RULE = re.compile(r"^Error in rule (\w+):")
_LOG_LINE = re.compile(r"^\s+log: (\S+)")
_GENOME = re.compile(r"genome=([^,\s]+)")
_NUMBER = re.compile(r"\d[\d,.]*")

# Which levels reach the screen at each verbosity. Files always get everything.
_SHOWN = {
    "quiet": {"ERROR"},
    "normal": {"OK", "WARN", "ERROR"},
    "verbose": {"DEBUG", "INFO", "OK", "WARN", "ERROR"},
}

# Symbol and word per level, and the colour that goes with them. Written as
# escapes: tracked files stay ASCII (tests/unit/test_ascii_only.py).
_SYMBOL = {
    "OK": "\u2714",
    "WARN": "\u26a0",
    "ERROR": "\u2716",
    "INFO": "\u00b7",
    "DEBUG": "\u00b7",
}
_WORD = {"OK": "ok", "WARN": "WARN", "ERROR": "ERROR", "INFO": "info", "DEBUG": "debug"}
_STYLE = {
    "OK": "green",
    "WARN": "yellow",
    "ERROR": "bold red",
    "INFO": "",
    "DEBUG": "dim",
}


@dataclass(frozen=True)
class Event:
    """One line of the contract: `time level step genome | message`."""

    time: str
    level: str
    step: str
    genome: str
    message: str


@dataclass(frozen=True)
class Failure:
    """A failed job as Snakemake reports it."""

    rule: str
    genome: str
    log: str


def parse_event(line: str) -> Event | None:
    """The Event on `line`, or None if it is not our line contract."""
    match = _EVENT.match(line)
    return Event(*match.groups()) if match else None


def parse_progress(line: str) -> tuple[int, int] | None:
    """(done, total) from Snakemake's "N of M steps (P%) done", else None."""
    match = _PROGRESS.match(line)
    return (int(match.group(1)), int(match.group(2))) if match else None


def shown(event: Event, verbosity: str) -> bool:
    """Whether `event` goes on screen at `verbosity`."""
    return event.level in _SHOWN[verbosity]


def elapsed(seconds: float) -> str:
    """A short duration: 59s, 2m05s, 1h02m."""
    s = int(seconds)
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m{s % 60:02d}s"
    return f"{s // 3600}h{s % 3600 // 60:02d}m"


def _step_of(rule: str) -> str:
    """Scripts name their step after the rule without `_setup`."""
    return rule.removesuffix("_setup")


@dataclass
class Tally:
    """What streamed past: progress, warnings, errors and failed jobs."""

    progress: tuple[int, int] | None = None
    warnings: list[Event] = field(default_factory=list)
    errors: dict[tuple[str, str], str] = field(default_factory=dict)
    failures: list[Failure] = field(default_factory=list)
    _open: dict[str, str] | None = None

    def feed(self, line: str) -> None:
        """Take one line of the merged stream into account."""
        if self._open is not None:
            if line.strip() and line.startswith((" ", "\t")):
                self._read_error_detail(line)
                return
            self._close_error_block()

        if event := parse_event(line):
            self._record(event)
        elif progress := parse_progress(line):
            self.progress = progress
        elif match := _ERROR_IN_RULE.match(line):
            self._open = {"rule": match.group(1), "genome": "all", "log": ""}

    def _record(self, event: Event) -> None:
        """Keep a warning for the summary, and the last error of each (step, genome)."""
        if event.level == "WARN":
            self.warnings.append(event)
        elif event.level == "ERROR":
            self.errors[(event.step, event.genome)] = event.message

    def finish(self) -> None:
        """Close an error block the stream ended inside."""
        if self._open is not None:
            self._close_error_block()

    def _read_error_detail(self, line: str) -> None:
        assert self._open is not None
        if match := _LOG_LINE.match(line):
            self._open["log"] = match.group(1)
            # Our logs are LOG_DIR/<step>/<genome>.log, so the path names the
            # genome; Snakemake 9 prints no wildcards line in this block.
            if self._open["genome"] == "all":
                self._open["genome"] = Path(match.group(1)).stem
        elif "wildcards:" in line and (match := _GENOME.search(line)):
            self._open["genome"] = match.group(1)

    def _close_error_block(self) -> None:
        assert self._open is not None
        failure = Failure(**self._open)
        # Snakemake prints each error block twice: when the job fails, and again
        # in its closing summary.
        if failure not in self.failures:
            self.failures.append(failure)
        self._open = None

    def last_error(self, rule: str, genome: str) -> str | None:
        """The last ERROR line the failed job printed, if any."""
        return self.errors.get((_step_of(rule), genome)) or self.errors.get(
            (rule, genome)
        )

    def warning_groups(self) -> list[tuple[str, str, int]]:
        """(step, message pattern, count), numbers folded to N, most frequent first."""
        counts: dict[tuple[str, str], int] = {}
        for event in self.warnings:
            key = (event.step, _NUMBER.sub("N", event.message))
            counts[key] = counts.get(key, 0) + 1
        ranked = sorted(counts.items(), key=lambda item: -item[1])
        return [(step, pattern, n) for (step, pattern), n in ranked]


def summary_lines(
    tally: Tally,
    status: str,
    elapsed_s: float,
    run_log: str,
    warnings_file: str | None = None,
) -> list[str]:
    """The closing summary, as plain lines (the Screen styles them)."""
    tally.finish()
    lines = [f"status    {status}", f"time      {elapsed(elapsed_s)}"]

    groups = tally.warning_groups()
    listed = f"   (all of them: {warnings_file})" if warnings_file else ""
    lines.append(f"warnings  {len(tally.warnings)}{listed}")
    lines.extend(f"  {step}: {pattern} ({n})" for step, pattern, n in groups[:10])
    if len(groups) > 10:
        lines.append(f"  ... {len(groups) - 10} more kinds in the run log")

    lines.append(f"failed    {len(tally.failures)}")
    for failure in tally.failures:
        why = tally.last_error(failure.rule, failure.genome) or "see its log"
        lines.append(f"  {failure.rule} {failure.genome}: {why}")
        if failure.log:
            lines.append(f"    log: {failure.log}")

    lines.append(f"run log   {run_log}")
    if status == "interrupted":
        lines.append("Rerun the same command to continue: finished work is kept.")
    elif tally.failures:
        lines.append(
            "Fix the cause, then rerun the same command: finished work is kept."
        )
    return lines


class Screen:
    """Draws events, a progress bar and the summary on the terminal.

    Args:
        verbosity: "quiet", "normal" or "verbose".
        file: Where to draw. Defaults to stdout.
    """

    def __init__(self, verbosity: str, file: IO[str] | None = None) -> None:
        self.verbosity = verbosity
        self.plain = bool(os.environ.get("NO_COLOR"))
        self.console = Console(
            file=file or sys.stdout,
            highlight=False,
            soft_wrap=True,
            no_color=self.plain,
            force_terminal=False if self.plain else None,
        )
        unicode_ok = (self.console.encoding or "").lower().startswith("utf")
        self.marks = _SYMBOL if unicode_ok else _WORD
        self._bar: Progress | None = None
        self._task: TaskID | None = None
        self._last_plain_percent = -10

    def _print(self, text: str, style: str = "") -> None:
        target = self._bar.console if self._bar else self.console
        target.print(text, style=None if self.plain else (style or None), markup=False)

    def show(self, event: Event | None) -> None:
        """Draw one event, if the verbosity shows it."""
        if event is None or not shown(event, self.verbosity):
            return
        genome = event.genome.replace("_", " ")
        mark = self.marks[event.level]
        self._print(
            f"{event.time} {mark} {event.step:<16} {genome:<24} {event.message}",
            _STYLE[event.level],
        )

    def raw(self, line: str) -> None:
        """Draw a line that is not ours (Snakemake, tools): only at verbose."""
        if self.verbosity == "verbose":
            self._print(line, "dim")

    def heading(self, title: str) -> None:
        """A section rule, e.g. "Run" or "Summary"."""
        if self.plain or not self.console.is_terminal:
            self._print(f"== {title} ==")
        else:
            (self._bar.console if self._bar else self.console).rule(
                f"[bold]{title}", align="left", style="cyan"
            )

    def lines(self, lines: list[str], style: str = "") -> None:
        """Plain lines, e.g. the banner or the summary."""
        for line in lines:
            self._print(line, style)

    def progress(self, done: int, total: int) -> None:
        """Advance the bar (on a terminal) or print every 10% (elsewhere)."""
        if self.console.is_terminal and not self.plain:
            if self._bar is None:
                self._bar = Progress(
                    TextColumn("  "),
                    BarColumn(bar_width=40),
                    TextColumn("{task.completed}/{task.total} steps"),
                    TextColumn("{task.percentage:>3.0f}%"),
                    TimeElapsedColumn(),
                    console=self.console,
                )
                self._bar.start()
                self._task = self._bar.add_task("run", total=total)
            assert self._task is not None
            self._bar.update(self._task, completed=done, total=total)
            return
        percent = 100 * done // max(total, 1)
        if percent >= self._last_plain_percent + 10 or done == total:
            self._last_plain_percent = percent
            self._print(f"  progress  {done} of {total} steps ({percent}%)")

    def stop(self) -> None:
        """Take the bar down before the summary."""
        if self._bar is not None:
            self._bar.stop()
            self._bar = None


class ScreenHandler(logging.Handler):
    """Route the launcher's own log messages (checks, guard) like a job's lines.

    Every line goes into the run log, and onto the screen as its verbosity allows.
    """

    def __init__(self, screen: Screen, run_log: Path) -> None:
        super().__init__(level=logging.DEBUG)
        self.screen = screen
        self.run_log = run_log
        self.setFormatter(LineFormatter("launcher", None))

    def emit(self, record: logging.LogRecord) -> None:
        """Append the formatted record to the run log and show each of its lines."""
        text = self.format(record)
        self.run_log.parent.mkdir(parents=True, exist_ok=True)
        with self.run_log.open("a", encoding="utf-8") as out:
            out.write(text + "\n")
        for line in text.splitlines():
            self.screen.show(parse_event(line))


def stream(cmd: list[str], screen: Screen, tally: Tally, run_log: Path) -> int:
    """Run `cmd`, route every output line, and return its exit code.

    Ctrl-C reaches Snakemake too (same process group): it stops its jobs and
    marks unfinished outputs as incomplete. Meanwhile we keep reading, because
    its last lines belong in the run log and a full pipe would stall it, and then
    report 130 (interrupted).
    """
    run_log.parent.mkdir(parents=True, exist_ok=True)
    interrupted = False
    with run_log.open("a", encoding="utf-8") as record:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            errors="replace",
        )
        assert proc.stdout is not None
        try:
            while True:
                try:
                    for raw in proc.stdout:
                        _route(raw, screen, tally, record)
                    break
                except KeyboardInterrupt:
                    interrupted = True
            code = proc.wait()
        finally:
            screen.stop()
    return 130 if interrupted else code


def _route(raw: str, screen: Screen, tally: Tally, record: IO[str]) -> None:
    """One line of the stream: into the run log, the tally and (maybe) the screen."""
    line = raw.rstrip("\n")
    record.write(raw)
    tally.feed(line)
    if event := parse_event(line):
        screen.show(event)
    elif progress := parse_progress(line):
        screen.progress(*progress)
    else:
        screen.raw(line)
