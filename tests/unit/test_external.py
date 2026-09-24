"""Unit tests for workflow/scripts/external.py, the one way to run a tool.

Five scripts each had their own copy of "run a tool; on failure dump the last
3,000 characters of its stderr and exit". The copies wrote raw text into the
stream the launcher parses and exited with a bare message. The shared version
keeps the stderr in the job log and stops with a PipelineError that says so.
"""

import logging
import sys

import pytest

import external
from log import PipelineError


def test_success_returns_the_output() -> None:
    result = external.run_tool([sys.executable, "-c", "print('hello')"])
    assert result.stdout.strip() == "hello"


def test_failure_names_the_tool_and_its_exit_code(caplog) -> None:
    script = "import sys; sys.stderr.write('bad input on line 3'); sys.exit(4)"
    with caplog.at_level(logging.INFO), pytest.raises(PipelineError) as caught:
        external.run_tool([sys.executable, "-c", script])
    assert "exit code 4" in str(caught.value)
    assert "job log" in str(caught.value)
    assert "bad input on line 3" in caplog.text  # kept for the job log


def test_missing_tool_says_how_to_get_it() -> None:
    with pytest.raises(PipelineError) as caught:
        external.run_tool(["retroseek-no-such-tool", "--version"])
    assert "retroseek-no-such-tool" in str(caught.value)
    assert "conda" in str(caught.value)


def test_stdout_can_go_to_a_file(tmp_path) -> None:
    target = tmp_path / "out.txt"
    with target.open("w") as fh:
        external.run_tool([sys.executable, "-c", "print('into the file')"], stdout=fh)
    assert target.read_text().strip() == "into the file"
