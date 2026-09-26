# =============================================================================
# external.py
# =============================================================================
# The one way scripts run an external tool (mafft, iqtree, gappa, blastx, ...).
#
# On failure the tool's own error output goes into the job log (as INFO, so it
# stays off the screen at normal verbosity) and the script stops with a
# PipelineError naming the tool and its exit code. That replaces five private
# copies that wrote the tail of the tool's stderr straight into the stream the
# launcher parses (ADR-021).
# =============================================================================

"""Run an external tool, stopping with a PipelineError when it fails."""

from __future__ import annotations

import logging
import subprocess
from typing import IO, Any

from log import PipelineError

logger = logging.getLogger(__name__)

# How much of a failing tool's stderr to keep; long enough for a stack trace.
_STDERR_TAIL = 3000


def run_tool(
    cmd: list[str], stdout: IO[Any] | None = None
) -> subprocess.CompletedProcess[str]:
    """Run `cmd` and return its result, or stop with a PipelineError.

    :param cmd: the command, tool first.
    :param stdout: a file to write the tool's output to; captured when None.
    """
    try:
        result = subprocess.run(
            cmd,
            stdout=stdout or subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except FileNotFoundError as error:
        raise PipelineError(
            f"{cmd[0]} is not installed",
            hint="activate the RetroSeek conda environment, or run make env-update",
        ) from error
    if result.returncode != 0:
        logger.info(
            "%s error output (last %d characters):\n%s",
            cmd[0],
            _STDERR_TAIL,
            (result.stderr or "")[-_STDERR_TAIL:],
        )
        raise PipelineError(
            f"{cmd[0]} failed with exit code {result.returncode}",
            hint="its error output is in the job log",
        )
    return result
