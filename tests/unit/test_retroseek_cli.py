# ruff: noqa: PLC0415
# Per-test imports of the module under test are deliberate (conftest
# defaults stub must take effect first).

"""Unit tests for the ``RetroSeek.py`` CLI wrapper around Snakemake.

Focus on ``run_snakemake_rule`` - a thin wrapper that shells out to the
``snakemake`` binary. Regression guards for three prior bugs:

1. ``subprocess.run`` was called without ``capture_output=True`` but the
   error branch logged ``result.stderr`` - which is ``None``. Failures
   therefore looked like empty log lines.
2. An ``except subprocess.CalledProcessError`` branch existed without
   ``check=True``, meaning it could never fire. Dead code that masked the
   real error handling.
3. On success the function returned ``None`` implicitly - fine, but the
   tests pin that contract so a future refactor doesn't break callers.
"""

from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

import pytest


class TestRunSnakemakeRule:
    """Behaviour of the Snakemake-invocation helper."""

    def test_successful_run_does_not_exit(self) -> None:
        """returncode == 0 returns normally; no sys.exit is called."""
        from RetroSeek import run_snakemake_rule

        completed = MagicMock(returncode=0, stderr=None, stdout=None)
        with patch("RetroSeek.subprocess.run", return_value=completed) as run:
            run_snakemake_rule(rule="probe_extractor", num_cores=2, display_info=True)
        run.assert_called_once()
        cmd = run.call_args.args[0]
        assert cmd[:2] == ["snakemake", "probe_extractor"]
        assert "--cores" in cmd
        assert "2" in cmd
        assert "--rerun-incomplete" in cmd

    def test_quiet_flag_appended_when_display_info_false(self) -> None:
        """display_info=False adds ``-q`` to the snakemake invocation."""
        from RetroSeek import run_snakemake_rule

        completed = MagicMock(returncode=0)
        with patch("RetroSeek.subprocess.run", return_value=completed) as run:
            run_snakemake_rule(rule="probe_extractor", num_cores=1, display_info=False)
        assert "-q" in run.call_args.args[0]

    def test_passthrough_flags_are_forwarded(self) -> None:
        """Extra flags are appended verbatim to the snakemake command."""
        from RetroSeek import run_snakemake_rule

        completed = MagicMock(returncode=0)
        passthrough = ["--keep-going", "--latency-wait", "60"]
        with patch("RetroSeek.subprocess.run", return_value=completed) as run:
            run_snakemake_rule(
                rule="ltr_harvester",
                num_cores=4,
                display_info=True,
                snakemake_flags=passthrough,
            )
        for flag in passthrough:
            assert flag in run.call_args.args[0]

    def test_nonzero_exit_code_propagates_via_sys_exit(self) -> None:
        """Nonzero Snakemake exit code is propagated via sys.exit."""
        from RetroSeek import run_snakemake_rule

        completed = MagicMock(returncode=2, stderr=None)
        with (
            patch("RetroSeek.subprocess.run", return_value=completed),
            pytest.raises(SystemExit) as excinfo,
        ):
            run_snakemake_rule(rule="ltr_digester", num_cores=1, display_info=True)
        assert excinfo.value.code == 2

    def test_nonzero_exit_does_not_raise_on_none_stderr(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Failure path must not raise when stderr wasn't captured.

        Previously the code did ``logging.error(result.stderr)`` with
        stderr = None because capture_output wasn't set. The fix is to not
        reference ``result.stderr`` at all - log the rule + exit code and
        rely on Snakemake's own terminal output.
        """
        from RetroSeek import run_snakemake_rule

        completed = MagicMock(returncode=3, stderr=None)
        with (
            caplog.at_level(logging.ERROR),
            patch("RetroSeek.subprocess.run", return_value=completed),
            pytest.raises(SystemExit),
        ):
            run_snakemake_rule(rule="ranges_analysis", num_cores=1, display_info=True)

        # We expect at least one ERROR record mentioning the rule name or
        # exit code. We do NOT want ``None`` to appear as a bare log line.
        error_messages = [
            rec.getMessage() for rec in caplog.records if rec.levelno == logging.ERROR
        ]
        assert error_messages, "expected at least one ERROR log on failure"
        assert any("ranges_analysis" in m for m in error_messages), (
            f"expected the failed rule name in error messages, got: {error_messages}"
        )
        assert not any(m.strip() == "None" for m in error_messages), (
            "bare 'None' log line means the code is logging an un-captured stderr"
        )

    def test_missing_snakemake_binary_is_handled(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A missing ``snakemake`` executable must surface as a clear error.

        ``FileNotFoundError`` from ``subprocess.run`` should be logged with
        the rule name + original exception, then sys.exit(1).
        """
        from RetroSeek import run_snakemake_rule

        with (
            caplog.at_level(logging.ERROR),
            patch(
                "RetroSeek.subprocess.run", side_effect=FileNotFoundError("snakemake")
            ),
            pytest.raises(SystemExit) as excinfo,
        ):
            run_snakemake_rule(rule="pair_detector", num_cores=1, display_info=True)
        assert excinfo.value.code == 1
        joined = " ".join(rec.getMessage() for rec in caplog.records)
        assert "pair_detector" in joined


class TestSnakemakeOptions:
    """How the launcher turns its own options into Snakemake's (ADR-020)."""

    def test_keep_going_is_the_default(self) -> None:
        from RetroSeek import snakemake_options

        assert snakemake_options(["-n"], stop_on_error=False) == ["--keep-going", "-n"]

    def test_stop_on_error_drops_keep_going(self) -> None:
        from RetroSeek import snakemake_options

        assert snakemake_options(["-n"], stop_on_error=True) == ["-n"]

    def test_user_options_keep_their_order_so_forcerun_stays_last(self) -> None:
        from RetroSeek import snakemake_options

        user = ["--configfile", "/c.yaml", "--forcerun", "rule_a", "rule_b"]
        assert snakemake_options(user, stop_on_error=False)[-5:] == user

    @pytest.mark.parametrize("flag", ["-n", "--dry-run", "--dryrun"])
    def test_dry_run_is_recognised(self, flag: str) -> None:
        from RetroSeek import is_dry_run

        assert is_dry_run(["--configfile", "/c.yaml", flag])

    def test_no_dry_run_without_the_flag(self) -> None:
        from RetroSeek import is_dry_run

        assert not is_dry_run(["--configfile", "/c.yaml"])

    @pytest.mark.parametrize("flag", ["--unlock", "--cleanup-metadata", "--cm"])
    def test_maintenance_modes_skip_the_guard(self, flag: str) -> None:
        """--unlock and --cleanup-metadata run nothing; a guard dry run with them
        would itself unlock or clean."""
        from RetroSeek import is_maintenance

        assert is_maintenance([flag, "/some/file"])

    def test_a_normal_run_is_not_maintenance(self) -> None:
        from RetroSeek import is_maintenance

        assert not is_maintenance(["--configfile", "/c.yaml", "-n"])


class TestCaptureDryRun:
    """The guard's dry run must ask Snakemake the same question as the real run."""

    def test_dry_run_matches_the_real_run_options(self) -> None:
        from RetroSeek import capture_dry_run

        done = MagicMock(returncode=0, stdout="Job stats:\n", stderr="")
        with patch("RetroSeek.subprocess.run", return_value=done) as run:
            capture_dry_run(
                ["taxonomy_classify"], ["--keep-going", "--configfile", "/c.yaml"]
            )
        cmd = run.call_args.args[0]
        assert cmd[:2] == ["snakemake", "taxonomy_classify"]
        assert "-n" in cmd
        # Without it, incomplete files from an interrupted run make the dry run
        # fail where the real run would simply redo them.
        assert "--rerun-incomplete" in cmd
        assert cmd[-3:] == ["--keep-going", "--configfile", "/c.yaml"]

    def test_a_failing_dry_run_stops_with_its_code(self) -> None:
        from RetroSeek import capture_dry_run

        failed = MagicMock(returncode=3, stdout="", stderr="DAG error")
        with (
            patch("RetroSeek.subprocess.run", return_value=failed),
            pytest.raises(SystemExit) as caught,
        ):
            capture_dry_run(["x"], [])
        assert caught.value.code == 3
