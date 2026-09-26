# ruff: noqa: PLC0415
# Per-test imports of the module under test are deliberate (conftest
# defaults stub must take effect first).

"""Unit tests for the pieces of ``RetroSeek.py`` that build Snakemake's command.

The launcher runs ONE Snakemake per run (ADR-020). These tests pin how that
command is put together: targets first, the launcher's options, then the user's
options in their order (``--forcerun`` swallows what follows), and a guard dry
run that asks exactly what the real run will.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest


class TestSnakemakeCommand:
    """The one Snakemake command line, shared by the real run and the guard."""

    def test_targets_cores_and_rerun_incomplete(self) -> None:
        from RetroSeek import snakemake_command

        cmd = snakemake_command(["probe_extractor"], cores=2, options=[])
        assert cmd == [
            "snakemake",
            "probe_extractor",
            "--cores",
            "2",
            "--rerun-incomplete",
        ]

    def test_options_come_last_in_order(self) -> None:
        """`--forcerun RULE...` swallows what follows, so nothing may follow it."""
        from RetroSeek import snakemake_command

        options = ["--keep-going", "--forcerun", "rule_a"]
        cmd = snakemake_command(["a", "b"], cores=4, options=options)
        assert cmd[:3] == ["snakemake", "a", "b"]
        assert cmd[-3:] == options

    def test_no_quiet_flag_whatever_the_verbosity(self) -> None:
        """Snakemake always speaks in full: the launcher filters, the run log keeps all."""
        from RetroSeek import snakemake_command

        assert "-q" not in snakemake_command(["a"], cores=1, options=[])


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
        """--unlock and --cleanup-metadata skip the guard dry run.

        They run nothing, and a guard dry run with them would itself unlock or
        clean.
        """
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

    def test_a_failing_dry_run_returns_its_code_and_output(self) -> None:
        from RetroSeek import capture_dry_run

        failed = MagicMock(returncode=3, stdout="", stderr="DAG error")
        with patch("RetroSeek.subprocess.run", return_value=failed):
            code, output = capture_dry_run(["x"], [])
        assert code == 3
        assert "DAG error" in output
