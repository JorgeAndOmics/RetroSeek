"""Path-resolution invariants for ``defaults.py``.

Guards one prior defect:

``PATH_DICT["CONFIG_DIR"]`` was derived from ``DATA_DIR`` (i.e. from
``root.data_root_folder``), but the two files read from it -
``schema.yaml`` (``validator.py``) and ``erv_class.tsv``
(``Snakefile`` rule ``taxonomy_reference``) - are repo-shipped SOURCE
artifacts versioned alongside the code, not pipeline outputs.

With the committed ``config.yaml`` the defect is invisible, because
``data_root_folder: 'data'`` already resolves inside the repo. It only
bites when the data root points elsewhere (the dev dataset on
``/mnt/v``): ``CONFIG_DIR`` then names a directory that the
mkdir-at-import loop creates but never populates, and validation dies
with ``FileNotFoundError: .../data/config/schema.yaml``.

``defaults`` is stubbed by ``tests/conftest.py`` for every other test, so
these cases run the real module in a subprocess to exercise the genuine
import-time behaviour.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"

# Emitted by the probe below; keep in sync with the keys asserted on.
_PROBE = (
    "import json, defaults; "
    "print(json.dumps({k: str(v) for k, v in defaults.PATH_DICT.items()}))"
)


def _path_dict_with_data_root(data_root: Path) -> dict[str, str]:
    """Import the real ``defaults`` with roots pointed at ``data_root``.

    Returns ``PATH_DICT`` as plain strings. Runs out-of-process because
    importing ``defaults`` has side effects (it mkdirs every path) and
    because the conftest stub owns the in-process ``defaults`` name.
    """
    config = yaml.safe_load((REPO_ROOT / "data" / "config" / "config.yaml").read_text())
    config["root"] = {
        "db_root_folder": str(data_root / "species"),
        "data_root_folder": str(data_root / "data"),
        "results_root_folder": str(data_root / "results"),
        "logs_root_folder": str(data_root / "logs"),
    }
    config_path = data_root / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))

    completed = subprocess.run(
        [sys.executable, "-c", _PROBE],
        capture_output=True,
        text=True,
        check=True,
        env={
            "PATH": "/usr/bin:/bin",
            "PYTHONPATH": str(SCRIPTS_DIR),
            "RETROSEEK_CONFIG": str(config_path),
        },
    )
    return json.loads(completed.stdout)


class TestConfigDirResolution:
    """``CONFIG_DIR`` must follow the repo, not the data root."""

    def test_config_dir_is_repo_relative_when_data_root_is_external(
        self, tmp_path: Path
    ) -> None:
        """An external data root must not move CONFIG_DIR away from the repo."""
        path_dict = _path_dict_with_data_root(tmp_path)

        assert Path(path_dict["CONFIG_DIR"]) == REPO_ROOT / "data" / "config"

    def test_repo_source_artifacts_exist_at_config_dir(self, tmp_path: Path) -> None:
        """The two files read from CONFIG_DIR must actually be there.

        This is the assertion that would have caught the original failure:
        ``validator.py`` reads ``schema.yaml`` and the ``taxonomy_reference``
        rule copies ``erv_class.tsv``, both from ``CONFIG_DIR``.
        """
        config_dir = Path(_path_dict_with_data_root(tmp_path)["CONFIG_DIR"])

        assert (config_dir / "schema.yaml").is_file()
        assert (config_dir / "erv_class.tsv").is_file()

    def test_output_dirs_still_follow_the_data_root(self, tmp_path: Path) -> None:
        """Pipeline OUTPUT paths must keep tracking the data root.

        Pins the distinction the fix rests on, so a later change cannot
        make everything repo-relative.
        """
        path_dict = _path_dict_with_data_root(tmp_path)

        assert Path(path_dict["PICKLE_DIR"]).is_relative_to(tmp_path)
        assert Path(path_dict["LOG_DIR"]).is_relative_to(tmp_path)
