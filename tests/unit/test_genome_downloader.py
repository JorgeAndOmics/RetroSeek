"""genome_downloader.sh, run with a stand-in `datasets` (no network).

The script must never report success when it could not record which assembly it
downloaded: the record is the only trace of the choice it made.
"""

import os
import stat
import subprocess
import sys
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "workflow"
    / "scripts"
    / "blast_search"
    / "genome_downloader.sh"
)

# A stand-in for `datasets download genome accession ... --filename F`: writes a zip
# holding one FASTA where the real archive keeps it.
FAKE_DATASETS = f"""#!{sys.executable}
import sys, zipfile
out = sys.argv[sys.argv.index("--filename") + 1]
with zipfile.ZipFile(out, "w") as z:
    z.writestr("ncbi_dataset/data/GCF_1/genome.fna", ">chr1\\nACGT\\n")
"""


def _run(tmp_path: Path, log: Path) -> subprocess.CompletedProcess[str]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "datasets"
    fake.write_text(FAKE_DATASETS)
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    out_dir = tmp_path / "genomes"
    out_dir.mkdir()
    env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}", "NCBI_API_KEY": "x"}
    return subprocess.run(
        ["bash", str(SCRIPT), "GCF_000001.1", str(out_dir), str(log)],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def test_a_download_is_recorded(tmp_path: Path) -> None:
    log = tmp_path / "download_log.log"
    result = _run(tmp_path, log)
    assert result.returncode == 0, result.stderr
    assert log.read_text() == "GCF_000001.1\tGCF_000001.1\tDirect_Accession\n"
    assert (tmp_path / "genomes" / "GCF_000001.1.fa").read_text() == ">chr1\nACGT\n"
    assert " OK genome_downloader GCF_000001.1 | " in result.stderr


def test_an_unwritable_record_fails_the_job(tmp_path: Path) -> None:
    """The old bug: the record path was a folder; the script said OK anyway."""
    log = tmp_path / "download_log.log"
    log.mkdir()
    result = _run(tmp_path, log)
    assert result.returncode != 0
    assert " ERROR genome_downloader " in result.stderr
    assert " OK " not in result.stderr


def test_the_api_key_reaches_datasets_as_one_option(tmp_path: Path) -> None:
    """The key goes to `datasets` as `--api-key <key>`: two words, not one, not split."""
    log = tmp_path / "download_log.log"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "datasets"
    fake.write_text(
        FAKE_DATASETS
        + "open(sys.argv[0] + '.args', 'w').write('\\n'.join(sys.argv[1:]))\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    out_dir = tmp_path / "genomes"
    out_dir.mkdir()
    env = {
        **os.environ,
        "PATH": f"{bin_dir}:{os.environ['PATH']}",
        "NCBI_API_KEY": "abc 123",
    }
    subprocess.run(
        ["bash", str(SCRIPT), "GCF_000001.1", str(out_dir), str(log)],
        capture_output=True,
        text=True,
        env=env,
        check=True,
    )
    args = (bin_dir / "datasets.args").read_text().split("\n")
    assert args[args.index("--api-key") + 1] == "abc 123"
