# =============================================================================
# pfam_release.py
# =============================================================================
# Record which Pfam release the domain scan used (ADR-019).
#
# Why this exists: Pfam-A.hmm carries no release number, and the pipeline once
# downloaded `current_release`, whatever that was on the day. A year later nobody
# could say which release produced a set of domain calls. The download is now
# pinned (config `input.pfam_release`), and this script writes a small record
# beside the curated subset:
#
#   pinned release: 38.2
#   Pfam release       : 38.2          <- EBI's own Pfam.version for that release
#   Pfam-A families    : 30134
#   ...
#   local Pfam-A.hmm.gz matches the pinned release: yes
#
# The last line compares the checksum the downloader saved (md5sum.txt, fetched
# with the file) with the pinned release's checksum list. "no" means the file on
# disk came from a different release: the record then says so rather than lying.
# The accession check in subset_pfam.py stays the hard gate; a mismatch here is a
# warning about provenance, not an error.
#
# Deliberately NOT an input of LTRdigest: a new file there would make every
# LTRdigest output look stale (about a day per genome to redo).
# =============================================================================

from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

PFAM_RELEASES = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases"
HMM_ARCHIVE = "Pfam-A.hmm.gz"


def release_url(release: str) -> str:
    """The EBI folder holding one Pfam release, e.g. `.../releases/Pfam38.2`."""
    return f"{PFAM_RELEASES}/Pfam{release}"


def checksum_for(md5_text: str, filename: str) -> str | None:
    """The md5 listed for `filename` in an `md5sum`-style text, or None."""
    for line in md5_text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == filename:
            return parts[0]
    return None


def release_record(
    release: str, version_text: str, local_md5_text: str, release_md5_text: str
) -> tuple[str, bool]:
    """Build the record text, and whether the local download is the pinned release.

    :param release: the pinned release, e.g. "38.2".
    :param version_text: EBI's Pfam.version for that release.
    :param local_md5_text: md5sum.txt saved by the downloader next to Pfam-A.hmm.
    :param release_md5_text: the pinned release's md5_checksums.
    """
    local = checksum_for(local_md5_text, HMM_ARCHIVE)
    matches = local is not None and local == checksum_for(release_md5_text, HMM_ARCHIVE)
    answer = "yes" if matches else "no"
    text = (
        f"pinned release: {release}\n"
        f"{version_text.rstrip()}\n"
        f"local {HMM_ARCHIVE} matches the pinned release: {answer}\n"
    )
    return text, matches


def _fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as response:
        data: bytes = response.read()
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description="Record the Pfam release in use.")
    parser.add_argument("--release", required=True, help="pinned release, e.g. 38.2")
    parser.add_argument(
        "--local-md5", type=Path, required=True, help="md5sum.txt from the downloader"
    )
    parser.add_argument("--out", type=Path, required=True, help="record to write")
    args = parser.parse_args()

    base = release_url(args.release)
    version_text = gzip.decompress(_fetch(f"{base}/Pfam.version.gz")).decode()
    release_md5_text = _fetch(f"{base}/md5_checksums").decode()
    local_md5_text = args.local_md5.read_text(encoding="utf-8")

    text, matches = release_record(
        args.release, version_text, local_md5_text, release_md5_text
    )
    args.out.write_text(text, encoding="utf-8")
    if not matches:
        print(
            f"WARNING: the local {HMM_ARCHIVE} is not Pfam {args.release}, the "
            "release in input.pfam_release. Domain calls will come from another "
            "release. Fix: ./RetroSeek --download-hmm --forcerun pfam_hmm_downloader",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
