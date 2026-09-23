"""Unit tests for workflow/scripts/domains/pfam_release.py.

The release record answers one question a year later: which Pfam release produced
these domain calls? The tests pin how the answer is built, without the network:
the checksum of the local download is compared with the pinned release's own
checksum list.
"""

import pfam_release

LOCAL_MD5 = (
    "7ab3c4e215d0daaea3004e37c4e24f8a  Pfam-A.hmm.gz\n"
    "6c7beb9426eb0b1075e778c20b842b24  Pfam.version.gz\n"
)
VERSION = "Pfam release       : 38.2\nPfam-A families    : 30134\n"


def test_checksum_for_finds_the_named_file() -> None:
    assert (
        pfam_release.checksum_for(LOCAL_MD5, "Pfam-A.hmm.gz")
        == "7ab3c4e215d0daaea3004e37c4e24f8a"
    )


def test_checksum_for_returns_none_when_the_file_is_not_listed() -> None:
    assert pfam_release.checksum_for(LOCAL_MD5, "Pfam-B.hmm.gz") is None


def test_record_says_yes_when_the_local_download_is_the_pinned_release() -> None:
    text, matches = pfam_release.release_record("38.2", VERSION, LOCAL_MD5, LOCAL_MD5)
    assert matches
    assert "Pfam release       : 38.2" in text
    assert "pinned release: 38.2" in text
    assert "matches the pinned release: yes" in text


def test_record_says_no_when_the_download_came_from_another_release() -> None:
    other = LOCAL_MD5.replace("7ab3c4e2", "00000000")
    text, matches = pfam_release.release_record("38.2", VERSION, other, LOCAL_MD5)
    assert not matches
    assert "matches the pinned release: no" in text


def test_record_says_no_when_the_local_checksum_is_missing() -> None:
    text, matches = pfam_release.release_record("38.2", VERSION, "", LOCAL_MD5)
    assert not matches
    assert "matches the pinned release: no" in text


def test_release_url_points_at_the_pinned_release_folder() -> None:
    assert pfam_release.release_url("38.2") == (
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2"
    )
