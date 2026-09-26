"""Pipeline-wide utility helpers.

Small, dependency-light functions shared across ``workflow/scripts``:
dill-based pickle round-trips, directory listing / cleanup, random
identifier generation, and the incomplete-record filter used during
GenBank enrichment.

All path-accepting helpers accept either :class:`pathlib.Path` or a
plain ``str``; Snakefile ``params:`` blocks hand in strings, and
interactive callers tend to use ``Path`` objects. Both work.
"""

from __future__ import annotations

import logging
import random
import string
from pathlib import Path
from typing import Any

import dill

logger = logging.getLogger(__name__)


def pickler(
    data: Any, output_directory_path: str | Path, output_file_name: str
) -> None:
    """Write ``data`` to ``<output_directory_path>/<output_file_name>`` via dill.

    The directory is expected to exist; the function does not create it.
    dill is used (not stdlib ``pickle``) because some pipeline objects
    embed BioPython records that stdlib pickle cannot serialise.

    Args:
        data: Any dill-serialisable object.
        output_directory_path: Directory to write the pickle file into.
        output_file_name: Basename of the pickle file (including extension).
    """
    target = Path(output_directory_path) / output_file_name
    with target.open("wb") as handle:
        dill.dump(data, handle)


def unpickler(input_directory_path: str | Path, input_file_name: str) -> Any:
    """Load a dill pickle from ``<input_directory_path>/<input_file_name>``.

    Args:
        input_directory_path: Directory containing the pickle.
        input_file_name: Basename of the pickle file.

    Returns:
        The deserialised payload.

    Raises:
        Exception: If the file cannot be opened or dill fails to load it; the
            message names the file, and the original error is chained.
    """
    file_path = Path(input_directory_path) / input_file_name
    try:
        with file_path.open("rb") as handle:
            return dill.load(handle)
    except Exception as exc:
        logger.error(f"Failed to unpickle {file_path}")
        raise Exception(f"Failed to unpickle {file_path}") from exc


def incomplete_dict_cleaner(object_dict: dict[str, Any]) -> dict[str, Any]:
    """Drop dict entries whose value fails ``.is_complete()``.

    Used to purge records that lack a GenBank / alignment / HSP payload
    after remote enrichment.
    """
    logger.debug("Cleaning up incomplete objects...")

    return {key: value for key, value in object_dict.items() if value.is_complete()}


def random_string_generator(length: int) -> str:
    """Generate an ``A-Z0-9`` random string of the given length.

    Not cryptographically strong - callers use this for short, human-
    distinguishable identifiers (typically 6 characters).
    """
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=length))
