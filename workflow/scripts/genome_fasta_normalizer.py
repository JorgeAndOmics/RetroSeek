"""Normalize genome FASTA filenames to canonical ``{genome}.fa`` via symlink.

Why this exists
---------------
Genome FASTAs arrive at RetroSeek with several plausible extensions:
``.fa`` (legacy), ``.fna`` (NCBI Datasets — the common case), ``.fasta``
(manual download convention), and ``.ffn`` (older GenBank exports).
Every downstream rule (BLAST DB build, suffixerator, LTRharvest,
LTR_retriever) hard-codes ``.fa`` as the input filename. This rule
canonicalises whatever shape upstream supplied into a ``{genome}.fa``
view of the file — without copying the (potentially gigabyte-scale)
genome.

Strategy: symlink. Repeat invocations are idempotent — a correct
symlink is left in place; a stale one is replaced atomically. A real
``.fa`` file is never overwritten.

Ambiguity policy
----------------
Extension preference, in order: ``fa`` > ``fna`` > ``fasta`` > ``ffn``.

If ``.fa`` is absent and *exactly one* of the other three is present,
that file becomes the symlink target. If two or more non-``.fa`` variants
coexist, the script refuses with a ``RuntimeError`` listing the
candidates — the user must disambiguate by removing duplicates.
Silent preference would mask the case where a working `.fna` from
NCBI Datasets accidentally coexists with an older `.fasta` from a
prior pipeline that pointed at a different genome.

CLI
---
::

    python genome_fasta_normalizer.py \\
        --genome-name <str> \\
        --species-dir <path>     # SPECIES_DB
        --output <path>          # SPECIES_DB / {genome}.fa
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

EXT_PREFERENCE: tuple[str, ...] = ("fa", "fna", "fasta", "ffn")


def pick_canonical_source(species_dir: Path, genome: str) -> Path:
    """Return the FASTA file that should back ``{genome}.fa``.

    Walks ``EXT_PREFERENCE`` in order. ``.fa`` always wins when present.
    For the non-``.fa`` extensions, refuses if more than one is present.

    Raises
    ------
    FileNotFoundError
        No file with any of the four extensions exists for ``genome``.
    RuntimeError
        Two or more non-``.fa`` variants coexist (user must disambiguate).
    """
    fa_path = species_dir / f"{genome}.fa"
    if fa_path.exists():
        return fa_path

    candidates = [
        species_dir / f"{genome}.{ext}"
        for ext in EXT_PREFERENCE[1:]
        if (species_dir / f"{genome}.{ext}").exists()
    ]
    if not candidates:
        tried = ", ".join(f"{genome}.{ext}" for ext in EXT_PREFERENCE)
        raise FileNotFoundError(
            f"No genome FASTA found for {genome!r} in {species_dir}. "
            f"Looked for: {tried}"
        )
    if len(candidates) > 1:
        names = ", ".join(p.name for p in candidates)
        raise RuntimeError(
            f"Ambiguous genome FASTA for {genome!r}: multiple variants present "
            f"({names}). Remove duplicates so exactly one .fa/.fna/.fasta/.ffn "
            f"file remains, or place the canonical .fa explicitly."
        )
    return candidates[0]


def _validate_fasta_first_byte(path: Path) -> None:
    """Cheap sanity check: a FASTA must start with ``>`` after whitespace.

    Catches misnamed binary blobs (gzip, tarballs) before downstream
    tools choke on them with cryptic errors. Reads at most the first
    16 bytes — fast even on huge files.
    """
    target = path.resolve() if path.is_symlink() else path
    with target.open("rb") as handle:
        head = handle.read(16)
    stripped = head.lstrip()
    if not stripped or stripped[0:1] != b">":
        raise RuntimeError(
            f"{path} does not look like a FASTA file (first non-whitespace "
            f"byte must be '>'). Got: {head!r}"
        )


def normalize(species_dir: Path, genome: str, output: Path) -> Path:
    """Make ``output`` (= ``{genome}.fa``) point at the canonical source.

    - If ``output`` is already a regular file (not a symlink), leave it.
    - If ``output`` is a correct symlink, leave it.
    - Otherwise replace the symlink (or create a new one) pointing at the
      file ``pick_canonical_source`` selected.

    Validates the FASTA shape before creating any symlink so a misnamed
    non-FASTA fails loudly here rather than tens of CPU-minutes later.
    """
    source = pick_canonical_source(species_dir, genome)

    # Real-file fast path: an existing regular file at output is the canonical source.
    if output.is_file() and not output.is_symlink() and source == output:
        _validate_fasta_first_byte(output)
        return output

    _validate_fasta_first_byte(source)

    # Replace any prior link / file (other than the validated real-file case above).
    if output.is_symlink() or output.exists():
        if output.is_symlink():
            try:
                if output.resolve() == source.resolve() and output.exists():
                    return output  # already correct
            except (OSError, RuntimeError):
                # Stale symlink whose target can't be resolved.
                pass
        output.unlink()

    output.parent.mkdir(parents=True, exist_ok=True)
    output.symlink_to(source.resolve())
    return output


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--genome-name", required=True)
    parser.add_argument("--species-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    canonical = normalize(args.species_dir, args.genome_name, args.output)
    logging.info("Canonical FASTA for %s: %s", args.genome_name, canonical)
    return 0


if __name__ == "__main__":
    sys.exit(main())
