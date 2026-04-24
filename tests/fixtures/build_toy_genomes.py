"""Build toy genomes for RetroSeek iteration.

Generates a small, deterministic set of synthetic ERV-containing genomes for
rapid pipeline testing. Each "species" gets a multi-scaffold FASTA with:

- Random background DNA (stable, seeded).
- Multiple planted LTR retrotransposons (matched direct repeats flanking
  an internal region seeded with real retroviral ORF fragments).
- Probe-matchable regions so tBLASTn produces non-trivial hits.

Output goes to ``/mnt/v/databases/toy-genomes/<species>/<species>.fa`` by
default, alongside a small probe CSV at ``.../toy_probes.csv`` suitable for
``config.input.probe_csv``.

Run inside the enERVate conda env::

    conda activate enERVate
    python tests/fixtures/build_toy_genomes.py

Pass ``--dest /somewhere/else`` to target a different directory, ``--force``
to overwrite existing outputs, ``--seed N`` to change the RNG seed.

The script is intentionally dependency-light (biopython only) and
reproducible — same seed, same output bytes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Real short retroviral ORF fragments (amino acid strings).
# Sources: trimmed fragments of published MLV, HERV-K, HIV-1 proteins —
# enough to produce non-trivial tBLASTn hits when back-translated.
PROBE_PROTEINS: dict[str, str] = {
    "GAG": (
        # Fragment of MLV Gag (p30/capsid-like region).
        "PKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEE"
        "MMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTA"
    ),
    "POL": (
        # Reverse-transcriptase-like fragment (RVT motifs).
        "FPISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYN"
        "TPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAY"
    ),
    "ENV": (
        # Transmembrane / fusion loop fragment.
        "QLNRLVSGKRCCFYADHTGVVRDSMAKLREQLEKRTQELLHTPYEPFRVDCGGDDWPI"
    ),
    "PRO": (
        # Aspartic protease active-site motif (DTG/LDTG).
        "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQ"
    ),
}

# Back-translation table: pick one codon per AA. Deterministic.
AA_TO_CODON: dict[str, str] = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
    "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
    "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
    "*": "TAA",
}


def back_translate(protein: str) -> str:
    """Back-translate a protein sequence to DNA using the canonical codon table."""
    return "".join(AA_TO_CODON[aa] for aa in protein if aa in AA_TO_CODON)


def random_dna(length: int, rng: random.Random) -> str:
    """Return a random DNA string of the given length."""
    return "".join(rng.choices("ACGT", k=length))


def make_ltr_pair(length: int, rng: random.Random) -> str:
    """Return a direct-repeat sequence of the given length.

    LTRharvest looks for direct repeats flanking an internal region. Both
    flanks of a planted ERV use the same string so the repeat is detectable.
    """
    return random_dna(length, rng)


@dataclass(frozen=True)
class ERV:
    """A planted ERV: two LTRs bracketing an internal region."""

    ltr: str
    internal: str

    @property
    def sequence(self) -> str:
        """Return the full ERV sequence (LTR + internal + LTR)."""
        return self.ltr + self.internal + self.ltr

    def __len__(self) -> int:
        """Return the total length (both LTRs + internal)."""
        return 2 * len(self.ltr) + len(self.internal)


def build_erv(
    rng: random.Random,
    ltr_length: int,
    probe_order: list[str],
    padding_per_gap: int,
) -> ERV:
    """Construct one ERV with back-translated retroviral ORFs inside.

    Args:
        rng: Seeded RNG.
        ltr_length: Length of each flanking LTR (bp).
        probe_order: Names of probes (keys of PROBE_PROTEINS) to embed,
            in 5'→3' order inside the internal region.
        padding_per_gap: Random DNA filler length between ORFs.

    Returns:
        ERV with direct-repeat LTRs and an internal region containing the
        back-translated ORFs separated by random filler.
    """
    ltr = make_ltr_pair(ltr_length, rng)
    parts: list[str] = [random_dna(padding_per_gap, rng)]
    for probe in probe_order:
        parts.append(back_translate(PROBE_PROTEINS[probe]))
        parts.append(random_dna(padding_per_gap, rng))
    return ERV(ltr=ltr, internal="".join(parts))


def build_scaffold(
    rng: random.Random,
    length: int,
    ervs: list[ERV],
) -> str:
    """Place ERVs into a random-DNA scaffold at quasi-uniform positions.

    ERVs are inserted in 5'→3' order with random DNA filler between them.
    The resulting scaffold length is approximately but not exactly ``length``.
    """
    total_erv_length = sum(len(erv) for erv in ervs)
    if total_erv_length >= length:
        msg = f"ERVs total {total_erv_length} bp exceed scaffold length {length}"
        raise ValueError(msg)
    filler_total = length - total_erv_length
    num_gaps = len(ervs) + 1
    base_gap = filler_total // num_gaps
    extra = filler_total - base_gap * num_gaps
    # Distribute "extra" bp across the first gaps so the total length matches.
    gaps = [base_gap + (1 if i < extra else 0) for i in range(num_gaps)]

    parts: list[str] = [random_dna(gaps[0], rng)]
    for i, erv in enumerate(ervs, start=1):
        parts.append(erv.sequence)
        parts.append(random_dna(gaps[i], rng))
    return "".join(parts)


@dataclass(frozen=True)
class SpeciesSpec:
    """Top-level plan for one toy species."""

    name: str                       # Folder / FASTA stem. Snake_case.
    scaffold_count: int             # How many chromosomes/scaffolds.
    scaffold_length: int            # bp per scaffold (approximate).
    ervs_per_scaffold: int          # Number of planted ERVs per scaffold.
    ltr_length: int                 # bp per LTR (direct repeat).
    probe_order: list[str]          # Probes (keys of PROBE_PROTEINS) per ERV.
    padding_per_gap: int            # bp of random filler between ORFs.


def build_species(spec: SpeciesSpec, rng: random.Random) -> list[SeqRecord]:
    """Return a list of SeqRecords (one per scaffold) for the species."""
    records: list[SeqRecord] = []
    for i in range(1, spec.scaffold_count + 1):
        ervs = [
            build_erv(
                rng=rng,
                ltr_length=spec.ltr_length,
                probe_order=spec.probe_order,
                padding_per_gap=spec.padding_per_gap,
            )
            for _ in range(spec.ervs_per_scaffold)
        ]
        sequence = build_scaffold(
            rng=rng,
            length=spec.scaffold_length,
            ervs=ervs,
        )
        seq_id = f"{spec.name}_scaffold_{i}"
        records.append(
            SeqRecord(
                Seq(sequence),
                id=seq_id,
                name=seq_id,
                description=f"toy scaffold {i} for {spec.name}",
            )
        )
    return records


def write_fasta(records: list[SeqRecord], destination: Path) -> int:
    """Write records to a FASTA file, following RetroSeek's SPECIES_DB layout.

    RetroSeek's convention (matching bat-experimental) is FASTA flat at
    ``SPECIES_DB/{genome}.fa`` with per-genome subdirectories reserved for
    tool-produced files (BLAST DB index, suffix-array index). Callers should
    pass ``destination`` as the flat path.

    Returns bytes written.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as handle:
        SeqIO.write(records, handle, "fasta")
    return destination.stat().st_size


def write_probe_csv(destination: Path) -> None:
    """Write a minimal probe CSV in the canonical schema.

    Canonical columns are ``Label, Name, Abbreviation, Probe, Accession``.
    Each row declares one probe; the ``Probe`` column is the AA fragment
    that was embedded in the toy genomes.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Label", "Name", "Abbreviation", "Probe", "Accession"])
        for probe, protein in PROBE_PROTEINS.items():
            writer.writerow(
                [
                    f"Toy_{probe}",                   # Label
                    f"Toy retrovirus {probe}",        # Name
                    probe,                            # Abbreviation
                    protein,                          # Probe AA sequence
                    f"TOY_{probe}_001",               # Accession (synthetic)
                ]
            )


def fasta_checksum(path: Path) -> str:
    """SHA256 of the FASTA bytes, for determinism verification."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def default_species_specs() -> list[SpeciesSpec]:
    """Three toy species spanning easy → moderately complex cases."""
    return [
        SpeciesSpec(
            name="Toyus_simplex",
            scaffold_count=2,
            scaffold_length=30_000,
            ervs_per_scaffold=2,
            ltr_length=200,
            probe_order=["GAG", "POL", "ENV"],
            padding_per_gap=150,
        ),
        SpeciesSpec(
            name="Toyus_complex",
            scaffold_count=3,
            scaffold_length=45_000,
            ervs_per_scaffold=3,
            ltr_length=250,
            probe_order=["GAG", "POL", "PRO", "ENV"],
            padding_per_gap=200,
        ),
        SpeciesSpec(
            name="Toyus_sparsus",
            scaffold_count=2,
            scaffold_length=60_000,
            ervs_per_scaffold=1,
            ltr_length=180,
            probe_order=["POL", "ENV"],
            padding_per_gap=300,
        ),
    ]


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path("/mnt/v/databases/toy-genomes"),
        help="Output root directory (default: /mnt/v/databases/toy-genomes).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1337,
        help="RNG seed. Same seed ⇒ byte-identical output (default: 1337).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing species FASTA files.",
    )
    args = parser.parse_args(argv)

    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    specs = default_species_specs()
    print(f"Building {len(specs)} toy species into {dest} (seed={args.seed})")

    manifest_rows: list[tuple[str, int, int, str]] = []
    for spec in specs:
        # Flat FASTA at dest/{name}.fa — matches RetroSeek SPECIES_DB layout.
        fasta_path = dest / f"{spec.name}.fa"
        if fasta_path.exists() and not args.force:
            print(f"  skip {spec.name}: {fasta_path} exists (pass --force to overwrite)")
            continue
        records = build_species(spec, rng)
        size = write_fasta(records, fasta_path)
        checksum = fasta_checksum(fasta_path)
        total_bp = sum(len(r.seq) for r in records)
        manifest_rows.append((spec.name, total_bp, size, checksum))
        print(
            f"  wrote {spec.name}: "
            f"{len(records)} scaffolds, {total_bp:,} bp, "
            f"{size:,} bytes — sha256 {checksum[:12]}…"
        )

    # Probe CSV always regenerated; trivially small.
    probe_csv = dest / "toy_probes.csv"
    write_probe_csv(probe_csv)
    print(f"  wrote probe CSV: {probe_csv}")

    # Manifest for quick diffing across regenerations.
    manifest_path = dest / "MANIFEST.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["species", "total_bp", "bytes", "sha256"])
        writer.writerows(manifest_rows)
    print(f"  wrote manifest: {manifest_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
