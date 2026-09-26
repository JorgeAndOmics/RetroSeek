# =============================================================================
# stages.py
# =============================================================================
# The launcher's stages, in one table (ADR-020).
#
# Every stage flag of `./RetroSeek` is one row below: its phase, the Snakemake
# targets it asks for, its help text, the heavy rules it is allowed to run and the
# external tools it needs. The command-line parser, the run order, the heavy-rule
# guard and the tool check are all built from this table, so they cannot drift
# apart the way the two hand-written parsers once did.
#
# Pure data plus the parser: this module imports nothing from the pipeline (not
# `defaults`, which creates directories when imported), so `./RetroSeek -h`
# works without a config.
# =============================================================================

"""The launcher's stages in one table, and the command-line parser built from it."""

from __future__ import annotations

import argparse
from dataclasses import dataclass


@dataclass(frozen=True)
class Stage:
    """One launcher stage.

    Attributes:
        flag: The command-line flag, e.g. "--classify".
        phase: The help section it is listed under (see PHASES).
        targets: Snakemake target rules the stage asks for.
        help: One line for `./RetroSeek -h`.
        heavy: Heavy rules this stage exists to run. Any other stage that
            would run one of them is stopped by the guard (see guard.py).
        tools: Executables the stage calls, checked before a run starts.
    """

    flag: str
    phase: str
    targets: tuple[str, ...]
    help: str
    heavy: tuple[str, ...] = ()
    tools: tuple[str, ...] = ()

    @property
    def dest(self) -> str:
        """The attribute argparse stores the flag under ("--ltr-domains" -> "ltr_domains")."""
        return self.flag.lstrip("-").replace("-", "_")


PHASES: dict[str, str] = {
    "Setup": "Fetch and prepare inputs. Run once per study.",
    "Indexing": "Build the per-genome search indexes.",
    "Discovery": "The heavy searches: about a day per genome for LTRdigest.",
    "Analysis": "Everything that reads the discovery results.",
    "Figures": "Stage-wide figure PDFs.",
}

STAGES: tuple[Stage, ...] = (
    # --- Setup
    Stage(
        "--download-genomes",
        "Setup",
        ("genome_downloader",),
        "Download the genomes listed under `species:` (skips files already present).",
        heavy=("genome_downloader_setup", "genome_downloader"),
        tools=("datasets", "jq", "unzip"),
    ),
    Stage(
        "--download-hmm",
        "Setup",
        ("pfam_hmm_downloader",),
        "Download the pinned Pfam release (input.pfam_release).",
        heavy=("pfam_hmm_downloader",),
        tools=("wget",),
    ),
    Stage(
        "--build-reference",
        "Setup",
        ("taxonomy_reference_trees",),
        "Build the classification reference: reference proteins, NCBI taxonomy "
        "and per-gene trees (network; built once).",
        tools=("mafft", "iqtree", "raxml-ng"),
    ),
    Stage(
        "--probe-extractor",
        "Setup",
        ("probe_extractor",),
        "Fetch the probe sequences named in the probe CSV from NCBI.",
        heavy=("probe_extractor",),
    ),
    # --- Indexing
    Stage(
        "--blast-dbs",
        "Indexing",
        ("blast_db_generator",),
        "Build a BLAST database per genome.",
        tools=("makeblastdb",),
    ),
    Stage(
        "--suffix-arrays",
        "Indexing",
        ("ltr_index_generator",),
        "Build a GenomeTools suffix array per genome.",
        heavy=("ltr_index_generator_setup", "ltr_index_generator"),
        tools=("gt",),
    ),
    # --- Discovery
    Stage(
        "--ltr-candidates",
        "Discovery",
        ("ltr_harvester",),
        "Find LTR element candidates with LTRharvest.",
        heavy=("ltr_harvester_setup", "ltr_harvester"),
        tools=("gt",),
    ),
    Stage(
        "--ltr-domains",
        "Discovery",
        ("ltr_digester",),
        "LTR element annotation with LTRdigest (protein domains and polypurine "
        "tracts inside each element).",
        heavy=("ltr_digester_setup", "ltr_digester"),
        tools=("gt",),
    ),
    Stage(
        "--blast",
        "Discovery",
        ("blast_pkl2parquet",),
        "Search every probe against every genome with tBLASTn.",
        # blast_pkl2parquet is the checkpoint that decides which genomes have
        # hits. Until it has run, Snakemake cannot list the jobs after it, so a
        # rerun would hide heavy jobs from the guard: it is guarded like one.
        heavy=("full_genome_blaster_setup", "full_genome_blaster", "blast_pkl2parquet"),
        tools=("tblastn",),
    ),
    # --- Analysis
    Stage(
        "--ranges-analysis",
        "Analysis",
        ("ranges_analysis",),
        "Join tBLASTn hits and LTR elements into element hits and orphans.",
        tools=("Rscript",),
    ),
    Stage(
        "--domain-scan",
        "Analysis",
        ("domain_scanner",),
        "Scan element and orphan loci for curated Pfam domains (hmmsearch).",
        tools=("hmmsearch", "Rscript"),
    ),
    Stage(
        "--classify",
        "Analysis",
        (
            "taxonomy_classify",
            "taxonomy_orphans",
            "taxonomy_plot_generator",
            "loss_analysis",
        ),
        "Classify each locus to a viral genus, plus the taxonomy and loss figures.",
        tools=(
            "blastx",
            "makeblastdb",
            "mafft",
            "epa-ng",
            "gappa",
            "hmmsearch",
            "Rscript",
        ),
    ),
    Stage(
        "--segment",
        "Analysis",
        ("taxonomy_segments",),
        "Split the catalog by taxon (classification.segment_rank), with a PDF each.",
        tools=("Rscript",),
    ),
    Stage(
        "--solo-ltr-detector",
        "Analysis",
        ("solo_ltr_detector",),
        "Find solo LTRs from the arms of ERV-bearing elements, with the evidence tree.",
        tools=("blastn", "makeblastdb", "mafft", "iqtree", "Rscript"),
    ),
    Stage(
        "--hotspot-detection",
        "Analysis",
        ("hotspot_detector",),
        "Find genomic windows enriched in ERV integrations.",
        tools=("Rscript",),
    ),
    Stage(
        "--pair-detection",
        "Analysis",
        ("pair_detector",),
        "Find nearby probe pairs (parameters.probe_to_pair).",
        tools=("Rscript",),
    ),
    Stage(
        "--placement-trees",
        "Analysis",
        ("placement_trees",),
        "Placement figures: heat-trees, placement uncertainty, co-phylogeny.",
        tools=("gappa",),
    ),
    # --- Figures
    Stage(
        "--generate-global-plots",
        "Figures",
        ("plot_generator", "stage_plot_generator", "erv_like_plot_generator"),
        "The homology, integration and structure PDFs.",
        tools=("Rscript",),
    ),
)

# `--downstream`: everything after the heavy discovery searches.
DOWNSTREAM: tuple[str, ...] = tuple(
    s.flag for s in STAGES if s.phase in ("Analysis", "Figures")
)

HEAVY_RULES: frozenset[str] = frozenset(rule for s in STAGES for rule in s.heavy)


def build_parser() -> argparse.ArgumentParser:
    """The `./RetroSeek` parser: one section per phase, then presets and options.

    Options it does not know are handed to Snakemake unchanged, which is why
    abbreviations are off: `--class` must not quietly mean `--classify`.
    """
    parser = argparse.ArgumentParser(
        prog="RetroSeek",
        description="RetroSeek: directed ERV detection and analysis.",
        allow_abbrev=False,
        epilog=(
            "Stages run as one Snakemake workflow, so their order on the command "
            "line does not matter. Any other option goes to Snakemake, e.g. -n "
            "(dry run), --configfile FILE, --forcerun RULE (put it last)."
        ),
    )
    for phase, description in PHASES.items():
        group = parser.add_argument_group(phase, description)
        for stage in STAGES:
            if stage.phase == phase:
                group.add_argument(stage.flag, action="store_true", help=stage.help)

    presets = parser.add_argument_group("Presets")
    presets.add_argument(
        "--downstream",
        action="store_true",
        help="Every Analysis and Figures stage: " + ", ".join(DOWNSTREAM) + ".",
    )

    options = parser.add_argument_group("Run options")
    options.add_argument(
        "-skp",
        "--skip-validation",
        action="store_true",
        help="Skip the slow checks: NCBI probe lookups and the prompts. The fast "
        "checks (config, tools, Pfam) always run.",
    )
    options.add_argument(
        "--allow-heavy",
        action="store_true",
        help="Let a heavy rule run even though its own stage was not requested.",
    )
    options.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop at the first failed job. By default the other jobs finish.",
    )
    options.add_argument(
        "--verbosity",
        choices=("quiet", "normal", "verbose"),
        default=None,
        help="How much the terminal shows for this run (overrides display.verbosity).",
    )
    options.add_argument(
        "--config-help",
        nargs="?",
        const=None,
        default=argparse.SUPPRESS,
        metavar="KEY",
        help="Print the documentation of one config field (or list them all) and exit.",
    )
    return parser


def selected(args: argparse.Namespace) -> list[Stage]:
    """The stages the command line asks for, in table order, presets expanded."""
    flags = {s.flag for s in STAGES if getattr(args, s.dest, False)}
    if getattr(args, "downstream", False):
        flags.update(DOWNSTREAM)
    return [s for s in STAGES if s.flag in flags]


def targets(stages: list[Stage]) -> list[str]:
    """Snakemake targets of `stages`, without repeats, in table order."""
    seen: list[str] = []
    for stage in stages:
        seen.extend(t for t in stage.targets if t not in seen)
    return seen


def allowed_heavy(stages: list[Stage]) -> set[str]:
    """Heavy rules that `stages` exist to run, so the guard lets them through."""
    return {rule for stage in stages for rule in stage.heavy}


def tools(stages: list[Stage]) -> list[str]:
    """Executables `stages` need, without repeats, in table order."""
    seen: list[str] = []
    for stage in stages:
        seen.extend(t for t in stage.tools if t not in seen)
    return seen
