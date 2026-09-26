"""DAG-contract regression tests for the Snakefile.

Goal: pin structural invariants of the workflow (rule presence, wildcard
constraints) without requiring a snakemake binary on PATH for every
contributor. We parse the Snakefile with regex; that's brittle for
fancy syntax but the existing rules are all top-level and conventionally
formatted, so it works.

Tests for rules that don't exist yet (e.g. ``genome_fasta_normalizer_setup``)
are marked ``xfail(strict=True)`` - they flip to PASS when the rule lands
in its target phase, so anyone removing a planned rule prematurely will
trip a failure.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import stages

RULE_DEF_RE = re.compile(r"^rule\s+(\w+)\s*:", re.MULTILINE)
WILDCARD_BLOCK_RE = re.compile(
    r"wildcard_constraints\s*:\s*\n\s*genome\s*=\s*(.+?)$",
    re.MULTILINE,
)


def _read_snakefile(project_root: Path) -> str:
    return (project_root / "workflow" / "Snakefile").read_text()


def _all_rules(project_root: Path) -> set[str]:
    return set(RULE_DEF_RE.findall(_read_snakefile(project_root)))


# ---------------------------------------------------------------------
# Rule presence
# ---------------------------------------------------------------------
@pytest.mark.parametrize("rule_name", ["ltr_harvester_setup", "species_tree_layout"])
def test_existing_rule_present(project_root: Path, rule_name: str) -> None:
    """Rules the rest of the workflow depends on must remain in the Snakefile."""
    rules = _all_rules(project_root)
    assert rule_name in rules, f"rule {rule_name!r} missing from Snakefile"


@pytest.mark.parametrize(
    "rule_name",
    [
        "ltr_retriever_prefilter_setup",
        "ltr_retriever_setup",
        "solo_ltr_integrator_setup",
        "solo_intact_ratio_aggregate",
    ],
)
def test_ltr_retriever_rule_is_gone(project_root: Path, rule_name: str) -> None:
    """The LTR_retriever workstream is archived (ADR-017).

    Asserted as absence rather than deleted outright, so that a partial revert
    that reinstates one rule without the rest is caught.
    """
    assert rule_name not in _all_rules(project_root)


def test_ltr_harvest_still_writes_the_scn(project_root: Path) -> None:
    """The .scn was LTR_retriever's input and now has no consumer, but it stays.

    Removing an output from ltr_harvester_setup changes that rule's signature and
    re-fires a 24-hour stage for no gain.
    """
    text = (project_root / "workflow" / "Snakefile").read_text()
    harvester = text[
        text.index("rule ltr_harvester_setup") : text.index("rule ltr_digester_setup")
    ]
    assert "LTR_SCN_DIR" in harvester


def test_genome_fasta_normalizer_rule_present(project_root: Path) -> None:
    """``genome_fasta_normalizer_setup`` rule must be in the workflow."""
    assert "genome_fasta_normalizer_setup" in _all_rules(project_root)


# ---------------------------------------------------------------------
# Native solo-LTR detection (ADR-017)
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "rule_name",
    [
        "solo_bait_builder_setup",
        "solo_blaster_setup",
        "solo_finder_setup",
        "solo_annotator_setup",
        "solo_tree_setup",
        "solo_tree_views_setup",
        "solo_plot_generator_setup",
        "solo_plot_summary",
        "solo_ltr_detector",
    ],
)
def test_solo_ltr_rule_present(project_root: Path, rule_name: str) -> None:
    """Every stage of the solo-LTR chain must stay in the workflow."""
    assert rule_name in _all_rules(project_root)


def test_solo_thresholds_come_from_config_not_literals(project_root: Path) -> None:
    """Every acceptance threshold must be read from config, never hard-coded.

    The prototype had these as literals; moving them into config.yaml is what lets
    a threshold be changed in one place and stay consistent across the scripts.
    """
    text = (project_root / "workflow" / "Snakefile").read_text()
    finder = text[
        text.index("rule solo_finder_setup") : text.index("rule solo_finder:")
    ]
    for key in (
        "min_identity",
        "min_coverage",
        "max_coverage",
        "min_alignment_length",
        "min_hit_length",
        "merge_gap",
        "orphan_pad",
    ):
        assert f"_SOLO['{key}']" in finder, f"{key} is not read from config"


# ---------------------------------------------------------------------
# Domain evidence (ADR-015)
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "rule_name",
    ["pfam_subset_builder", "domain_scanner_setup", "domain_scanner"],
)
def test_domain_evidence_rule_present(project_root: Path, rule_name: str) -> None:
    """The domain scan rule must stay in the workflow.

    Without it every locus falls back to domain_source=not_scanned.
    """
    assert rule_name in _all_rules(project_root)


def test_domain_scan_does_not_rerun_ltrdigest(project_root: Path) -> None:
    """The domain scan reads finished LTRdigest tracks, never reruns gt ltrdigest.

    ADR-015 keeps LTRdigest untouched: rerunning it costs ~24 h per genome.
    """
    text = (project_root / "workflow" / "Snakefile").read_text(encoding="utf-8")
    start = text.index("rule domain_scanner_setup:")
    end = text.index("rule domain_scanner:")
    assert "gt ltrdigest" not in text[start:end]


def test_classification_consumes_the_domain_scan(project_root: Path) -> None:
    """Both classification tiers must be wired to the domain scan.

    Otherwise orphans silently regress to the hard-coded non_domain default
    this ADR removed.
    """
    text = (project_root / "workflow" / "Snakefile").read_text(encoding="utf-8")
    assert text.count("--domains-parquet") == 2
    assert text.count("--domains-scanned") == 2


def _rule_block(project_root: Path, start: str, end: str) -> str:
    text = (project_root / "workflow" / "Snakefile").read_text(encoding="utf-8")
    return text[text.index(start) : text.index(end)]


def test_pfam_download_is_pinned_to_a_release(project_root: Path) -> None:
    """The Pfam download must name a fixed release, not `current_release`.

    `current_release` changes under our feet: an older download once lacked 6
    accessions of the curated table and stopped the domain scan (ADR-019).
    """
    downloader = _rule_block(
        project_root, "rule pfam_hmm_downloader:", "rule pfam_release_recorder:"
    )
    assert "Pfam/current_release" not in downloader
    assert "releases/Pfam" in downloader


def test_pfam_downloader_outputs_are_unchanged(project_root: Path) -> None:
    """The Pfam downloader keeps exactly its two outputs.

    A new output would rerun the downloader, and a new Pfam-A.hmm makes every
    LTRdigest output stale (about a day per genome).
    """
    downloader = _rule_block(
        project_root, "rule pfam_hmm_downloader:", "rule pfam_release_recorder:"
    )
    outputs = downloader[downloader.index("output:") : downloader.index("params:")]
    assert outputs.count("os.path.join") == 2
    assert "'md5sum.txt'" in outputs
    assert "'Pfam-A.hmm'" in outputs


def test_pfam_release_record_never_feeds_ltrdigest(project_root: Path) -> None:
    digester = _rule_block(
        project_root, "rule ltr_digester_setup:", "rule ltr_digester:"
    )
    assert "pfam_release_recorder" not in digester
    assert "Pfam.version" not in digester


def test_pfam_release_record_is_made_with_the_domain_scan(project_root: Path) -> None:
    """The record rides on the aggregate rules only, so adding it reruns nothing."""
    scanner = _rule_block(
        project_root,
        "rule domain_scanner:",
        "# ----------------------------- Taxonomic",
    )
    assert "pfam_release_recorder" in scanner


def test_retired_domain_regex_config_is_gone(project_root: Path) -> None:
    """The retired domain regex settings stay out of the config and schema.

    `config.domains` and `parameters.hit_domain_mode` were retired by ADR-015;
    a reintroduced block would silently resurrect the name-regex labelling.
    """
    config = (project_root / "data" / "config" / "config.yaml").read_text(
        encoding="utf-8"
    )
    schema = (project_root / "data" / "config" / "schema.yaml").read_text(
        encoding="utf-8"
    )
    assert "\ndomains:" not in config
    assert "hit_domain_mode" not in config
    assert "domains_schema" not in schema


def test_curated_class_table_is_committed(project_root: Path) -> None:
    """The class table is the only human artifact in the feature; it must ship."""
    tsv = project_root / "data" / "config" / "pfam_domain_classes.tsv"
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == ["pfam_acc", "pfam_name", "class"]
    classes = {line.split("\t")[2] for line in lines[1:]}
    assert classes <= {
        "retroviral_diagnostic",
        "retroelement_shared",
        "non_ltr",
        "dna_transposon",
        "other",
    }
    # The families the retired regex could not see must now be present.
    names = {line.split("\t")[1] for line in lines[1:]}
    assert {"rve", "RVP", "IN_DBD_C", "GP41"} <= names


def test_ruleorder_normalizer_wins_over_downloader(project_root: Path) -> None:
    """The normalizer must take precedence when both can produce {genome}.fa."""
    text = _read_snakefile(project_root)
    assert "ruleorder: genome_fasta_normalizer_setup > genome_downloader_setup" in text


def test_blast_db_generator_no_longer_renames_inline(project_root: Path) -> None:
    """The inline `find / parallel mv` workaround must be gone - normalizer owns it."""
    text = _read_snakefile(project_root)
    assert "parallel 'mv" not in text
    assert "-iname '*.fna'" not in text


# ---------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------
def test_genome_wildcard_constraint_pinned_to_species_list(
    project_root: Path,
) -> None:
    """``wildcard_constraints.genome`` must be a regex pinned to SPECIES.

    Protects the recent fix (commit f1a844b) that prevents Snakemake's
    default ``.+`` matcher from greedily absorbing suffixes like
    ``_retroviral`` or ``_full`` into ``{genome}``.
    """
    text = _read_snakefile(project_root)
    match = WILDCARD_BLOCK_RE.search(text)
    assert match is not None, (
        "wildcard_constraints block for {genome} missing - "
        "see commit f1a844b for why this matters"
    )
    expr = match.group(1)
    # The expression should reference SPECIES and the no-match fallback.
    assert "SPECIES" in expr, "wildcard regex no longer pinned to SPECIES list"


# ---------------------------------------------------------------------
# Phase-2/3 contract checks (xfail until those phases land)
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Retired erv_like assembly tier - its producer outputs must be GONE; the
# erv-like plot panel now reads the genus-founded taxonomy loci table.
# ---------------------------------------------------------------------
def test_erv_like_producer_tier_retired(project_root: Path) -> None:
    """``ranges_analysis`` must no longer build the probe-label erv_like tier."""
    text = _read_snakefile(project_root)
    assert "TRACK_ERV_LIKE_DIR" not in text
    assert "erv_like_tracks=" not in text
    assert "{genome}.erv_like_loci.parquet" not in text
    assert "{genome}.erv_like_members.parquet" not in text
    assert "--erv_like_ranges" not in text
    # The producer module is deleted too.
    assert not (
        project_root / "workflow" / "scripts" / "range_analysis" / "erv_assembly.R"
    ).exists()


def test_erv_like_plot_generator_reads_genus_loci(project_root: Path) -> None:
    """The structure plot panel survives, repointed at the taxonomy loci table."""
    text = _read_snakefile(project_root)
    rules = _all_rules(project_root)
    assert "erv_like_plot_generator_setup" in rules
    assert "erv_like_plot_generator" in rules
    assert "STRUCTURE_PLOT_DIR" in text  # panel renamed erv-like -> structure
    # Its input is now the genus-founded loci table, not ranges_analysis tables.
    assert "TAXONOMY_TABLES_PARQUET_DIR" in text


def _stage(flag: str) -> stages.Stage:
    return next(s for s in stages.STAGES if s.flag == flag)


def test_generate_global_plots_includes_erv_like_panel() -> None:
    """--generate-global-plots must drive the structure panel alongside ranges."""
    assert "erv_like_plot_generator" in _stage("--generate-global-plots").targets


def test_stage_pdfs_declared(project_root: Path) -> None:
    """Every plotting stage declares its one stage PDF (docs/visual_style.md).

    Declared, not left to params: an undeclared figure is never rebuilt when
    stale and a half-finished run looks complete.
    """
    text = _read_snakefile(project_root)
    for pdf in (
        "'homology.pdf'",
        "'integration.pdf'",
        "'taxonomy.pdf'",
        "'structure.pdf'",
        "'loss.pdf'",
        "'all_species.solo_ltr.pdf'",
        "'{genome}.hotspots.pdf'",
    ):
        assert pdf in text, f"stage PDF {pdf} not declared in the Snakefile"


def test_provirus_tables_declared(project_root: Path) -> None:
    """The overlap / LTR-interaction tables behind the integration pages exist."""
    text = _read_snakefile(project_root)
    for table in (
        "provirus_overlap",
        "ltr_interaction",
        "probe_domain_overlap",
        "reduction_coverage",
    ):
        assert f"{{genome}}.{table}.parquet" in text, f"table {table!r} not declared"


def test_original_candidate_reduced_exports_removed(project_root: Path) -> None:
    """Only the valid tier keeps a reduced export; original/candidate dropped."""
    text = _read_snakefile(project_root)
    assert "original_tracks_reduced" not in text
    assert "candidate_tracks_reduced" not in text
    assert "--original_ranges_reduced" not in text
    assert "--candidate_ranges_reduced" not in text
    # The valid reduced track and its CLI flag must remain.
    assert "element_hits_tracks_reduced" in text
    assert "--element_hits_ranges_reduced" in text


def test_erv_like_config_block_removed(project_root: Path) -> None:
    """The retired erv_like assembly tier leaves no config/schema surface."""
    config = (project_root / "data" / "config" / "config.yaml").read_text()
    schema = (project_root / "data" / "config" / "schema.yaml").read_text()
    assert "erv_like:" not in config
    assert "max_join_distance:" not in config
    assert "erv_like_schema:" not in schema
    # hotspot input enum no longer offers the retired tier.
    assert "erv_like" not in schema


# ---------------------------------------------------------------------
# Taxonomic-classification stage (reference build + per-locus genus calls)
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "rule_name",
    [
        "taxonomy_reference",
        "taxonomy_reference_trees_setup",
        "taxonomy_reference_trees",
        "taxonomy_classify_setup",
        "taxonomy_classify",
        "taxonomy_plot_generator_setup",
        "taxonomy_plot_generator",
    ],
)
def test_taxonomy_rule_present(project_root: Path, rule_name: str) -> None:
    """Every taxonomic-classification rule must be in the workflow."""
    assert rule_name in _all_rules(project_root)


def test_taxonomy_classify_declares_outputs(project_root: Path) -> None:
    """The classify rule emits the genus-founded loci tables + IGV tracks."""
    text = _read_snakefile(project_root)
    assert "{genome}.loci.parquet" in text
    assert "TRACK_TAXONOMY_DIR" in text
    assert "TAXONOMY_TABLES_PARQUET_DIR" in text


def test_taxonomy_reference_seed_from_config(project_root: Path) -> None:
    """Placement trees must be seeded from parameters.seed (reproducibility)."""
    text = _read_snakefile(project_root)
    assert "_TAX_SEED" in text
    assert "config['parameters'].get('seed'" in text


def test_classify_cli_flags_present() -> None:
    """RetroSeek CLI exposes --build-reference and --classify."""
    assert _stage("--build-reference").targets == ("taxonomy_reference_trees",)
    assert "taxonomy_classify" in _stage("--classify").targets


def test_classification_config_and_schema(project_root: Path) -> None:
    """The classification config block + its Yamale sub-schema must be present."""
    config = (project_root / "data" / "config" / "config.yaml").read_text()
    schema = (project_root / "data" / "config" / "schema.yaml").read_text()
    assert "classification:" in config
    assert "placement_genes:" in config
    assert "classification_schema:" in schema


def test_curated_erv_class_committed(project_root: Path) -> None:
    """The one curated reference piece (erv_class.tsv) is tracked under data/config."""
    erv_class = project_root / "data" / "config" / "erv_class.tsv"
    assert erv_class.exists()
    assert "erv_class" in erv_class.read_text()


def test_solo_blaster_threads_are_configurable(project_root: Path) -> None:
    """A 102-genome run needs several blastn jobs side by side (2026-09-23)."""
    text = _read_snakefile(project_root)
    block = text.split("rule solo_blaster_setup:", 1)[1].split("\nrule ", 1)[0]
    assert "_SOLO.get('blast_threads') or workflow.cores" in block


def test_heavy_rules_keep_their_commands(project_root: Path) -> None:
    """Heavy rules (and the blast checkpoint) keep their shell commands unchanged.

    A changed command reruns them on every genome, about a day each for
    LTRdigest (ADR-020/021). The per-job log layout is therefore added to every
    rule except these.
    """
    text = _read_snakefile(project_root)
    for rule in sorted(stages.HEAVY_RULES):
        start = re.search(rf"^(?:rule|checkpoint) {rule}:", text, flags=re.MULTILINE)
        assert start, rule
        following = re.search(
            r"^(?:rule|checkpoint|ruleorder) ", text[start.end() :], flags=re.MULTILINE
        )
        block = text[
            start.start() : start.end()
            + (following.start() if following else len(text))
        ]
        assert "job_log(" not in block, f"{rule} gained a job log"
