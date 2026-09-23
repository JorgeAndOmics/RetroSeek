"""Unit tests for workflow/scripts/stages.py, the launcher's stage table.

The table drives the parser, the run order, the heavy-rule guard and the tool
check. These tests pin what would silently break if a row were edited: every
target must exist in the Snakefile, phases must follow the order the pipeline
runs in, the preset must hold exactly the downstream stages, and unknown or
abbreviated flags must reach Snakemake untouched.
"""

import re
from pathlib import Path

import pytest

import stages


def _snakefile_rules(project_root: Path) -> set[str]:
    text = (project_root / "workflow" / "Snakefile").read_text(encoding="utf-8")
    return set(re.findall(r"^(?:rule|checkpoint) (\w+):", text, flags=re.MULTILINE))


def test_flags_are_unique() -> None:
    flags = [s.flag for s in stages.STAGES]
    assert len(flags) == len(set(flags))


def test_stages_are_listed_in_pipeline_order() -> None:
    """The table is read top to bottom as the order the phases run in."""
    order = list(stages.PHASES)
    positions = [order.index(s.phase) for s in stages.STAGES]
    assert positions == sorted(positions)


def test_every_target_is_a_snakefile_rule(project_root: Path) -> None:
    rules = _snakefile_rules(project_root)
    for stage in stages.STAGES:
        for target in stage.targets:
            assert target in rules, f"{stage.flag}: no rule {target}"


def test_heavy_rules_exist_and_match_the_protected_list(project_root: Path) -> None:
    """The heavy steps named in the run protocol, each owned by exactly one stage."""
    rules = _snakefile_rules(project_root)
    assert rules >= stages.HEAVY_RULES
    assert {
        "genome_downloader_setup",
        "genome_downloader",
        "pfam_hmm_downloader",
        "probe_extractor",
        "ltr_index_generator_setup",
        "ltr_index_generator",
        "ltr_harvester_setup",
        "ltr_harvester",
        "ltr_digester_setup",
        "ltr_digester",
        "full_genome_blaster_setup",
        "full_genome_blaster",
    } == stages.HEAVY_RULES
    owners = [rule for s in stages.STAGES for rule in s.heavy]
    assert len(owners) == len(set(owners))


def test_domain_scan_is_its_own_stage() -> None:
    """The element and orphan domain scan runs after ranges, apart from LTRdigest."""
    scan = next(s for s in stages.STAGES if s.flag == "--domain-scan")
    assert scan.targets == ("domain_scanner",)
    assert scan.phase == "Analysis"
    assert not scan.heavy


def test_help_lists_the_phases_in_order() -> None:
    text = stages.build_parser().format_help()
    headings = [*stages.PHASES, "Presets", "Run options"]
    positions = [text.index(f"\n{h}:") for h in headings]
    assert positions == sorted(positions)


def test_abbreviations_go_to_snakemake_not_to_a_stage() -> None:
    args, unknown = stages.build_parser().parse_known_args(["--class"])
    assert not args.classify
    assert unknown == ["--class"]


def test_snakemake_options_pass_through() -> None:
    argv = ["--classify", "-n", "--configfile", "/x.yaml", "--forcerun", "r"]
    args, unknown = stages.build_parser().parse_known_args(argv)
    assert args.classify
    assert unknown == ["-n", "--configfile", "/x.yaml", "--forcerun", "r"]


def _select(*argv: str) -> list[str]:
    args, _ = stages.build_parser().parse_known_args(list(argv))
    return [s.flag for s in stages.selected(args)]


def test_selection_follows_table_order_not_command_line_order() -> None:
    assert _select("--segment", "--classify", "--ranges-analysis") == [
        "--ranges-analysis",
        "--classify",
        "--segment",
    ]


def test_downstream_is_every_analysis_and_figures_stage_but_circle_plots() -> None:
    chosen = _select("--downstream")
    assert "--ranges-analysis" in chosen
    assert "--domain-scan" in chosen
    assert "--solo-ltr-detector" in chosen
    assert "--generate-global-plots" in chosen
    assert "--generate-circle-plots" not in chosen
    assert not {"--ltr-domains", "--blast", "--download-hmm"} & set(chosen)


def test_nothing_selected_without_stage_flags() -> None:
    assert _select("-skp", "-n") == []


def test_targets_and_tools_drop_repeats() -> None:
    args, _ = stages.build_parser().parse_known_args(["--classify", "--domain-scan"])
    chosen = stages.selected(args)
    targets = stages.targets(chosen)
    assert targets[0] == "domain_scanner"
    assert len(targets) == len(set(targets))
    tools = stages.tools(chosen)
    assert tools.count("hmmsearch") == 1


def test_allowed_heavy_comes_only_from_requested_stages() -> None:
    args, _ = stages.build_parser().parse_known_args(["--ltr-domains", "--classify"])
    assert stages.allowed_heavy(stages.selected(args)) == {
        "ltr_digester_setup",
        "ltr_digester",
    }


@pytest.mark.parametrize("flag", ["--config-help"])
def test_config_help_is_absent_unless_asked(flag: str) -> None:
    args, _ = stages.build_parser().parse_known_args(["--classify"])
    assert not hasattr(args, flag.lstrip("-").replace("-", "_"))
