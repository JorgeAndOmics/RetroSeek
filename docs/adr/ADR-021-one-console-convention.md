# ADR-021: One console convention, drawn by the launcher

- **Status**: Accepted
- **Date**: 2026-09-23
- **Deciders**: Jorge González García
- **Builds on**: [ADR-020](ADR-020-one-stage-table-one-workflow.md) (one Snakemake run per launch), [ADR-018](ADR-018-one-visual-system.md) (the same idea for figures)

## Context

What the pipeline printed had no owner:

- **Seven Python output styles**: a coloured logger (`colored_logging`, fixed at
  DEBUG, with no level in its files), private `basicConfig` setups in two
  formats, bare `print`, `print("WARNING: ...")` to stdout, raw
  `sys.stderr.write` of tool output, and progress bars.
- **R**: nine copies of the same `log_section` helper, bare `message()` calls,
  warnings written as messages, and `options(warn = -1)` switching all warnings
  off in three scripts.
- **Three log layouts**: flat files overwritten by parallel jobs
  (`full_genome_blaster.txt`), per-rule files that captured only stderr, and none
  at all for most R rules. `validator.log` was announced but never written.
- **Switches instead of levels**: `display_snakemake_info`,
  `display_operation_info` and `display_requests_warning`, the last of which could
  hide a failed NCBI fetch.
- **Silent failures**: helper functions returning `None` where a stop was due, and
  five private copies of a tool runner that dumped a cut-off stderr tail.

## Decision

- **A line contract.** Every script writes `HH:MM:SS LEVEL step genome | message`
  to stderr and to its job log, through `log.py` (Python), `utils/log.R` (R) or a
  small `say` function (Bash). Levels are DEBUG, INFO, OK, WARN and ERROR.
- **The launcher is the only renderer.** It reads the one stream Snakemake and
  every job write into, keeps all of it in `LOG_DIR/runs/<time>.log`, counts
  warnings and failures, and draws the screen with rich: a banner, a progress bar,
  coloured lines, a summary. rich is used in `console.py` and nowhere else.
- **One log layout.** Each non-heavy rule has `log: job_log('<step>', '{genome}')`,
  `LOG_DIR/<step>/<genome>.log`; scripts read their step and genome from that
  path. Heavy rules stay byte-identical, so their scripts build the same path
  themselves.
- **One setting.** `display.verbosity` (`quiet` | `normal` | `verbose`), or
  `--verbosity` for one run, replaces the three switches and the `logging:` colour
  block. Retired keys stop the run with a message naming their replacement.
- **Errors carry their fix.** `PipelineError(message, hint)` (Python) and
  `abort_hint(message, hint)` (R), raised where the cause is known, become one
  ERROR line ending in `Fix: ...`. `run_main` turns any other exception into one
  line on screen and its traceback in the job log; in R it also turns every
  `warning()` into a counted WARN line. Tools run through `external.run_tool`.
- **Guards.** `tests/unit/test_console_style.py` rejects private logging setups,
  `print`, copied helpers and switched-off warnings.

## Consequences

- Positive:
  - One readable screen for a run of a hundred genomes, with the failures, their
    fixes and their logs gathered in the summary.
  - Every job has its own complete log; nothing is overwritten or lost.
  - A failed NCBI fetch or tool can no longer pass silently.
  - One dependency less (`coloredlogs`); `rich` was already installed.
- Negative:
  - Adding `log:` to the non-heavy rules changed their commands, so the next run
    recomputes the downstream stages once.
  - The launcher reads four patterns from Snakemake's text: the job table and
    the reasons (the guard), the progress line and the error block (the screen).
    A Snakemake that changes them stops the bar or empties the failure list; it
    never stops the run, and the unit tests notice.
  - Configs with the old display or logging keys must be edited once.
- Neutral:
  - The heavy rules keep writing their tools' own output as before.

## Alternatives considered

- **A Snakemake logger plugin**: the official way to restyle Snakemake's output,
  but an installable package to maintain for what a few dozen lines do here.
- **rich or cli in every script**: colour decisions in forty places, and escape
  codes in files whenever a script guesses wrong about the terminal.
- **Everything on screen**: with 64 parallel jobs the screen becomes unreadable;
  the logs keep everything instead.

## Revisit trigger

- Snakemake changing its "steps done" or "Error in rule" lines.
- A cluster executor where jobs' stderr no longer reaches the launcher.

## References

- [console_style.md](../console_style.md)
- `workflow/scripts/log.py`, `workflow/scripts/utils/log.R`,
  `workflow/scripts/console.py`, `workflow/scripts/external.py`
