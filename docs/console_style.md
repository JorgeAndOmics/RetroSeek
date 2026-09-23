# Console style

What RetroSeek writes to the terminal and to its logs follows one convention, as
its figures follow one visual system ([visual_style.md](visual_style.md)). This
page explains the convention and the reasons for it. Its executable form is
`workflow/scripts/log.py` and `workflow/scripts/utils/log.R` (what scripts write),
`workflow/scripts/external.py` (running tools) and `workflow/scripts/console.py`
(what the launcher draws). The decision is recorded in ADR-021.

## Four rules

1. **Scripts write plain lines; only the launcher draws.** A script never colours
   anything. It writes one line per message, and the launcher, which reads every
   job's output, adds the colour, the symbols and the progress bar.
2. **One line, one shape.** Every message from Python, R or Bash is

   ```
   14:02:40 WARN solo_finder Mus_musculus | 3 baits under 300 bp
   time     level step        genome         message
   ```

   The step is the rule name without `_setup`; the genome is `all` for a job over
   every genome. With dozens of jobs running at once, the step and genome are
   what make interleaved lines readable.
3. **A level means one thing.**
   - `ERROR`: the job cannot produce a correct output. It stops, non-zero, and
     says what to do.
   - `WARN`: the output exists but a person should look. It says what, where, what
     it means and what to do; if it cannot, it is not a warning. A clean run of
     the model genomes shows no warnings, so a warning always deserves a look.
   - `OK`: the one headline a job prints when it finishes, with its key number.
   - `INFO` and `DEBUG`: detail, kept in the job log.
4. **Every line is kept.** The screen shows what the verbosity asks for; the logs
   always hold everything.

## Where things go

| What | Where |
|---|---|
| A job's messages and its tools' error output | `LOG_DIR/<step>/<genome>.log`, one file per job, never shared |
| Everything the run printed, ours and Snakemake's | `LOG_DIR/runs/<time>.log`, written only by the launcher |
| Every warning of the run | `LOG_DIR/runs/<time>.warnings.txt` |
| The screen | the banner, the checks, the lines the verbosity shows, the progress bar, the summary |

A job log opens with a `started:` line (the exact command) and closes with
`finished: done` or `finished: failed` and the time taken.

## Verbosity

`display.verbosity` in the config, or `--verbosity` for one run:

| Verbosity | The screen shows |
|---|---|
| `quiet` | banner, progress bar, errors, summary |
| `normal` | also one `OK` line per finished job, and every warning |
| `verbose` | also every message of every job, and Snakemake's own output |

## Colour

Colour carries meaning, never decoration, and never alone: each level also has a
symbol (or word, where the terminal cannot show symbols).

| Level | Symbol | Colour |
|---|---|---|
| OK | a tick | green |
| WARN | a warning sign | yellow |
| ERROR | a cross | bold red |
| INFO, DEBUG | a dot | plain, dim |

The colours are the terminal's own, so they follow its light or dark theme. There
is no colour when the output is not a terminal, and the standard `NO_COLOR=1` and
`FORCE_COLOR=1` switch it off or on. Files never contain colour codes, so there is
no need to pipe a run through `tee`: the run log already holds it.

## Many jobs, one screen

Snakemake runs many jobs at once, all writing to one stream. Linux delivers a
write of up to 4 KB to a pipe in one piece, and each message is written in one
call, so lines interleave but never tear. The launcher is the only reader and the
only process that draws, so nothing competes for the screen, and each job writes
only its own log file.

## Words

- Say what happened in the reader's terms: "3 baits under 300 bp", not a variable
  name or a Python repr.
- Numbers with thousands separators: `26,499`.
- A warning or error says what to do. A known error is raised with its fix,
  `PipelineError("what went wrong", hint="what to do")` in Python or
  `abort_hint("what went wrong", "what to do")` in R, and shows as one line ending
  in `Fix: ...`.
- No dash as punctuation, no arrows; ASCII in the source.

## Adding a script

1. Python: `job_logging(args.log, "<step>")` after parsing the arguments (with a
   `--log` option), then `run_main(main)`. R: `log_job(args$log, "<step>")`, then
   `run_main(main)` inside the script's if-main guard.
2. Give its rule `log: job_log('<step>', '{genome}')` and pass `--log {log:q}`.
3. Log with the standard loggers (Python) or `log_info` / `log_warn` / `log_ok`
   (R). End with one `OK` headline.
4. Run tools through `external.run_tool(cmd)`: a failure then stops the job with
   one line and keeps the tool's error output in the job log.
5. `make check` runs the guards: `tests/unit/test_console_style.py` (no private
   logging setup, no `print`, no copied helpers, no warnings switched off),
   `tests/unit/test_log.py` and `workflow/tests/testthat/test-log.R` (the line
   contract), `tests/unit/test_console.py` (the launcher's parsing and summary).
