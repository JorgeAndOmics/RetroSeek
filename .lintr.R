# lintr configuration for the R sources (read by `make lint-r`).
#
# lintr's defaults, with three changes:
# - line length 88, the same limit ruff uses for Python (pyproject.toml);
# - cyclomatic complexity at most 8 per function, the house limit that qlty
#   applies to Python (.qlty/qlty.toml);
# - object_usage_linter off: the scripts source() their helpers and name data
#   frame columns bare (dplyr, ggplot2), so it reports thousands of "no visible
#   binding" and "no visible global function" lines that are not defects.
linters <- linters_with_defaults(
  line_length_linter = line_length_linter(88L),
  cyclocomp_linter = cyclocomp_linter(complexity_limit = 8L),
  object_usage_linter = NULL
)
