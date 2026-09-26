# testthat coverage for workflow/scripts/utils/log.R
#
# R scripts write the same line as Python (ADR-021):
#   HH:MM:SS LEVEL step genome | message
# The launcher parses these lines from many jobs at once, so the shape, the
# per-line prefix and the verbosity threshold are pinned here.
#
# Run via: make test-r

suppressMessages(library(testthat))

source(file.path("..", "..", "scripts", "utils", "log.R"))

.LINE <- "^[0-9]{2}:[0-9]{2}:[0-9]{2} (DEBUG|INFO|OK|WARN|ERROR) \\S+ \\S+ \\| .*$"

.capture_stderr <- function(expr) {
  sink_file <- tempfile()
  con <- file(sink_file, open = "wt")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    close(con)
  })
  force(expr)
  sink(type = "message")
  close(con)
  on.exit()
  readLines(sink_file)
}

test_that("a line has the contract's shape", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  log_setup("solo_plots", "Mus_musculus")
  out <- .capture_stderr(log_warn("%d panels empty", 3L))
  expect_length(out, 1)
  expect_match(out, .LINE)
  expect_match(out, " WARN solo_plots Mus_musculus \\| 3 panels empty$")
})

test_that("no genome reads 'all'", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  log_setup("hotspot")
  out <- .capture_stderr(log_ok("done"))
  expect_match(out, " OK hotspot all \\| done$")
})

test_that("every physical line carries the prefix", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "verbose")
  log_setup("s", "g")
  out <- .capture_stderr(log_info("one\ntwo"))
  expect_length(out, 2)
  expect_true(all(grepl(.LINE, out)))
})

test_that("a literal percent sign survives without arguments", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "verbose")
  log_setup("s", "g")
  out <- .capture_stderr(log_info("100% done"))
  expect_match(out, "\\| 100% done$")
})

test_that("normal verbosity keeps INFO off the console but in the file", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- file.path(tempfile(), "rule", "g.log")
  log_setup("rule", "g", job_log)
  out <- .capture_stderr({
    log_info("detail")
    log_ok("headline")
  })
  expect_false(any(grepl("detail", out)))
  expect_true(any(grepl("OK rule g \\| headline", out)))
  written <- readLines(job_log)
  expect_true(any(grepl("INFO rule g \\| detail", written)))
  expect_true(any(grepl("OK rule g \\| headline", written)))
})

test_that("the job log is appended to, not replaced", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- tempfile(fileext = ".log")
  writeLines("tool output first", job_log)
  log_setup("rule", "g", job_log)
  .capture_stderr(log_info("after"))
  expect_equal(readLines(job_log)[1], "tool output first")
})

test_that("verbosity thresholds match the Python side", {
  expect_equal(log_console_threshold("quiet"), 25)
  expect_equal(log_console_threshold("normal"), 25)
  expect_equal(log_console_threshold("verbose"), 10)
  expect_equal(log_console_threshold(""), 20)
  expect_error(log_console_threshold("loud"), "loud")
})

test_that("log_section reports the elapsed time since setup", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "verbose")
  log_setup("s", "g")
  out <- .capture_stderr(log_section("Phase 1: reading"))
  expect_match(out, "INFO s g \\| Phase 1: reading \\([0-9.]+ s\\)$")
})

test_that("log_job reads step and genome from the log path", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- file.path(tempfile(), "solo_plot_summary", "all.log")
  out <- .capture_stderr({
    log_job(job_log, "fallback")
    log_warn("w")
  })
  expect_true(any(grepl("WARN solo_plot_summary all \\| w$", out)))
  expect_true(any(grepl("started:", readLines(job_log))))
})

test_that("run_main turns abort_hint into one line with its fix and status 1", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- file.path(tempfile(), "hotspot_detector", "Mus_musculus.log")
  status <- NULL
  out <- .capture_stderr({
    log_job(job_log, "x")
    status <- run_main(
      function() abort_hint("no windows passed", "lower hotspot.min_hits"),
      quit_on_error = FALSE
    )
  })
  expect_equal(status, 1L)
  errors <- grep(" ERROR ", out, value = TRUE)
  expect_length(errors, 1)
  expect_match(errors, "\\| no windows passed\\. Fix: lower hotspot.min_hits$")
  expect_true(any(grepl("finished: failed", readLines(job_log))))
})

test_that("run_main routes warnings through log_warn and keeps going", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- file.path(tempfile(), "s", "g.log")
  status <- NULL
  out <- .capture_stderr({
    log_job(job_log, "x")
    status <- run_main(function() {
      warning("3 panels empty")
      log_ok("done anyway")
    }, quit_on_error = FALSE)
  })
  expect_equal(status, 0L)
  expect_true(any(grepl("WARN s g \\| 3 panels empty$", out)))
  expect_true(any(grepl("OK s g \\| done anyway$", out)))
  expect_true(any(grepl("finished: done", readLines(job_log))))
})

test_that("an unexpected error is one screen line, the call stack in the log", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  job_log <- file.path(tempfile(), "s", "g.log")
  status <- NULL
  out <- .capture_stderr({
    log_job(job_log, "x")
    status <- run_main(function() stop("object 'x' not found"), quit_on_error = FALSE)
  })
  expect_equal(status, 1L)
  errors <- grep(" ERROR ", out, value = TRUE)
  expect_length(errors, 1)
  expect_match(errors, "object 'x' not found")
  expect_true(any(grepl("call stack", readLines(job_log))))
})

test_that("a multi-line warning becomes one line, so it counts once", {
  withr::local_envvar(RETROSEEK_VERBOSITY = "normal")
  log_setup("s", "g")
  out <- .capture_stderr(run_main(function() {
    warning("Each of the 2 combined objects has sequence levels not in the other:
  - in 'x': chr9
  Make sure to always combine objects based on the same reference")
  }, quit_on_error = FALSE))
  warns <- grep(" WARN ", out, value = TRUE)
  expect_length(warns, 1)
  expect_match(warns, "not in the other: - in 'x': chr9 Make sure")
})
