# =============================================================================
# log.R
# =============================================================================
# The console line contract for every R script (ADR-021), the twin of log.py.
#
# Each message is one plain line on stderr, and in the job's log file:
#
#   14:02:40 WARN solo_plots Mus_musculus | 3 panels empty
#   time     level step      genome         message
#
# Scripts never colour anything: the launcher parses these lines from every job
# and draws them (console.py). Levels are DEBUG, INFO, OK (the one headline a job
# prints when it finishes), WARN and ERROR. RETROSEEK_VERBOSITY, set by the
# launcher from `display.verbosity`, decides which reach stderr:
#   quiet / normal -> OK and up; verbose -> everything; unset -> INFO and up.
#
# Usage: a script sources this file once near its top, calls log_setup() with
# its step name, genome and log path, then writes through log_section(),
# log_info(), log_warn() and log_ok() as it goes.
#
# Messages use sprintf() when arguments follow, so "100% done" alone is safe.
# Each line is written with one cat() call so parallel jobs never tear it.
# =============================================================================

.LOG_LEVELS <- c(DEBUG = 10, INFO = 20, OK = 25, WARN = 30, ERROR = 40)

.log_state <- new.env(parent = emptyenv())
.log_state$step <- "R"
.log_state$genome <- "all"
.log_state$file <- NULL
.log_state$console <- 20
.log_state$file_level <- 20
.log_state$t0 <- Sys.time()

#' The lowest level written to stderr at a verbosity ("" means run by hand).
log_console_threshold <- function(verbosity) {
  if (is.null(verbosity) || is.na(verbosity) || verbosity == "") {
    return(20)
  }
  switch(verbosity,
    quiet = 25,
    normal = 25,
    verbose = 10,
    stop(sprintf("verbosity '%s' is not one of quiet, normal, verbose", verbosity))
  )
}

#' Set the step, genome and job log for every later message.
#'
#' @param step the rule or script name shown on each line (the rule name
#'   without `_setup`, so the launcher can match a failed job to its message).
#' @param genome the genome the job works on; NULL means all genomes.
#' @param log_file the job log (Snakemake's `log:`); appended to, folder created.
log_setup <- function(step, genome = NULL, log_file = NULL) {
  verbosity <- Sys.getenv("RETROSEEK_VERBOSITY", unset = "")
  .log_state$step <- step
  .log_state$genome <- if (is.null(genome) || !nzchar(genome)) "all" else genome
  .log_state$console <- log_console_threshold(verbosity)
  .log_state$file_level <- if (verbosity == "verbose") 10 else 20
  .log_state$t0 <- Sys.time()
  .log_state$file <- NULL
  if (!is.null(log_file) && nzchar(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    .log_state$file <- log_file
  }
  invisible(NULL)
}

.log_emit <- function(level, fmt, ...) {
  message_text <- if (...length() > 0) sprintf(fmt, ...) else fmt
  pieces <- strsplit(message_text, "\n", fixed = TRUE)[[1]]
  if (length(pieces) == 0) pieces <- ""
  prefix <- sprintf(
    "%s %s %s %s | ",
    format(Sys.time(), "%H:%M:%S"), level, .log_state$step, .log_state$genome
  )
  text <- paste0(paste0(prefix, pieces, collapse = "\n"), "\n")
  rank <- .LOG_LEVELS[[level]]
  if (rank >= .log_state$console) cat(text, file = stderr())
  if (!is.null(.log_state$file) && rank >= .log_state$file_level) {
    cat(text, file = .log_state$file, append = TRUE)
  }
  invisible(NULL)
}

log_debug <- function(fmt, ...) .log_emit("DEBUG", fmt, ...)
log_info <- function(fmt, ...) .log_emit("INFO", fmt, ...)
log_ok <- function(fmt, ...) .log_emit("OK", fmt, ...)
log_warn <- function(fmt, ...) .log_emit("WARN", fmt, ...)
log_error <- function(fmt, ...) .log_emit("ERROR", fmt, ...)

#' An INFO line marking a phase, with the seconds since log_setup().
log_section <- function(name) {
  seconds <- as.numeric(difftime(Sys.time(), .log_state$t0, units = "secs"))
  .log_emit("INFO", "%s (%.1f s)", name, seconds)
}

#' Set up logging for a pipeline job from its Snakemake log path.
#'
#' Log paths are LOG_DIR/<step>/<genome>.log (`all.log` for jobs over every
#' genome), so the path names the job even when one script serves two rules.
#' Without a path (run by hand) lines go to stderr under `step`.
log_job <- function(log_file, step) {
  if (is.null(log_file) || !nzchar(log_file)) {
    log_setup(step)
  } else {
    genome <- tools::file_path_sans_ext(basename(log_file))
    log_setup(basename(dirname(log_file)), if (genome == "all") NULL else genome,
              log_file)
  }
  log_info("started: %s", paste(commandArgs(), collapse = " "))
}

#' Stop with a failure the script understands: what went wrong, and the fix.
abort_hint <- function(message, hint = "") {
  stop(structure(
    class = c("pipeline_error", "error", "condition"),
    list(message = message, call = NULL, hint = hint)
  ))
}

.pipeline_message <- function(cond) {
  text <- sub("\\.$", "", conditionMessage(cond))
  if (nzchar(cond$hint)) sprintf("%s. Fix: %s", text, cond$hint) else paste0(text, ".")
}

#' Run a script's main function so every ending is recorded the same way.
#'
#' Warnings become WARN lines (counted by the launcher) and the script goes on.
#' An abort_hint() becomes one ERROR line with its fix; any other error one ERROR
#' line, with the call stack in the job log. The job log ends with "finished:
#' done" or "finished: failed". On failure R quits with status 1, unless
#' `quit_on_error = FALSE` (tests), which returns the status instead.
run_main <- function(fn, quit_on_error = TRUE) {
  started <- Sys.time()
  stack <- NULL
  status <- tryCatch(
    withCallingHandlers(
      {
        fn()
        0L
      },
      warning = function(w) {
        # One warning, one line: package warnings often wrap over several lines,
        # which would otherwise count as several warnings in the summary.
        log_warn("%s", gsub("\\s*\n\\s*", " ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        # Captured here, while the failing calls are still on the stack.
        if (!inherits(e, "pipeline_error")) stack <<- sys.calls()
      }
    ),
    pipeline_error = function(e) {
      log_error("%s", .pipeline_message(e))
      1L
    },
    error = function(e) {
      calls <- vapply(
        stack, function(call) paste(deparse(call, nlines = 1L), collapse = ""), ""
      )
      log_info("call stack of the error below:\n%s", paste(calls, collapse = "\n"))
      log_error("unexpected error: %s (call stack in the job log)", conditionMessage(e))
      1L
    }
  )
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  log_info("finished: %s in %.1f s", if (status == 0L) "done" else "failed", seconds)
  if (status != 0L && quit_on_error) quit(save = "no", status = status)
  invisible(status)
}
