# Shared test helpers. testthat sources every helper-*.R file before the tests.

# The lines an expression writes to stderr. utils/log.R writes its console lines
# there with cat(), which a message sink catches (a bare message() condition
# would be intercepted by testthat instead).
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
