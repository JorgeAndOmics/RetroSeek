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


# The words of a page's subtitle. A species-led subtitle is a plotmath call
# (italic species, then the sentence: add_titles()); its string pieces, joined,
# read as the plain subtitle would. A plain subtitle comes back unchanged.
.subtitle_text <- function(p) {
  s <- p$labels$subtitle
  if (is.character(s)) return(s)
  code <- paste(deparse(s, width.cutoff = 500L), collapse = "")
  pieces <- regmatches(code, gregexpr('"[^"]*"', code))[[1]]
  paste(gsub('^"|"$', "", pieces), collapse = "")
}
