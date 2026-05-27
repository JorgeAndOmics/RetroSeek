# testthat scaffold for workflow/scripts/circle_plot_generator.R
#
# circle_plot_generator is an EXPERIMENTAL rule (ggbio Circos-style rendering —
# visual output) — see ROADMAP.md. Real unit tests are deferred until the
# karyotype builder is extracted into a pure function. Skipped (not fake-passing)
# so the coverage gap is honest rather than hidden behind expect_true(TRUE).

suppressMessages({
  library(testthat)
})

test_that("circle_plot_generator unit tests are deferred (experimental)", {
  skip("experimental rule — unit tests deferred; see ROADMAP.md")
})
