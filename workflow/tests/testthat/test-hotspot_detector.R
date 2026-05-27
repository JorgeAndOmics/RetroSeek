# testthat scaffold for workflow/scripts/hotspot_detector.R
#
# hotspot_detector is an EXPERIMENTAL rule (regioneR permutation testing —
# stochastic and slow). Real unit tests are deferred until the window-tiling
# helper is extracted into a pure function. Skipped (not fake-passing) so the
# coverage gap is honest rather than hidden behind expect_true(TRUE).

suppressMessages({
  library(testthat)
})

test_that("hotspot_detector unit tests are deferred (experimental)", {
  skip("experimental rule — unit tests deferred")
})
