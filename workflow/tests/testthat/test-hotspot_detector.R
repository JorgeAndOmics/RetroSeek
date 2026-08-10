# testthat scaffold for workflow/scripts/hotspot_detector.R
#
# hotspot_detector.R is now an orchestrator only - every pure transform lives in
# workflow/scripts/hotspot/*.R and is unit-tested directly:
#   * test-hotspot_io.R          - config reader + FASTA/GFF loaders
#   * test-hotspot_masking.R     - N-mask, effective_bp, scaffold pooling
#   * test-hotspot_windowing.R   - tiling, counting, window-table assembly
#   * test-hotspot_models.R      - NB GLM fit + scoring (deterministic)
#   * test-hotspot_postprocess.R - select / merge / filter / id-assignment
#
# The end-to-end orchestrator (deterministic NB GLM) is exercised at the
# integration layer by the Snakemake dry-run (`make test-snakemake`) and the
# real-genome sanity check, so there is nothing left to fake-pass here.

suppressMessages({
  library(testthat)
})

test_that("hotspot pure transforms are covered by the hotspot_analysis module tests", {
  skip("orchestrator is integration-level; pure logic covered by test-hotspot_* module tests")
})
