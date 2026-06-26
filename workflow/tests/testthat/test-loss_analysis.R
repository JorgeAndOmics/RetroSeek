# testthat coverage for loss_analysis.R pure aggregation functions.
# Sources only the script (the CLI/plot block is guarded by sys.nframe()), so no
# plot environment is needed.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "loss_analysis.R"))


test_that("build_loss_funnel orders stages and computes step retention", {
  counts <- tribble(
    ~genome, ~metric,                ~value,
    "g1",    "raw_blast_hits",        1000,
    "g1",    "filtered_blast_hits",    600,
    "g1",    "first_reduced_ranges",   500,
    "g1",    "candidate_ranges",       300,
    "g1",    "valid_ranges",            60,
    "g1",    "global_reduced_ranges",  250,
    "g1",    "unanchored_fragments",   180,
    "g1",    "fragments_recovered",     45
  )
  f <- build_loss_funnel(counts)

  # stages come back in funnel order
  expect_equal(f$metric[f$genome == "g1"][1], "raw_blast_hits")
  expect_true(all(diff(f$stage_order[f$genome == "g1"]) > 0))

  # step retention is value / parent value
  filt <- f %>% filter(metric == "filtered_blast_hits")
  expect_equal(filt$step_retained, 0.6)            # 600 / 1000

  # candidate descends from FIRST reduction (anchored spine), not global reduction
  cand <- f %>% filter(metric == "candidate_ranges")
  expect_equal(cand$branch, "main")
  expect_equal(cand$step_retained, 300 / 500)      # vs first_reduced parent -> <= 1
  expect_true(cand$step_retained <= 1)

  # frac_of_input is value / raw hits
  expect_equal((f %>% filter(metric == "valid_ranges"))$frac_of_input, 0.06)

  # global reduction heads the fragments branch (sibling of candidate)
  glob <- f %>% filter(metric == "global_reduced_ranges")
  expect_equal(glob$branch, "fragments")
  expect_equal(glob$step_retained, 250 / 500)      # vs first_reduced parent

  # fragments are anchored to the global-reduced parent
  frag <- f %>% filter(metric == "unanchored_fragments")
  expect_equal(frag$branch, "fragments")
  expect_equal(frag$step_retained, 180 / 250)
  rec <- f %>% filter(metric == "fragments_recovered")
  expect_equal(rec$step_retained, 45 / 180)
})

test_that("build_loss_funnel ignores unknown metrics and tolerates empties", {
  counts <- tribble(
    ~genome, ~metric,          ~value,
    "g1",    "raw_blast_hits",   100,
    "g1",    "some_other_thing",  50
  )
  f <- build_loss_funnel(counts)
  expect_false("some_other_thing" %in% f$metric)
  expect_equal(nrow(build_loss_funnel(counts[0, ])), 0L)
})

test_that("build_loss_funnel guards against zero/absent parents", {
  counts <- tribble(
    ~genome, ~metric,               ~value,
    "g1",    "filtered_blast_hits",   10   # no raw_blast_hits present
  )
  f <- build_loss_funnel(counts)
  expect_true(is.na((f %>% filter(metric == "filtered_blast_hits"))$step_retained))
  expect_true(is.na((f %>% filter(metric == "filtered_blast_hits"))$frac_of_input))
})

test_that("pick_novel_candidates keeps only zero-blastx-hit loci", {
  loci <- tribble(
    ~id,  ~genus_call,        ~n_blastx_hits,
    "L0", "Gammaretrovirus",  "5",
    "L1", "UNCLASSIFIED",     "0",
    "L2", "UNCLASSIFIED",     "0",
    "L3", "Lentivirus",       "2"
  )
  novel <- pick_novel_candidates(loci)
  expect_equal(novel$id, c("L1", "L2"))
})

test_that("pick_novel_candidates tolerates missing column / empty frame", {
  expect_equal(nrow(pick_novel_candidates(tibble())), 0L)
  expect_equal(nrow(pick_novel_candidates(tibble(id = "L0"))), 0L)
})


test_that("novel_burden_table computes per-genome novel count + fraction", {
  funnel <- build_loss_funnel(tribble(
    ~genome, ~metric,              ~value,
    "g1",    "loci_total",          100,
    "g1",    "loci_no_blastx_hit",   25,
    "g2",    "loci_total",           40,
    "g2",    "loci_no_blastx_hit",    0
  ))
  bt <- novel_burden_table(funnel)
  g1 <- bt[bt$genome == "g1", ]
  expect_equal(g1$n_novel, 25)
  expect_equal(g1$n_total, 100)
  expect_equal(g1$frac, 0.25)
  expect_equal(bt$frac[bt$genome == "g2"], 0)
})

test_that("novel_burden_table guards zero/absent totals and empty funnel", {
  expect_equal(nrow(novel_burden_table(build_loss_funnel(tibble()))), 0L)
  # total present but zero -> frac 0, not NaN
  funnel <- build_loss_funnel(tribble(
    ~genome, ~metric,              ~value,
    "g1",    "loci_total",           0,
    "g1",    "loci_no_blastx_hit",   0
  ))
  expect_equal(novel_burden_table(funnel)$frac, 0)
})
