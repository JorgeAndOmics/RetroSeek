# testthat coverage for workflow/scripts/hotspot/chromosomes.R
#
# Two checks the window test cannot make: whether a whole chromosome carries
# more than its share of loci (each chromosome has its own baseline in the
# window model, so chromosome-wide excess is invisible there), and randomised
# p-values for the Q-Q page (a discrete count gives tied p-values, which bend a
# Q-Q plot below the line even for a calibrated model).

suppressMessages(library(testthat))

source(file.path("..", "..", "scripts", "hotspot", "chromosomes.R"))

test_that("randomised p-values from a correct model are uniform", {
  withr::local_seed(11)
  mu <- rep(c(0.2, 1, 3), each = 2000)
  counts <- stats::rnbinom(length(mu), mu = mu, size = 1.5)
  p <- randomised_nb_pvalues(counts, mu, theta = 1.5, seed = 1L)
  expect_true(all(p > 0 & p <= 1))
  expect_gt(stats::ks.test(p, "punif")$p.value, 0.01)
})

test_that("randomised p-values are reproducible for a seed", {
  p1 <- randomised_nb_pvalues(c(0, 1, 5), c(1, 1, 1), theta = 2, seed = 3L)
  p2 <- randomised_nb_pvalues(c(0, 1, 5), c(1, 1, 1), theta = 2, seed = 3L)
  expect_identical(p1, p2)
})

.windows <- function() {
  # Twelve chromosomes of equal size; chrD holds four times its share.
  chroms <- c("chrA", "chrB", "chrC", "chrD", sprintf("chr%02d", 1:8))
  d <- data.frame(chrom = rep(chroms, each = 10), effective_bp = 5e5, count = 1L,
                  label = "Betaretrovirus")
  d$count[d$chrom == "chrD"] <- 4L
  d
}

test_that("the chromosome test flags the chromosome with excess loci", {
  out <- chromosome_rate_test(.windows())
  d <- out[out$chrom == "chrD", ]
  expect_equal(d$observed, 40L)
  expect_equal(d$expected, 12.5)
  expect_gt(d$rate_ratio, 2)
  expect_lt(d$qval, 0.05)
  expect_true(all(out$qval[out$chrom != "chrD"] > 0.05))
})

test_that("a label with one chromosome or no loci is left untested", {
  one <- .windows()[.windows()$chrom == "chrA", ]
  expect_equal(nrow(chromosome_rate_test(one)), 0L)
  empty <- transform(.windows(), count = 0L)
  expect_equal(nrow(chromosome_rate_test(empty)), 0L)
})

test_that("one locus on a chromosome expecting a tiny share is not called", {
  # 1 observed against 0.06 expected: the exact chance of one or more is ~0.06,
  # which a normal approximation would turn into p ~ 1e-4.
  d <- .windows()
  d$count <- 0L
  d$count[seq(1, nrow(d), by = 2)] <- 1L  # 5 loci on every chromosome
  short <- d$chrom == "chr02"
  d$effective_bp[short] <- 5e3            # a hundredth of the others' size
  d$count[short] <- c(1L, rep(0L, 9))
  out <- chromosome_rate_test(d)
  expect_lt(out$expected[out$chrom == "chr02"], 0.1)
  expect_gt(out$pval[out$chrom == "chr02"], 0.01)
})
