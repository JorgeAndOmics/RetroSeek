# -----------------------------------------------------------------------------
# Module hotspot/chromosomes.R
# -----------------------------------------------------------------------------
# Two checks the per-window negative binomial test cannot make on its own.
#
#   chromosome_rate_test()   does a whole chromosome carry more loci than its
#                            share of the genome? The window model gives each
#                            chromosome its own baseline, so excess spread over
#                            a whole chromosome (an enriched Y, for one) is
#                            absorbed there and never called.
#   randomised_nb_pvalues()  p-values for the Q-Q page. A count is discrete, so
#                            many windows share one p-value (every empty window
#                            has p = 1) and a calibrated model still bends a Q-Q
#                            plot below the line. Spreading each tie evenly over
#                            its step gives p-values that are uniform when the
#                            model is right (Dunn and Smyth 1996).
# Measured on the model 5 in Workbench measurements/2026-10-02/hotspot_calibration.md.
# -----------------------------------------------------------------------------

#' Randomised upper-tail p-values for counts under a negative binomial.
#'
#' P(X > x) + U * P(X = x) with U uniform on (0, 1), from a fixed seed so the
#' page draws the same way every run.
#' @param count integer counts; @param mu fitted means; @param theta NB size
#' @return numeric p-values in (0, 1], uniform when the model fits
randomised_nb_pvalues <- function(count, mu, theta, seed = 1L) {
  u <- withr::with_seed(seed, stats::runif(length(count)))
  above <- stats::pnbinom(count, mu = mu, size = theta, lower.tail = FALSE)
  above + u * stats::dnbinom(count, mu = mu, size = theta)
}

# The median of a chi-square with one degree of freedom, to turn a median
# Pearson contribution into a dispersion estimate.
.CHISQ1_MEDIAN <- stats::qchisq(0.5, df = 1)

#' Per label and chromosome, observed loci against the chromosome's share of
#' the genome's loci, with a one-sided test for excess.
#'
#' Expected = all of the label's loci x the chromosome's share of effective bp.
#' Loci cluster, so chromosome counts vary more than Poisson; the dispersion is
#' the median Pearson contribution over the chromosomes (robust to the one
#' chromosome being tested, which an overall estimate would let hide itself),
#' never below 1. p from z = (O - E) / sqrt(dispersion x E); BH within a label.
#' @param window_df window table: label, chrom, effective_bp, count
#' @return one row per label and chromosome: label, chrom, observed, expected,
#'   rate_ratio, pval, qval. A label with fewer than two chromosomes or no
#'   loci is left out.
chromosome_rate_test <- function(window_df) {
  per_label <- split(window_df, window_df$label)
  rows <- lapply(names(per_label),
                 function(label) .label_chromosome_test(per_label[[label]], label))
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(label = character(), chrom = character(), observed = integer(),
                      expected = numeric(), rate_ratio = numeric(), pval = numeric(),
                      qval = numeric()))
  }
  out
}

.label_chromosome_test <- function(d, label) {
  observed <- tapply(as.integer(d$count), d$chrom, sum)
  bp <- tapply(as.numeric(d$effective_bp), d$chrom, sum)
  keep <- bp > 0
  observed <- observed[keep]
  bp <- bp[keep]
  if (length(observed) < 2L || sum(observed) == 0L) return(NULL)
  expected <- sum(observed) * bp / sum(bp)
  pearson <- (observed - expected)^2 / expected
  dispersion <- max(1, stats::median(pearson) / .CHISQ1_MEDIAN)
  z <- (observed - expected) / sqrt(dispersion * expected)
  pval <- stats::pnorm(z, lower.tail = FALSE)
  data.frame(label = label, chrom = names(observed),
             observed = as.integer(observed), expected = as.numeric(expected),
             rate_ratio = as.numeric(observed / expected), pval = as.numeric(pval),
             qval = stats::p.adjust(pval, method = "BH"), row.names = NULL)
}
