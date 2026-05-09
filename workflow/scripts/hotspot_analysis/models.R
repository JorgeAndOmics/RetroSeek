# -----------------------------------------------------------------------------
# hotspot_analysis / models.R
# -----------------------------------------------------------------------------
# Per-label enrichment models. Two complementary paths:
#
#   * `fit_nb_model()` + `score_windows_nb()` — Negative-Binomial GLM on
#     per-window counts, mask-aware via offset, optionally chromosome-
#     stratified. Convergence ladder: NB (`MASS::glm.nb`) -> Poisson
#     (`stats::glm`, family = poisson) -> skip with explicit status. Catches
#     the `theta.ml: iteration limit reached` warning that `glm.nb` emits
#     while still returning a "model" — that warning is treated as a fit
#     failure because the resulting theta is unreliable.
#
#   * `score_windows_perm()` — opt-in permutation pass via
#     `regioneR::randomizeRegions`, seeded from the global
#     `parameters.seed`. Streams a per-window >= counter instead of
#     materialising the full null matrix, so memory stays O(W).
#
# Pure module: callers pass tibbles / GRanges in, get tibbles out. No I/O.

suppressMessages({
  library(MASS)
  library(stats)
  library(tibble)
  library(dplyr)
  library(GenomicRanges)
  library(regioneR)
})


# Threshold under which a window's effective_bp is too small to trust as a
# rate denominator. 10% of the configured window_size, per the Plan agent's
# recommendation: smaller-than-this windows let one stray hit produce a huge
# rate that biases theta.
.MIN_EFFECTIVE_FRACTION <- 0.1


# Internal: run `MASS::glm.nb` with explicit handling of the silent
# "iteration limit reached" warning that the inner `theta.ml` estimator
# emits while still returning a model object whose theta is unreliable.
.fit_glm_nb_safely <- function(formula, data) {
  caught_theta_warning <- FALSE
  model <- withCallingHandlers(
    tryCatch(
      MASS::glm.nb(formula, data = data),
      error = function(e) NULL
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("iteration limit reached", msg, ignore.case = TRUE) ||
          grepl("alternation limit reached", msg, ignore.case = TRUE)) {
        caught_theta_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  list(model = model, theta_warning = caught_theta_warning)
}


#' Fit a Negative-Binomial GLM on per-window counts for one label group.
#'
#' Pre-filters windows whose `effective_bp / window_size < 0.1` (rate-
#' denominator too small) and refuses to fit if fewer than `min_nonzero`
#' non-zero windows remain — theta is unidentifiable below that.
#'
#' Convergence ladder:
#'   1. `MASS::glm.nb` succeeds and produces a stable theta -> status = "ok".
#'   2. `glm.nb` fails or theta-MLE warning fires -> Poisson fallback;
#'      status = "poisson_fallback". The dispersion bias makes Poisson
#'      anti-conservative, so the orchestrator logs a warning.
#'   3. Poisson also fails -> status = "failed".
#'   4. Pre-filter leaves < min_nonzero non-zero windows -> status =
#'      "insufficient_data".
#'
#' @return list with:
#'   $model           — fitted glm/glm.nb object or NULL
#'   $family          — "nb" | "poisson" | NA
#'   $theta           — NB dispersion or NA
#'   $status          — one of: ok | poisson_fallback | failed | insufficient_data
#'   $fit_data        — the (filtered) tibble actually fit, for diagnostics
#'   $strata          — TRUE if chromosome was a covariate
fit_nb_model <- function(window_df,
                         window_size,
                         strata_by_chromosome = TRUE,
                         min_nonzero = 50L) {
  fit_data <- dplyr::filter(
    window_df,
    .data$effective_bp > 0L,
    .data$effective_bp / as.integer(window_size) >= .MIN_EFFECTIVE_FRACTION
  )
  n_nonzero <- sum(fit_data$count > 0L)
  if (n_nonzero < as.integer(min_nonzero)) {
    return(list(
      model = NULL, family = NA_character_, theta = NA_real_,
      status = "insufficient_data", fit_data = fit_data,
      strata = isTRUE(strata_by_chromosome)
    ))
  }

  # Drop chrom_stratum levels with no variability (a stratum that exists
  # in fit_data but has only zeros provides no information and risks
  # rank-deficient fits).
  formula_str <- if (isTRUE(strata_by_chromosome) &&
                     length(unique(fit_data$chrom_stratum)) > 1L) {
    "count ~ chrom_stratum + offset(log(effective_bp))"
  } else {
    "count ~ 1 + offset(log(effective_bp))"
  }
  formula <- stats::as.formula(formula_str)

  nb_attempt <- .fit_glm_nb_safely(formula, fit_data)
  if (!is.null(nb_attempt$model) && !nb_attempt$theta_warning) {
    return(list(
      model = nb_attempt$model, family = "nb",
      theta = nb_attempt$model$theta,
      status = "ok", fit_data = fit_data,
      strata = grepl("chrom_stratum", formula_str, fixed = TRUE)
    ))
  }

  pois_model <- tryCatch(
    suppressWarnings(stats::glm(formula, data = fit_data, family = stats::poisson())),
    error = function(e) NULL
  )
  if (!is.null(pois_model)) {
    return(list(
      model = pois_model, family = "poisson", theta = NA_real_,
      status = "poisson_fallback", fit_data = fit_data,
      strata = grepl("chrom_stratum", formula_str, fixed = TRUE)
    ))
  }

  list(
    model = NULL, family = NA_character_, theta = NA_real_,
    status = "failed", fit_data = fit_data,
    strata = grepl("chrom_stratum", formula_str, fixed = TRUE)
  )
}


#' Predict a per-window mean given a fit and a tibble of windows.
#'
#' `predict.glm` with `type = "response"` returns mu (counts per the supplied
#' offset). When `newdata$chrom_stratum` contains levels unseen at fit time
#' (windows we filtered out before fitting), prediction errors. We coalesce
#' those to the most-populous fitted level so prediction yields a sensible
#' baseline rate rather than failing — flagged as a small caveat in the
#' manifest. When `model = NULL` (insufficient_data / failed), returns NA mu.
.predict_mu <- function(window_df, fit) {
  if (is.null(fit$model)) {
    return(rep(NA_real_, nrow(window_df)))
  }
  newdata <- window_df
  if (isTRUE(fit$strata)) {
    fitted_levels <- unique(fit$fit_data$chrom_stratum)
    fallback <- names(sort(table(fit$fit_data$chrom_stratum), decreasing = TRUE))[1]
    miss <- !(newdata$chrom_stratum %in% fitted_levels)
    newdata$chrom_stratum[miss] <- fallback
  }
  # log(0) is undefined; clamp effective_bp >= 1 in the prediction frame so
  # the offset is finite. Such rows already have NA p-value rules applied
  # downstream (count must be > 0 for a hotspot anyway).
  newdata$effective_bp <- pmax(newdata$effective_bp, 1L)
  as.numeric(stats::predict(fit$model, newdata = newdata, type = "response"))
}


#' Score per-window enrichment given a fit. Adds mu_hat, pval_nb, qval_nb.
#'
#' p-value uses the upper tail of the fitted family at observed count:
#'   NB:      pnbinom(count - 1, mu = mu_hat, size = theta, lower.tail = FALSE)
#'   Poisson: ppois (count - 1, lambda = mu_hat,             lower.tail = FALSE)
#'
#' BH adjustment is computed within the per-label call (one label = one
#' family of tests). Windows with effective_bp == 0 receive NA p-values
#' (no rate denominator, no signal).
#'
#' Note: per-window p-values from a model fit on those same windows is
#' in-sample. Extreme-count windows pull the fit toward themselves, mildly
#' deflating their own p-values. Acknowledged caveat — see plan.
score_windows_nb <- function(window_df, fit) {
  mu <- .predict_mu(window_df, fit)
  count <- as.integer(window_df$count)
  pval <- rep(NA_real_, nrow(window_df))
  callable <- !is.na(mu) & window_df$effective_bp > 0L
  if (any(callable)) {
    if (identical(fit$family, "nb") && !is.na(fit$theta)) {
      pval[callable] <- stats::pnbinom(
        q = count[callable] - 1L,
        mu = mu[callable],
        size = fit$theta,
        lower.tail = FALSE
      )
    } else if (identical(fit$family, "poisson")) {
      pval[callable] <- stats::ppois(
        q = count[callable] - 1L,
        lambda = mu[callable],
        lower.tail = FALSE
      )
    }
  }
  qval <- rep(NA_real_, nrow(window_df))
  if (any(!is.na(pval))) {
    qval[!is.na(pval)] <- stats::p.adjust(pval[!is.na(pval)], method = "BH")
  }
  dplyr::mutate(
    window_df,
    mu_nb   = mu,
    pval_nb = pval,
    qval_nb = qval
  )
}


#' Permutation-based scoring (opt-in via `hotspot_validate_permutation`).
#'
#' Streams a per-window >=-count tally across `n_perm` randomizations rather
#' than materialising the full null matrix; memory stays O(W). The seed is
#' set immediately before the loop so each call is reproducible given the
#' global `parameters.seed`.
#'
#' Adds `pval_perm` and `qval_perm` columns.
#'
#' @param window_df    tibble produced by `assemble_window_table()`
#' @param windows      GRanges of windows aligned with `window_df` rows
#' @param hits         GRanges of label-specific hits to randomize
#' @param genome_gr    GRanges describing the genome (per-chromosome bounds)
#' @param mask         GRanges to exclude from randomization starts
#' @param n_perm       integer >= 1
#' @param seed         integer; consumed by `set.seed`
score_windows_perm <- function(window_df, windows, hits, genome_gr,
                               mask, n_perm, seed) {
  observed <- as.integer(window_df$count)
  n_windows <- length(windows)
  geq_counter <- integer(n_windows)

  randomize_ervs <- function(A) {
    regioneR::randomizeRegions(
      A,
      genome         = genome_gr,
      mask           = mask,
      per.chromosome = TRUE
    )
  }

  set.seed(as.integer(seed))
  n_perm <- as.integer(n_perm)
  perm_totals <- integer(n_perm)
  for (i in seq_len(n_perm)) {
    rand <- randomize_ervs(hits)
    null_counts <- as.integer(GenomicRanges::countOverlaps(windows, rand))
    geq_counter <- geq_counter + as.integer(null_counts >= observed)
    perm_totals[i] <- sum(null_counts)
  }

  # Add-one smoothing so no p-value is exactly 0 (which would log-transform
  # to +Inf in Q-Q and Manhattan plots).
  pval <- (geq_counter + 1L) / (n_perm + 1L)
  qval <- stats::p.adjust(pval, method = "BH")

  result <- dplyr::mutate(window_df, pval_perm = pval, qval_perm = qval)
  # Stash per-perm total counts as an attribute so the orchestrator can
  # render the histogram / density plots without re-running permutations.
  attr(result, "perm_totals") <- perm_totals
  result
}
