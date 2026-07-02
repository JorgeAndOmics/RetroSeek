# -----------------------------------------------------------------------------
# hotspot / models.R
# -----------------------------------------------------------------------------
# Per-label enrichment model: a Negative-Binomial GLM on per-window counts,
# mask-aware via an offset and optionally chromosome-stratified.
#
#   * `fit_nb_model()`     — fit `MASS::glm.nb` on per-window counts with an
#     `offset(log(effective_bp))` and an optional `chrom_stratum` covariate.
#     Convergence ladder: NB (`MASS::glm.nb`) succeeds with a stable theta ->
#     "ok"; otherwise -> "failed" (NA p-values + a logged warning in the
#     orchestrator). We deliberately do NOT fall back to Poisson: a Poisson
#     fit ignores overdispersion and is anti-conservative, so emitting NA is
#     the honest, reliability-first choice.
#   * `score_windows_nb()` — per-window upper-tail p-value from the fitted NB
#     plus BH adjustment.
#
# The model is fully deterministic — no RNG is involved — so results are
# reproducible by construction. The global `parameters.seed` is recorded in
# the run manifest for provenance but the core result does not depend on it.
#
# Pure module: callers pass tibbles / GRanges in, get tibbles out. No I/O.

suppressMessages({
  library(MASS)
  library(stats)
  library(tibble)
  library(dplyr)
  library(GenomicRanges)
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
#'   2. `glm.nb` fails or the theta-MLE warning fires -> status = "failed";
#'      NA p-values are emitted downstream and the orchestrator logs a warning.
#'      (No Poisson fallback: it would be anti-conservative under the
#'      overdispersion that NB exists to model.)
#'   3. Pre-filter leaves < min_nonzero non-zero windows -> status =
#'      "insufficient_data".
#'
#' @return list with:
#'   $model           — fitted glm.nb object or NULL
#'   $family          — "nb" | NA
#'   $theta           — NB dispersion or NA
#'   $status          — one of: ok | failed | insufficient_data
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


#' Score per-window enrichment given a fit. Adds mu_nb, pval_nb, qval_nb.
#'
#' p-value uses the upper tail of the fitted NB at observed count:
#'   pnbinom(count - 1, mu = mu_hat, size = theta, lower.tail = FALSE)
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
  if (any(callable) && identical(fit$family, "nb") && !is.na(fit$theta)) {
    pval[callable] <- stats::pnbinom(
      q = count[callable] - 1L,
      mu = mu[callable],
      size = fit$theta,
      lower.tail = FALSE
    )
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
