# =============================================================================
# balance.R — covariate balance assessment and plots
# =============================================================================

#' Weighted F-statistics for T_j ~ X per treatment (internal)
#'
#' Fits `lm(T_j ~ ., data = X_df, weights = w)` for each treatment column and
#' returns the F-statistic before (unweighted) and after (GPS-weighted) testing
#' H0: T_j is linearly uncorrelated with all covariates.  Under perfect balance
#' the weighted F converges to 1 in expectation regardless of n and p.
#'
#' Categorical covariates must already be expanded to dummy columns before
#' calling this function (as done by `assess_balance()`).
#'
#' @param T_mat  n × d numeric treatment matrix.
#' @param X_mat  n × p numeric covariate matrix (factors already expanded).
#' @param w_valid  n-length weight vector; NAs treated as 0.
#' @return data.frame with one row per treatment: `treatment`, `f_pre`,
#'   `f_post`, `p_pre`, `p_post`, `df1`, `df2`.
#' @keywords internal
.compute_fstats <- function(T_mat, X_mat, w_valid) {
  X_df  <- as.data.frame(X_mat)
  n_trt <- ncol(T_mat)
  trt_nm <- if (!is.null(colnames(T_mat))) colnames(T_mat) else
              paste0("T", seq_len(n_trt))

  do.call(rbind, lapply(seq_len(n_trt), function(j) {
    df_j <- cbind(.T = T_mat[, j], X_df)

    get_fstat <- function(wts = NULL) {
      tryCatch({
        fit <- if (is.null(wts)) stats::lm(.T ~ ., data = df_j) else
                 stats::lm(.T ~ ., data = df_j, weights = wts)
        fs  <- summary(fit)$fstatistic   # [F, df1, df2]
        c(f   = unname(fs[1L]),
          df1 = unname(fs[2L]),
          df2 = unname(fs[3L]))
      }, error = function(e)
        c(f = NA_real_, df1 = NA_real_, df2 = NA_real_))
    }

    pre  <- get_fstat()
    post <- get_fstat(w_valid)

    data.frame(
      treatment = trt_nm[j],
      f_pre     = pre["f"],
      f_post    = post["f"],
      p_pre     = stats::pf(pre["f"],  pre["df1"],  pre["df2"],  lower.tail = FALSE),
      p_post    = stats::pf(post["f"], post["df1"], post["df2"], lower.tail = FALSE),
      df1       = pre["df1"],
      df2       = pre["df2"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}

#' Vectorised weighted correlation between two numeric matrices
#' @keywords internal
.weighted_cor_matrix <- function(X, Y, w = NULL) {
  if (is.null(w)) return(stats::cor(X, Y))
  w  <- w / sum(w)
  XY <- cbind(X, Y)
  p  <- ncol(X); q <- ncol(Y)
  mu    <- colSums(XY * w)
  XY_c  <- sweep(XY, 2L, mu)
  sw    <- sqrt(w)
  Sigma <- crossprod(XY_c * sw)
  stds  <- sqrt(diag(Sigma))
  Cor   <- Sigma / outer(stds, stds)
  Cor[seq_len(p), p + seq_len(q), drop = FALSE]
}

#' Assess covariate balance before and after GPS weighting
#'
#' Computes weighted and unweighted correlations between covariates and
#' treatments, summarising balance improvement.
#'
#' @param gps Object of class `"gps_weights"`
#' @param covariates Matrix or data frame of covariates (same rows as treatments)
#' @param treatments Optional treatment matrix; defaults to `gps$data$y`
#' @return Object of class `"balance_assessment"`
#' @export
assess_balance <- function(gps, covariates, treatments = NULL) {
  if (is.null(treatments)) treatments <- gps$data$y
  weights <- gps$weights

  if (is.vector(covariates)) covariates <- matrix(covariates, ncol = 1L)
  if (is.data.frame(covariates)) {
    fac <- sapply(covariates, is.factor)
    if (any(fac)) {
      fac_mats   <- lapply(names(covariates)[fac], function(v)
                      stats::model.matrix(~ . - 1, data = covariates[v]))
      covariates <- cbind(as.matrix(covariates[!fac]), do.call(cbind, fac_mats))
    } else {
      covariates <- as.matrix(covariates)
    }
  }

  # Use only non-NA-weighted observations for the "after" correlation
  w_valid <- weights
  w_valid[is.na(w_valid)] <- 0

  corr_before <- .weighted_cor_matrix(covariates, treatments)
  corr_after  <- .weighted_cor_matrix(covariates, treatments, w = w_valid)

  mk <- function(corr) list(
    correlations  = corr,
    max_abs_corr  = max(abs(corr)),
    mean_abs_corr = mean(abs(corr))
  )
  before      <- mk(corr_before)
  after       <- mk(corr_after)
  pct_red     <- function(a, b) (a - b) / max(a, .Machine$double.eps) * 100
  improvement <- list(
    corr_reduction     = pct_red(before$mean_abs_corr, after$mean_abs_corr),
    max_corr_reduction = pct_red(before$max_abs_corr,  after$max_abs_corr)
  )

  # F-statistic balance: lm(T_j ~ X, weights=w) per treatment.
  # Covariates matrix has already been expanded (factors → dummies) above.
  fstats <- tryCatch(
    .compute_fstats(treatments, covariates, w_valid),
    error = function(e) NULL
  )

  structure(
    list(before = before, after = after, improvement = improvement,
         fstats          = fstats,
         n_covariates    = ncol(covariates), n_treatments = ncol(treatments),
         covariate_names = colnames(covariates),
         treatment_names = colnames(treatments)),
    class = "balance_assessment"
  )
}

#' Summarise covariate balance using MSMD and weighted F-statistics
#'
#' Returns a compact scalar summary of balance: the maximum absolute
#' treatment-covariate correlation (MSMD) before and after GPS weighting,
#' per-treatment max |r|, and — when available — the weighted F-statistics
#' from `lm(T_j ~ X, weights = w)` following Brown et al. (2021).
#'
#' The F-statistic complement handles categorical covariates correctly and
#' captures any linear T–X association; F_post ≈ 1 indicates good balance.
#'
#' @param balance Object of class `"balance_assessment"`
#' @return A list with elements `msmd_before`, `msmd_after`,
#'   `per_treatment_before`, `per_treatment_after`, `reduction_pct`,
#'   and (if F-statistics were computed) `fstat_pre`, `fstat_post`,
#'   `fstats_per_treatment`.
#' @export
summarise_balance <- function(balance) {
  cor_b <- abs(balance$before$correlations)
  cor_a <- abs(balance$after$correlations)

  msmd_before <- max(cor_b, na.rm = TRUE)
  msmd_after  <- max(cor_a, na.rm = TRUE)

  per_trt_before <- apply(cor_b, 2, max, na.rm = TRUE)
  per_trt_after  <- apply(cor_a, 2, max, na.rm = TRUE)

  reduction_pct <- (msmd_before - msmd_after) /
                   max(msmd_before, .Machine$double.eps) * 100

  out <- list(
    msmd_before          = msmd_before,
    msmd_after           = msmd_after,
    per_treatment_before = per_trt_before,
    per_treatment_after  = per_trt_after,
    reduction_pct        = reduction_pct
  )

  # Append F-statistic summaries when available (backward compatible)
  fs <- balance$fstats
  if (!is.null(fs)) {
    out$fstat_pre            <- max(fs$f_pre,  na.rm = TRUE)
    out$fstat_post           <- max(fs$f_post, na.rm = TRUE)
    out$fstats_per_treatment <- fs
  }

  out
}

#' Plot covariate balance (ggplot2 + patchwork)
#'
#' Four-panel figure showing before vs after correlations, per-variable
#' improvement, and marginal correlation distributions.  When `summary_only =
#' TRUE` (recommended for d > 4), plots per-treatment MSMD bars instead of the
#' full scatter panel.
#'
#' @param balance Object of class `"balance_assessment"`
#' @param max_vars Maximum number of variables to display (default: 20)
#' @param summary_only When `TRUE`, show per-treatment MSMD bar chart only.
#'   Default: `FALSE`; auto-set to `TRUE` when d > 4.
#' @return A patchwork or ggplot2 object (printed and returned invisibly)
#' @export
plot_balance <- function(balance, max_vars = 20L, summary_only = FALSE) {
  if (!summary_only && balance$n_treatments > 4L) {
    message("d = ", balance$n_treatments,
            " treatments: switching to summary_only = TRUE for readability. ",
            "Pass summary_only = FALSE to override.")
    summary_only <- TRUE
  }

  if (summary_only) {
    sm  <- summarise_balance(balance)
    trt <- names(sm$per_treatment_before) %||%
           paste0("T", seq_along(sm$per_treatment_before))
    df_sm <- data.frame(
      treatment = rep(trt, 2),
      max_r     = c(sm$per_treatment_before, sm$per_treatment_after),
      timing    = rep(c("Before", "After"), each = length(trt))
    )
    df_sm$timing <- factor(df_sm$timing, levels = c("Before", "After"))

    p <- ggplot2::ggplot(df_sm, ggplot2::aes(x = treatment, y = max_r, fill = timing)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::scale_fill_manual(values = c("Before" = "#fc8d59", "After" = "#91bfdb"),
                                  name = NULL) +
      ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey40") +
      ggplot2::labs(
        title = paste0("Per-treatment max |r|  \u2014  MSMD reduction: ",
                       round(sm$reduction_pct, 1), "%"),
        x = "Treatment", y = "Max |r| with covariates"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    print(p)
    return(invisible(p))
  }
  cb <- as.vector(balance$before$correlations)
  ca <- as.vector(balance$after$correlations)

  if (length(cb) > max_vars) {
    idx <- order(abs(cb), decreasing = TRUE)[seq_len(max_vars)]
    cb <- cb[idx]; ca <- ca[idx]
  }

  df_sc  <- data.frame(before = cb, after = ca)
  df_imp <- data.frame(
    var      = seq_along(cb),
    delta    = abs(cb) - abs(ca),
    improved = (abs(cb) - abs(ca)) > 0
  )

  p1 <- ggplot2::ggplot(df_sc, ggplot2::aes(x = before, y = after)) +
    ggplot2::geom_point(alpha = 0.6, colour = "#2c7bb6") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted") +
    ggplot2::labs(title = "Before vs After Weighting",
                  x = "Correlation (before)", y = "Correlation (after)") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df_imp, ggplot2::aes(x = var, y = delta, fill = improved)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1a9641", "FALSE" = "#d7191c"),
      labels = c("TRUE" = "Improved", "FALSE" = "Worsened"), name = NULL
    ) +
    ggplot2::labs(title = "Improvement per Variable",
                  x = "Variable index", y = "Reduction in |r|") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  p3 <- ggplot2::ggplot(data.frame(r = cb), ggplot2::aes(x = r)) +
    ggplot2::geom_histogram(bins = 20L, fill = "#fc8d59", colour = "white") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Correlations Before", x = "Correlation", y = "Count") +
    ggplot2::theme_minimal()

  p4 <- ggplot2::ggplot(data.frame(r = ca), ggplot2::aes(x = r)) +
    ggplot2::geom_histogram(bins = 20L, fill = "#91bfdb", colour = "white") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Correlations After", x = "Correlation", y = "Count") +
    ggplot2::theme_minimal()

  combined <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = paste0("Covariate Balance \u2014 mean |r| reduction: ",
                     round(balance$improvement$corr_reduction, 1), "%"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  print(combined)
  invisible(combined)
}
