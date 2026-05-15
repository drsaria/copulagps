# =============================================================================
# estimate_gps.R — main GPS entry point; auto-dispatches on d
# =============================================================================

#' Estimate generalized propensity scores using copulas (d = 2 or d > 2)
#'
#' For d = 2, fits a bivariate copula (BiCopSelect).  For d > 2, fits a vine
#' copula (RVineStructureSelect) with optional parallel marginal fitting and
#' a configurable vine-refitting schedule.  All downstream functions
#' (`estimate_ate`, `assess_balance`, `run_chains`, etc.) work unchanged for
#' any d.
#'
#' @param data List from [prepare_data()], or a hand-built list with `$y`
#'   (n × d treatment matrix) and `$x` (list of d design matrices).
#' @param margins Character vector of length d specifying marginal
#'   distributions.  Options: `"normal"`, `"lognormal"`, `"t"`, `"gamma"`,
#'   `"invgauss"`.
#' @param params Optional list with pre-specified `$mean` and `$sd`.
#' @param copula_type d = 2 only. One of `"clayton"`, `"gumbel"`, `"frank"`,
#'   `"normal"`, `"t"`, `"auto"` (default).
#' @param vine_type d > 2 only. `"auto"` (R-vine, default), `"C"` (C-vine),
#'   or `"D"` (D-vine).
#' @param auto_adjust d = 2 only. Switch to residual copula when strong
#'   confounding is detected.  Default: `TRUE`.
#' @param stabilization Weight stabilisation: `"none"`, `"std"`, `"trim"`,
#'   or `"both"`.  Default: `"std"`.
#' @param trim_quantiles Quantile bounds for weight winsorization.
#'   Default: `c(0.01, 0.99)`.
#' @param trim_method Overlap trimming method applied to treatments.
#'   `"quantile"` (default, weight-space only), `"hull"` (convex hull, d = 2),
#'   `"mahalanobis"` (chi-sq distance, any d).
#' @param hull_trim Deprecated shortcut for d = 2 hull trimming.  Prefer
#'   `trim_method = "hull"`.  Default: `FALSE`.
#' @param mahal_quantile Chi-squared quantile for Mahalanobis trimming.
#'   Default: `0.975`.
#' @param df Degrees of freedom for Student-t marginals (`NULL` = estimate).
#' @param lod Optional LOD thresholds (length d, or scalar recycled).
#' @param n_impute_iter Maximum imputation iterations.  Default: `30L`.
#' @param conv_tol Convergence tolerance for imputation.  Default: `1e-4`.
#' @param vine_refit_every d > 2. Refit the vine every k imputation iterations.
#'   `1L` (default) = every iteration.  `Inf` = first iteration only
#'   (equivalent to `freeze_structure = TRUE`).
#' @param freeze_structure d > 2. When `TRUE`, freezes the vine tree structure
#'   after the first fit; subsequent iterations re-estimate pair copula
#'   parameters only.  Default: `FALSE`.
#' @param impute_method d > 2. `"gaussian_approx"` (default, fast) or
#'   `"vine_exact"` (slower but exact conditional imputation).
#' @param parallel_margins d > 2. Parallelise marginal fitting via
#'   `future_lapply`.  `NULL` = auto (enabled for d >= 4).
#' @param verbose Print informational messages.  Default: `TRUE`.
#'
#' @return Object of class `"gps_weights"`.
#' @export
estimate_gps <- function(
    data,
    margins,
    params            = NULL,
    # d = 2 copula options
    copula_type       = "auto",
    auto_adjust       = TRUE,
    # d > 2 vine options
    vine_type         = "auto",
    vine_refit_every  = 1L,
    freeze_structure  = FALSE,
    impute_method     = "gaussian_approx",
    parallel_margins  = NULL,
    # shared options
    stabilization     = "std",
    trim_quantiles    = c(0.01, 0.99),
    trim_method       = "quantile",
    hull_trim         = FALSE,       # legacy d=2 shortcut
    mahal_quantile    = 0.975,
    df                = NULL,
    lod               = NULL,
    n_impute_iter     = 30L,
    conv_tol          = 1e-4,
    verbose           = TRUE
) {
  if (!all(c("y", "x") %in% names(data)))
    stop("`data` must have components 'y' and 'x'; use prepare_data() to build it")

  y <- if (is.matrix(data$y)) data$y else as.matrix(data$y)
  d <- ncol(y)

  if (length(margins) != d)
    stop("'margins' must have length d = ", d)

  valid_margins <- c("normal", "lognormal", "t", "gamma", "invgauss")
  bad <- setdiff(margins, valid_margins)
  if (length(bad))
    stop("Invalid margins: ", paste(bad, collapse = ", "))

  valid_stab <- c("none", "std", "trim", "both")
  if (!stabilization %in% valid_stab)
    stop("'stabilization' must be one of: ", paste(valid_stab, collapse = ", "))

  # Back-compat: hull_trim = TRUE → trim_method = "hull"
  if (isTRUE(hull_trim) && trim_method == "quantile")
    trim_method <- "hull"

  # Hull trimming only valid for d = 2
  if (trim_method == "hull" && d > 2) {
    warning("trim_method = 'hull' is only valid for d = 2; switching to 'mahalanobis'")
    trim_method <- "mahalanobis"
  }

  # ---- Dispatch ---------------------------------------------------------------
  if (d == 2) {
    return(.estimate_gps_bivariate(
      data = data, y = y, d = d, margins = margins, params = params,
      copula_type = copula_type, auto_adjust = auto_adjust,
      stabilization = stabilization, trim_quantiles = trim_quantiles,
      trim_method = trim_method, mahal_quantile = mahal_quantile,
      df = df, lod = lod, n_impute_iter = n_impute_iter,
      conv_tol = conv_tol, verbose = verbose
    ))
  }

  .estimate_gps_vine(
    data = data, y = y, d = d, margins = margins, params = params,
    vine_type = vine_type, vine_refit_every = vine_refit_every,
    freeze_structure = freeze_structure, impute_method = impute_method,
    parallel_margins = parallel_margins,
    stabilization = stabilization, trim_quantiles = trim_quantiles,
    trim_method = trim_method, mahal_quantile = mahal_quantile,
    df = df, lod = lod, n_impute_iter = n_impute_iter,
    conv_tol = conv_tol, verbose = verbose
  )
}

# =============================================================================
# Bivariate (d = 2) implementation  — unchanged from original
# =============================================================================

#' @keywords internal
.estimate_gps_bivariate <- function(
    data, y, d, margins, params, copula_type, auto_adjust,
    stabilization, trim_quantiles, trim_method, mahal_quantile,
    df, lod, n_impute_iter, conv_tol, verbose
) {
  n <- nrow(y)

  valid_cop <- c("clayton", "gumbel", "frank", "normal", "t", "auto")
  if (!copula_type %in% valid_cop)
    stop("'copula_type' must be one of: ", paste(valid_cop, collapse = ", "))

  if (n < 50) warning("n < 50: copula fitting may be unstable")

  # LOD initialisation
  cens_mask <- matrix(FALSE, nrow = n, ncol = d)
  if (!is.null(lod)) {
    lod <- rep_len(as.numeric(lod), d)
    for (j in seq_len(d)) cens_mask[, j] <- y[, j] <= lod[j]
    n_cens <- colSums(cens_mask)
    if (verbose && any(n_cens > 0))
      message("LOD censoring: ",
              paste(paste0("T", seq_len(d), " = ", n_cens, "/", n,
                           " (", round(n_cens / n * 100, 1), "%)"),
                    collapse = ", "))
    y <- .init_censored(y, cens_mask, lod)
  }
  do_impute <- !is.null(lod) && any(cens_mask)
  n_iter    <- if (do_impute) max(1L, as.integer(n_impute_iter)) else 1L

  imp <- .run_imputation_loop(
    y = y, data = data, margins = margins, copula_type = copula_type,
    cens_mask = cens_mask, lod = lod, params = params,
    df = df, n_iter = n_iter, conv_tol = conv_tol, d = d, verbose = verbose
  )
  y             <- imp$y
  marginal_fits <- imp$marginal_fits
  n_iter_done   <- imp$n_iter_done
  y_imputed     <- if (do_impute) y else NULL

  use_residual_copula <- FALSE
  marginal_r2         <- rep(NA_real_, d)
  if (auto_adjust && is.null(params)) {
    adj                 <- .should_use_residual_copula(y, marginal_fits, copula_type, d, verbose)
    use_residual_copula <- adj$use
    marginal_r2         <- adj$r2
  }

  eps       <- .Machine$double.eps
  null_fits <- .fit_null_marginals(y, margins, marginal_fits, lod, eps)

  if (use_residual_copula) {
    resid_mat    <- sapply(seq_len(d), function(j) y[, j] - marginal_fits[[j]]$mean)
    u_emp        <- VineCopula::pobs(resid_mat)
    u_parametric <- sapply(seq_len(d), function(j) {
      std_r <- resid_mat[, j] / pmax(marginal_fits[[j]]$sd, eps)
      pmax(pmin(stats::pnorm(std_r), 1 - eps), eps)
    })
  } else {
    u_emp        <- sapply(seq_len(d), function(j) {
      df_null <- if (!is.na(null_fits[[j]]$df)) null_fits[[j]]$df else 10
      pmax(pmin(.pdist_margin(y[, j], margins[j], null_fits[[j]]$mean,
                               null_fits[[j]]$sd, df = df_null), 1 - eps), eps)
    })
    u_parametric <- sapply(seq_len(d), function(j) {
      df_j <- if (!is.na(marginal_fits[[j]]$df)) marginal_fits[[j]]$df else 10
      pmax(pmin(.pdist_margin(y[, j], margins[j], marginal_fits[[j]]$mean,
                               marginal_fits[[j]]$sd, df = df_j), 1 - eps), eps)
    })
  }

  fam_set <- .get_copula_family(copula_type)
  copula_emp <- tryCatch(
    VineCopula::BiCopSelect(u_emp[, 1], u_emp[, 2], familyset = fam_set),
    error = function(e) { warning("Empirical copula failed: ", e$message)
                           list(family = 0L, par = 0, par2 = 0) }
  )
  copula_par <- tryCatch(
    VineCopula::BiCopSelect(u_parametric[, 1], u_parametric[, 2], familyset = fam_set),
    error = function(e) { warning("Parametric copula failed: ", e$message)
                           list(family = 0L, par = 0, par2 = 0) }
  )

  dens_emp <- VineCopula::BiCopPDF(u_emp[, 1], u_emp[, 2],
                                    copula_emp$family, copula_emp$par, copula_emp$par2)
  dens_par <- VineCopula::BiCopPDF(u_parametric[, 1], u_parametric[, 2],
                                    copula_par$family, copula_par$par, copula_par$par2)

  log_marg_ratio <- .compute_log_marg_ratio(y, margins, marginal_fits, lod, lod, d, n, eps,
                                             null_fits = null_fits)
  marg_ratio <- exp(rowSums(log_marg_ratio))
  marg_ratio[!is.finite(marg_ratio) | marg_ratio <= 0] <- 1

  log_w_cop <- log(pmax(dens_emp, eps)) - log(pmax(dens_par, eps))
  log_w     <- log_w_cop + rowSums(log_marg_ratio)

  weights <- .finalise_weights(log_w, n, stabilization, trim_quantiles,
                                trim_method, mahal_quantile, y, verbose)

  w_valid <- weights[!is.na(weights)]
  ess     <- min(sum(w_valid)^2 / sum(w_valid^2), length(w_valid))

  overlap_mask <- attr(weights, "overlap_mask") %||% rep(TRUE, n)
  attr(weights, "overlap_mask") <- NULL

  structure(
    list(
      weights           = weights,
      marginal_fits     = marginal_fits,
      copula_empirical  = copula_emp,
      copula_parametric = copula_par,
      u_empirical       = u_emp,
      u_parametric      = u_parametric,
      data              = data,
      y_imputed         = y_imputed,
      hull_membership   = overlap_mask,
      diagnostics = list(
        effective_sample_size    = ess,
        weight_range             = range(w_valid),
        weight_cv                = stats::sd(w_valid) / mean(w_valid),
        params_estimated         = is.null(params),
        copula_method            = "bivariate",
        empirical_copula_family  = .family_name(copula_emp$family),
        parametric_copula_family = .family_name(copula_par$family),
        residual_copula_used     = use_residual_copula,
        marginal_r_squared       = marginal_r2,
        marg_ratio_range         = range(marg_ratio),
        n_censored               = colSums(cens_mask),
        n_overlap_trimmed        = sum(!overlap_mask),
        imputation_iters         = if (do_impute) n_iter_done else 0L
      ),
      specs = list(
        n = n, d = d, margins = margins, copula_type = copula_type,
        vine_type = NULL, stabilization = stabilization,
        trim_quantiles = trim_quantiles, trim_method = trim_method,
        auto_adjust = auto_adjust, df_specified = df,
        lod = lod, n_impute_iter = n_impute_iter, conv_tol = conv_tol
      )
    ),
    class = "gps_weights"
  )
}

# =============================================================================
# Shared helpers used by both d=2 and d>2 paths
# =============================================================================

#' Compute log marginal density ratios log[f_j(T_j)] - log[f_j(T_j|X_j)]
#' @keywords internal
# Fit intercept-only (null, covariate-free) marginal parameters for each treatment.
# Used both for the marginal log-ratio and for computing null-CDF-based copula uniforms.
.fit_null_marginals <- function(y, margins, marginal_fits, lod_vec, eps) {
  d <- ncol(y)
  lapply(seq_len(d), function(j) {
    yj         <- y[, j]
    use_tobit0 <- !is.null(lod_vec) && margins[j] %in% c("normal", "lognormal") &&
                  requireNamespace("AER", quietly = TRUE)
    if (use_tobit0) {
      if (margins[j] == "lognormal") {
        ly      <- log(pmax(yj, eps)); lod_log <- log(max(lod_vec[j], eps))
        fit0    <- tryCatch(AER::tobit(ly ~ 1, left = lod_log), error = function(e) NULL)
        if (!is.null(fit0))
          return(list(mean = unname(coef(fit0)["(Intercept)"]),
                      sd = fit0$scale, df = NA_real_))
        ly_obs <- ly[ly > lod_log]
        return(list(mean = mean(ly_obs), sd = stats::sd(ly_obs), df = NA_real_))
      } else {
        fit0 <- tryCatch(AER::tobit(yj ~ 1, left = lod_vec[j]), error = function(e) NULL)
        if (!is.null(fit0))
          return(list(mean = unname(coef(fit0)["(Intercept)"]),
                      sd = fit0$scale, df = NA_real_))
      }
    }
    switch(margins[j],
      normal    = list(mean = mean(yj),      sd = stats::sd(yj),  df = NA_real_),
      lognormal = { ly <- log(pmax(yj, eps))
                    list(mean = mean(ly),     sd = stats::sd(ly),  df = NA_real_) },
      t         = { df_j <- if (!is.na(marginal_fits[[j]]$df)) marginal_fits[[j]]$df else 10
                    list(mean = mean(yj),     sd = stats::mad(yj), df = df_j)   },
      gamma     = list(mean = mean(yj),      sd = sqrt(var(yj)),  df = NA_real_),
      invgauss  = list(mean = mean(yj),      sd = sqrt(var(yj)),  df = NA_real_)
    )
  })
}

.compute_log_marg_ratio <- function(y, margins, marginal_fits, lod, lod_vec,
                                     d, n, eps, null_fits = NULL) {
  if (is.null(null_fits))
    null_fits <- .fit_null_marginals(y, margins, marginal_fits, lod_vec, eps)

  vapply(seq_len(d), function(j) {
    df_null <- if (!is.na(null_fits[[j]]$df))      null_fits[[j]]$df     else 10
    df_cond <- if (!is.na(marginal_fits[[j]]$df))  marginal_fits[[j]]$df else 10
    lf_marg <- log(pmax(.ddist_margin(y[, j], margins[j], null_fits[[j]]$mean,
                                       null_fits[[j]]$sd, df = df_null), eps))
    lf_cond <- log(pmax(.ddist_margin(y[, j], margins[j], marginal_fits[[j]]$mean,
                                       marginal_fits[[j]]$sd, df = df_cond), eps))
    lf_marg - lf_cond
  }, numeric(n))
}

#' Apply positivity cap, stabilisation, and overlap trimming to log-weights
#'
#' Returns the final weight vector with NA for trimmed observations.
#' Attaches an "overlap_mask" attribute (logical, TRUE = kept).
#' @keywords internal
.finalise_weights <- function(log_w, n, stabilization, trim_quantiles,
                               trim_method, mahal_quantile, y, verbose) {
  eps <- .Machine$double.eps

  log_cap  <- 0.25 * log(n)
  n_capped <- sum(log_w > log_cap, na.rm = TRUE)
  if (n_capped > 0 && verbose)
    message(n_capped, " weight(s) exceeded positivity cap [log(sqrt(n)) = ",
            round(log_cap, 2), "] — possible overlap violation")
  log_w <- pmin(log_w, log_cap)

  weights     <- exp(log_w)
  finite_mask <- is.finite(weights) & weights > 0
  if (!all(finite_mask)) {
    warning(sum(!finite_mask), "/", n, " non-finite weights replaced with 1")
    weights[!finite_mask] <- 1
  }

  if (stabilization %in% c("std", "both")) {
    mw <- mean(weights, na.rm = TRUE)
    if (mw > 0) weights <- weights / mw
  }
  if (stabilization %in% c("trim", "both")) {
    q       <- stats::quantile(weights, probs = trim_quantiles, na.rm = TRUE)
    weights <- pmin(pmax(weights, q[1]), q[2])
  }

  # Overlap trimming (treatment-space)
  overlap_mask <- rep(TRUE, n)
  if (trim_method == "hull") {
    overlap_mask <- .check_convex_hull(y)
  } else if (trim_method == "mahalanobis") {
    overlap_mask <- .check_mahalanobis_trim(y, quantile = mahal_quantile)
  }

  n_trimmed <- sum(!overlap_mask)
  if (n_trimmed > 0) {
    weights[!overlap_mask] <- NA_real_
    pct <- round(100 * n_trimmed / n, 1)
    if (verbose)
      message(trim_method, " trimming: ", n_trimmed, " obs excluded (", pct, "%)")
    if (n_trimmed / n > 0.05)
      warning(trim_method, " trimming removed > 5% of observations (",
              pct, "%) — consider checking treatment overlap")
  }

  attr(weights, "overlap_mask") <- overlap_mask
  weights
}
