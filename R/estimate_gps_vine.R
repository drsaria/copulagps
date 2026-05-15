# =============================================================================
# estimate_gps_vine.R — vine copula GPS engine for d > 2
# =============================================================================

#' Estimate GPS via vine copula (d > 2)
#'
#' Called automatically by [estimate_gps()] when d > 2.  Fits marginal models,
#' runs the iterative vine imputation loop (if LOD censoring is present), then
#' fits a vine copula for both the empirical and parametric uniforms to compute
#' GPS weights in log space.
#'
#' @inheritParams estimate_gps
#' @param y n × d treatment matrix (already coerced to matrix)
#' @param d number of treatments (>= 3)
#' @return Object of class `"gps_weights"`
#' @keywords internal
.estimate_gps_vine <- function(
    data, y, d, margins, params,
    vine_type, vine_refit_every, freeze_structure, impute_method,
    parallel_margins,
    stabilization, trim_quantiles, trim_method, mahal_quantile,
    df, lod, n_impute_iter, conv_tol, verbose
) {
  n   <- nrow(y)
  eps <- .Machine$double.eps

  if (n < 50)
    warning("n < 50: vine copula fitting may be unstable with d = ", d, " treatments")

  # ---- LOD initialisation -----------------------------------------------------
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

  # ---- Iterative imputation loop ----------------------------------------------
  imp <- .run_imputation_loop_vine(
    y = y, data = data, margins = margins, vine_type = vine_type,
    cens_mask = cens_mask, lod = lod, params = params,
    df = df, n_iter = n_iter, conv_tol = conv_tol,
    vine_refit_every = vine_refit_every, freeze_structure = freeze_structure,
    impute_method = impute_method, parallel_margins = parallel_margins,
    d = d, verbose = verbose
  )
  y             <- imp$y
  marginal_fits <- imp$marginal_fits
  n_iter_done   <- imp$n_iter_done
  vine_imp      <- imp$vine          # vine fit from imputation loop (reused if no LOD)
  y_imputed     <- if (do_impute) y else NULL

  # ---- GPS vine fit on final imputed data -------------------------------------
  # Null (covariate-free) marginal fits — used for u_emp and log-marginal-ratio
  null_fits <- .fit_null_marginals(y, margins, marginal_fits, lod, eps)

  # Empirical uniforms via null marginal CDFs (avoids X-contamination from pobs on imputed y)
  u_emp <- sapply(seq_len(d), function(k) {
    df_null <- if (!is.na(null_fits[[k]]$df)) null_fits[[k]]$df else 10
    pmax(pmin(.pdist_margin(y[, k], margins[k], null_fits[[k]]$mean,
                             null_fits[[k]]$sd, df = df_null), 1 - eps), eps)
  })

  # Parametric uniforms via fitted marginal CDFs
  u_parametric <- sapply(seq_len(d), function(k) {
    df_k <- if (!is.na(marginal_fits[[k]]$df)) marginal_fits[[k]]$df else 10
    pmax(pmin(.pdist_margin(y[, k], margins[k], marginal_fits[[k]]$mean,
                             marginal_fits[[k]]$sd, df = df_k), 1 - eps), eps)
  })

  fam_set       <- c(1L, 3L:5L)   # Gaussian, Clayton, Gumbel, Frank
  vine_type_code <- .vine_type_code(vine_type)

  if (verbose) message("Fitting vine copula for GPS weight computation...")

  vine_emp <- tryCatch(
    VineCopula::RVineStructureSelect(u_emp, familyset = fam_set,
                                      type = vine_type_code),
    error = function(e) {
      warning("Vine fit (empirical) failed: ", e$message, " — using Gaussian vine")
      VineCopula::RVineStructureSelect(u_emp, familyset = 1L, type = 0L)
    }
  )

  vine_par <- tryCatch(
    VineCopula::RVineCopSelect(u_parametric, Matrix = vine_emp$Matrix,
                                familyset = fam_set),
    error = function(e) {
      warning("Vine fit (parametric) failed: ", e$message,
              " — using empirical vine structure with Gaussian families")
      VineCopula::RVineStructureSelect(u_parametric, familyset = 1L, type = 0L)
    }
  )

  # ---- Densities + log-weights ------------------------------------------------
  dens_emp <- tryCatch(
    VineCopula::RVinePDF(u_emp, vine_emp),
    error = function(e) { warning("Empirical vine PDF failed: ", e$message); rep(1, n) }
  )
  dens_par <- tryCatch(
    VineCopula::RVinePDF(u_parametric, vine_par),
    error = function(e) { warning("Parametric vine PDF failed: ", e$message); rep(1, n) }
  )

  dens_emp <- pmax(dens_emp, eps)
  dens_par <- pmax(dens_par, eps)

  log_marg_ratio <- .compute_log_marg_ratio(y, margins, marginal_fits,
                                              lod, lod, d, n, eps, null_fits = null_fits)
  log_w_cop <- log(dens_emp) - log(dens_par)
  log_w     <- log_w_cop + rowSums(log_marg_ratio)

  # ---- Finalise weights -------------------------------------------------------
  weights <- .finalise_weights(log_w, n, stabilization, trim_quantiles,
                                trim_method, mahal_quantile, y, verbose)

  w_valid <- weights[!is.na(weights)]
  ess     <- if (length(w_valid) > 0)
               min(sum(w_valid)^2 / sum(w_valid^2), length(w_valid))
             else 0

  overlap_mask <- attr(weights, "overlap_mask") %||% rep(TRUE, n)
  attr(weights, "overlap_mask") <- NULL

  # ---- Vine Kendall tau summary -----------------------------------------------
  vine_tau_mat <- tryCatch(VineCopula::RVinePar2Tau(vine_emp),
                            error = function(e) matrix(NA_real_, d, d))
  diag(vine_tau_mat) <- 1

  # ---- Return -----------------------------------------------------------------
  marg_ratio <- exp(rowSums(log_marg_ratio))
  marg_ratio[!is.finite(marg_ratio) | marg_ratio <= 0] <- 1

  structure(
    list(
      weights           = weights,
      marginal_fits     = marginal_fits,
      vine_empirical    = vine_emp,
      vine_parametric   = vine_par,
      vine_tau_matrix   = vine_tau_mat,
      u_empirical       = u_emp,
      u_parametric      = u_parametric,
      data              = data,
      y_imputed         = y_imputed,
      hull_membership   = overlap_mask,
      diagnostics = list(
        effective_sample_size = ess,
        weight_range          = range(w_valid),
        weight_cv             = if (length(w_valid) > 1)
                                  stats::sd(w_valid) / mean(w_valid) else NA_real_,
        params_estimated      = is.null(params),
        copula_method         = "vine",
        vine_type             = vine_type,
        vine_families_emp     = vine_emp$family,
        marg_ratio_range      = range(marg_ratio),
        n_censored            = colSums(cens_mask),
        n_overlap_trimmed     = sum(!overlap_mask),
        imputation_iters      = if (do_impute) n_iter_done else 0L
      ),
      specs = list(
        n = n, d = d, margins = margins, copula_type = NULL,
        vine_type = vine_type, vine_refit_every = vine_refit_every,
        freeze_structure = freeze_structure, impute_method = impute_method,
        stabilization = stabilization, trim_quantiles = trim_quantiles,
        trim_method = trim_method, mahal_quantile = mahal_quantile,
        df_specified = df, lod = lod,
        n_impute_iter = n_impute_iter, conv_tol = conv_tol
      )
    ),
    class = "gps_weights"
  )
}
