# =============================================================================
# impute.R — LOD censoring initialisation and copula-based imputation
# =============================================================================

#' Impute censored treatment values using the fitted copula conditional distribution
#'
#' Draws from P(T_j | T_j < LOD_j, T_{-j}, X) by inverting the copula
#' h-function sequentially (j = 1 then j = 2).  After imputing T_1, the
#' updated value is immediately used as the conditioning partner when
#' imputing T_2 (Algorithm A2).
#'
#' @param y n x d treatment matrix (current values; censored cells contain
#'   Beta-initialised or previous-iteration imputed values)
#' @param cens_mask n x d logical matrix; TRUE = censored
#' @param lod numeric vector of length d
#' @param mfits list of marginal fit objects from `.fit_marginal()`
#' @param margins character vector of length d
#' @param cop BiCop object (family, par, par2)
#' @param d number of treatments (must be 2)
#' @return Updated n x d matrix with imputed values for censored cells
#' @keywords internal
.impute_censored_copula <- function(y, cens_mask, lod, mfits, margins, cop, d) {
  y_work <- y
  eps    <- .Machine$double.eps

  for (j in seq_len(d)) {
    idx <- which(cens_mask[, j])
    if (length(idx) == 0) next

    j2   <- ifelse(j == 1L, 2L, 1L)
    df_j  <- if (!is.na(mfits[[j]]$df))  mfits[[j]]$df  else 10
    df_j2 <- if (!is.na(mfits[[j2]]$df)) mfits[[j2]]$df else 10

    u_lod <- pmax(pmin(
      .pdist_margin(rep(lod[j], length(idx)), margins[j],
                    mfits[[j]]$mean[idx], mfits[[j]]$sd, df = df_j),
      1 - eps), eps)

    u_partner <- pmax(pmin(
      .pdist_margin(y_work[idx, j2], margins[j2],
                    mfits[[j2]]$mean[idx], mfits[[j2]]$sd, df = df_j2),
      1 - eps), eps)

    h_upper <- if (j == 1L)
      VineCopula::BiCopHfunc1(u_lod, u_partner, cop$family, cop$par, cop$par2)
    else
      VineCopula::BiCopHfunc2(u_partner, u_lod, cop$family, cop$par, cop$par2)
    h_upper <- pmax(h_upper, eps)

    v     <- pmax(pmin(stats::runif(length(idx)) * h_upper, 1 - eps), eps)
    u_imp <- if (j == 1L)
      VineCopula::BiCopHinv1(v, u_partner, cop$family, cop$par, cop$par2)
    else
      VineCopula::BiCopHinv2(u_partner, v, cop$family, cop$par, cop$par2)
    u_imp <- pmax(pmin(u_imp, 1 - eps), eps)

    t_imp <- .qdist_margin(u_imp, margins[j],
                            mfits[[j]]$mean[idx], mfits[[j]]$sd, df = df_j)
    y_work[idx, j] <- pmin(pmax(t_imp, eps), lod[j] - eps)
  }
  y_work
}

#' Run the iterative copula imputation loop
#'
#' Fits marginals and re-imputes censored values until the maximum relative
#' change in imputed cells drops below `conv_tol` or `n_iter` is reached.
#'
#' @param y n x d treatment matrix (censored cells initialised with Beta draws)
#' @param data internal data list (for design matrices)
#' @param margins character vector of length d
#' @param copula_type copula family string
#' @param cens_mask n x d logical censoring mask
#' @param lod numeric vector of length d
#' @param params optional user-supplied parameters list
#' @param df degrees of freedom for t margin
#' @param n_iter maximum iterations
#' @param conv_tol convergence tolerance
#' @param d number of treatments
#' @param verbose print iteration messages
#' @return List: y (final imputed matrix), marginal_fits, n_iter_done
#' @keywords internal
.run_imputation_loop <- function(y, data, margins, copula_type, cens_mask, lod,
                                  params, df, n_iter, conv_tol, d, verbose) {
  n_iter_done <- 1L
  cop_prev    <- NULL

  for (iter in seq_len(n_iter)) {
    if (is.null(params)) {
      marginal_fits <- lapply(seq_len(d), function(j) {
        lod_j <- if (!is.null(lod)) lod[j] else NULL
        tryCatch(
          .fit_marginal(y[, j], data$x[[j]], margins[j], df = df, lod = lod_j),
          error = function(e) stop("Marginal ", j, " (", margins[j], "): ", e$message)
        )
      })
    } else {
      n <- nrow(y)
      marginal_fits <- .process_user_params(params, y, margins, n, d)
    }

    if (iter < n_iter) {
      fam_set <- .get_copula_family(copula_type)
      eps     <- .Machine$double.eps
      u_cond  <- sapply(seq_len(d), function(j) {
        df_j <- if (!is.na(marginal_fits[[j]]$df)) marginal_fits[[j]]$df else 10
        pmax(pmin(.pdist_margin(y[, j], margins[j], marginal_fits[[j]]$mean,
                                marginal_fits[[j]]$sd, df = df_j),
                  1 - eps), eps)
      })
      cop_imp <- tryCatch(
        VineCopula::BiCopSelect(u_cond[, 1], u_cond[, 2], familyset = fam_set),
        error = function(e) {
          warning("Imputation copula failed (iter ", iter, "): ", e$message,
                  " — using Gaussian")
          list(family = 1L, par = 0, par2 = 0)
        }
      )
      y_new <- .impute_censored_copula(y, cens_mask, lod, marginal_fits,
                                        margins, cop_imp, d)

      # Convergence: relative change in copula dependence parameter.
      # Checking per-observation imputed-value changes is incorrect for a
      # stochastic algorithm (fresh random draws change O(sigma_T) each
      # iteration, never converging).  The marginal parameters decouple from
      # imputed values via the Tobit likelihood and are consistent from
      # iteration 1; only alpha needs monitoring.
      if (!is.null(cop_prev)) {
        par_chg <- abs(cop_imp$par - cop_prev$par) /
                     max(abs(cop_prev$par), eps)
        if (verbose)
          message("Imputation iter ", iter, ": copula par change = ",
                  formatC(par_chg, format = "e", digits = 2),
                  "  (alpha = ", round(cop_imp$par, 4), ")")
        if (par_chg < conv_tol) {
          y           <- y_new
          n_iter_done <- iter
          if (verbose) message("Converged at iteration ", iter,
                               " (copula par = ", round(cop_imp$par, 4), ")")
          break
        }
      } else {
        if (verbose) message("Imputation iteration ", iter, "/", n_iter - 1L)
      }

      y           <- y_new
      cop_prev    <- cop_imp
      n_iter_done <- iter
    }
  }

  list(y = y, marginal_fits = marginal_fits, n_iter_done = n_iter_done)
}

#' Initialise censored cells with Beta(2,1) draws scaled to [0, LOD)
#'
#' @param y n x d treatment matrix
#' @param cens_mask n x d logical censoring mask
#' @param lod numeric vector of length d
#' @return y with censored cells replaced by Beta(2,1) * LOD initialisations
#' @keywords internal
.init_censored <- function(y, cens_mask, lod) {
  d <- ncol(y)
  for (j in seq_len(d)) {
    nc <- sum(cens_mask[, j])
    if (nc > 0)
      y[cens_mask[, j], j] <- stats::rbeta(nc, 2, 1) * lod[j]
  }
  y
}
