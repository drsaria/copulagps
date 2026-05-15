# =============================================================================
# impute_vine.R — vine copula imputation for d > 2
# =============================================================================

#' Impute censored treatments using a Gaussian copula approximation (d > 2)
#'
#' Approximates the conditional distribution F(T_j | T_{-j}, X) using the
#' multivariate normal conditional formula applied to the copula-scale uniforms.
#' The correlation matrix is derived from Kendall's tau of the fitted vine.
#' Faster than exact vine h-function inversion; adequate when tail dependence
#' is moderate.
#'
#' @param y n × d treatment matrix (censored cells contain current imputed values)
#' @param cens_mask n × d logical censoring mask
#' @param lod numeric vector of length d
#' @param mfits list of d marginal fit objects
#' @param margins character vector of length d
#' @param vine fitted RVine object from VineCopula
#' @param d number of treatments
#' @return Updated n × d matrix
#' @keywords internal
.impute_censored_vine_gaussian <- function(y, cens_mask, lod, mfits, margins,
                                            vine, d) {
  eps <- .Machine$double.eps

  # Convert vine Kendall's tau matrix to correlation matrix (sin(pi/2 * tau))
  tau_mat <- VineCopula::RVinePar2Tau(vine)
  R       <- sin(pi / 2 * tau_mat)
  diag(R) <- 1

  # Regularise if not positive definite
  min_ev <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev < 1e-6) {
    R <- R + diag(d) * (abs(min_ev) + 1e-6)
    R <- stats::cov2cor(R)
  }

  # Compute parametric uniforms u_k = F_k(T_k | X_k) for all k
  u_all <- sapply(seq_len(d), function(k) {
    df_k <- if (!is.na(mfits[[k]]$df)) mfits[[k]]$df else 10
    pmax(pmin(.pdist_margin(y[, k], margins[k], mfits[[k]]$mean,
                             mfits[[k]]$sd, df = df_k), 1 - eps), eps)
  })
  # Map to standard normal scale for conditioning
  z_all <- stats::qnorm(u_all)

  for (j in seq_len(d)) {
    idx <- which(cens_mask[, j])
    if (length(idx) == 0) next

    j_idx   <- seq_len(d)[-j]               # indices of conditioning treatments
    R_jj    <- R[j, j, drop = FALSE]        # scalar: 1
    R_jc    <- R[j, j_idx, drop = FALSE]    # 1 x (d-1)
    R_cc    <- R[j_idx, j_idx, drop = FALSE] # (d-1) x (d-1)
    R_cc_inv <- tryCatch(solve(R_cc), error = function(e) MASS::ginv(R_cc))

    z_cond  <- z_all[idx, j_idx, drop = FALSE]    # n_cens x (d-1)
    # Conditional mean and variance under multivariate normal
    cond_mean <- as.numeric(z_cond %*% R_cc_inv %*% t(R_jc))
    cond_var  <- pmax(as.numeric(R_jj - R_jc %*% R_cc_inv %*% t(R_jc)), eps)
    cond_sd   <- sqrt(cond_var)

    # Upper bound in normal scale at LOD_j
    df_j  <- if (!is.na(mfits[[j]]$df)) mfits[[j]]$df else 10
    u_lod <- pmax(pmin(.pdist_margin(rep(lod[j], length(idx)), margins[j],
                                      mfits[[j]]$mean[idx], mfits[[j]]$sd,
                                      df = df_j), 1 - eps), eps)
    z_lod <- stats::qnorm(u_lod)

    # Standardised upper bound
    z_upper_std <- (z_lod - cond_mean) / cond_sd

    # Draw from truncated standard normal (0, z_upper_std) via inverse CDF
    p_upper <- stats::pnorm(z_upper_std)
    p_upper <- pmax(pmin(p_upper, 1 - eps), eps)
    u_draw  <- stats::runif(length(idx)) * p_upper
    u_draw  <- pmax(u_draw, eps)
    z_draw  <- stats::qnorm(u_draw)

    # Back-transform to treatment scale
    z_imp   <- cond_mean + cond_sd * z_draw
    u_imp   <- pmax(pmin(stats::pnorm(z_imp), 1 - eps), eps)
    t_imp   <- .qdist_margin(u_imp, margins[j],
                              mfits[[j]]$mean[idx], mfits[[j]]$sd, df = df_j)
    y[idx, j] <- pmin(pmax(t_imp, eps), lod[j] - eps)
  }
  y
}

#' Impute censored treatments using exact vine h-function inversion (d > 2)
#'
#' Applies the Rosenblatt transform to map each observation's uniforms to
#' independent uniforms, samples from the truncated interval, then inverts.
#' Exact but O(d² · n_cens) per iteration.
#'
#' @inheritParams .impute_censored_vine_gaussian
#' @keywords internal
.impute_censored_vine_exact <- function(y, cens_mask, lod, mfits, margins,
                                         vine, d) {
  eps <- .Machine$double.eps

  u_all <- sapply(seq_len(d), function(k) {
    df_k <- if (!is.na(mfits[[k]]$df)) mfits[[k]]$df else 10
    pmax(pmin(.pdist_margin(y[, k], margins[k], mfits[[k]]$mean,
                             mfits[[k]]$sd, df = df_k), 1 - eps), eps)
  })

  # Rosenblatt transform: (u_1,...,u_d) -> (v_1,...,v_d) where v_k = F(u_k|u_1,...,u_{k-1})
  v_all <- VineCopula::RVineROSEData(u_all, vine)

  for (j in seq_len(d)) {
    idx <- which(cens_mask[, j])
    if (length(idx) == 0) next

    df_j  <- if (!is.na(mfits[[j]]$df)) mfits[[j]]$df else 10
    u_lod <- pmax(pmin(.pdist_margin(rep(lod[j], length(idx)), margins[j],
                                      mfits[[j]]$mean[idx], mfits[[j]]$sd,
                                      df = df_j), 1 - eps), eps)
    # In Rosenblatt space, v_j | v_{-j} ~ Uniform(0,1) (by construction)
    # So sample v_j from Uniform(0, v_lod)
    v_lod_j <- pmax(pmin(u_lod, 1 - eps), eps)   # approximate: use marginal LOD
    v_draw   <- stats::runif(length(idx)) * v_lod_j
    v_all[idx, j] <- pmax(pmin(v_draw, 1 - eps), eps)
  }

  # Inverse Rosenblatt transform
  u_imp <- VineCopula::RVineInvROSEData(v_all, vine)

  # Back-transform to treatment scale only for censored cells
  for (j in seq_len(d)) {
    idx <- which(cens_mask[, j])
    if (length(idx) == 0) next
    df_j <- if (!is.na(mfits[[j]]$df)) mfits[[j]]$df else 10
    t_imp <- .qdist_margin(pmax(pmin(u_imp[idx, j], 1 - eps), eps),
                            margins[j], mfits[[j]]$mean[idx],
                            mfits[[j]]$sd, df = df_j)
    y[idx, j] <- pmin(pmax(t_imp, eps), lod[j] - eps)
  }
  y
}

#' Run the iterative vine imputation loop (d > 2)
#'
#' @param y n × d treatment matrix (censored cells initialised)
#' @param data internal data list
#' @param margins character vector of length d
#' @param vine_type vine type string ("auto", "C", "D")
#' @param cens_mask n × d logical censoring mask
#' @param lod numeric vector of length d
#' @param params optional user-supplied parameters list
#' @param df degrees of freedom for t margin
#' @param n_iter maximum imputation iterations
#' @param conv_tol convergence tolerance
#' @param vine_refit_every refit vine every k iterations
#' @param freeze_structure freeze vine tree structure after first fit
#' @param impute_method "gaussian_approx" or "vine_exact"
#' @param parallel_margins use future_lapply for marginal fitting
#' @param d number of treatments
#' @param verbose print progress messages
#' @return List: y, marginal_fits, vine, n_iter_done
#' @keywords internal
.run_imputation_loop_vine <- function(y, data, margins, vine_type, cens_mask,
                                       lod, params, df, n_iter, conv_tol,
                                       vine_refit_every, freeze_structure,
                                       impute_method, parallel_margins, d,
                                       verbose) {
  n_iter_done      <- 1L
  y_prev           <- NULL
  vine_obj         <- NULL
  vine_matrix_fixed <- NULL
  fam_set          <- c(1L, 3L:5L)   # Gaussian, Clayton, Gumbel, Frank

  vine_refit_every <- if (isTRUE(is.infinite(vine_refit_every))) .Machine$integer.max
                      else as.integer(vine_refit_every)

  for (iter in seq_len(n_iter)) {
    # Marginal fitting: parallel for d >= 4 when enabled
    if (is.null(params)) {
      do_par <- if (is.null(parallel_margins)) (d >= 4L) else isTRUE(parallel_margins)
      fit_fn <- function(j) {
        lod_j <- if (!is.null(lod)) lod[j] else NULL
        tryCatch(
          .fit_marginal(y[, j], data$x[[j]], margins[j], df = df, lod = lod_j),
          error = function(e) stop("Marginal ", j, " (", margins[j], "): ", e$message)
        )
      }
      marginal_fits <- if (do_par)
        future.apply::future_lapply(seq_len(d), fit_fn, future.seed = TRUE)
      else
        lapply(seq_len(d), fit_fn)
    } else {
      marginal_fits <- .process_user_params(params, y, margins, nrow(y), d)
    }

    if (iter < n_iter) {
      eps    <- .Machine$double.eps
      u_cond <- sapply(seq_len(d), function(k) {
        df_k <- if (!is.na(marginal_fits[[k]]$df)) marginal_fits[[k]]$df else 10
        pmax(pmin(.pdist_margin(y[, k], margins[k], marginal_fits[[k]]$mean,
                                 marginal_fits[[k]]$sd, df = df_k),
                  1 - eps), eps)
      })

      # Vine fitting / refitting schedule
      should_refit <- (iter == 1L) ||
                      (!freeze_structure && (iter %% vine_refit_every == 0L))

      if (should_refit && is.null(vine_matrix_fixed)) {
        # Full structure selection (first time or every k iters)
        vine_obj <- tryCatch(
          VineCopula::RVineStructureSelect(u_cond, familyset = fam_set,
                                            type = .vine_type_code(vine_type)),
          error = function(e) {
            warning("Vine structure selection failed (iter ", iter, "): ", e$message,
                    " — using Gaussian vine")
            VineCopula::RVineStructureSelect(u_cond, familyset = 1L, type = 0L)
          }
        )
        if (freeze_structure) vine_matrix_fixed <- vine_obj$Matrix

      } else if (!is.null(vine_matrix_fixed)) {
        # Re-estimate pair copula parameters with fixed tree structure
        vine_obj <- tryCatch(
          VineCopula::RVineCopSelect(u_cond, Matrix = vine_matrix_fixed,
                                      familyset = fam_set),
          error = function(e) {
            warning("Vine parameter update failed (iter ", iter, "): ", e$message)
            vine_obj  # keep previous
          }
        )
      } else if (!should_refit && !is.null(vine_obj)) {
        # Keep vine from previous iteration (periodic refitting not triggered)
        vine_obj <- vine_obj
      }

      # Imputation
      impute_fn <- if (impute_method == "vine_exact")
        .impute_censored_vine_exact
      else
        .impute_censored_vine_gaussian

      y_new <- tryCatch(
        impute_fn(y, cens_mask, lod, marginal_fits, margins, vine_obj, d),
        error = function(e) {
          warning("Imputation failed (iter ", iter, "): ", e$message,
                  " — using Gaussian approx fallback")
          .impute_censored_vine_gaussian(y, cens_mask, lod, marginal_fits,
                                          margins, vine_obj, d)
        }
      )

      # Convergence check
      if (!is.null(y_prev)) {
        lod_mat <- matrix(rep(lod, nrow(y)), nrow = nrow(y), byrow = TRUE)
        max_chg <- max(abs(y_new[cens_mask] - y_prev[cens_mask]) /
                         pmax(lod_mat[cens_mask], .Machine$double.eps),
                       na.rm = TRUE)
        if (verbose)
          message("Vine imputation iter ", iter, ": max rel. change = ",
                  formatC(max_chg, format = "e", digits = 2))
        if (max_chg < conv_tol) {
          y           <- y_new
          n_iter_done <- iter
          if (verbose) message("Converged at iteration ", iter)
          break
        }
      } else {
        if (verbose) message("Vine imputation iteration ", iter, "/", n_iter - 1L)
      }

      y           <- y_new
      y_prev      <- y
      n_iter_done <- iter
    }
  }

  list(y = y, marginal_fits = marginal_fits, vine = vine_obj,
       n_iter_done = n_iter_done)
}
