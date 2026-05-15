# =============================================================================
# drf.R — Dose-response function estimation and plotting
# =============================================================================

#' Estimate the dose-response function from a fitted ATE object
#'
#' Evaluates the GPS-weighted G-computation estimator
#' \eqn{\hat\mu(\mathbf{t}) = n^{-1}\sum_i \hat g(\mathbf{t}, \mathbf{X}_i)}
#' over a grid of treatment values using the fitted outcome model stored in
#' an `"gps_ate"` object.
#'
#' @param ate_obj Object of class `"gps_ate"` returned by [estimate_ate()].
#'   Must have been fit with `model_type` in `c("lm", "glm", "gam")`.
#'   BART objects are not supported (no closed-form prediction at new treatment
#'   values).
#' @param grid Optional.  For `vary = "marginal"`: a named list with one
#'   numeric vector per treatment, e.g.
#'   `list(T1 = seq(0, 3, length.out = 50))`.  For `vary = "joint"` (d = 2
#'   only): a named list with vectors for both treatments.  If `NULL` (default)
#'   the grid spans the 5th–95th percentile of each observed treatment.
#' @param vary `"marginal"` (default) or `"joint"`.
#'   \describe{
#'     \item{`"marginal"`}{Vary each treatment one at a time, holding all
#'       other treatments at their GPS-weighted median.  Returns one data
#'       frame per treatment.}
#'     \item{`"joint"`}{Vary all treatments simultaneously over a 2D grid
#'       (d = 2 only).  Returns a single data frame with columns
#'       `T1`, `T2`, `mu_hat`.  Standard errors are not computed for the
#'       joint grid (use marginal slices for inference).}
#'   }
#' @param n_grid Number of grid points per treatment when `grid = NULL`.
#'   Default: 50.
#' @param alpha Significance level for confidence intervals.  Default: 0.05.
#' @return Object of class `"gps_drf"`, a list with elements:
#'   \describe{
#'     \item{`slices`}{For `vary = "marginal"`: a named list of data frames,
#'       one per treatment, each with columns `treatment`, `t`, `mu_hat`,
#'       `se`, `lower`, `upper`.}
#'     \item{`surface`}{For `vary = "joint"`: a data frame with columns
#'       `T1`, `T2`, `mu_hat`.  `NULL` for marginal mode.}
#'     \item{`vary`}{The `vary` argument used.}
#'     \item{`trt_names`}{Character vector of treatment column names.}
#'     \item{`model_type`}{Outcome model type from the ATE object.}
#'   }
#' @export
drf <- function(ate_obj, grid = NULL, vary = c("marginal", "joint"),
                n_grid = 50L, alpha = 0.05) {
  if (!inherits(ate_obj, "gps_ate"))
    stop("'ate_obj' must be a \"gps_ate\" object from estimate_ate()")
  if (is.null(ate_obj$outcome_model))
    stop("No outcome model found; drf() does not support model_type = \"bart\"")
  if (identical(ate_obj$estimand, "OR"))
    stop("drf() is not applicable when estimand = \"OR\"; use estimand = \"ATE\"")

  vary <- match.arg(vary)

  gps        <- ate_obj$gps
  d          <- gps$specs$d
  treatments <- gps$data$y
  weights    <- gps$weights
  active     <- !is.na(weights)
  w_active   <- weights[active]
  covariates <- if (!is.null(ate_obj$gps$data$x)) ate_obj$gps$data$x else NULL
  trt_names  <- paste0("T", seq_len(d))
  model      <- ate_obj$outcome_model
  z_crit     <- stats::qnorm(1 - alpha / 2)

  if (vary == "joint" && d != 2L)
    stop("vary = \"joint\" is only available for d = 2 treatments")

  # GPS-weighted median for each treatment (used to hold fixed in marginal mode)
  trt_medians <- apply(treatments[active, , drop = FALSE], 2, function(t)
    stats::weighted.mean(t, w_active))

  # Build a prediction data frame at a fixed treatment value `treat`
  make_pred_df_full <- function(treat) {
    df_t <- as.data.frame(matrix(rep(treat, nrow(treatments)),
                                 nrow = nrow(treatments), byrow = TRUE))
    colnames(df_t) <- trt_names
    if (!is.null(covariates)) {
      cov_df <- as.data.frame(covariates)
      cov_names <- colnames(covariates) %||% paste0("covar_", seq_len(ncol(covariates)))
      colnames(cov_df) <- cov_names
      df_t <- cbind(df_t, cov_df)
    }
    df_t
  }

  # G-computation at a single treatment vector: average predicted outcome over X
  gcomp_at <- function(treat) {
    pd   <- make_pred_df_full(treat)
    pred <- predict(model, newdata = pd, type = "response")
    stats::weighted.mean(pred[active], w_active)
  }

  # Delta-method SE for a single treatment vector (lm/glm only)
  gcomp_se <- function(treat) {
    if (is.null(ate_obj$vcov_matrix)) return(NA_real_)
    pd   <- make_pred_df_full(treat)
    coef_est <- coef(model)
    keep     <- !is.na(coef_est)
    logistic <- identical(ate_obj$model_type, "glm") &&
                !is.null(ate_obj$family) &&
                identical(ate_obj$family$family, "binomial")
    mm <- tryCatch(
      stats::model.matrix(model, data = pd)[active, keep, drop = FALSE],
      error = function(e) NULL
    )
    if (is.null(mm)) return(NA_real_)
    if (logistic) {
      eta  <- as.numeric(mm %*% coef_est[keep])
      mu_p <- exp(eta) / (1 + exp(eta))^2
      grad <- apply(mu_p * mm, 2, stats::weighted.mean, w = w_active)
    } else {
      grad <- apply(mm, 2, stats::weighted.mean, w = w_active)
    }
    V <- ate_obj$vcov_matrix[names(grad), names(grad), drop = FALSE]
    as.numeric(sqrt(max(0, t(grad) %*% V %*% grad)))
  }

  # ---- Build auto grid --------------------------------------------------------
  build_auto_grid <- function(j) {
    x <- treatments[, j]
    seq(stats::quantile(x, 0.05), stats::quantile(x, 0.95), length.out = n_grid)
  }

  # ---- Marginal mode ----------------------------------------------------------
  if (vary == "marginal") {
    slices <- vector("list", d)
    names(slices) <- trt_names

    for (j in seq_len(d)) {
      t_grid <- if (!is.null(grid[[trt_names[j]]])) grid[[trt_names[j]]]
                else build_auto_grid(j)

      mu_hat <- numeric(length(t_grid))
      se_hat <- numeric(length(t_grid))

      for (k in seq_along(t_grid)) {
        treat      <- trt_medians
        treat[j]   <- t_grid[k]
        mu_hat[k]  <- gcomp_at(treat)
        se_hat[k]  <- gcomp_se(treat)
      }

      slices[[j]] <- data.frame(
        treatment = trt_names[j],
        t         = t_grid,
        mu_hat    = mu_hat,
        se        = se_hat,
        lower     = mu_hat - z_crit * se_hat,
        upper     = mu_hat + z_crit * se_hat
      )
    }

    return(structure(
      list(slices = slices, surface = NULL, vary = "marginal",
           trt_names = trt_names, model_type = ate_obj$model_type,
           alpha = alpha),
      class = "gps_drf"
    ))
  }

  # ---- Joint mode (d = 2) -----------------------------------------------------
  t1_grid <- if (!is.null(grid[[trt_names[1L]]])) grid[[trt_names[1L]]]
             else build_auto_grid(1L)
  t2_grid <- if (!is.null(grid[[trt_names[2L]]])) grid[[trt_names[2L]]]
             else build_auto_grid(2L)

  grid_df <- expand.grid(T1 = t1_grid, T2 = t2_grid)
  colnames(grid_df) <- trt_names

  mu_hat <- apply(grid_df, 1, function(treat) gcomp_at(as.numeric(treat)))
  grid_df$mu_hat <- mu_hat

  structure(
    list(slices = NULL, surface = grid_df, vary = "joint",
         trt_names = trt_names, model_type = ate_obj$model_type,
         alpha = alpha),
    class = "gps_drf"
  )
}
