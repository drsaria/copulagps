# =============================================================================
# copula_utils.R — copula family helpers, auto-adjustment, convex hull trimming
# =============================================================================

#' Map a copula type string to VineCopula integer family code(s)
#' @keywords internal
.get_copula_family <- function(copula_type) {
  switch(copula_type,
    clayton = 3L,
    gumbel  = 4L,
    frank   = 5L,
    normal  = 1L,
    t       = 2L,
    c(1L, 3L:5L)   # "auto" or unrecognised: let BiCopSelect choose
  )
}

#' Auto-detect whether to use residual copula based on R² and Kendall-tau shift
#'
#' Returns `TRUE` if any marginal R² > 0.3 AND the Kendall-tau shift between
#' the raw and residual copulas exceeds 0.05.
#'
#' @param y n x d treatment matrix
#' @param marginal_fits list of marginal fit objects
#' @param copula_type copula type string
#' @param d number of treatments
#' @param verbose print message if residual copula selected
#' @return logical scalar
#' @keywords internal
.should_use_residual_copula <- function(y, marginal_fits, copula_type, d, verbose) {
  marginal_r2 <- sapply(seq_len(d), function(j) {
    fitted_mean <- marginal_fits[[j]]$mean
    if (is.null(fitted_mean)) return(NA_real_)
    y_ref  <- y[, j]
    resids  <- y_ref - fitted_mean
    ss_res  <- sum(resids^2, na.rm = TRUE)
    ss_tot  <- sum((y_ref - mean(y_ref, na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot < .Machine$double.eps) return(NA_real_)
    max(0, 1 - ss_res / ss_tot)
  })

  if (!any(marginal_r2 > 0.3, na.rm = TRUE))
    return(list(use = FALSE, r2 = marginal_r2))

  fam_set   <- .get_copula_family(copula_type)
  u_raw     <- VineCopula::pobs(y)
  raw_cop   <- VineCopula::BiCopSelect(u_raw[, 1], u_raw[, 2], familyset = fam_set)

  resid_mat <- sapply(seq_len(d), function(j) y[, j] - marginal_fits[[j]]$mean)
  u_resid   <- VineCopula::pobs(resid_mat)
  resid_cop <- VineCopula::BiCopSelect(u_resid[, 1], u_resid[, 2], familyset = fam_set)

  raw_tau   <- VineCopula::BiCopPar2Tau(raw_cop$family,   raw_cop$par,   raw_cop$par2)
  resid_tau <- VineCopula::BiCopPar2Tau(resid_cop$family, resid_cop$par, resid_cop$par2)
  delta_tau <- abs(raw_tau - resid_tau)

  use <- delta_tau > 0.05
  if (use && verbose)
    message("Auto-adjustment: residual copula selected (tau shift = ",
            round(delta_tau, 3), ", max R\u00b2 = ",
            round(max(marginal_r2, na.rm = TRUE), 3), ")")

  list(use = use, r2 = marginal_r2)
}

# =============================================================================
# Convex hull overlap trimming (d = 2)
# =============================================================================

#' Test whether 2-D points lie inside a convex hull
#'
#' Uses the cross-product half-plane test, which works for both CW and CCW
#' vertex orderings returned by `grDevices::chull()`.  All hull vertices are
#' included (boundary points return `TRUE`).
#'
#' @param pts n x 2 matrix of query points
#' @param hull_pts k x 2 matrix of hull vertices **in the order returned by
#'   `grDevices::chull()`**
#' @return Logical vector of length n; `TRUE` = inside or on boundary
#' @keywords internal
.in_convex_hull_2d <- function(pts, hull_pts) {
  k <- nrow(hull_pts)
  n <- nrow(pts)
  # Cross-product for each edge; columns = edges, rows = query points
  cp <- matrix(NA_real_, n, k)
  for (i in seq_len(k)) {
    j  <- if (i < k) i + 1L else 1L
    dx <- hull_pts[j, 1] - hull_pts[i, 1]
    dy <- hull_pts[j, 2] - hull_pts[i, 2]
    cp[, i] <- dx * (pts[, 2] - hull_pts[i, 2]) -
                dy * (pts[, 1] - hull_pts[i, 1])
  }
  # Inside iff all cross products >= 0 (CCW hull) or all <= 0 (CW hull)
  rowSums(cp >= 0) == k | rowSums(cp <= 0) == k
}

#' Identify observations inside the convex hull of the joint treatment distribution
#'
#' For d = 2, the convex hull is computed on the observed treatment values and
#' each observation is tested for membership.  Observations outside the hull
#' are in a region of the joint treatment space with no support in the data —
#' their GPS weights would extrapolate beyond any observed unit.
#'
#' For d > 2 use `check_hull_overlap()` with `method = "mahalanobis"` instead;
#' convex hulls in high dimensions become geometrically loose and expensive to
#' compute.
#'
#' @param y n x 2 (or n x d) treatment matrix
#' @return Logical vector of length n: `TRUE` = inside hull (keep),
#'   `FALSE` = outside (candidate for trimming)
#' @keywords internal
.check_convex_hull <- function(y) {
  if (ncol(y) != 2)
    stop(".check_convex_hull() supports only d = 2. Use Mahalanobis trimming for d > 2.")

  hull_idx <- grDevices::chull(y[, 1], y[, 2])
  hull_pts  <- y[hull_idx, , drop = FALSE]
  .in_convex_hull_2d(y, hull_pts)
}

#' Check and report convex hull overlap for a fitted GPS object (d = 2 only)
#'
#' Deprecated in favour of [check_overlap()], which works for any d.
#'
#' @param gps Object of class `"gps_weights"`
#' @return Invisible logical vector of length n (`TRUE` = inside hull)
#' @export
check_hull_overlap <- function(gps) {
  check_overlap(gps, method = "hull")
}

# =============================================================================
# Mahalanobis distance trimming (any d)
# =============================================================================

#' Flag observations outside the joint treatment distribution via Mahalanobis distance
#'
#' Computes the squared Mahalanobis distance of each observation's treatment
#' vector from the sample centroid.  Observations exceeding the chi-squared
#' quantile at `df = d` degrees of freedom are candidates for trimming.
#'
#' This generalises convex hull trimming to d > 2: the ellipsoidal decision
#' boundary adapts to the correlation structure of the treatments and has a
#' direct probabilistic interpretation (under multivariate normality,
#' `quantile = 0.975` retains 97.5% of the marginal distribution).
#'
#' @param y n × d treatment matrix
#' @param quantile Chi-squared quantile.  Default: `0.975`.
#' @return Logical vector of length n; `TRUE` = inside ellipsoid (keep)
#' @keywords internal
.check_mahalanobis_trim <- function(y, quantile = 0.975) {
  d      <- ncol(y)
  mu     <- colMeans(y)
  S      <- stats::cov(y)
  # Use pseudoinverse-safe solve via tryCatch
  S_inv  <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  md_sq  <- mahalanobis(y, center = mu, cov = S_inv, inverted = TRUE)
  cutoff <- stats::qchisq(quantile, df = d)
  md_sq <= cutoff
}

# =============================================================================
# Unified overlap check (any d)
# =============================================================================

#' Check treatment overlap for a fitted GPS object
#'
#' Identifies observations in sparse regions of the joint treatment distribution.
#' For d = 2, defaults to the convex hull method; for d > 2, defaults to
#' Mahalanobis distance.  Use `method = "mahalanobis"` to override for d = 2.
#'
#' @param gps Object of class `"gps_weights"`
#' @param method `"auto"` (default), `"hull"` (d = 2 only), or `"mahalanobis"`
#' @param mahal_quantile Chi-squared quantile for Mahalanobis method.
#'   Default: `0.975`.
#' @return Invisible logical vector of length n (`TRUE` = inside overlap region)
#' @export
check_overlap <- function(gps, method = "auto", mahal_quantile = 0.975) {
  y <- gps$data$y
  d <- ncol(y)
  n <- nrow(y)

  if (method == "auto") method <- if (d == 2) "hull" else "mahalanobis"

  if (method == "hull" && d != 2)
    stop("method = 'hull' requires d = 2; use method = 'mahalanobis' for d > 2")

  mask <- if (method == "hull")
    .check_convex_hull(y)
  else
    .check_mahalanobis_trim(y, quantile = mahal_quantile)

  n_out <- sum(!mask)
  cat(sprintf(
    "Overlap check (%s)  |  n = %d  |  outside: %d (%.1f%%)\n",
    method, n, n_out, 100 * n_out / n
  ))
  if (n_out > 0)
    cat(sprintf("  Tip: refit with trim_method = '%s' to exclude these observations.\n",
                method))
  invisible(mask)
}

#' Map vine type string to VineCopula integer type code
#' @keywords internal
.vine_type_code <- function(vine_type) {
  switch(vine_type, "auto" = 0L, "C" = 1L, "D" = 2L, 0L)
}
