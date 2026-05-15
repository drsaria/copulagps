# =============================================================================
# simulate.R — data simulation utilities
# =============================================================================

#' Simulate bivariate treatment data with a known copula structure
#'
#' Generates n observations from a bivariate copula with Normal marginals,
#' optionally with covariate effects on treatment means.
#'
#' @param n Sample size
#' @param copula_family Copula family code: 1 = Gaussian, 3 = Clayton,
#'   4 = Gumbel, 5 = Frank.  Default: `3` (Clayton).
#' @param copula_par Copula dependence parameter.  Default: `2`.
#' @param marginal_params List of 2 sublists, each with `$mean` and `$sd`
#'   (scalar).  Default: `list(list(mean=2, sd=1), list(mean=-1, sd=1.5))`.
#' @param covariates Optional n × p covariate matrix.
#' @param treatment_effects Optional list of 2 numeric vectors (length p each)
#'   giving covariate-to-treatment effects on the mean.
#' @param lod Optional numeric vector of length 2 for left-censoring thresholds.
#'   When supplied, treatment values below the LOD are set to the LOD.
#' @return List with `$y` (n × 2 treatment matrix), `$x` (list of design
#'   matrices), `$u_true`, `$covariates`, `$true_copula`, `$marginal_params`,
#'   and optionally `$lod`.
#' @export
simulate_treatments <- function(
    n,
    copula_family   = 3L,
    copula_par      = 2,
    marginal_params = list(list(mean = 2, sd = 1), list(mean = -1, sd = 1.5)),
    covariates      = NULL,
    treatment_effects = NULL,
    lod             = NULL
) {
  if (length(marginal_params) != 2)
    stop("marginal_params must have length 2 (d = 2 only in this version)")

  u_data     <- VineCopula::BiCopSim(n, family = copula_family, par = copula_par)
  treatments <- matrix(0, n, 2L)

  for (j in seq_len(2L)) {
    p  <- marginal_params[[j]]
    mu <- p$mean
    if (!is.null(covariates) && !is.null(treatment_effects) &&
        length(treatment_effects) >= j) {
      eff <- treatment_effects[[j]]
      if (length(eff) == ncol(covariates))
        mu <- mu + as.numeric(covariates %*% eff)
    }
    treatments[, j] <- stats::qnorm(u_data[, j], mean = mu, sd = p$sd)
  }

  x_matrices <- if (is.null(covariates))
    lapply(seq_len(2L), function(j) matrix(1, n, 1L))
  else
    lapply(seq_len(2L), function(j) cbind(1, covariates))

  out <- list(
    y               = treatments,
    x               = x_matrices,
    u_true          = u_data,
    covariates      = covariates,
    true_copula     = list(family = copula_family, par = copula_par),
    marginal_params = marginal_params
  )

  if (!is.null(lod)) {
    lod <- rep_len(as.numeric(lod), 2L)
    for (j in seq_len(2L))
      out$y[out$y[, j] < lod[j], j] <- lod[j]
    out$lod <- lod
  }

  out
}
