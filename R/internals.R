# =============================================================================
# internals.R — core distribution functions, marginal fitting, helpers
# =============================================================================

# ---- Null-coalescing operator ------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Human-readable copula family name --------------------------------------
.family_name <- function(code) {
  lut <- c("0" = "Independence", "1" = "Gaussian", "2" = "Student-t",
           "3" = "Clayton", "4" = "Gumbel", "5" = "Frank")
  nm  <- lut[as.character(code)]
  if (is.na(nm)) paste("Family", code) else unname(nm)
}

# ---- Student-t df MLE -------------------------------------------------------
.estimate_t_df <- function(residuals) {
  scaled <- residuals / max(stats::mad(residuals), .Machine$double.eps)
  tryCatch({
    fit <- MASS::fitdistr(scaled, densfun = "t")
    max(3, fit$estimate[["df"]])
  }, error = function(e) {
    warning("t-df MLE failed (", conditionMessage(e), "). Using df = 10.")
    10
  })
}

# ---- Marginal distribution functions ----------------------------------------

#' @keywords internal
.ddist_margin <- function(y, margin, mean, sd, df = 10) {
  switch(margin,
    normal    = stats::dnorm(y, mean = mean, sd = sd),
    lognormal = stats::dlnorm(y, meanlog = mean, sdlog = sd),
    t         = stats::dt((y - mean) / sd, df = df) / sd,
    gamma     = stats::dgamma(y, shape = (mean / sd)^2, rate = mean / sd^2),
    invgauss  = statmod::dinvgauss(y, mean = mean, shape = 1 / (sd^2))
  )
}

#' @keywords internal
.pdist_margin <- function(y, margin, mean, sd, df = 10) {
  switch(margin,
    normal    = stats::pnorm(y, mean = mean, sd = sd),
    lognormal = stats::plnorm(y, meanlog = mean, sdlog = sd),
    t         = stats::pt((y - mean) / sd, df = df),
    gamma     = stats::pgamma(y, shape = (mean / sd)^2, rate = mean / sd^2),
    invgauss  = statmod::pinvgauss(y, mean = mean, shape = 1 / (sd^2))
  )
}

#' @keywords internal
.qdist_margin <- function(p, margin, mean, sd, df = 10) {
  switch(margin,
    normal    = stats::qnorm(p, mean = mean, sd = sd),
    lognormal = stats::qlnorm(p, meanlog = mean, sdlog = sd),
    t         = mean + sd * stats::qt(p, df = df),
    gamma     = {
      shape <- (mean / sd)^2; rate <- mean / sd^2
      stats::qgamma(p, shape = shape, rate = rate)
    },
    invgauss  = {
      disp <- sd^2 / mean^3
      statmod::qinvgauss(p, mean = mean, dispersion = disp)
    }
  )
}

# ---- Marginal regression model fitting --------------------------------------

#' Fit a marginal regression model for one treatment dimension
#'
#' For `normal` and `lognormal` margins with LOD supplied, uses Tobit MLE
#' (`AER::tobit`) which is censoring-aware and gives consistent estimates even
#' at high censoring rates.
#'
#' @param y Numeric response vector (length n)
#' @param x Design matrix (no intercept column)
#' @param margin Distribution family: `"normal"`, `"lognormal"`, `"t"`,
#'   `"gamma"`, `"invgauss"`
#' @param df Degrees of freedom for t margin. `NULL` estimates via MLE.
#' @param lod Left-censoring threshold. When non-`NULL` and margin is
#'   `"normal"` or `"lognormal"`, activates Tobit MLE.
#' @return Named list: mean, sd, df, model, coefficients, vcov
#' @keywords internal
.fit_marginal <- function(y, x, margin, df = NULL, lod = NULL) {
  if (margin %in% c("gamma", "invgauss") && any(y <= 0, na.rm = TRUE))
    stop(margin, " requires all y > 0")
  if (margin == "lognormal" && any(y <= 0, na.rm = TRUE))
    warning("lognormal requires y > 0; values <= 0 will be clamped")

  use_tobit <- !is.null(lod) && margin %in% c("normal", "lognormal") &&
               requireNamespace("AER", quietly = TRUE)

  .strip_scale <- function(fit) {
    cf  <- coef(fit);  cf  <- cf[names(cf)  != "Log(scale)"]
    vc  <- vcov(fit);  keep <- rownames(vc) != "Log(scale)"
    vc  <- vc[keep, keep, drop = FALSE]
    list(coef = cf, vcov = vc)
  }

  switch(margin,

    normal = {
      if (use_tobit) {
        fit <- tryCatch(
          AER::tobit(y ~ x - 1, left = lod),
          error = function(e) {
            warning("Tobit failed for normal margin (", conditionMessage(e),
                    "); falling back to OLS.")
            NULL
          }
        )
        if (!is.null(fit)) {
          cs <- .strip_scale(fit)
          return(list(mean = predict(fit, type = "linear"), sd = fit$scale,
                      df = NA_real_, model = fit,
                      coefficients = cs$coef, vcov = cs$vcov))
        }
      }
      fit <- stats::lm(y ~ x - 1)
      list(mean = predict(fit), sd = stats::sigma(fit), df = NA_real_,
           model = fit, coefficients = coef(fit), vcov = vcov(fit))
    },

    lognormal = {
      y_log <- log(pmax(y, .Machine$double.eps))
      if (use_tobit) {
        lod_log <- log(max(lod, .Machine$double.eps))
        fit <- tryCatch(
          AER::tobit(y_log ~ x - 1, left = lod_log),
          error = function(e) {
            warning("Tobit failed for lognormal margin (", conditionMessage(e),
                    "); falling back to OLS.")
            NULL
          }
        )
        if (!is.null(fit)) {
          cs <- .strip_scale(fit)
          return(list(mean = predict(fit, type = "linear"), sd = fit$scale,
                      df = NA_real_, model = fit,
                      coefficients = cs$coef, vcov = cs$vcov))
        }
      }
      fit <- stats::lm(y_log ~ x - 1)
      list(mean = predict(fit), sd = stats::sigma(fit), df = NA_real_,
           model = fit, coefficients = coef(fit), vcov = vcov(fit))
    },

    t = {
      fit    <- stats::lm(y ~ x - 1)
      resids <- stats::residuals(fit)
      df_est <- if (is.null(df)) .estimate_t_df(resids) else max(3, as.numeric(df))
      list(mean = predict(fit), sd = stats::mad(resids), df = df_est,
           model = fit, coefficients = coef(fit), vcov = vcov(fit))
    },

    gamma = {
      fit <- stats::glm(y ~ x - 1, family = stats::Gamma(link = "log"))
      mu  <- exp(predict(fit))
      phi <- summary(fit)$dispersion
      list(mean = mu, sd = sqrt(mu^2 * phi), df = NA_real_,
           model = fit, coefficients = coef(fit), vcov = vcov(fit),
           dispersion = phi)
    },

    invgauss = {
      fit    <- suppressWarnings(
        stats::glm(y ~ x - 1, family = statmod::inverse.gaussian(link = "log"))
      )
      mu     <- exp(predict(fit))
      lambda <- 1 / summary(fit)$dispersion
      list(mean = mu, sd = sqrt(mu^3 / lambda), df = NA_real_,
           model = fit, coefficients = coef(fit), vcov = vcov(fit),
           dispersion = summary(fit)$dispersion)
    }
  )
}

# ---- User-supplied parameter validation -------------------------------------
#' @keywords internal
.process_user_params <- function(params_user, y, margins, n, d) {
  if (!all(c("mean", "sd") %in% names(params_user)))
    stop("params must contain 'mean' and 'sd' components")
  if (length(params_user$mean) != d) stop("params$mean must have length d = ", d)
  if (length(params_user$sd)   != d) stop("params$sd must have length d = ", d)

  lapply(seq_len(d), function(j) {
    expand <- function(v, nm) {
      if (length(v) == 1) return(rep(v, n))
      if (length(v) != n) stop(nm, "[[", j, "]] must have length 1 or n (", n, ")")
      v
    }
    mean_j <- expand(params_user$mean[[j]], "params$mean")
    sd_j   <- expand(params_user$sd[[j]],   "params$sd")

    m <- margins[j]
    if (m %in% c("gamma", "invgauss", "lognormal") && any(sd_j   <= 0))
      stop("sd must be positive for ", m)
    if (m %in% c("gamma", "invgauss") && any(mean_j <= 0))
      stop("mean must be positive for ", m)

    list(mean = mean_j, sd = sd_j, df = NA_real_,
         coefficients = NA, vcov = NA, model = NULL)
  })
}
