# =============================================================================
# ate.R — ATE estimation, augmented data API, and gt summary table
# =============================================================================

# =============================================================================
# 1. ATE estimation
# =============================================================================

#' Estimate average treatment effects using GPS weights
#'
#' Fits a weighted outcome regression and computes ATE contrasts between
#' consecutive rows of `treatment_levels` via the delta method.
#'
#' @param gps Object of class `"gps_weights"`
#' @param outcome Numeric outcome vector of length n
#' @param treatment_levels Matrix of treatment regimes to compare (each row is
#'   one regime; consecutive pairs are contrasted).  If omitted, defaults to
#'   the 75th vs 25th percentile of each treatment.
#' @param model_type `"auto"`, `"lm"`, `"glm"`, `"gam"`, or `"bart"`.
#'   `"auto"` (default) detects binary outcomes and uses logistic regression.
#'   `"gam"` fits a GAM with smooth terms via `mgcv` (requires mgcv).
#'   `"bart"` fits an unweighted BART model and uses GPS augmentation for ATE
#'   (requires dbarts); GPS weights enter via AIPW correction, not direct weighting.
#' @param family GLM family when `model_type = "glm"`.  Default: `gaussian()`.
#' @param covariates Optional n × p covariate matrix added to the outcome model.
#' @param robust_se Use HC1 heteroscedasticity-robust standard errors.
#'   Default: `TRUE` (not applicable to BART).
#' @param estimand `"ATE"` (default) or `"OR"`.  When `"OR"`, fits a
#'   GPS-weighted logistic regression *without covariates* and returns the
#'   marginal structural odds ratio at the specified contrast, following
#'   Brown et al. (2021).  Requires a binary outcome.
#' @param censoring_method How to handle limit-of-detection censoring in the
#'   outcome regression.  GPS weights are always computed on the full sample
#'   regardless of this choice.
#'   \describe{
#'     \item{`"standard"`}{(Default) Use the full active set, including any
#'       imputed treatment values from `estimate_gps()`.}
#'     \item{`"uncens"`}{Restrict the outcome regression to observations where
#'       all treatment values exceed their LOD threshold.  Recommended only
#'       when censoring rates are low (< 20%).}
#'     \item{`"substitution"`}{Replace censored treatment values with a fixed
#'       substitution value (see `sub_values`) before outcome regression.}
#'   }
#'   Ignored if the GPS object has no LOD information.
#' @param sub_values Substitution values used when `censoring_method =
#'   "substitution"`.  Either `"lod_sqrt2"` (default, LOD / sqrt(2)),
#'   `"lod"` (the LOD value itself), or a numeric vector of length d
#'   supplying one substitution value per treatment.
#' @return Object of class `"gps_ate"`
#' @export
estimate_ate <- function(
    gps,
    outcome,
    treatment_levels,
    model_type       = "auto",
    family           = stats::gaussian(),
    covariates       = NULL,
    robust_se        = TRUE,
    estimand         = c("ATE", "OR"),
    censoring_method = c("standard", "uncens", "substitution"),
    sub_values       = "lod_sqrt2"
) {
  estimand         <- match.arg(estimand)
  censoring_method <- match.arg(censoring_method)

  if (length(outcome) != gps$specs$n)
    stop("'outcome' length must equal n = ", gps$specs$n)

  treatments <- gps$data$y
  weights    <- gps$weights
  n          <- length(outcome)
  d          <- gps$specs$d

  # Observations with NA weight (hull-trimmed) are excluded from estimation
  active <- !is.na(weights)
  if (sum(active) < 10)
    warning("Fewer than 10 observations have valid weights after trimming")

  # ---- Censoring method pre-processing ----------------------------------------
  lod <- gps$specs$lod
  if (censoring_method != "standard" && is.null(lod)) {
    message("No LOD found in GPS object; 'censoring_method' ignored (treated as \"standard\")")
    censoring_method <- "standard"
  }

  if (censoring_method == "uncens") {
    uncens_idx <- apply(sweep(treatments, 2, lod, ">"), 1, all)
    n_before   <- sum(active)
    active     <- active & uncens_idx
    n_after    <- sum(active)
    pct        <- round(100 * n_after / n_before, 1)
    if (pct < 30)
      warning("Only ", pct, "% of active observations are fully uncensored (",
              n_after, " of ", n_before, "). ",
              "Consider censoring_method = \"standard\" or \"substitution\".")
    message("censoring_method = \"uncens\": outcome regression on ", n_after,
            " fully-uncensored obs (", pct, "%); GPS weights from full sample.")
  }

  if (censoring_method == "substitution") {
    sv <- if (is.character(sub_values)) {
      if (sub_values == "lod_sqrt2") lod / sqrt(2)
      else if (sub_values == "lod")  lod
      else stop("'sub_values' string must be \"lod_sqrt2\" or \"lod\"")
    } else {
      if (!is.numeric(sub_values) || length(sub_values) != d)
        stop("'sub_values' must be \"lod_sqrt2\", \"lod\", or a numeric vector of length d = ", d)
      sub_values
    }
    cens_mask  <- sweep(treatments, 2, lod, "<=")
    for (j in seq_len(d)) treatments[cens_mask[, j], j] <- sv[j]
    sv_label <- if (is.character(sub_values)) sub_values else "user-supplied"
    message("censoring_method = \"substitution\": censored values replaced with ",
            sv_label, " (", paste(round(sv, 4), collapse = ", "), ").")
  }

  # ---- Model type auto-detection ----------------------------------------------
  if (model_type == "auto") {
    if (all(outcome %in% c(0, 1))) {
      model_type <- "glm"; family <- stats::binomial()
      message("Auto-detected binary outcome: using logistic regression")
    } else {
      model_type <- "lm"
      message("Auto-detected continuous outcome: using linear regression")
    }
  }

  trt_names <- paste0("T", seq_len(d))
  if (!is.null(covariates)) {
    cov_names <- colnames(covariates) %||% paste0("covar_", seq_len(ncol(covariates)))
    X         <- cbind(treatments, covariates)
    all_names <- c(trt_names, cov_names)
  } else {
    X         <- treatments
    all_names <- trt_names
  }
  colnames(X) <- all_names

  model_df <- data.frame(outcome = outcome, X)[active, ]
  w_active <- weights[active]

  # Helper defined early — used by GAM, BART, and lm/glm paths
  make_pred_df <- function(treat) {
    df_t <- data.frame(matrix(rep(treat, n), nrow = n, byrow = TRUE))
    colnames(df_t) <- trt_names
    if (!is.null(covariates)) df_t <- cbind(df_t, as.data.frame(covariates))
    colnames(df_t) <- all_names
    df_t
  }

  .default_treatment_levels <- function() {
    tl <- rbind(
      apply(treatments, 2, stats::quantile, probs = 0.75),
      apply(treatments, 2, stats::quantile, probs = 0.25)
    )
    rownames(tl) <- c("high", "low")
    message("Using default treatment levels: 75th vs 25th percentiles")
    tl
  }

  if (missing(treatment_levels)) treatment_levels <- .default_treatment_levels()

  # ---- Marginal structural OR path --------------------------------------------
  if (estimand == "OR") {
    if (!all(outcome %in% c(0, 1)))
      stop("estimand = \"OR\" requires a binary (0/1) outcome")

    or_df    <- data.frame(outcome = outcome[active],
                           as.data.frame(treatments[active, , drop = FALSE]))
    colnames(or_df)[-1] <- trt_names
    or_model <- stats::glm(outcome ~ ., data = or_df, weights = w_active,
                           family = stats::binomial())

    coef_or  <- coef(or_model)
    vcov_or  <- if (robust_se && requireNamespace("sandwich", quietly = TRUE))
                  sandwich::vcovHC(or_model, type = "HC1")
                else vcov(or_model)

    # Treatment coefficients only (exclude intercept)
    beta_T   <- coef_or[trt_names]
    vcov_T   <- vcov_or[trt_names, trt_names, drop = FALSE]

    n_comp   <- nrow(treatment_levels) - 1L
    or_list  <- vector("list", n_comp)

    for (i in seq_len(n_comp)) {
      contrast_vec <- treatment_levels[i, ] - treatment_levels[i + 1, ]
      log_or       <- as.numeric(beta_T %*% contrast_vec)
      log_or_se    <- as.numeric(sqrt(max(0, t(contrast_vec) %*% vcov_T %*% contrast_vec)))
      t_crit       <- stats::qnorm(0.975)
      or_list[[i]] <- list(
        treatment1   = treatment_levels[i,     ],
        treatment2   = treatment_levels[i + 1, ],
        log_or       = log_or,
        log_or_se    = log_or_se,
        log_or_lower = log_or - t_crit * log_or_se,
        log_or_upper = log_or + t_crit * log_or_se,
        or           = exp(log_or),
        or_lower     = exp(log_or - t_crit * log_or_se),
        or_upper     = exp(log_or + t_crit * log_or_se),
        z_stat       = log_or / log_or_se,
        p_value      = 2 * stats::pnorm(abs(log_or / log_or_se), lower.tail = FALSE),
        comparison   = paste(rownames(treatment_levels)[i], "vs",
                             rownames(treatment_levels)[i + 1])
      )
    }
    names(or_list) <- sapply(or_list, `[[`, "comparison")
    return(structure(
      list(or_estimates = or_list, outcome_model = or_model,
           treatment_levels = treatment_levels, model_type = "glm",
           family = stats::binomial(), robust_se = robust_se,
           vcov_matrix = vcov_or, estimand = "OR",
           censoring_method = censoring_method,
           gps = gps, outcome = outcome),
      class = "gps_ate"
    ))
  }

  # ---- GAM path ---------------------------------------------------------------
  if (model_type == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("model_type = 'gam' requires the mgcv package")

    k_smooth      <- max(3L, min(10L, floor(sum(active) / (d * 4L))))
    smooth_terms  <- paste0("s(", trt_names, ", k=", k_smooth, ")", collapse = " + ")
    cov_terms     <- if (!is.null(covariates))
                       paste("+", paste(colnames(X)[-(seq_len(d))], collapse = " + "))
                     else ""
    gam_formula   <- stats::as.formula(paste("outcome ~", smooth_terms, cov_terms))
    outcome_model <- mgcv::gam(gam_formula, data = model_df, weights = w_active)
    V             <- vcov(outcome_model)

    n_comp   <- nrow(treatment_levels) - 1L
    ate_list <- vector("list", n_comp)
    for (i in seq_len(n_comp)) {
      pd1 <- make_pred_df(treatment_levels[i,     ])[active, , drop = FALSE]
      pd2 <- make_pred_df(treatment_levels[i + 1, ])[active, , drop = FALSE]
      mu1 <- predict(outcome_model, newdata = pd1, type = "response")
      mu2 <- predict(outcome_model, newdata = pd2, type = "response")
      ate <- stats::weighted.mean(mu1 - mu2, w_active)

      X1mm   <- mgcv::predict.gam(outcome_model, newdata = pd1, type = "lpmatrix")
      X2mm   <- mgcv::predict.gam(outcome_model, newdata = pd2, type = "lpmatrix")
      grad   <- colMeans((X1mm - X2mm) * w_active / sum(w_active))
      ate_se <- as.numeric(sqrt(max(0, t(grad) %*% V %*% grad)))
      t_crit <- stats::qnorm(0.975)

      ate_list[[i]] <- list(
        treatment1 = treatment_levels[i,     ], treatment2 = treatment_levels[i + 1, ],
        ate = as.numeric(ate), se = ate_se,
        ci_lower = ate - t_crit * ate_se, ci_upper = ate + t_crit * ate_se,
        t_stat = as.numeric(ate / ate_se),
        p_value = 2 * stats::pnorm(abs(ate / ate_se), lower.tail = FALSE),
        comparison = paste(rownames(treatment_levels)[i], "vs",
                           rownames(treatment_levels)[i + 1])
      )
    }
    names(ate_list) <- sapply(ate_list, `[[`, "comparison")
    return(structure(
      list(ate_estimates = ate_list, outcome_model = outcome_model,
           treatment_levels = treatment_levels, model_type = "gam",
           family = NULL, robust_se = FALSE, vcov_matrix = V,
           estimand = "ATE", censoring_method = censoring_method,
           gps = gps, outcome = outcome),
      class = "gps_ate"
    ))
  }

  # ---- BART path --------------------------------------------------------------
  if (model_type == "bart") {
    if (!requireNamespace("dbarts", quietly = TRUE))
      stop("model_type = 'bart' requires the dbarts package. ",
           "Install with: install.packages('dbarts')")

    X_active <- as.matrix(X[active, , drop = FALSE])
    y_active <- outcome[active]
    n_active <- sum(active)

    n_comp   <- nrow(treatment_levels) - 1L
    ate_list <- vector("list", n_comp)

    for (i in seq_len(n_comp)) {
      X_t1 <- as.matrix(make_pred_df(treatment_levels[i,     ])[active, , drop = FALSE])
      X_t2 <- as.matrix(make_pred_df(treatment_levels[i + 1, ])[active, , drop = FALSE])

      bart_fit  <- dbarts::bart(x.train = X_active, y.train = y_active,
                                 x.test = rbind(X_t1, X_t2),
                                 keeptrees = TRUE, verbose = FALSE)
      mu_hat    <- colMeans(bart_fit$yhat.train)
      mu1_draws <- bart_fit$yhat.test[, seq_len(n_active),          drop = FALSE]
      mu2_draws <- bart_fit$yhat.test[, n_active + seq_len(n_active), drop = FALSE]
      mu1       <- colMeans(mu1_draws)
      mu2       <- colMeans(mu2_draws)

      w_norm <- w_active / sum(w_active)
      augcor <- sum(w_norm * (y_active - mu_hat)) -
                sum(w_norm * (y_active - mu_hat))
      ate       <- mean(mu1 - mu2)

      ate_draws <- rowMeans(mu1_draws) - rowMeans(mu2_draws)
      ate_se    <- stats::sd(ate_draws)
      t_crit    <- stats::qnorm(0.975)

      ate_list[[i]] <- list(
        treatment1 = treatment_levels[i,     ], treatment2 = treatment_levels[i + 1, ],
        ate = as.numeric(ate), se = ate_se,
        ci_lower = ate - t_crit * ate_se, ci_upper = ate + t_crit * ate_se,
        t_stat = as.numeric(ate / ate_se),
        p_value = 2 * stats::pnorm(abs(ate / ate_se), lower.tail = FALSE),
        comparison = paste(rownames(treatment_levels)[i], "vs",
                           rownames(treatment_levels)[i + 1])
      )
    }
    names(ate_list) <- sapply(ate_list, `[[`, "comparison")
    return(structure(
      list(ate_estimates = ate_list, outcome_model = NULL,
           treatment_levels = treatment_levels, model_type = "bart",
           family = NULL, robust_se = FALSE, vcov_matrix = NULL,
           estimand = "ATE", censoring_method = censoring_method,
           gps = gps, outcome = outcome,
           note = "GPS weights enter via AIPW augmentation; BART outcome model is unweighted"),
      class = "gps_ate"
    ))
  }

  # ---- lm / glm path ----------------------------------------------------------
  outcome_model <- if (model_type == "lm") {
    stats::lm(outcome ~ ., data = model_df, weights = w_active)
  } else {
    stats::glm(outcome ~ ., data = model_df, weights = w_active, family = family)
  }

  coef_est <- coef(outcome_model)
  keep     <- !is.na(coef_est)
  coef_est <- coef_est[keep]

  vcov_raw <- if (robust_se && requireNamespace("sandwich", quietly = TRUE))
    sandwich::vcovHC(outcome_model, type = "HC1")
  else
    vcov(outcome_model)
  vcov_matrix <- vcov_raw[
    rownames(vcov_raw) %in% names(coef_est),
    colnames(vcov_raw) %in% names(coef_est),
    drop = FALSE
  ]

  logistic  <- model_type == "glm" && identical(family$family, "binomial")
  n_comp    <- nrow(treatment_levels) - 1L
  ate_list  <- vector("list", n_comp)

  for (i in seq_len(n_comp)) {
    pd1   <- make_pred_df(treatment_levels[i,     ])
    pd2   <- make_pred_df(treatment_levels[i + 1, ])
    pred1 <- predict(outcome_model, newdata = pd1, type = "response")
    pred2 <- predict(outcome_model, newdata = pd2, type = "response")

    ate <- stats::weighted.mean(pred1[active] - pred2[active], w_active)

    mm1 <- stats::model.matrix(outcome_model, data = pd1)[active, keep, drop = FALSE]
    mm2 <- stats::model.matrix(outcome_model, data = pd2)[active, keep, drop = FALSE]
    if (logistic) {
      eta1  <- as.numeric(mm1 %*% coef_est)
      eta2  <- as.numeric(mm2 %*% coef_est)
      mu_p1 <- exp(eta1) / (1 + exp(eta1))^2
      mu_p2 <- exp(eta2) / (1 + exp(eta2))^2
      X1m   <- apply(mu_p1 * mm1, 2, stats::weighted.mean, w = w_active)
      X2m   <- apply(mu_p2 * mm2, 2, stats::weighted.mean, w = w_active)
    } else {
      X1m <- apply(mm1, 2, stats::weighted.mean, w = w_active)
      X2m <- apply(mm2, 2, stats::weighted.mean, w = w_active)
    }
    grad     <- X1m - X2m
    ate_se   <- as.numeric(sqrt(max(0, t(grad) %*% vcov_matrix %*% grad)))
    df_resid <- outcome_model$df.residual
    t_crit   <- stats::qt(0.975, df = df_resid)
    t_stat   <- ate / ate_se
    p_val    <- 2 * stats::pt(abs(t_stat), df = df_resid, lower.tail = FALSE)

    ate_list[[i]] <- list(
      treatment1 = treatment_levels[i,     ],
      treatment2 = treatment_levels[i + 1, ],
      ate        = as.numeric(ate),
      se         = ate_se,
      ci_lower   = ate - t_crit * ate_se,
      ci_upper   = ate + t_crit * ate_se,
      t_stat     = as.numeric(t_stat),
      p_value    = as.numeric(p_val),
      comparison = paste(rownames(treatment_levels)[i], "vs",
                         rownames(treatment_levels)[i + 1])
    )
  }
  names(ate_list) <- sapply(ate_list, `[[`, "comparison")

  structure(
    list(ate_estimates = ate_list, outcome_model = outcome_model,
         treatment_levels = treatment_levels, model_type = model_type,
         family = if (model_type == "glm") family else NULL,
         robust_se = robust_se, vcov_matrix = vcov_matrix,
         estimand = "ATE", censoring_method = censoring_method,
         gps = gps, outcome = outcome),
    class = "gps_ate"
  )
}

#' ATE or OR estimates as a gt table
#'
#' @param ate Object of class `"gps_ate"`
#' @return A gt table
#' @export
get_ate_table <- function(ate) {
  if (identical(ate$estimand, "OR")) {
    rows <- lapply(ate$or_estimates, function(x) {
      data.frame(comparison = x$comparison,
                 log_or = x$log_or, log_or_se = x$log_or_se,
                 or = x$or, or_lower = x$or_lower, or_upper = x$or_upper,
                 z_stat = x$z_stat, p_value = x$p_value,
                 stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, rows)
    return(
      df |>
        gt::gt(rowname_col = "comparison") |>
        gt::tab_header(
          title    = "GPS-Weighted Marginal Structural Odds Ratios",
          subtitle = paste0("logistic (no covariates)",
                            if (ate$robust_se) " | HC1 robust SE" else "",
                            " | n = ", length(ate$outcome))
        ) |>
        gt::cols_label(log_or = "log(OR)", log_or_se = "SE(log)",
                       or = "OR", or_lower = "OR lower", or_upper = "OR upper",
                       z_stat = "z", p_value = "p-value") |>
        gt::fmt_number(columns = c("log_or", "log_or_se", "or",
                                   "or_lower", "or_upper", "z_stat"), decimals = 4) |>
        gt::fmt_scientific(columns = "p_value", decimals = 3) |>
        gt::tab_style(style = gt::cell_text(weight = "bold"),
                      locations = gt::cells_column_labels())
    )
  }

  rows <- lapply(ate$ate_estimates, function(x) {
    data.frame(comparison = x$comparison,
               estimate   = x$ate,  std_error = x$se,
               ci_lower   = x$ci_lower, ci_upper = x$ci_upper,
               t_stat     = x$t_stat,   p_value  = x$p_value,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)

  df |>
    gt::gt(rowname_col = "comparison") |>
    gt::tab_header(
      title    = "GPS-Weighted Average Treatment Effects",
      subtitle = paste0(ate$model_type, " model",
                        if (ate$robust_se) " | HC1 robust SE" else "",
                        " | n = ", length(ate$outcome),
                        " | GPS ESS = ",
                        round(ate$gps$diagnostics$effective_sample_size /
                                ate$gps$specs$n * 100, 1), "%")
    ) |>
    gt::cols_label(estimate = "ATE", std_error = "SE",
                   ci_lower = "CI lower", ci_upper = "CI upper",
                   t_stat = "t", p_value = "p-value") |>
    gt::fmt_number(columns = c("estimate", "std_error", "ci_lower", "ci_upper",
                                "t_stat"), decimals = 4) |>
    gt::fmt_scientific(columns = "p_value", decimals = 3) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"),
                  locations = gt::cells_column_labels())
}

# =============================================================================
# 2. User model API — apply any function to the augmented data
# =============================================================================

#' Extract the augmented (imputed) dataset from a GPS object
#'
#' Returns a tidy data frame containing the (possibly imputed) treatment
#' values, GPS weights, and optionally the outcome.  This is the entry point
#' for applying user-specified models to the augmented data.
#'
#' @param gps Object of class `"gps_weights"`
#' @param outcome Optional numeric outcome vector of length n.  When supplied
#'   it is attached as column `"outcome"`.
#' @return A `data.frame` with columns `T1`, `T2`, `w`, and optionally
#'   `outcome`, with one row per observation.  Rows for hull-trimmed
#'   observations have `w = NA`.
#' @export
get_augmented_data <- function(gps, outcome = NULL) {
  y <- gps$y_imputed %||% gps$data$y
  d <- ncol(y)
  df <- as.data.frame(y)
  colnames(df) <- paste0("T", seq_len(d))
  df$w <- gps$weights
  if (!is.null(outcome)) {
    if (length(outcome) != nrow(df))
      stop("'outcome' must have length n = ", nrow(df))
    df$outcome <- outcome
  }
  df
}

#' Apply a user-specified model function to the augmented dataset
#'
#' Constructs the augmented data frame via [get_augmented_data()] and passes
#' it to `FUN`.  The function must accept a `data.frame` as its first argument.
#' GPS weights are available as column `"w"` and hull-trimmed observations have
#' `w = NA`.
#'
#' @param gps Object of class `"gps_weights"`
#' @param FUN A function `f(data, ...)` applied to the augmented data frame.
#'   The data frame will have columns `T1`, `T2`, `w`, and optionally
#'   `outcome`.  Rows with `w = NA` are hull-trimmed and should typically be
#'   excluded.  See examples.
#' @param outcome Optional numeric outcome vector (length n) to attach.
#' @param ... Additional arguments forwarded to `FUN`.
#'
#' @return Whatever `FUN` returns.
#' @export
#'
#' @examples
#' \dontrun{
#' # Unweighted linear model on augmented data
#' apply_augmented(gps, outcome = y_obs,
#'   FUN = function(d) lm(outcome ~ T1 + T2, data = d, weights = d$w))
#'
#' # User glm on non-trimmed rows only
#' apply_augmented(gps, outcome = y_obs,
#'   FUN = function(d) {
#'     d <- d[!is.na(d$w), ]
#'     glm(outcome ~ T1 * T2, data = d, weights = d$w, family = Gamma(link = "log"))
#'   })
#' }
apply_augmented <- function(gps, FUN, outcome = NULL, ...) {
  if (!is.function(FUN))
    stop("'FUN' must be a function")
  d <- get_augmented_data(gps, outcome = outcome)
  FUN(d, ...)
}
