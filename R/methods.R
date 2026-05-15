# =============================================================================
# methods.R — GPS method dispatcher: copulagps, AIPW, and user-defined methods
# =============================================================================

#' Estimate GPS weights using a specified method
#'
#' A unified interface for computing GPS weights and estimating ATEs under
#' multiple methods.  Built-in options are `"copulagps"` (the default copula
#' approach) and `"aipw"` (augmented IPW / doubly-robust estimator).  You can
#' also register up to two user-defined methods via `user_methods`.
#'
#' @param data Prepared data list from [prepare_data()].
#' @param margins Character vector of length d (required for `"copulagps"`).
#' @param outcome Optional numeric outcome vector.  Required for `"aipw"` and
#'   to return ATE estimates for any method.
#' @param method Character scalar or function: `"copulagps"`, `"aipw"`, or
#'   the name of a method registered in `user_methods`.  Alternatively, supply
#'   a function directly — it will be called as `method(data, ...)`.
#' @param user_methods Named list of up to 2 user-defined method functions.
#'   Each function must accept `(data, outcome, ...)` and return a list with
#'   at minimum `$weights` (numeric vector, length n).  The name becomes the
#'   method label.
#' @param treatment_levels Treatment regime matrix for ATE contrasts (passed to
#'   [estimate_ate()]).
#' @param ... Additional arguments forwarded to [estimate_gps()] for the
#'   `"copulagps"` method, or to the user method function.
#'
#' @return Object of class `"gps_method_result"` with components:
#'   \describe{
#'     \item{`method`}{Method name used}
#'     \item{`gps`}{`gps_weights` object (or equivalent list for custom methods)}
#'     \item{`ate`}{`gps_ate` object, or `NULL` if no outcome supplied}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' d <- prepare_data(df, treatments = c("T1","T2"), covariates = c("age","sex"),
#'                   outcome = "Y")
#'
#' # Built-in copulagps
#' res <- fit_gps_method(d, margins = c("normal","normal"),
#'                       outcome = d$outcome, method = "copulagps")
#'
#' # Doubly-robust AIPW
#' res <- fit_gps_method(d, margins = c("normal","normal"),
#'                       outcome = d$outcome, method = "aipw")
#'
#' # Custom method (e.g., CBGPS-style balancing)
#' my_method <- function(data, outcome, ...) {
#'   w <- rep(1 / nrow(data$y), nrow(data$y))  # uniform weights placeholder
#'   list(weights = w, gps = NULL)
#' }
#' res <- fit_gps_method(d, margins = c("normal","normal"),
#'                       outcome = d$outcome, method = "my_method",
#'                       user_methods = list(my_method = my_method))
#' }
fit_gps_method <- function(data, margins = NULL, outcome = NULL,
                             method = "copulagps", user_methods = NULL,
                             treatment_levels, ...) {
  if (length(user_methods) > 2)
    stop("At most 2 user-defined methods are supported")

  # ---- Resolve method --------------------------------------------------------
  if (is.function(method)) {
    method_fn   <- method
    method_name <- "user_fn"
  } else if (is.character(method)) {
    method_name <- method
    method_fn   <- NULL
  } else {
    stop("'method' must be a character string or a function")
  }

  gps_obj <- NULL
  ate_obj <- NULL

  # ---- Dispatch --------------------------------------------------------------
  if (identical(method_name, "copulagps")) {
    if (is.null(margins))
      stop("'margins' is required for method = 'copulagps'")
    gps_obj <- estimate_gps(data, margins = margins, ...)
    if (!is.null(outcome)) {
      ate_obj <- if (missing(treatment_levels))
        estimate_ate(gps_obj, outcome = outcome)
      else
        estimate_ate(gps_obj, outcome = outcome,
                     treatment_levels = treatment_levels)
    }

  } else if (identical(method_name, "aipw")) {
    if (is.null(margins))
      stop("'margins' is required for method = 'aipw' (needs GPS weights)")
    if (is.null(outcome))
      stop("'outcome' is required for method = 'aipw'")
    gps_obj <- estimate_gps(data, margins = margins, ...)
    ate_obj <- .estimate_aipw(gps_obj, outcome = outcome,
                               treatment_levels = if (missing(treatment_levels))
                                 NULL else treatment_levels)

  } else if (!is.null(user_methods) && method_name %in% names(user_methods)) {
    fn      <- user_methods[[method_name]]
    result  <- fn(data, outcome = outcome, ...)
    if (!is.list(result) || !"weights" %in% names(result))
      stop("User method '", method_name, "' must return a list with a 'weights' element")
    gps_obj <- result$gps   # may be NULL
    # Wrap result as a minimal gps_weights-like object for downstream use
    if (is.null(gps_obj)) {
      gps_obj <- structure(
        list(weights = result$weights, data = data,
             y_imputed = NULL, hull_membership = rep(TRUE, nrow(data$y)),
             diagnostics = list(effective_sample_size = NA, weight_cv = NA,
                                n_hull_trimmed = 0L),
             specs = list(n = nrow(data$y), d = ncol(data$y),
                          margins = margins)),
        class = "gps_weights"
      )
    }
    ate_obj <- if (!is.null(outcome))
      tryCatch(
        if (missing(treatment_levels))
          estimate_ate(gps_obj, outcome = outcome)
        else
          estimate_ate(gps_obj, outcome = outcome,
                       treatment_levels = treatment_levels),
        error = function(e) NULL
      )

  } else if (!is.null(method_fn)) {
    result  <- method_fn(data, outcome = outcome, ...)
    if (!is.list(result) || !"weights" %in% names(result))
      stop("Method function must return a list with a 'weights' element")
    gps_obj <- result$gps %||% structure(
      list(weights = result$weights, data = data, y_imputed = NULL,
           hull_membership = rep(TRUE, nrow(data$y)),
           diagnostics = list(effective_sample_size = NA, weight_cv = NA,
                              n_hull_trimmed = 0L),
           specs = list(n = nrow(data$y), d = ncol(data$y), margins = margins)),
      class = "gps_weights"
    )
    ate_obj <- NULL

  } else {
    stop("Unknown method: '", method_name, "'. ",
         "Built-in: 'copulagps', 'aipw'. ",
         "Register custom methods via user_methods = list(<name> = <fn>).")
  }

  structure(
    list(method = method_name, gps = gps_obj, ate = ate_obj),
    class = "gps_method_result"
  )
}

# =============================================================================
# AIPW (doubly-robust) estimator for continuous treatment
# =============================================================================

#' Augmented IPW estimator for continuous treatment
#'
#' Combines the GPS-weighted outcome estimator with an outcome regression
#' (plug-in) correction for double robustness.
#'
#' For treatment contrast (t1, t2):
#'   ATE_AIPW = E[μ(t1,X) - μ(t2,X)]
#'            + E[W_i * (Y_i - μ(T_i,X))]
#' where W_i are the GPS weights and μ(t,X) = E[Y|T=t,X] from an unweighted
#' outcome regression.
#'
#' @keywords internal
.estimate_aipw <- function(gps, outcome, treatment_levels = NULL, covariates = NULL) {
  treatments <- gps$data$y
  weights    <- gps$weights
  n          <- length(outcome)
  d          <- gps$specs$d
  active     <- !is.na(weights)

  trt_names <- paste0("T", seq_len(d))
  X         <- as.data.frame(treatments)
  colnames(X) <- trt_names
  if (!is.null(covariates)) {
    cov_names   <- colnames(covariates) %||% paste0("covar_", seq_len(ncol(covariates)))
    X           <- cbind(X, as.data.frame(covariates))
    colnames(X) <- c(trt_names, cov_names)
  }

  # Unweighted outcome regression (step 1 of AIPW)
  model_df   <- data.frame(outcome = outcome, X)
  mu_model   <- stats::lm(outcome ~ ., data = model_df)
  mu_fitted  <- predict(mu_model)   # E[Y | T_i, X_i] at observed treatments

  if (is.null(treatment_levels)) {
    treatment_levels <- rbind(
      apply(treatments, 2, stats::quantile, probs = 0.75),
      apply(treatments, 2, stats::quantile, probs = 0.25)
    )
    rownames(treatment_levels) <- c("high", "low")
  }

  make_pred_df <- function(treat) {
    df_t <- as.data.frame(matrix(rep(treat, n), nrow = n, byrow = TRUE))
    colnames(df_t) <- trt_names
    if (!is.null(covariates)) df_t <- cbind(df_t, as.data.frame(covariates))
    colnames(df_t) <- colnames(X)
    df_t
  }

  n_comp   <- nrow(treatment_levels) - 1L
  ate_list <- vector("list", n_comp)

  for (i in seq_len(n_comp)) {
    pd1    <- make_pred_df(treatment_levels[i,     ])
    pd2    <- make_pred_df(treatment_levels[i + 1, ])
    mu1    <- predict(mu_model, newdata = pd1)
    mu2    <- predict(mu_model, newdata = pd2)

    # Plug-in term: E[mu(t1,X) - mu(t2,X)]
    plugin <- mean(mu1 - mu2)

    # Augmentation correction: E[W * (Y - mu(T,X))] — uses observed weights
    w_std  <- weights[active] / sum(weights[active]) * sum(active)
    aug    <- mean(w_std * (outcome[active] - mu_fitted[active]))

    ate    <- plugin + aug
    # Bootstrap SE placeholder (delta method for AIPW requires more structure)
    se_plugin <- stats::sd(mu1 - mu2) / sqrt(n)

    ate_list[[i]] <- list(
      treatment1 = treatment_levels[i,     ],
      treatment2 = treatment_levels[i + 1, ],
      ate        = as.numeric(ate),
      se         = se_plugin,
      ci_lower   = ate - 1.96 * se_plugin,
      ci_upper   = ate + 1.96 * se_plugin,
      t_stat     = ate / se_plugin,
      p_value    = 2 * stats::pnorm(abs(ate / se_plugin), lower.tail = FALSE),
      comparison = paste(rownames(treatment_levels)[i], "vs",
                         rownames(treatment_levels)[i + 1]),
      note       = "AIPW: SE is plug-in approximation; use bootstrap for inference"
    )
  }
  names(ate_list) <- sapply(ate_list, `[[`, "comparison")

  structure(
    list(ate_estimates = ate_list, outcome_model = mu_model,
         treatment_levels = treatment_levels,
         model_type = "aipw", family = NULL, robust_se = FALSE,
         vcov_matrix = NULL, gps = gps, outcome = outcome),
    class = "gps_ate"
  )
}

#' @export
print.gps_method_result <- function(x, ...) {
  cat("GPS Method Result\n")
  cat("=================\n")
  cat("Method:", x$method, "\n")
  if (!is.null(x$gps) && !is.null(x$gps$diagnostics$effective_sample_size))
    cat("ESS:", round(x$gps$diagnostics$effective_sample_size, 1), "\n")
  if (!is.null(x$ate)) print(x$ate)
  invisible(x)
}
