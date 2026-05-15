# =============================================================================
# multi_chain.R — multiple augmentation chains and Rubin's rules combination
# =============================================================================

#' Run M independent augmentation chains
#'
#' Runs `estimate_gps()` M times from different random seeds, varying the
#' stochastic initialisation and copula-imputation draws each time.  This is
#' the recommended approach for propagating imputation uncertainty: M
#' independent chains give proper between-imputation variance without requiring
#' MCMC mixing diagnostics.
#'
#' For each chain you can either:
#' - Collect the raw `gps_weights` objects (default, for later analysis), or
#' - Immediately estimate ATEs and store results (set `outcome`), or
#' - Apply any custom function via `user_fn` (overrides ATE estimation).
#'
#' @param data Prepared data list from [prepare_data()] (or a hand-built list
#'   with `$y` and `$x`).
#' @param margins Character vector of length d for marginal distributions.
#' @param M Integer number of chains.  Minimum 2 (Rubin's rules require M ≥ 2).
#'   Default: `5L`.
#' @param seed Integer base seed.  Chain m uses `seed + m` so results are
#'   exactly reproducible.  `NULL` (default) does not set seeds.
#' @param outcome Optional numeric outcome vector.  When supplied and
#'   `user_fn = NULL`, [estimate_ate()] is called on each chain.
#' @param user_fn Optional function `f(gps)` applied to each chain's
#'   `gps_weights` object.  Must return a list or atomic vector that can be
#'   combined across chains.  When supplied, `outcome` is still passed to
#'   [get_augmented_data()] and available inside `user_fn` via
#'   `get_augmented_data(gps, outcome)`.  See examples.
#' @param verbose Print chain-level progress messages.  Default: `TRUE`.
#' @param ... Additional arguments forwarded to [estimate_gps()].
#'
#' @return Object of class `"gps_chains"` with components:
#'   \describe{
#'     \item{`chains`}{List of M results (GPS objects, ATE objects, or
#'       user-function outputs)}
#'     \item{`M`}{Number of chains}
#'     \item{`mode`}{"gps", "ate", or "user"}
#'     \item{`seeds`}{Seeds used (or NULL)}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' d    <- prepare_data(df, treatments = c("T1","T2"), covariates = c("age","sex"),
#'                      outcome = "Y", lod = c(0.5, 0.5))
#'
#' # Collect GPS objects; combine ATEs later
#' chains <- run_chains(d, margins = c("lognormal","lognormal"), M = 10, seed = 42)
#'
#' # Estimate ATE in each chain, then pool with Rubin's rules
#' chains <- run_chains(d, margins = c("lognormal","lognormal"), M = 10,
#'                      seed = 42, outcome = d$outcome)
#' pooled <- combine_rubin(chains)
#'
#' # Custom user function: return coefficients of a weighted lm
#' chains <- run_chains(d, margins = c("lognormal","lognormal"), M = 10,
#'                      seed = 42, outcome = d$outcome,
#'                      user_fn = function(gps) {
#'                        aug  <- get_augmented_data(gps, outcome = d$outcome)
#'                        fit  <- lm(outcome ~ T1 + T2, data = aug, weights = aug$w)
#'                        list(estimate = coef(fit)["T1"], se = sqrt(vcov(fit)["T1","T1"]))
#'                      })
#' pooled <- combine_rubin(chains)
#' }
run_chains <- function(data, margins, M = 1L, seed = NULL, outcome = NULL,
                        user_fn = NULL, verbose = TRUE, ...) {
  M <- as.integer(M)
  if (M < 2L) stop("M must be >= 2 for Rubin's rules to be applicable")
  if (M > 5L)
    message("M = ", M, ": running many chains is expensive for d > 2. ",
            "Consider M = 3-5 and parallel execution via future::plan('multisession').")

  mode <- if (!is.null(user_fn)) "user" else if (!is.null(outcome)) "ate" else "gps"

  chains <- vector("list", M)
  seeds  <- if (!is.null(seed)) seed + seq_len(M) else NULL

  for (m in seq_len(M)) {
    if (!is.null(seed)) set.seed(seeds[m])
    if (verbose) message("Chain ", m, "/", M, " ...")

    gps_m <- tryCatch(
      estimate_gps(data, margins = margins, verbose = FALSE, ...),
      error = function(e) {
        warning("Chain ", m, " failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(gps_m)) next

    chains[[m]] <- switch(mode,
      gps  = gps_m,
      ate  = tryCatch(
               estimate_ate(gps_m, outcome = outcome),
               error = function(e) {
                 warning("ATE failed in chain ", m, ": ", conditionMessage(e))
                 NULL
               }),
      user = tryCatch(
               user_fn(gps_m),
               error = function(e) {
                 warning("user_fn failed in chain ", m, ": ", conditionMessage(e))
                 NULL
               })
    )
  }

  n_ok <- sum(!vapply(chains, is.null, logical(1)))
  if (verbose) message("Completed ", n_ok, "/", M, " chains successfully")
  if (n_ok < 2L) stop("Fewer than 2 chains completed — cannot apply Rubin's rules")

  structure(
    list(chains = chains, M = M, mode = mode, seeds = seeds, outcome = outcome),
    class = "gps_chains"
  )
}

# =============================================================================
# Rubin's rules combination
# =============================================================================

#' Combine M chain estimates using Rubin's rules
#'
#' Pools point estimates and standard errors across M chains following Rubin
#' (1987).  The total variance accounts for both within-chain variance (model
#' uncertainty) and between-chain variance (imputation uncertainty).
#'
#' **Supported inputs:**
#' - `chains` with `mode = "ate"`: extracts ATE and SE from each
#'   `"gps_ate"` object automatically.
#' - `chains` with `mode = "user"`: expects each chain result to be a list
#'   with named elements `$estimate` (numeric scalar or vector) and `$se`
#'   (same shape).  Multiple named estimands are pooled independently.
#' - `chains` with `mode = "gps"`: not directly combinable by Rubin's rules;
#'   a warning is issued and GPS diagnostics are averaged instead.
#'
#' @param chains Object of class `"gps_chains"` from [run_chains()].
#' @param comparison Integer or character scalar selecting which ATE comparison
#'   to pool when `mode = "ate"` and multiple contrasts exist.  Default: `1L`
#'   (first contrast).
#'
#' @return Object of class `"rubin_result"` with components:
#'   \describe{
#'     \item{`estimate`}{Pooled point estimate(s)}
#'     \item{`se`}{Total standard error(s)}
#'     \item{`ci_lower`, `ci_upper`}{95% CI based on Barnard-Rubin df}
#'     \item{`within_var`}{Within-chain variance (Ū)}
#'     \item{`between_var`}{Between-chain variance (B)}
#'     \item{`total_var`}{Total variance (T = Ū + (1 + 1/M)B)}
#'     \item{`fmi`}{Fraction of missing information}
#'     \item{`df`}{Barnard-Rubin degrees of freedom}
#'     \item{`M`}{Number of chains used}
#'   }
#' @export
combine_rubin <- function(chains, comparison = 1L) {
  if (!inherits(chains, "gps_chains"))
    stop("'chains' must be a 'gps_chains' object from run_chains()")

  ok <- !vapply(chains$chains, is.null, logical(1))
  M  <- sum(ok)
  if (M < 2L) stop("Need at least 2 valid chains for Rubin's rules")

  valid_chains <- chains$chains[ok]
  mode         <- chains$mode

  if (mode == "gps") {
    warning("mode = 'gps': Rubin's rules require scalar estimates. ",
            "Returning averaged GPS diagnostics instead.")
    return(.pool_gps_diagnostics(valid_chains, M))
  }

  # Extract (estimate, se) pairs from each chain
  ests <- ses <- NULL

  if (mode == "ate") {
    ests <- vapply(valid_chains, function(ch) {
      if (is.character(comparison))
        ch$ate_estimates[[comparison]]$ate
      else
        ch$ate_estimates[[comparison]]$ate
    }, numeric(1))
    ses  <- vapply(valid_chains, function(ch) {
      if (is.character(comparison))
        ch$ate_estimates[[comparison]]$se
      else
        ch$ate_estimates[[comparison]]$se
    }, numeric(1))
    estimand_name <- if (is.character(comparison)) comparison
                     else names(valid_chains[[1]]$ate_estimates)[comparison]
  } else {
    # user mode: expect list(estimate = ..., se = ...)
    first <- valid_chains[[1]]
    if (!all(c("estimate", "se") %in% names(first)))
      stop("user_fn results must be named lists with 'estimate' and 'se' elements")
    if (length(first$estimate) == 1L) {
      ests <- vapply(valid_chains, function(ch) as.numeric(ch$estimate), numeric(1))
      ses  <- vapply(valid_chains, function(ch) as.numeric(ch$se),       numeric(1))
      estimand_name <- "user_estimate"
    } else {
      # Multiple named estimands — pool each independently
      est_names <- names(first$estimate)
      result_list <- lapply(seq_along(est_names), function(k) {
        ests_k <- vapply(valid_chains, function(ch) as.numeric(ch$estimate[k]), numeric(1))
        ses_k  <- vapply(valid_chains, function(ch) as.numeric(ch$se[k]),       numeric(1))
        .rubin_combine_scalar(ests_k, ses_k, M, label = est_names[k])
      })
      names(result_list) <- est_names
      return(structure(list(estimates = result_list, M = M, mode = mode),
                       class = "rubin_result"))
    }
  }

  .rubin_combine_scalar(ests, ses, M, label = estimand_name)
}

#' Apply Rubin's rules to a single scalar estimand
#' @keywords internal
.rubin_combine_scalar <- function(ests, ses, M, label = "estimate") {
  Q_bar   <- mean(ests)
  U_bar   <- mean(ses^2)                        # within-imputation variance
  B       <- stats::var(ests)                   # between-imputation variance
  T_total <- U_bar + (1 + 1 / M) * B           # total variance
  se_pool <- sqrt(T_total)

  # Barnard-Rubin degrees of freedom
  r   <- (1 + 1 / M) * B / U_bar
  df  <- (M - 1) * (1 + 1 / r)^2
  t_c <- stats::qt(0.975, df = df)

  fmi <- (r + 2 / (df + 3)) / (r + 1)           # fraction of missing information

  structure(
    list(
      estimand   = label,
      estimate   = Q_bar,
      se         = se_pool,
      ci_lower   = Q_bar - t_c * se_pool,
      ci_upper   = Q_bar + t_c * se_pool,
      within_var = U_bar,
      between_var = B,
      total_var  = T_total,
      fmi        = fmi,
      df         = df,
      M          = M,
      chain_estimates = ests,
      chain_ses       = ses
    ),
    class = "rubin_result"
  )
}

#' Average GPS diagnostics across chains (fallback for mode = "gps")
#' @keywords internal
.pool_gps_diagnostics <- function(valid_chains, M) {
  diag_list <- lapply(valid_chains, `[[`, "diagnostics")
  avg_ess   <- mean(vapply(diag_list, `[[`, numeric(1), "effective_sample_size"))
  avg_cv    <- mean(vapply(diag_list, `[[`, numeric(1), "weight_cv"))
  cat(sprintf("GPS diagnostics (averaged over %d chains):\n  ESS: %.1f  CV: %.3f\n",
              M, avg_ess, avg_cv))
  invisible(list(mean_ess = avg_ess, mean_cv = avg_cv, M = M))
}

# =============================================================================
# Print methods
# =============================================================================

#' @export
print.gps_chains <- function(x, ...) {
  ok  <- sum(!vapply(x$chains, is.null, logical(1)))
  cat("GPS Augmentation Chains\n")
  cat("=======================\n")
  cat("Chains: ", ok, "/", x$M, " completed | mode: ", x$mode, "\n", sep = "")
  if (!is.null(x$seeds)) cat("Seeds: ", x$seeds[1], "\u2026", tail(x$seeds, 1), "\n")
  invisible(x)
}

#' @export
print.rubin_result <- function(x, ...) {
  if (!is.null(x$estimates)) {
    cat("Rubin's Rules \u2014 Multiple Estimands\n")
    cat("=====================================\n")
    for (nm in names(x$estimates)) {
      r <- x$estimates[[nm]]
      cat(sprintf("  %s: %.4f \u00b1 %.4f  [%.4f, %.4f]  FMI = %.3f\n",
                  nm, r$estimate, r$se, r$ci_lower, r$ci_upper, r$fmi))
    }
    return(invisible(x))
  }
  cat("Rubin's Rules Result\n")
  cat("====================\n")
  cat("Estimand:", x$estimand, "| M =", x$M, "\n")
  cat(sprintf("Estimate:  %.4f\nSE:        %.4f\n95%% CI:   [%.4f, %.4f]\n",
              x$estimate, x$se, x$ci_lower, x$ci_upper))
  cat(sprintf("FMI:       %.3f  (fraction of missing information)\n", x$fmi))
  cat(sprintf("Between-chain SD: %.4f | Within-chain SE: %.4f\n",
              sqrt(x$between_var), sqrt(x$within_var)))
  invisible(x)
}
