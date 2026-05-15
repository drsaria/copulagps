# =============================================================================
# generics.R — S3 print / summary methods
# =============================================================================

#' @export
print.gps_weights <- function(x, ...) {
  d    <- x$diagnostics
  sp   <- x$specs
  epct <- d$effective_sample_size / sp$n * 100
  status <- if (epct > 90 && d$weight_cv < 0.5) "EXCELLENT"
             else if (epct > 70 && d$weight_cv < 1.0) "GOOD"
             else "NEEDS ATTENTION"

  cat("Copula-Based GPS Weights\n")
  cat("========================\n")
  cat("n:", sp$n, "| d:", sp$d, "| margins:", paste(sp$margins, collapse = ", "), "\n")
  cat("Copula:", sp$copula_type, "| stabilization:", sp$stabilization, "\n")
  if (identical(sp$trim_method, "hull"))
    cat("Hull-trimmed:", d$n_overlap_trimmed, "obs\n")
  cat("Auto-adjust:", sp$auto_adjust)
  if (isTRUE(d$residual_copula_used)) cat("  [residual copula active]")
  cat("\nParameters:", if (isTRUE(d$params_estimated)) "estimated" else "user-provided", "\n\n")

  cat("Weight Summary:\n"); print(summary(x$weights[!is.na(x$weights)])); cat("\n")

  cat("Diagnostics:\n")
  cat("  ESS:        ", round(d$effective_sample_size, 1),
      " (", round(epct, 1), "%)\n", sep = "")
  cat("  Weight CV:  ", round(d$weight_cv, 3), "\n")
  cat("  Range:      [", round(d$weight_range[1], 3), ", ",
      round(d$weight_range[2], 3), "]\n", sep = "")
  cat("  Emp. cop.:  ", d$empirical_copula_family,  "\n")
  cat("  Par. cop.:  ", d$parametric_copula_family, "\n")
  if (any(!is.na(d$marginal_r_squared)))
    cat("  Marginal R\u00b2:", paste(round(d$marginal_r_squared, 3), collapse = ", "), "\n")
  cat("  Status:     ", status, "\n")
  invisible(x)
}

#' @export
summary.gps_weights <- function(object, ...) {
  cat("Copula-Based GPS Weights \u2014 Detailed Summary\n")
  cat("===============================================\n\n")
  print(object)
  if (object$diagnostics$params_estimated) {
    cat("\nMarginal Coefficients:\n")
    print(get_coef_table(object, as_gt = FALSE))
  }
  cat("\nCopula Information:\n")
  cp <- get_copula_params(object)
  for (nm in c("empirical", "parametric")) {
    label <- if (nm == "empirical") "Empirical" else "Parametric"
    cat(label, "copula:\n")
    cat("  Family:       ", cp[[nm]]$family_name, "\n")
    cat("  Parameter 1:  ", round(cp[[nm]]$parameter1, 4), "\n")
    if (!is.na(cp[[nm]]$parameter2) && cp[[nm]]$parameter2 != 0)
      cat("  Parameter 2:  ", round(cp[[nm]]$parameter2, 4), "\n")
    cat("  Kendall's \u03c4: ", round(cp[[nm]]$tau, 4), "\n")
  }
  invisible(object)
}

#' @export
print.balance_assessment <- function(x, ...) {
  rating <- if (x$improvement$corr_reduction > 20) "EXCELLENT"
             else if (x$improvement$corr_reduction > 10) "GOOD"
             else if (x$improvement$corr_reduction > 5)  "MODEST"
             else "LIMITED"
  cat("Covariate Balance Assessment\n")
  cat("============================\n")
  cat("Covariates:", x$n_covariates, "| Treatments:", x$n_treatments, "\n\n")
  cat("Before:  mean|r| =", round(x$before$mean_abs_corr, 4),
      " max|r| =", round(x$before$max_abs_corr, 4), "\n")
  cat("After:   mean|r| =", round(x$after$mean_abs_corr, 4),
      " max|r| =", round(x$after$max_abs_corr, 4), "\n\n")
  cat("Mean |r| reduction:", round(x$improvement$corr_reduction, 1), "%\n")
  cat("Max  |r| reduction:", round(x$improvement$max_corr_reduction, 1), "%\n")
  cat("Assessment:", rating, "\n")
  invisible(x)
}

#' @export
print.gps_ate <- function(x, ...) {
  cat("GPS-Weighted ATE Estimates\n")
  cat("==========================\n")
  cat("Model:", x$model_type)
  if (!is.null(x$family)) cat(" (", x$family$family, ")", sep = "")
  cat(" | robust SE:", x$robust_se, "| n:", length(x$outcome), "\n")
  if (!is.null(x$gps$diagnostics$effective_sample_size))
    cat("GPS ESS:",
        round(x$gps$diagnostics$effective_sample_size / x$gps$specs$n * 100, 1), "%\n")
  cat("\n")
  for (a in x$ate_estimates) {
    cat(a$comparison, ":\n")
    cat("  ATE:", round(a$ate, 4), " SE:", round(a$se, 4), "\n")
    cat("  95% CI: [", round(a$ci_lower, 4), ", ", round(a$ci_upper, 4), "]\n",
        sep = "")
    if (!is.null(a$note)) cat("  Note:", a$note, "\n")
    cat("\n")
  }
  invisible(x)
}

#' @export
summary.gps_ate <- function(object, ...) {
  if (!is.null(object$gps$diagnostics$effective_sample_size)) {
    epct <- object$gps$diagnostics$effective_sample_size / object$gps$specs$n * 100
    cat("GPS ESS:", round(epct, 1), "%  CV:",
        round(object$gps$diagnostics$weight_cv, 3), "\n")
    cat("GPS quality:", if (epct > 80) "GOOD" else "SUBOPTIMAL", "\n\n")
  }
  cat("Outcome Model:\n"); print(summary(object$outcome_model)); cat("\n")
  print(object)
  invisible(object)
}

# =============================================================================
# gps_drf methods
# =============================================================================

#' @export
print.gps_drf <- function(x, ...) {
  cat("GPS Dose-Response Function\n")
  cat("  Treatments:", paste(x$trt_names, collapse = ", "), "\n")
  cat("  Mode:", x$vary, "\n")
  cat("  Outcome model:", x$model_type, "\n")
  if (x$vary == "marginal") {
    for (nm in names(x$slices)) {
      sl <- x$slices[[nm]]
      cat("  ", nm, ": grid [", round(min(sl$t), 3), ",",
          round(max(sl$t), 3), "] (", nrow(sl), " points)\n", sep = "")
    }
  } else {
    cat("  Joint grid:", nrow(x$surface), "points\n")
  }
  invisible(x)
}

#' Plot a GPS dose-response function
#'
#' @param x Object of class `"gps_drf"` from [drf()].
#' @param ... Additional arguments (unused).
#' @return A `ggplot2` object.
#' @export
plot.gps_drf <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot.gps_drf requires the ggplot2 package")

  if (x$vary == "marginal") {
    df <- do.call(rbind, x$slices)
    has_se <- !all(is.na(df$se))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$mu_hat)) +
      ggplot2::geom_line(linewidth = 0.8, colour = "#2166ac") +
      ggplot2::facet_wrap(~ treatment, scales = "free_x") +
      ggplot2::labs(x = "Treatment value", y = expression(hat(mu)(t)),
                    title = "GPS-Weighted Dose-Response Function") +
      ggplot2::theme_bw()
    if (has_se)
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.2, fill = "#2166ac")
    return(p)
  }

  # Joint surface (d = 2)
  ggplot2::ggplot(x$surface,
                  ggplot2::aes(x = .data[[x$trt_names[1]]],
                               y = .data[[x$trt_names[2]]],
                               fill = .data$mu_hat)) +
    ggplot2::geom_tile() +
    ggplot2::geom_contour(ggplot2::aes(z = .data$mu_hat),
                          colour = "white", alpha = 0.5) +
    ggplot2::scale_fill_distiller(palette = "RdYlBu", name = expression(hat(mu))) +
    ggplot2::labs(x = x$trt_names[1], y = x$trt_names[2],
                  title = "GPS-Weighted Dose-Response Surface") +
    ggplot2::theme_bw()
}
