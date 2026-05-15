# =============================================================================
# plots.R — diagnostic plots for marginal fit and GPS weights
# =============================================================================

#' PIT uniformity plots for marginal distribution fit assessment
#'
#' Plots probability-integral-transform (PIT) values against Uniform(0,1)
#' quantiles.  Points on the diagonal indicate a good marginal fit.
#'
#' @param gps Object of class `"gps_weights"`
#' @return A patchwork object (printed and returned invisibly)
#' @export
check_marginal_fit <- function(gps) {
  y       <- gps$data$y
  margins <- gps$specs$margins
  fits    <- gps$marginal_fits

  plots <- lapply(seq_len(ncol(y)), function(j) {
    y_j   <- y[, j]
    fit_j <- fits[[j]]
    df_j  <- if (!is.na(fit_j$df)) fit_j$df else 10

    pit <- tryCatch(
      .pdist_margin(y_j, margins[j], fit_j$mean, fit_j$sd, df = df_j),
      error = function(e) {
        warning("PIT failed for treatment ", j, ": ", e$message)
        as.numeric(stats::pnorm(scale(y_j)))
      }
    )
    pit   <- pmax(pmin(pit, 1 - 1e-10), 1e-10)
    n_j   <- length(pit)
    r2_lb <- if (!is.na(gps$diagnostics$marginal_r_squared[j]))
               paste0("R\u00b2 = ", round(gps$diagnostics$marginal_r_squared[j], 3))
             else ""

    df_plot <- data.frame(empirical   = sort(pit),
                          theoretical = seq_len(n_j) / (n_j + 1L))
    ggplot2::ggplot(df_plot, ggplot2::aes(x = theoretical, y = empirical)) +
      ggplot2::geom_point(alpha = 0.4, size = 1.2, colour = "#2c7bb6") +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red",
                           linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(title    = paste0("T", j, " \u2014 ", margins[j]),
                    subtitle = r2_lb,
                    x = "Theoretical Uniform Quantile", y = "PIT Quantile") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  })

  combined <- Reduce(`+`, plots) +
    patchwork::plot_layout(ncol = min(2L, length(plots)))
  print(combined)
  invisible(combined)
}

#' GPS weight diagnostic plots
#'
#' Four-panel figure: weight histogram, weights vs observation index,
#' normal Q-Q of weights, and a text summary panel.
#'
#' @param gps Object of class `"gps_weights"`
#' @return A patchwork object (printed and returned invisibly)
#' @export
plot_gps_weights <- function(gps) {
  w    <- gps$weights[!is.na(gps$weights)]
  n    <- gps$specs$n
  diag <- gps$diagnostics
  df_w <- data.frame(weight = w, index = seq_along(w))

  p1 <- ggplot2::ggplot(df_w, ggplot2::aes(x = weight)) +
    ggplot2::geom_histogram(bins = 30L, fill = "#a6cee3", colour = "white") +
    ggplot2::geom_vline(xintercept = mean(w), colour = "red",
                        linetype = "dashed", linewidth = 0.9) +
    ggplot2::annotate("text", x = mean(w), y = Inf,
                      label = paste0("Mean = ", round(mean(w), 3)),
                      hjust = -0.1, vjust = 2, size = 3.2, colour = "red") +
    ggplot2::labs(title = "Weight Distribution", x = "Weight", y = "Count") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df_w, ggplot2::aes(x = index, y = weight)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.9, colour = "#1f78b4") +
    ggplot2::geom_hline(yintercept = mean(w), colour = "red", linetype = "dashed") +
    ggplot2::labs(title = "Weights vs Observation Index",
                  x = "Observation", y = "Weight") +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(df_w, ggplot2::aes(sample = weight)) +
    ggplot2::geom_qq(alpha = 0.4, size = 0.9, colour = "#1f78b4") +
    ggplot2::geom_qq_line(colour = "red", linetype = "dashed") +
    ggplot2::labs(title = "Weight Q-Q (vs Normal)",
                  x = "Theoretical", y = "Sample") +
    ggplot2::theme_minimal()

  ess_pct   <- round(diag$effective_sample_size / n * 100, 1)
  n_trimmed <- diag$n_hull_trimmed %||% 0L
  stats_txt <- paste0(
    "ESS: ", round(diag$effective_sample_size, 1), " (", ess_pct, "%)\n",
    "CV:  ", round(diag$weight_cv, 3), "\n",
    "Min: ", round(diag$weight_range[1], 3), "\n",
    "Max: ", round(diag$weight_range[2], 3), "\n",
    "n:   ", n,
    if (n_trimmed > 0) paste0("\nHull-trimmed: ", n_trimmed) else ""
  )
  p4 <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = stats_txt,
                      size = 4, hjust = 0.5, vjust = 0.5, family = "mono") +
    ggplot2::labs(title = "Weight Statistics") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  combined <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = "GPS Weight Diagnostics",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  print(combined)
  invisible(combined)
}
