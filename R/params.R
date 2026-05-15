# =============================================================================
# params.R — parameter extraction and gt summary tables
# =============================================================================

#' Extract marginal distribution parameters from a GPS object
#'
#' @param gps Object of class `"gps_weights"`
#' @return Named list of per-treatment parameter summaries
#' @export
get_marginal_params <- function(gps) {
  d       <- gps$specs$d
  margins <- gps$specs$margins
  result  <- vector("list", d)
  names(result) <- paste0("T", seq_len(d), "_", margins)

  for (j in seq_len(d)) {
    fit <- gps$marginal_fits[[j]]
    result[[j]] <- list(
      distribution    = margins[j],
      coefficients    = fit$coefficients,
      vcov            = fit$vcov,
      fitted_values   = fit$mean,
      scale_parameter = fit$sd,
      df              = fit$df,
      scale_summary   = if (length(unique(fit$sd)) > 1) summary(fit$sd) else fit$sd[1],
      model_object    = fit$model,
      r_squared       = gps$diagnostics$marginal_r_squared[j]
    )
  }
  result
}

#' Marginal model coefficients as a gt table
#'
#' @param gps Object of class `"gps_weights"`
#' @param treatment_index Integer vector selecting treatments (default: all)
#' @param as_gt Return a gt table (default: `TRUE`); `FALSE` returns a
#'   plain `data.frame`
#' @return A gt table or `data.frame`
#' @export
get_coef_table <- function(gps, treatment_index = NULL, as_gt = TRUE) {
  d <- gps$specs$d
  if (is.null(treatment_index)) treatment_index <- seq_len(d)

  rows <- lapply(treatment_index, function(j) {
    fit <- gps$marginal_fits[[j]]
    if (!is.null(fit$model)) {
      s  <- summary(fit$model)
      ct <- if (!is.null(s$coefficients)) {
        s$coefficients
      } else if (!is.null(s$table)) {
        tbl <- s$table[!grepl("Log\\(scale\\)", rownames(s$table)), , drop = FALSE]
        colnames(tbl) <- c("Estimate", "Std. Error", "z value",
                           "Pr(>|z|)")[seq_len(ncol(tbl))]
        tbl
      } else NULL

      if (is.null(ct)) {
        return(data.frame(
          treatment = j, distribution = gps$specs$margins[j],
          coefficient = "(unavailable)",
          estimate = NA_real_, std_error = NA_real_,
          t_value = NA_real_, p_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        treatment    = j, distribution = gps$specs$margins[j],
        coefficient  = rownames(ct),
        estimate     = ct[, 1], std_error = ct[, 2],
        t_value      = ct[, 3], p_value   = ct[, 4],
        stringsAsFactors = FALSE, row.names = NULL
      )
    } else {
      data.frame(treatment = j, distribution = gps$specs$margins[j],
                 coefficient = "user_specified",
                 estimate = NA_real_, std_error = NA_real_,
                 t_value = NA_real_, p_value = NA_real_,
                 stringsAsFactors = FALSE)
    }
  })
  df <- do.call(rbind, rows)
  if (!as_gt) return(df)

  df |>
    gt::gt(rowname_col = "coefficient") |>
    gt::tab_header(
      title    = "Marginal Model Coefficients",
      subtitle = paste0(gps$specs$n, " observations | treatments: ",
                        paste(gps$specs$margins, collapse = ", "))
    ) |>
    gt::cols_label(treatment = "Trt", distribution = "Family",
                   estimate = "Estimate", std_error = "SE",
                   t_value = "t", p_value = "p-value") |>
    gt::fmt_number(columns = c("estimate", "std_error", "t_value"), decimals = 4) |>
    gt::fmt_scientific(columns = "p_value", decimals = 3) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"),
                  locations = gt::cells_column_labels()) |>
    gt::tab_options(table.font.size = "small")
}

#' Extract copula parameters as a list
#'
#' @param gps Object of class `"gps_weights"`
#' @return List with empirical and parametric copula details
#' @export
get_copula_params <- function(gps) {
  make_entry <- function(cop) {
    list(family      = cop$family,
         family_name = .family_name(cop$family),
         parameter1  = cop$par,
         parameter2  = cop$par2,
         tau         = if (cop$family == 0) 0
                       else VineCopula::BiCopPar2Tau(cop$family, cop$par, cop$par2))
  }
  list(method     = "bivariate",
       empirical  = make_entry(gps$copula_empirical),
       parametric = make_entry(gps$copula_parametric))
}

#' Copula parameters as a gt table
#'
#' @param gps Object of class `"gps_weights"`
#' @return A gt table
#' @export
get_copula_table <- function(gps) {
  cp  <- get_copula_params(gps)
  df <- data.frame(
    role        = c("Empirical", "Parametric"),
    family      = c(cp$empirical$family_name,  cp$parametric$family_name),
    parameter1  = c(cp$empirical$parameter1,   cp$parametric$parameter1),
    parameter2  = c(cp$empirical$parameter2,   cp$parametric$parameter2),
    kendall_tau = c(cp$empirical$tau,           cp$parametric$tau),
    stringsAsFactors = FALSE
  )
  df |>
    gt::gt(rowname_col = "role") |>
    gt::tab_header(
      title    = "Fitted Copula Parameters",
      subtitle = paste0("type: ", gps$specs$copula_type,
                        if (gps$diagnostics$residual_copula_used)
                          " [residual copula]" else "")
    ) |>
    gt::cols_label(family = "Family", parameter1 = "\u03b8\u2081",
                   parameter2 = "\u03b8\u2082", kendall_tau = "Kendall \u03c4") |>
    gt::fmt_number(columns = c("parameter1", "parameter2", "kendall_tau"),
                   decimals = 4) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"),
                  locations = gt::cells_column_labels())
}
