# =============================================================================
# preprocess.R — user-facing data preparation from a data frame
# =============================================================================

#' Prepare data for GPS estimation from a data frame
#'
#' Converts a standard data frame into the internal list format expected by
#' `estimate_gps()`, validating inputs and constructing design matrices along
#' the way.
#'
#' @param data A `data.frame` (or `tibble`) containing all variables.
#' @param treatments Character vector of column names for the treatment
#'   variables.  Must identify exactly `d` columns (currently d = 2 for the
#'   copula GPS estimator).
#' @param covariates Character vector of column names to use as covariates in
#'   the treatment models.  All selected columns are included in every treatment
#'   model unless per-treatment formulas are supplied via `formulas`.
#' @param outcome Optional character scalar naming the outcome column.  When
#'   supplied the outcome is extracted and stored in `$outcome`; it is **not**
#'   included in the treatment design matrices.
#' @param lod Optional numeric vector of limit-of-detection thresholds, one per
#'   treatment (recycled to length `d`).  Observations with treatment values
#'   below their respective LOD are flagged as left-censored.
#' @param formulas Optional named list of length `d` containing one-sided R
#'   formulae (e.g. `~ age + sex`) specifying the covariates for each treatment
#'   model individually.  Names should match `treatments`.  When `NULL` (the
#'   default) the same set of `covariates` is used for all treatments, with an
#'   intercept prepended.
#' @param scale_covariates Logical.  If `TRUE` (default), numeric covariates
#'   are mean-centred and scaled to unit SD before constructing design matrices.
#'   Factor/character columns are left as-is.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`y`}{n × d numeric treatment matrix}
#'     \item{`x`}{List of d design matrices (one per treatment), each n × p
#'       with an intercept column}
#'     \item{`outcome`}{Numeric vector of length n, or `NULL` if not supplied}
#'     \item{`covariates`}{The covariate data frame (numeric columns possibly
#'       scaled) used to build the design matrices}
#'     \item{`lod`}{Numeric vector of length d, or `NULL`}
#'     \item{`treatment_names`}{Character vector of treatment column names}
#'     \item{`covariate_names`}{Character vector of covariate column names}
#'   }
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   T1 = rnorm(200, 2), T2 = rnorm(200, -1),
#'   age = rnorm(200, 50, 10), sex = rbinom(200, 1, 0.5),
#'   Y   = rnorm(200)
#' )
#' d <- prepare_data(df,
#'   treatments = c("T1", "T2"),
#'   covariates = c("age", "sex"),
#'   outcome    = "Y"
#' )
#' str(d)
prepare_data <- function(data,
                          treatments,
                          covariates,
                          outcome          = NULL,
                          lod              = NULL,
                          formulas         = NULL,
                          scale_covariates = TRUE) {

  if (!is.data.frame(data))
    stop("`data` must be a data.frame or tibble")

  # ---- Validate column names --------------------------------------------------
  all_cols <- names(data)

  bad_trt <- setdiff(treatments, all_cols)
  if (length(bad_trt))
    stop("Treatment columns not found in data: ", paste(bad_trt, collapse = ", "))

  bad_cov <- setdiff(covariates, all_cols)
  if (length(bad_cov))
    stop("Covariate columns not found in data: ", paste(bad_cov, collapse = ", "))

  if (!is.null(outcome)) {
    if (!outcome %in% all_cols)
      stop("Outcome column '", outcome, "' not found in data")
    if (outcome %in% treatments)
      stop("Outcome column cannot also be a treatment")
  }

  overlap <- intersect(treatments, covariates)
  if (length(overlap))
    stop("Columns appear in both treatments and covariates: ",
         paste(overlap, collapse = ", "))

  d <- length(treatments)

  # ---- Extract treatment matrix -----------------------------------------------
  y <- as.matrix(data[, treatments, drop = FALSE])
  storage.mode(y) <- "double"

  # ---- Check for missing values -----------------------------------------------
  if (any(is.na(y)))
    stop("Treatment matrix contains NA values; impute or remove before calling prepare_data()")

  # ---- Extract and optionally scale covariates --------------------------------
  cov_df <- data[, covariates, drop = FALSE]

  if (scale_covariates) {
    num_cols <- sapply(cov_df, is.numeric)
    if (any(num_cols)) {
      cov_df[num_cols] <- lapply(cov_df[num_cols], function(v) {
        s <- stats::sd(v, na.rm = TRUE)
        if (is.na(s) || s < .Machine$double.eps) v - mean(v, na.rm = TRUE)
        else (v - mean(v, na.rm = TRUE)) / s
      })
    }
  }

  # ---- Build design matrices --------------------------------------------------
  if (!is.null(formulas)) {
    if (!is.list(formulas) || length(formulas) != d)
      stop("`formulas` must be a list of length d = ", d)
    x_list <- lapply(formulas, function(f) {
      stats::model.matrix(f, data = cov_df)
    })
  } else {
    cov_mat <- stats::model.matrix(~ ., data = cov_df)   # includes intercept
    x_list  <- replicate(d, cov_mat, simplify = FALSE)
  }

  # ---- LOD -------------------------------------------------------------------
  if (!is.null(lod)) {
    lod <- as.numeric(lod)
    if (length(lod) == 1L) lod <- rep(lod, d)
    if (length(lod) != d)
      stop("`lod` must have length 1 or d = ", d)
  }

  # ---- Outcome ---------------------------------------------------------------
  y_out <- if (!is.null(outcome)) as.numeric(data[[outcome]]) else NULL

  list(
    y               = y,
    x               = x_list,
    outcome         = y_out,
    covariates      = cov_df,
    lod             = lod,
    treatment_names = treatments,
    covariate_names = covariates
  )
}

#' Create design matrices from a covariate set (low-level helper)
#'
#' @param covariates Matrix or data frame of covariates
#' @param d Number of treatments (required when `formulas = NULL`)
#' @param formulas Optional list of d one-sided formulae, one per treatment
#' @param include_intercept Prepend intercept column (default: `TRUE`)
#' @return List of d design matrices
#' @export
create_design_matrix <- function(covariates, d = NULL, formulas = NULL,
                                  include_intercept = TRUE) {
  if (is.vector(covariates))    covariates <- matrix(covariates, ncol = 1L)
  if (is.data.frame(covariates)) covariates <- as.matrix(covariates)

  if (is.null(formulas)) {
    if (is.null(d) || !is.numeric(d) || d < 1L)
      stop("Specify d (number of treatments) when formulas is NULL")
    dm <- if (include_intercept) cbind(1, covariates) else covariates
    replicate(as.integer(d), dm, simplify = FALSE)
  } else {
    cov_df <- as.data.frame(covariates)
    if (is.null(colnames(cov_df)))
      colnames(cov_df) <- paste0("X", seq_len(ncol(cov_df)))
    lapply(formulas, function(f) stats::model.matrix(f, data = cov_df))
  }
}
