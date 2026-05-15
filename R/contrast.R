# =============================================================================
# contrast.R — treatment contrast helper for d >= 2
# =============================================================================

#' Build a treatment contrast matrix for ATE estimation
#'
#' Constructs a matrix of treatment regimes to pass to [estimate_ate()] as
#' `treatment_levels`.  Consecutive row pairs are contrasted.
#'
#' @param data Prepared data list from [prepare_data()] (must have `$y`).
#' @param type One of:
#'   \describe{
#'     \item{`"quantile"`}{All treatments simultaneously at `probs[1]` vs
#'       `probs[2]` quantile (default).  Returns a 2-row matrix.}
#'     \item{`"single_vary"`}{Vary each treatment one at a time between
#'       `probs[1]` and `probs[2]`; hold all others at `hold_at`.  Returns a
#'       2d-row matrix with one high/low pair per treatment.}
#'   }
#' @param vary Integer vector of treatment indices to vary when
#'   `type = "single_vary"`.  Default: all treatments.
#' @param probs Length-2 numeric vector `c(high_prob, low_prob)`.
#'   Default: `c(0.75, 0.25)`.
#' @param hold_at How to hold fixed treatments: `"median"` (default) or
#'   `"mean"`.
#' @return A numeric matrix with named rows and d columns (one per treatment).
#'   Column names match treatment column names from `data$y` or `T1, T2, ...`.
#' @export
make_contrast <- function(data, type = "quantile", vary = NULL,
                           probs = c(0.75, 0.25), hold_at = "median") {
  if (!all(c("y") %in% names(data)))
    stop("`data` must have a `$y` component; use prepare_data() to build it")

  y   <- if (is.matrix(data$y)) data$y else as.matrix(data$y)
  d   <- ncol(y)
  trt <- colnames(y) %||% paste0("T", seq_len(d))

  if (length(probs) != 2 || probs[1] <= probs[2])
    stop("`probs` must be a length-2 vector with probs[1] > probs[2] ",
         "(e.g. c(0.75, 0.25))")

  if (!hold_at %in% c("median", "mean"))
    stop("`hold_at` must be 'median' or 'mean'")

  if (type == "quantile") {
    high <- apply(y, 2, stats::quantile, probs = probs[1])
    low  <- apply(y, 2, stats::quantile, probs = probs[2])
    mat  <- rbind(high = high, low = low)
    colnames(mat) <- trt
    return(mat)
  }

  if (type == "single_vary") {
    if (is.null(vary)) vary <- seq_len(d)
    vary <- as.integer(vary)
    bad  <- vary[vary < 1L | vary > d]
    if (length(bad))
      stop("'vary' contains out-of-range treatment indices: ",
           paste(bad, collapse = ", "))

    hold <- if (hold_at == "median")
              apply(y, 2, stats::median)
            else
              colMeans(y)

    mat_list <- lapply(vary, function(j) {
      h <- rbind(
        replace(hold, j, stats::quantile(y[, j], probs = probs[1])),
        replace(hold, j, stats::quantile(y[, j], probs = probs[2]))
      )
      rownames(h) <- paste0(trt[j], c("_high", "_low"))
      h
    })
    mat <- do.call(rbind, mat_list)
    colnames(mat) <- trt
    return(mat)
  }

  stop("`type` must be 'quantile' or 'single_vary'; got: '", type, "'")
}
