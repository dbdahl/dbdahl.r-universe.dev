#' Univariate Slice Sampler using Doubling Method and Shrinkage
#'
#' This function implements a slice sampler using the doubling and shrinkage
#' methods.
#'
#' @param x A numeric scalar giving the current value.
#' @param target A function taking a numeric scalar and returning the target function evaluated at that value.
#' @param on_log_scale A logical scalar indicating whether \code{target} is on the log scale.
#' @param width A numeric scalar giving the estimated slice width.
#'
#' @return A numeric vector.
#' @export
#' @examples
#' n_samples <- 1000
#' sum <- 0.0
#' state <- 0.5
#' target <- function(x) if (0 < x && x < 1) x else 0
#' for (i in seq_len(n_samples)) {
#'   state <- univariate_slice_sampler_doubling(state, target, FALSE, 1)$value
#'   sum <- sum + state
#' }
#' sum / n_samples
#' 
#' @references 
#' Radford M Neal. (2003). Slice Sampling.
#' *The Annals of Statistics*. Vol 31, No 3, 705-767.
#' 
univariate_slice_sampler_doubling <- function(x, target, on_log_scale, width) {
  if (!is.numeric(x) || length(x) != 1) stop("'x' must be a scalar.")
  if (!is.function(target)) stop("'target' must be function.")
  if (!is.logical(on_log_scale)) stop("'on_log_scale' should be a logical.")
  if (!is.numeric(x) || length(width) != 1 || width <= 0.0) stop("'width' should be a positive scalar.")
  .Call(.univariate_slice_sampler_doubling, x, target, on_log_scale, width)
}
