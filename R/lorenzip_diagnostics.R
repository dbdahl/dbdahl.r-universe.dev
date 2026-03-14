#' Buffered Mean Diagnostics for Lorenz-IP Construction
#'
#' Computes per-step diagnostics for the buffered mean rule, either along a given
#' integer partition path or via Monte Carlo sampling.
#'
#' @param partition Integer vector of length `k` giving an integer partition.
#'   If provided, path diagnostics are computed. The integer partition is sorted
#'   to non-decreasing order internally.
#' @param n_samples Integer specifying the number of Monte Carlo samples. If provided
#'   (and `partition` is `NULL`), Monte Carlo diagnostics are computed.
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#' @inheritParams params
#'
#' @return A list with elements:
#' \describe{
#'   \item{steps}{A data frame with per-step bounds, buffered means, and flags,
#'     including `mu_star`, `mu_used`, `clamped`, and kernel-parameter diagnostics.
#'     For path diagnostics, columns include `log_skew`, `log_skew_eff`,
#'     `log_skew_adjusted` (and `tail_shape` for TaDPoLe); for Monte Carlo diagnostics,
#'     columns include the baseline `log_skew` (and `tail_shape` for TaDPoLe). For
#'     TiDaL, `log_skew` is the log-skew parameter in the paper. For TaDPoLe,
#'     `tail_shape` is the tail-shape parameter and is not adjusted by the
#'     buffered-mean rule.}
#'   \item{summary}{A list with `any_clamped` and `max_abs_delta` for path diagnostics,
#'     or `any_clamped_rate`, `mean_max_abs_delta`, `max_abs_delta`, and `n_samples`
#'     for Monte Carlo diagnostics.}
#'   \item{assumptions}{A list reporting strict-interior and first-step support
#'     diagnostics from the paper.}
#'   \item{sample_max_abs_delta}{For Monte Carlo diagnostics, a numeric vector of
#'     per-sample maxima.}
#' }
#'
#' @example man/examples/lorenzip_diagnostics.R
#' @export
lorenzip_diagnostics <- function(
  partition = NULL,
  n_samples = NULL,
  n_items,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  n_clusters = NULL
) {
  if (!is.null(partition) && !is.null(n_samples)) {
    stop("Provide either partition or n_samples, not both.", call. = FALSE)
  }
  if (is.null(partition) && is.null(n_samples)) {
    stop("Provide partition or n_samples.", call. = FALSE)
  }

  if (!is.null(partition)) {
    .Call(
      .lorenzip_diagnostics_path,
      partition,
      n_items,
      target,
      concentration,
      log_skew,
      tail_shape,
      buffer_c,
      n_clusters
    )
  } else {
    .Call(
      .lorenzip_diagnostics_mc,
      n_samples,
      n_items,
      target,
      concentration,
      log_skew,
      tail_shape,
      buffer_c,
      n_clusters
    )
  }
}
