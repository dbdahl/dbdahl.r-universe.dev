#' Monte Carlo Downsampling Diagnostics
#'
#' Approximates total variation and KL divergence between the size-`n` Lorenz
#' partition distribution on set partitions `p_n` and the downsampled law `q_n` obtained by deleting
#' one item from the size-`n + 1` model.
#'
#' @param k_weights Numeric vector of non-negative weights for the prior on `k`.
#'   Entry `j` is proportional to Pr(`k = j`). Weights are renormalized after
#'   truncation to `k <= n_items` (and `k <= n_items + 1` for the downsampled
#'   model). If `log_k_weights = TRUE`, interpret `k_weights` as log weights or
#'   log probabilities.
#' @param log_k_weights Logical; if `TRUE`, interpret `k_weights` as log weights
#'   or log probabilities.
#' @param n_samples Integer specifying the number of Monte Carlo samples used in
#'   the approximation.
#'
#' @inheritParams params
#'
#' @return A list with elements `tv` and `kl` containing the Monte Carlo estimates.
#'
#' @details
#' If `target`, `concentration`, `log_skew`, or `tail_shape` are supplied as
#' vectors of fixed length, they imply a fixed number of clusters `k`; in that
#' case `k_weights` must assign all mass to that `k`.
#'
#' @example man/examples/downsample.R
#' @export
downsample <- function(
  n_items,
  k_weights,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  n_samples = 1000L,
  log_k_weights = FALSE
) {
  .Call(
    .downsample,
    n_items,
    k_weights,
    target,
    concentration,
    log_skew,
    tail_shape,
    buffer_c,
    n_samples,
    log_k_weights
  )
}
