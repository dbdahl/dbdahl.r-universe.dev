#' Monte Carlo Co-clustering Probability
#'
#' Estimates the co-clustering probability c_n by Monte Carlo sampling from
#' the Lorenz-IP distribution. The estimate is the average of
#' sum_h X_h (X_h - 1) / (n_items (n_items - 1)) over sampled set partitions.
#'
#' @param n_samples Integer specifying the number of Monte Carlo samples.
#' @param n_items Integer specifying the total number of items in the set partition.
#' @inheritParams params
#'
#' @return A numeric scalar giving the Monte Carlo estimate.
#'
#' @example man/examples/cocluster_mc.R
#' @export
cocluster_mc <- function(
  n_samples,
  n_items,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  n_clusters = NULL
) {
  .Call(
    .cocluster_mc,
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
