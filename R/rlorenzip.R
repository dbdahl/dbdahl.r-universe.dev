#' Sample from the Lorenz Integer Partition Distribution
#'
#' Generates random samples from the Lorenz integer partition (Lorenz-IP) distribution.
#' This distribution is used to model integer partitions for cluster sizes,
#' where the cumulative ordered cluster sizes follow a specified Lorenz curve pattern.
#'
#' @param n Integer specifying the number of samples to generate.
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#'
#' @inheritParams params
#'
#' @return A matrix with `n` rows and `k` columns (even when `n = 1`), where
#'   each row represents a sampled integer partition. Each element indicates the
#'   number of items assigned to the corresponding cluster, in non-decreasing
#'   order (smallest to largest).
#'
#' @example man/examples/rlorenzip.R
#' @export
rlorenzip <- function(
  n,
  n_items,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  n_clusters = NULL
) {
  .Call(
    .rlorenzip,
    n,
    n_items,
    target,
    concentration,
    log_skew,
    tail_shape,
    buffer_c,
    n_clusters
  )
}

#' Sample from the Chinese Restaurant Process Integer Partition
#'
#' Generates random samples from the integer partition distribution of the Chinese Restaurant Process (CRP)
#'
#' @param n Integer specifying the number of samples to generate.
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#' @param n_clusters Integer specifying the number of clusters in the integer partition.
#' @param discount Numeric value in [0, 1) specifying the discount parameter for
#'   the CRP. Larger values favor more clusters (and relatively smaller cluster
#'   sizes); discount = 0 recovers the Ewens/CRP case.
#'
#' @return A matrix with `n` rows and `k` columns (even when `n = 1`), where
#'   each row represents a sampled integer partition. Each element indicates the
#'   number of items assigned to the corresponding cluster, in non-decreasing
#'   order (smallest to largest).
#'
#' @example man/examples/rcrpip.R
#' @export
rcrpip <- function(n, n_items, n_clusters, discount) {
  .Call(.rcrpip, n, n_items, n_clusters, discount)
}
