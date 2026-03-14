#' Probability Mass Function of the Lorenz Integer Partition Distribution
#'
#' Evaluates the probability mass function of the Lorenz integer partition (Lorenz-IP) distribution.
#' This distribution is used to model integer partitions for cluster sizes,
#' where the cumulative ordered cluster sizes follow a specified Lorenz curve pattern.
#'
#' @param x A matrix with `n` rows and `k` columns, where each row represents an
#'   integer partition. Each element indicates the number of items assigned to the corresponding
#'   cluster. Integer partitions are sorted internally to non-decreasing order; any ordering is accepted.
#'   A vector of length `k` is also permissible.
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#'
#' @inheritParams params
#'
#' @return A numeric vector of length `n`. Invalid partitions (wrong sum,
#'   non-positive entries, wrong length) return 0 (or `-Inf` when `log = TRUE`),
#'   following the convention of base R density functions such as [dpois()].
#'
#' @example man/examples/dlorenzip.R
#' @export
dlorenzip <- function(
  x,
  n_items,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  log = FALSE,
  n_clusters = NULL
) {
  .Call(
    .dlorenzip,
    x,
    n_items,
    target,
    concentration,
    log_skew,
    tail_shape,
    buffer_c,
    log,
    n_clusters
  )
}
