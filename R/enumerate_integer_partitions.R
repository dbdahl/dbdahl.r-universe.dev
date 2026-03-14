#' Enumerate All Integer Partitions
#'
#' Enumerates all possible integer partitions of `n_items` into exactly `n_clusters`
#' positive parts in non-decreasing order. This is useful for computing exact
#' probabilities over the entire support of the Lorenz-IP distribution.
#'
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`). Use `n_items = 0` only with
#'   `n_clusters = 0`, which returns a single empty integer partition.
#' @param n_clusters Integer specifying the number of clusters (parts) in each integer partition.
#'   Use `n_clusters = 0` only when `n_items = 0`.
#'
#' @return A matrix where each row represents a unique integer partition. The matrix
#'   has `n_clusters` columns, with elements in non-decreasing order within each row.
#'   The number of rows equals the number of integer partitions of `n_items` into exactly
#'   `n_clusters` positive parts.
#'
#' @details
#' An integer partition of \eqn{n} into \eqn{k} parts is a sequence
#' \eqn{(x_1, x_2, \ldots, x_k)} of positive integers satisfying:
#' \itemize{
#'   \item \eqn{x_1 + x_2 + \cdots + x_k = n} (parts sum to n)
#'   \item \eqn{x_1 \le x_2 \le \cdots \le x_k} (non-decreasing order)
#'   \item \eqn{x_i \ge 1} for all \eqn{i} (all parts are positive)
#' }
#'
#' The number of such integer partitions grows rapidly with `n_items` and `n_clusters`.
#' For large values, this function may consume significant memory and time.
#'
#' @example man/examples/enumerate_integer_partitions.R
#' @export
enumerate_integer_partitions <- function(n_items, n_clusters) {
  .Call(.enumerate_integer_partitions, n_items, n_clusters)
}
