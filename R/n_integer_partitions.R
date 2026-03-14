#' Count Integer Partitions
#'
#' Counts the number of integer partitions (ways to write a positive integer as
#' a sum of positive integers, where order does not matter). The function can
#' compute either the integer partition function p(n) or the number of integer
#' partitions into exactly k parts.
#'
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#' @param n_clusters Integer specifying the number of parts. If `NULL` (default),
#'   the function returns the total number of integer partitions (integer partition function p(n)).
#'   If specified, returns the number of integer partitions into exactly that many parts.
#'
#' @inheritParams params
#'
#' @return A numeric value. If `log = FALSE`, the count (as a whole number).
#'   If `log = TRUE`, the natural logarithm of the count.
#'
#' @details
#' An integer partition of n is a way of writing n as a sum of positive integers,
#' where the order of the summands does not matter. For example, the integer partitions
#' of 5 are: 5, 1:4, 2:3, 1:1:3, 1:2:2, 1:1:1:2, and 1:1:1:1:1, giving p(5) = 7.
#'
#' This function computes two related quantities:
#'
#' \describe{
#'   \item{Integer partition function p(n)}{When only `n_items` is specified.
#'     The total number of integer partitions of n.}
#'   \item{Integer partitions into k parts p(n, k)}{When `n_items` and `n_clusters` are specified.
#'     The number of ways to write n as a sum of exactly k positive integers.}
#' }
#'
#' Note: Integer partitions differ from set partitions (see \code{\link{n_partitions}}).
#' Integer partitions count ways to decompose a number; set partitions count ways
#' to group distinguishable items. For example, p(5) = 7 but B(5) = 52.
#'
#' @seealso \code{\link{enumerate_integer_partitions}} for listing all integer partitions,
#'   \code{\link{n_partitions}} for counting set partitions.
#'
#' @example man/examples/n_integer_partitions.R
#' @export
n_integer_partitions <- function(
  n_items,
  n_clusters = NULL,
  log = FALSE
) {
  .Call(.n_integer_partitions, n_items, n_clusters, log)
}
