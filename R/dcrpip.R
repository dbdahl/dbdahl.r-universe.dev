#' Probability Mass Function of the CRP Integer Partition Distribution
#'
#' Evaluates the probability mass function (PMF) for integer partitions under
#' the two-parameter Chinese Restaurant Process (CRP), conditional on the number
#' of clusters.
#'
#' @param x A matrix with `n` rows and `k` columns, where each row represents an
#'   integer partition. Each element indicates the number of items assigned to the
#'   corresponding cluster, in non-decreasing order (smallest to largest). For
#'   convenience, decreasing-order integer partitions are accepted and are sorted internally.
#'   A vector of length `k` is also permissible.
#' @param n_items Integer specifying the total number of items being partitioned
#'   (the integer partition sums to `n_items`).
#' @param n_clusters Integer specifying the number of clusters in the integer partition.
#' @param discount Numeric value in `[0, 1)` specifying the discount parameter.
#' @param log Logical; if `TRUE`, returns the natural logarithm. Default is `FALSE`.
#'
#' @return A numeric vector of length `n`.
#'
#' @details
#' For an integer partition with sizes `n_1, ..., n_k` and multiplicities `r_j` of equal sizes,
#' the conditional PMF is
#' \deqn{P(n_1, ..., n_k | n, k, \sigma) = \frac{n!}{\prod_i n_i! \prod_j r_j!}
#'   \frac{\prod_i (1-\sigma)_{n_i-1}}{S_\sigma(n, k)}}{P = [n! / (prod n_i! prod r_j!)] * prod (1-sigma)_{n_i-1} / S_sigma(n, k)}
#'
#' where `S_\sigma(n, k)` is the generalized Stirling number. Invalid integer partitions
#' return 0 (or `-Inf` if `log = TRUE`).
#'
#' @seealso [rcrpip()] for sampling from this distribution, [dcrpk()] for `P(K = k)`.
#'
#' @example man/examples/dcrpip.R
#' @export
dcrpip <- function(x, n_items, n_clusters, discount, log = FALSE) {
  .Call(.dcrpip, x, n_items, n_clusters, discount, log)
}
