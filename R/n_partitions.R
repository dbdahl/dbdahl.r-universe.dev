#' Count Set Partitions
#'
#' Counts the number of set partitions (ways to partition distinguishable items into
#' non-empty, unordered groups). The function can compute Bell numbers, Stirling numbers
#' of the second kind, or multinomial coefficients for specific integer partitions.
#'
#' @param n_clusters Integer specifying the number of clusters (groups). If `NULL`,
#'   the function returns the Bell number (total number of set partitions).
#' @param integer_partition Integer vector specifying the sizes of each cluster.
#'   If provided, `n_clusters` must also be specified and must equal the length of
#'   this vector. The function then returns the number of set partitions with
#'   exactly these cluster sizes.
#'
#' @inheritParams params
#'
#' @return A numeric value. If `log = FALSE`, the count (as a whole number).
#'   If `log = TRUE`, the natural logarithm of the count.
#'
#' @details
#' This function computes three related quantities depending on the arguments:
#'
#' \describe{
#'   \item{Bell number \eqn{B(n)}}{When only `n_items` is specified.
#'     The total number of ways to partition \eqn{n} distinguishable items into
#'     any number of non-empty groups.}
#'   \item{Stirling number of the second kind \eqn{S(n,k)}}{When `n_items` and
#'     `n_clusters` are specified. The number of ways to partition \eqn{n}
#'     distinguishable items into exactly \eqn{k} non-empty groups.}
#'   \item{Multinomial coefficient for set partitions}{When all three arguments
#'     are specified. The number of ways to partition \eqn{n} distinguishable items
#'     into groups of specified sizes. This equals
#'     \eqn{n! / (\prod_{i} m_i! \cdot \prod_{j} r_j!)}
#'     where \eqn{m_i} are the group sizes and \eqn{r_j} are the multiplicities
#'     of each distinct group size.}
#' }
#'
#' @note
#' Bell numbers and Stirling numbers grow very rapidly. For example, \eqn{B(25) \approx 4.6 \times 10^{18}}.
#' For very large values (beyond \eqn{10^{15}}), consider using `log = TRUE` for better precision.
#'
#' @example man/examples/n_partitions.R
#' @export
n_partitions <- function(
  n_items,
  n_clusters = NULL,
  integer_partition = NULL,
  log = FALSE
) {
  .Call(.n_partitions, n_items, n_clusters, integer_partition, log)
}
