#' Enumerate All Set Partitions
#'
#' Enumerates all possible set partitions of `n_items` labeled items into exactly
#' `n_clusters` non-empty, non-overlapping groups (clusters). This is useful for
#' computing exact probabilities over the entire support of distributions over
#' set partitions.
#'
#' @param n_items Integer specifying the total number of items to partition.
#'   Items are labeled 1, 2, ..., n_items. Use `n_items = 0` to return the single
#'   empty partition (with zero clusters).
#' @param n_clusters Integer specifying the number of clusters in each partition.
#'   If `NULL` (the default), all partitions for any number of clusters (1 to n_items)
#'   are returned, with the total count equal to the Bell number B(n_items). Use
#'   `n_clusters = 0` only when `n_items = 0`.
#'
#' @return A matrix in cluster label notation where each row represents a unique
#'   set partition. The matrix has `n_items` columns, where column \eqn{j} contains
#'   the cluster label (an integer from 1 to the number of clusters) for item \eqn{j}.
#'
#'   When `n_clusters` is specified, the number of rows equals the Stirling number
#'   of the second kind S(n_items, n_clusters). When `n_clusters` is `NULL`, the
#'   number of rows equals the Bell number B(n_items).
#'
#'   Cluster labels are assigned in order of first appearance: item 1 is always
#'   in cluster 1, and new cluster labels are introduced sequentially as items
#'   are assigned to previously unused clusters.
#'
#' @details
#' A set partition of \eqn{n} labeled items into \eqn{k} clusters is a collection of
#' \eqn{k} non-empty, pairwise disjoint subsets whose union is \eqn{\{1, 2, \ldots, n\}}.
#'
#' Unlike integer partitions (see \code{\link{enumerate_integer_partitions}}), set partitions
#' consider items as distinguishable. For example, partitioning \eqn{\{1, 2, 3\}} into 2 clusters
#' yields 3 distinct partitions, represented in cluster label notation as:
#' \itemize{
#'   \item \code{[1, 1, 2]} - items 1 and 2 in cluster 1, item 3 in cluster 2
#'   \item \code{[1, 2, 1]} - items 1 and 3 in cluster 1, item 2 in cluster 2
#'   \item \code{[1, 2, 2]} - item 1 in cluster 1, items 2 and 3 in cluster 2
#' }
#'
#' The number of such partitions grows very rapidly with `n_items` and `n_clusters`.
#' For large values, this function may consume significant memory and time.
#'
#' @seealso \code{\link{enumerate_integer_partitions}} for enumerating integer partitions,
#'   \code{\link{n_partitions}} for counting set partitions without enumeration.
#'
#' @example man/examples/enumerate_partitions.R
#' @export
enumerate_partitions <- function(n_items, n_clusters = NULL) {
  .Call(.enumerate_partitions, n_items, n_clusters)
}
