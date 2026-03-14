#' Splintered Clustering
#'
#' This function implements the splintered clustering algorithm.
#'
#' @param splinters A nonempty list containing lists.  Each list must be of
#'   length two containing, as its first element, a numeric vector named
#'   \dQuote{indices} and, as its second element, a numeric matrix named
#'   \dQuote{partitions} with number of rows equaling the length of
#'   \dQuote{indices} and at least one column.
#' @param maxNClusters The maximum number of clusters that can be considered by
#'   the optimization algorithm, which has important implications for the
#'   interpretability of the resulting clustering and can greatly influence the
#'   RAM needed for the optimization algorithm. If the supplied value is zero,
#'   the optimization is constrained by the maximum number of clusters among the
#'   clusterings in \code{splinters}.
#' @param nRuns The number of runs to try, although the actual number is a
#'   multiple of the number of cores specified by the \code{nCores} argument.
#' @param useVI Use the (generalized) variation of information loss?  If not,
#'   Binder loss is used.
#' @param a The argument \code{a} is a nonnegative scalar in \eqn{[0,2]} giving
#'   (for Binder loss) the cost of placing two items in separate clusters when in
#'   truth they belong to the same cluster. Without loss of generality, the cost
#'   (under Binder loss) of placing two items in the same cluster when in truth
#'   they belong to separate clusters is fixed \code{2-a}. For VI, \code{a} has a
#'   similar interpretation, although is not a unit cost. See Dahl, Johnson,
#'   Müller (2022).
#' @param nCores The number of CPU cores to use, i.e., the number of
#'   simultaneous runs at any given time. A value of zero indicates to use all
#'   cores on the system.
#'
#' @return A numeric vector giving an estimated clustering.
#'
#' @export
#' @examples
#' cat("To do.\n")
#'
splinclust <- function(splinters, maxNClusters=0, nRuns=16, useVI=TRUE, a=1.0, nCores=0) {
  out <- mk_splinter_ptrs(splinters)
  .Call(.splintered_clustering, out$ptrs, out$nItems, maxNClusters, nRuns, useVI, a, nCores)
}
