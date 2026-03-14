#' Update the Anchor Partition of Correlated Partitions
#'
#' This function updates the common anchor partitions in a collection of
#' shrinkage partition distributions.
#'
#' @param distr A specification of the partition distribution, i.e., an object
#'   of class \code{PartitionDistribution} as returned by, for example, a
#'   function such as \code{\link{CRPPartition}}.
#' @param correlated A list of lists, with each list containing two elements: i.
#'   a partition (in cluster label notation), and ii. a Shrinkage partition
#'   distribution.
#
#' @export
updateAnchorOfCorrelatedPartitions <- function(distr, correlated) {
  if ( ! is.list(correlated) || length(correlated) == 0 ) stop("'correlated' should be a non-empty list.")
  correlated <- lapply(correlated, function(x) {
    x$distr <- distrUnlock(distrClone(x$distr))
    x
  })
  nItems <- NULL
  partition <- NULL
  for ( x in correlated ) {
    if ( ! isTRUE(all.equal(names(x),c("partition","distr"))) ) stop("Each element of the 'correlated' list should be a list having elements 'partition' and 'distr'.")
    if ( ! inherits(x$distr, "ShrinkagePartition") ) stop("All elements of the 'correlated' list should be of class 'ShrinkagePartition'.")
    if ( is.null(partition) ) {
      partition <- x$distr$anchor
    } else {
      if ( ! identical(x$distr$anchor, partition) ) stop("The 'anchor' is not consistent among the elements of the 'correlated' list.")
    }
    if ( is.null(nItems) ) {
      nItems <- length(x$partition)
    }
    if ( ( length(x$partition) != nItems ) || ( x$distr$nItems != nItems ) ) stop("'Inconsistent number of items among the elements of 'list'.")
  }
  if ( distr$nItems != nItems ) stop("'Inconsistent number of items among the elements of 'list'.")
  for ( i in seq_len(nItems) ) {
    labels <- unique(partition[-i])
    labels <- c(labels, max(labels)+1)
    log_probs <- numeric(length(labels))
    for ( j in seq_along(labels) ) {
      candidate <- partition
      candidate[i] <- labels[j]
      candidatesCorrelated <- lapply(correlated, function(x) { x$distr$anchor <- candidate; x })
      log_probs[j] <- prPartition(distr, candidate, log=TRUE) + sum(sapply(candidatesCorrelated, function(x) prPartition(x$distr, x$partition, log=TRUE)))
    }
    m <- max(log_probs)
    weights <- exp(log_probs - m)
    partition[i] <- sample(labels, 1, prob=weights)
  }
  lapply(correlated, function(x) {
    x$distr$anchor <- canonicalForm(partition)
    distrLock(x$distr)
    x
  })
}
