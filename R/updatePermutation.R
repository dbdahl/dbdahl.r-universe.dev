#' Update the Permutation Parameter of a Partition Distribution
#'
#' This function updates the permutation parameter of a partition distribution (
#' e.g, \code{\link{CRPPartition}} and \code{\link{ShrinkagePartition}})
#' using a Metropolis proposal.  The prior is uniform on all permutations.
#'
#' @inheritParams nealAlgorithm3
#' @param k Among all \eqn{n} elements in the permutation of \code{distr}, the
#'   number of items to shuffle when forming the Metropolis proposal. The other
#'   elements stay fixed. This value must be at least \code{2} and defaults to
#'   all of the elements.
#'
#' @return A list containing the potentially updated partition distribution and
#'   a logical named \code{accepted} (indicating whether the proposal was
#'   accepted).
#'
#' @seealso \code{\link{nealAlgorithm3}}, \code{\link{nealAlgorithm8}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @export
#'
updatePermutation <- function(distr, partition, k=NULL) {
  if ( is.null(distr[['permutation']]) ) return(list(distr=distr, accepted=FALSE))
  n <- length(distr[['permutation']])
  if ( is.null(k) ) k <- n
  if ( ( n == 1 ) || ( k < 2 ) ) return(list(distr=distr, accepted=FALSE))
  if ( k > n ) k <- n
  proposalDistr <- distrClone(distr)
  if ( k == n ) {
    proposalDistr$.permutation <- sample(distr$.permutation)
  } else {
    activePositions <- sample(n, k)
    proposalDistr$.permutation[activePositions] <- sample(proposalDistr$.permutation[activePositions])
  }
  proposalDistr$permutation <- proposalDistr$.permutation + 1L
  distrLock(proposalDistr)
  logMHRatio <- prPartition(proposalDistr, partition, log=TRUE) - prPartition(distr, partition, log=TRUE)
  if ( log(runif(1)) < logMHRatio ) list(distr=proposalDistr, accepted=TRUE)
  else list(distr=distr, accepted=FALSE)
}
