#' Update the Shrinkage Parameter of a Partition Distribution
#'
#' This function updates the shrinkage parameter of a partition distribution (e.g,
#' \code{\link{ShrinkagePartition}} and \code{\link{LocationScalePartition}}) using a
#' Gaussian random walk.  The prior can be specified by the user.
#'
#' @inheritParams updatePermutation
#' @inheritParams samplePartition
#' @param w The initial size of the interval for the slice sampler.
#' @param logPrior A function of a scalar returning the log of the density of
#'   the prior on the shrinkage parameter.
#'
#' @return An updated partition distribution
#'
#' @seealso \code{\link{nealAlgorithm3}}, \code{\link{nealAlgorithm8}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @importFrom stats dbeta
#' @export
#'
updateShrinkage <- function(distr, partition, randomizeShrinkage=c("fixed","common","cluster","idiosyncratic")[1], w=1, logPrior=function(shrinkage) dbeta(shrinkage/5, 2, 2, log=TRUE)) {
  if ( is.null(distr[['shrinkage']]) ) return(distr)
  if ( randomizeShrinkage == "fixed" ) return(distr)
  if ( inherits(distr,"CenteredPartition") ) return(distr)
  if ( randomizeShrinkage == "common" ) {
    unique <- unique(distr[['shrinkage']])
    if ( ( length(unique) > 2 ) || ( length(unique) == 2 && ( ! 0.0 %in% unique ) ) || ( ( length(unique) == 1 ) && ( unique == 0.0 ) ) ) {
      stop("There should be exactly one non-zero unique shrinkage value, so 'common' is not an appropriate value for 'randomizeShrinkage'.")
    }
    distrUnlock(distr)
    zeros <- distr[['shrinkage']] == 0.0
    lf <- function(shrinkage)  {
      if ( shrinkage < 0.0 ) -Inf
      else {
        distr[['shrinkage']][!zeros] <- shrinkage
        distr[['shrinkage']][ zeros] <- 0.0
        prPartition(distr, partition, log=TRUE) + logPrior(shrinkage)
      }
    }
    shrinkage <- .Call(.slice_sampler, max(unique), lf, w, Inf, TRUE)
    distr[['shrinkage']][!zeros] <- shrinkage
    distr[['shrinkage']][ zeros] <- 0.0
    distrLock(distr)
    return(distr)
  }
  if ( randomizeShrinkage == "cluster" ) {
    distrUnlock(distr)
    for ( k in unique(distr[['baselinePartition']]) ) {
      items <- distr[['baselinePartition']] == k
      y <- distr[['shrinkage']][items]
      x <- y == 0.0
      if ( all(x) ) next
      if ( any(x) ) stop(paste0("Cluster ",k," in baseline partition has some items with zero shrinkage, but not all items."))
      currentShrinkage <- unique(y)
      lf <- function(shrinkage)  {
        if ( shrinkage < 0.0 ) -Inf
        else {
          distr[['shrinkage']][items] <- shrinkage
          prPartition(distr, partition, log=TRUE) + logPrior(shrinkage)
        }
      }
      distr[['shrinkage']][items] <- .Call(.slice_sampler, currentShrinkage, lf, w, Inf, TRUE)
    }
    distrLock(distr)
    return(distr)
  }
  if ( randomizeShrinkage == "idiosyncratic" ) {
    distrUnlock(distr)
    for ( i in seq_len(distr[['nItems']]) ) {
      if ( distr[['shrinkage']][i] == 0.0 ) next
      lf <- function(shrinkage)  {
        if ( shrinkage < 0.0 ) -Inf
        else {
          distr[['shrinkage']][i] <- shrinkage
          prPartition(distr, partition, log=TRUE) + logPrior(shrinkage)
        }
      }
      distr[['shrinkage']][i] <- .Call(.slice_sampler, distr[['shrinkage']][i], lf, w, Inf, TRUE)
    }
    distrLock(distr)
    return(distr)
  }
  stop("Unrecognized value for 'randomShrinkage'")
}
