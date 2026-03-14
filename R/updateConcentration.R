#' Update the Concentration Parameter of a Partition Distribution
#'
#' This function updates the concentration parameter of a partition distribution (e.g,
#' \code{\link{CRPPartition}} and \code{\link{ShrinkagePartition}}) using a
#' Gaussian random walk.  The prior can be specified by the user.
#'
#' @inheritParams updateShrinkage
#' @param rwsd Standard deviation of the Gaussian random walk proposal density
#'   for the Metropolis algorithm.
#' @param logPriorDensity A function taking a numeric scalar and returning the
#'   log of the prior density evaluated at the value.
#'
#' @return A list containing the (potentially updated) partition distribution,
#'   the (potentially new value for the) concentration, and a logical named
#'   \code{accepted} (indicating whether the proposal was accepted).
#'
#' @seealso \code{\link{nealAlgorithm3}}, \code{\link{nealAlgorithm8}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @importFrom stats dunif dnorm runif rnorm
#' @export
#'
updateConcentration <- function(distr, partition, rwsd=0.5,
                       logPriorDensity=function(w) dnorm(w, 1, 0.5, log=TRUE)) {
  UseMethod("updateConcentration")
}

#' @export
#'
updateConcentration.default <- function(distr, partition, rwsd=0.5,
                               logPriorDensity=function(w) dnorm(w, 1, 0.5, log=TRUE)) {
  if ( is.null(distr[['concentration']]) || rwsd <= 0 || inherits(distr,"CenteredPartition") ) {
    return(list(distr=distr, accepted=FALSE))
  }
  oldM <- distr[['concentration']]
  newM <- rnorm(1, mean=oldM, sd=rwsd)
  discount <- distr[['discount']]
  if ( is.null(discount) ) discount <- 0.0
  if  ( newM <= -discount ) return(list(distr=distr, concentration=oldM, accepted=FALSE))
  proposalDistr <- distrClone(distr)
  proposalDistr$concentration <- newM
  distrLock(proposalDistr)
  logMHRatio <- prPartition(proposalDistr, partition, log=TRUE) + logPriorDensity(newM) - prPartition(distr, partition, log=TRUE) - logPriorDensity(oldM)
  if ( log(runif(1)) < logMHRatio ) list(distr=proposalDistr, concentration=newM, accepted=TRUE)
  else list(distr=distr, concentration=oldM, accepted=FALSE)
}

#' @export
#'
updateConcentration.ShrinkagePartition <- function(distr, partition, rwsd=0.5,
                                          logPriorDensity=function(w) dnorm(w, 1, 0.5, log=TRUE)) {
  if ( is.null(distr$baseline[['concentration']]) || rwsd <= 0 || inherits(distr$baseline,"CenteredPartition") ) {
    return(list(distr=distr, accepted=FALSE))
  }
  oldM <- distr$baseline[['concentration']]
  newM <- rnorm(1, mean=oldM, sd=rwsd)
  discount <- distr$baseline[['discount']]
  if ( is.null(discount) ) discount <- 0.0
  if  ( newM <= -discount ) return(list(distr=distr, concentration=oldM, accepted=FALSE))
  proposalDistr <- distrClone(distr)
  proposalDistr$baseline <- distrClone(proposalDistr$baseline)
  proposalDistr$baseline$concentration <- newM
  distrLock(proposalDistr$baseline)
  distrLock(proposalDistr)
  logMHRatio <- prPartition(proposalDistr, partition, log=TRUE) + logPriorDensity(newM) - prPartition(distr, partition, log=TRUE) - logPriorDensity(oldM)
  if ( log(runif(1)) < logMHRatio ) list(distr=proposalDistr, concentration=newM, accepted=TRUE)
  else list(distr=distr, concentration=oldM, accepted=FALSE)
}

