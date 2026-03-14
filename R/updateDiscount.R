#' Update the Discount Parameter of a Partition Distribution
#'
#' This function updates the discount parameter of a partition distribution (e.g,
#' \code{\link{CRPPartition}} and \code{\link{ShrinkagePartition}}) using a
#' Gaussian random walk.  The prior can be specified by the user.
#'
#' @inheritParams updateConcentration
#'
#' @return A list containing the potentially updated partition distribution and
#'   a logical named \code{accepted} (indicating whether the proposal was
#'   accepted).
#'
#' @seealso \code{\link{nealAlgorithm3}}, \code{\link{nealAlgorithm8}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @importFrom stats dunif runif rnorm
#' @export
#'
updateDiscount <- function(distr, partition, rwsd=0.1,
                           logPriorDensity=function(w) dunif(w, 0, 1, log=TRUE)) {
  UseMethod("updateDiscount")
}

#' @export
#'
updateDiscount.default <- function(distr, partition, rwsd=0.5,
                                   logPriorDensity=function(w) dunif(w, 0, 15, log=TRUE)) {
  if ( is.null(distr[['discount']]) || rwsd <= 0 || inherits(distr,"CenteredPartition") ) {
    return(list(distr=distr, accepted=FALSE))
  }
  oldD <- distr[['discount']]
  newD <- rnorm(1, mean=oldD, sd=rwsd)
  if ( distr[['concentration']] <= -newD || newD < 0.0 || newD >= 1.0 ) return(list(distr=distr, accepted=FALSE))
  proposalDistr <- distrClone(distr)
  proposalDistr$discount <- newD
  distrLock(proposalDistr)
  logMHRatio <- prPartition(proposalDistr, partition, log=TRUE) + logPriorDensity(newD) - prPartition(distr, partition, log=TRUE) - logPriorDensity(oldD)
  if ( log(runif(1)) < logMHRatio ) list(distr=proposalDistr, accepted=TRUE)
  else list(distr=distr, accepted=FALSE)
}

#' @export
#'
updateDiscount.ShrinkagePartition <- function(distr, partition, rwsd=0.5,
                                       logPriorDensity=function(w) dunif(w, 0, 15, log=TRUE)) {
  if ( is.null(distr$baseline[['discount']]) || rwsd <= 0 || inherits(distr$baseline,"CenteredPartition") ) {
    return(list(distr=distr, accepted=FALSE))
  }
  oldD <- distr$baseline[['discount']]
  newD <- rnorm(1, mean=oldD, sd=rwsd)
  if ( distr$baseline[['concentration']] <= -newD || newD < 0.0 || newD >= 1.0 ) return(list(distr=distr, accepted=FALSE))
  proposalDistr <- distrClone(distr)
  proposalDistr$baseline <- distrClone(proposalDistr$baseline)
  proposalDistr$baseline$discount <- newD
  distrLock(proposalDistr$baseline)
  distrLock(proposalDistr)
  logMHRatio <- prPartition(proposalDistr, partition, log=TRUE) + logPriorDensity(newD) - prPartition(distr, partition, log=TRUE) - logPriorDensity(oldD)
  if ( log(runif(1)) < logMHRatio ) list(distr=proposalDistr, accepted=TRUE)
  else list(distr=distr, accepted=FALSE)
}

