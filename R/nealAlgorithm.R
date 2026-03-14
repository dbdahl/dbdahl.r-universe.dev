#' Neal's Algorithm 3 for a Partition
#'
#' This function performs Algorithm 3 of Neal (2000), which updates a partition
#' based on a partition prior distribution and a user-supplied function giving
#' an item's contribution to the log integrated likelihood function (in which
#' the cluster-specific model parameters have been integrated out).
#'
#' \code{logIntegratedLikelihoodItem} is a function giving the contribution of
#' item \eqn{i} to the log integrate likelihood, defined as \deqn{ p( y_i |
#' y_s ) = \int p( y_i | \theta ) p( \theta | y_s ) d \theta}. This is
#' typically available for conditionally conjugate models.
#'
#' @inheritParams samplePartition
#' @param partition An integer vector of cluster labels representing the current
#'   partition. Items \eqn{i} and \eqn{j} are in the same cluster if and only if
#'   \code{partition[i] == partition[j]}.
#' @param logIntegratedLikelihoodItem A function taking an index \eqn{i} (as a
#'   single integer) and a subset of integers \eqn{s} (as a numeric vector) and
#'   returning the natural logarithm of \eqn{p( y_i | y_s )}, i.e., that item's
#'   contribution to the log integrated likelihood given the subset \eqn{s} of
#'   integers.
#' @param mcmcTuning An optional list containing \code{nUpdates} (an integer
#'   giving the number of Gibbs scans before returning, defaulting to \eqn{1})
#'   and other variables specific to the particular prior partition
#'   distribution.
#'
#' @return A list containing \code{partition} (an integer vector giving the
#'   updated partition in cluster label notation) and possibly other items
#'   specific to a particular partition distribution.
#'
#' @seealso \code{\link{CRPPartition}}, \code{\link{ShrinkagePartition}},
#'   \code{\link{LocationScalePartition}}, \code{\link{CenteredPartition}},
#'   \code{\link{prPartition}}, \code{\link{samplePartition}},
#'   \code{\link{updatePermutation}} \code{\link{nealAlgorithm8}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @export
#'
nealAlgorithm3 <- function(partition, distr, logIntegratedLikelihoodItem=function(i, subset) 0, mcmcTuning=list()) {
  result <- nealAlgorithEngine(partition, logIntegratedLikelihoodItem, TRUE, distr, mcmcTuning)
  list(partition=result$partition)
}

#' Neal's Algorithm 8 for a Partition
#'
#' This function performs Algorithm 8 of Neal (2000), which updates a partition
#' based on a partition prior distribution and user-supplied functions to
#' evaluate an item's contribution to the loglikelihood sample a model parameter
#' for its prior distribution, and update a model parameter based on a subset.
#'
#' @inheritParams nealAlgorithm3
#' @param parameters A list whose \eqn{j}th element is the shared model
#'   parameter for all those items having cluster label \eqn{j} in the
#'   \code{partition} argument.
#' @param logLikelihood A function taking an index \eqn{i} (as a single
#'   integer), a cluster \eqn{label} (as a single integer), and an
#'   indication as to whether a new parameter should be sampled (as an single
#'   logical) and returning the natural logarithm of \eqn{p( y_i | parameter
#'   associated with "label" )}, i.e., that item's contribution to the log
#'   likelihood given the parameter associated with the cluster \dQuote{label}.
#' @param sampleParameter A function taking no arguments and returning a draw
#'   from the prior parameter distribution (i.e., the centering distribution).
#' @param updateParameter A function taking a \code{parameter} (a current value
#'   of a parameter) and a subset of integers \code{subset} (as a numeric
#'   vector) and returning an updated value for the parameter by some valid MCMC
#'   update scheme for the full conditional distribution of the parameter.
#'
#' @return A list containing \code{partition} (an integer vector giving the
#'   updated partition encoded using cluster labels), \code{parameters} (a list
#'   like the \code{parameters} argument, with updated values), and possibly
#'   other items specific to a particular partition distribution.
#'
#' @seealso \code{\link{CRPPartition}}, \code{\link{ShrinkagePartition}},
#'   \code{\link{LocationScalePartition}}, \code{\link{CenteredPartition}},
#'   \code{\link{prPartition}}, \code{\link{samplePartition}},
#'   \code{\link{updatePermutation}}, \code{\link{nealAlgorithm3}}
#'
#' @example man/examples/nealAlgorithm.R
#'
#' @export
#'
nealAlgorithm8 <- function(partition, parameters, distr, logLikelihood=function(i, parameter) 0, sampleParameter=function() 0, updateParameter=function(parameter, subset) 0, mcmcTuning=list()) {
  logLike <- function(i, label, isNew) {
    if ( isNew ) parameters[[label]] <<- sampleParameter()
    logLikelihood(i, parameters[[label]])
  }
  result <- nealAlgorithEngine(partition, logLike, FALSE, distr, mcmcTuning)
  # Standardize the order of the parameters
  parameters <- parameters[result$map]
  # Update parameters
  sapply(unique(result$partition), function(label) {
    subset <- which(result$partition==label)
    parameters[[label]] <<- if ( length(subset) == 0 ) NA
    else updateParameter(parameters[[label]], subset)
  })
  list(partition=result$partition, parameters=parameters)
}

nealAlgorithEngine <- function(partition, logLikelihood, isAlgorithm3, distr, mcmcTuning) {
  nUpdatesForPartition <- coerceInteger(getOr(mcmcTuning$nUpdatesForPartition, 1L))
  if ( nUpdatesForPartition < 1 ) stop("The number of updates should be at least one.")
  .Call(.nealAlgorithm, partition, logLikelihood, isAlgorithm3, nUpdatesForPartition, mkDistrPtr(distr))
}
