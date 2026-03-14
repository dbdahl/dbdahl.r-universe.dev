#' Probabilities for the Pegged Time Dependent Random Partition Distribution
#'
#' This function specifies the pegged time dependent random partition
#' distribution for given concentration, discount, and pegging
#' probability parameters and the given number of time slices.
#'
#' @inheritParams CRPPartition
#' @inheritParams ShrinkagePartition
#' @param nSlices The number of time slices to consider in {1, 2, ...}.
#' @param peggedProbability The probability that an item is not reallocated
#'   between time points.
#'
#' @return An object of class \code{PartitionDistribution} representing this
#'   partition distribution.
#'
#' @example man/examples/PeggedTimeDependentPartition.R
#' @noRd
#'
PeggedTimeDependentPartition <- function(nItems, concentration, discount, nSlices, peggedProbability, anchor=NULL) {
  nItems <- coerceNItems(nItems)
  md <- coerceConcentrationDiscount(concentration, discount)
  nSlices <- coerceInteger(nSlices)
  if ( nSlices < 1 ) stop("'nSlices' must be at least one.")
  peggedProbability <- coerceDouble(peggedProbability)
  if ( peggedProbability < 0 || peggedProbability > 1 ) stop("'peggedProbability' must be in [0,1].")
  if ( ! is.null(anchor) ) {
    anchor <- coerceAnchor(anchor)
    if ( length(anchor) != nItems ) stop("'anchor' implied an incompatible number of items.")
  }
  distrEnv("PeggedTimeDependentPartition", list(nItems=nItems, concentration=md$concentration, discount=md$discount, nSlices=nSlices, peggedProbability=peggedProbability, anchor=anchor))
}

#' @noRd
print.PeggedTimeDependentPartition <- function(x, ...) {
  cat("\nPegged time dependent partition distribution\n\n")
  print(distrAsList(x))
}

#' @importFrom stats rbinom
#' @noRd
#'
samplePartition.PeggedTimeDependentPartition <- function(distr, nSamples, randomizePermutation=FALSE, randomizeShrinkage=TRUE, max=4, shape1=1.0, shape2=1.0, nCores=0) {
  x <- array(0L, dim=c(nSamples, distr$nItems, distr$nSlices, 2))
  x[,,-1,1] <- rbinom(nSamples*distr$nItems*(distr$nSlices-1), 1, distr$peggedProbability)
  x[,,1,2] <- if ( is.null(distr$anchor) ) {
    samplePartition(distr=CRPPartition(distr$nItems, distr$concentration, distr$discount), nSamples=nSamples, randomizePermutation=randomizePermutation, max=max, shape1=shape1, shape2=shape2, nCores=nCores)
  } else {
    rep(distr$anchor, each=nSamples)
  }
  for ( t in seq_len(distr$nSlices-1)+1 ) {
    x[,,t,2] <- x[,,t-1,2,drop=FALSE]
    x[,,t,2][x[,,t,1,drop=FALSE]!=1] <- NA
    x[,,t,2] <- t(apply(x[,,t,2,drop=FALSE], 1, function(y) {
      tab <- if ( all(is.na(y)) ) 0
      else table(factor(y, levels=seq_len(max(y,na.rm=TRUE))))
      for ( i in which(is.na(y)) ) {
        y[i] <- sample(c(seq_along(tab),length(tab)+1), 1, prob=c(tab-ifelse(tab==0,0,distr$discount), distr$concentration+distr$discount*sum(tab!=0)))
        if ( y[i] > length(tab) ) tab[y[i]] <- 0
        tab[y[i]] <- tab[y[i]] + 1
      }
      canonicalForm(y)
    }))
  }
  x
}

# expectedRI <- function(x, y) {
#   s <- 0
#   for ( i in seq_len(nrow(x)) ) {
#     s <- s + salso::RI(x[i,],y[i,])
#   }
#   s/nrow(x)
# }

probPaired <- function(x, y, togetherAtTime1=TRUE, i=1, j=2) {
  w <- if ( is.null(togetherAtTime1) ) rep(TRUE,nrow(x))
  else if ( togetherAtTime1 ) x[,i] == x[,j]
  else if ( ! togetherAtTime1 ) x[,i] != x[,j]
  else stop("Error")
  mean(y[w,i] == y[w,j])
}
