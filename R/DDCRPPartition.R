#' Probabilities for the Distance-Dependent CRP (ddCRP) Partition Distribution
#'
#' This function specifies the Distance-Dependent CRP (ddCRP) partition
#' distribution for given similarity matrix and concentration.
#'
#' @inheritParams EPAPartition
#' @param similarity An n-by-n symmetric matrix giving the similarity between
#'   pairs of items.  Matrix entries must be strictly positive.
#'
#' @return An object of class \code{PartitionDistribution} representing this
#'   partition distribution.
#'
#' @example man/examples/DDCRPPartition.R
#' @export
#'
DDCRPPartition <- function(similarity, concentration) {
  similarity <- coerceSimilarity(similarity)
  md <- coerceConcentrationDiscount(concentration, 0.0)
  distrEnv("DDCRPPartition", list(nItems=nrow(similarity), similarity=similarity, concentration=md$concentration))
}

#' @export
print.DDCRPPartition <- function(x, ...) {
  cat("\nDDCRP partition distribution\n\n")
  print(distrAsList(x))
}

edgesToLabels <- function(edges) {
  n <- length(edges)
  labels <- integer(n)
  visited <- logical(n)
  nextAvailableLabel <- 0
  findLabelFor <- function(i) {
    if ( labels[i] == 0 ) {
      labels[i] <<- if ( visited[i] ) {
        nextAvailableLabel <<- nextAvailableLabel + 1
        nextAvailableLabel
      } else {
        visited[i] <<- TRUE
        findLabelFor(edges[i])
      }
    }
    labels[i]
  }
  sapply(seq_along(edges), findLabelFor)
}

#' @export
#'
samplePartition.DDCRPPartition <- function(distr, nSamples, randomizePermutation=FALSE, randomizeShrinkage=TRUE, max=4, shape1=1.0, shape2=1.0, nCores=0) {
  m <- distr$similarity
  diag(m) <- distr$concentration
  t(sapply(seq_len(nSamples), function(k) edgesToLabels(sapply(seq_len(distr$nItems), function(i) sample(distr$nItems, 1, prob=m[i,])))))
}
