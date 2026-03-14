#' Probabilities for the Ewens-Pitman Attaction (EPA) Partition Distribution
#'
#' This function specifies the Ewens-Pitman attraction partition distribution
#' given similarity matrix, permutation, concentration, and
#' discount parameters.
#'
#' @inheritParams CRPPartition
#' @inheritParams ShrinkagePartition
#' @param similarity An n-by-n symmetric matrix giving the similarity between
#'   pairs of items.  Matrix entries must be strictly positive.
#'
#' @return An object of class \code{PartitionDistribution} representing this
#'   partition distribution.
#'
#' @example man/examples/EPAPartition.R
#' @export
#'
EPAPartition <- function(similarity, permutation, concentration, discount=0) {
  similarity <- coerceSimilarity(similarity)
  nItems <- nrow(similarity)
  permutation <- coercePermutation(permutation, nItems)
  md <- coerceConcentrationDiscount(concentration, discount)
  distrEnv("EPAPartition", list(nItems=nItems, similarity=similarity, permutation=permutation, .permutation=permutation-1L, concentration=md$concentration, discount=md$discount))
}

#' @export
print.EPAPartition <- function(x, ...) {
  cat("\nEPA partition distribution\n\n")
  print(distrAsList(x))
}
