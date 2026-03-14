samplePartition <- function(distr, nSamples, randomizePermutation=FALSE) {
  UseMethod("samplePartition")
}

#' @export
samplePartition.default <- function(distr, nSamples, randomizePermutation=FALSE) {
  stop("Not yet implemented.")
}
