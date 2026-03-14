DDCRPPartition <- function(similarity, mass) {
  checkSimilarity(similarity)
  nItems <- nrow(similarity)
  if ( nItems < 1 ) stop("The number of rows in 'similarity' must be at least one.")
  checkMassDiscount(mass, 0.0)
  result <- list(nItems=nItems, similarity=similarity, mass=mass)
  class(result) <- c("DDCRPPartition", "partitionDistribution")
  result
}

checkSimilarity <- function(similarity) {
  if ( ! is.matrix(similarity) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
  if ( ! isSymmetric(similarity) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
  if ( any( similarity <= 0 ) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
}

checkMassDiscount <- function(mass, discount) {
  if ( ( discount < 0.0 ) || ( discount >= 1 ) ) stop("'discount' must be in [0,1).")
  if ( mass <= -discount ) stop("'mass' must be greater than -'discount'.")
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
samplePartition.DDCRPPartition <- function(distr, nSamples, randomizePermutation=FALSE) {
  m <- distr$similarity
  diag(m) <- distr$mass
  t(sapply(seq_len(nSamples), function(k) edgesToLabels(sapply(seq_len(distr$nItems), function(i) sample(distr$nItems, 1, prob=m[i,])))))
}
