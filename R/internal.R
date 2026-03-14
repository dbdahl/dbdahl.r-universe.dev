mk_splinter_ptrs <- function(splinters) {
  if ( ! is.list(splinters) ) stop("'splinters' should be a list.")
  if ( length(splinters) == 0 ) stop("'splinters' should be a non-empty list.")
  if ( ! all(sapply(splinters, \(x) is.list(x))) ) {
    stop("Each element of 'splinters' must be a list.")
  }
  if ( ! all(sapply(splinters, \(x) identical(names(x),c("indices","partitions")))) ) {
    stop("Each element of 'splinters' must be a list with names c('indices','partitions').")
  }
  if ( ! all(sapply(splinters, \(x) is.integer(x$indices))) ) {
    stop("'indices' must be a vector of integers'.")
  }
  if ( ! all(sapply(splinters, \(x) (is.integer(x$partitions) || is.double(x$partitions)) && is.matrix(x$partitions))) ) {
    stop("'partitions' must be a matrix of integers'.")
  }
  if ( ! all(sapply(splinters, \(x) length(x$indices) > 0)) ) {
    stop("The length of 'indices' not be greater than zero.")
  }
  indices <- Reduce(union, sapply(splinters, \(x) x$indices))
  if ( min(indices) != 1 ) stop(sprintf("The minimum index should be 1, but was found to be %s.", min(indices)))
  nItems <- max(indices)
  len1 <- length(indices)
  if ( len1 != nItems ) stop(sprintf("Not all indices are covered; only %s indices found but %s were expected.", len1, nItems))
  ptrs <- lapply(splinters, \(splinter) {
    partitions <- if ( is.character(splinter$partitions) ) {
      readRDS(splinter$partitions)$samples
    } else {
      splinter$partitions
    }
    if ( length(splinter$indices) != nrow(partitions) ) {
      stop("The number of rows in 'partitions' must equal the length of 'indices'.")
    }
    if ( ncol(partitions) == 0 ) {
      stop("The number of columns in 'partitions' must be greater than zero.")
    }
    .Call(.new_splinter, nItems, splinter$indices, partitions)
  })
  list(nItems = nItems, ptrs = ptrs)
}
