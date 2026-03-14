densityEstimate <- function(partition, parameters, distr, sampleParameter, density, grid) {
  if ( ! inherits(distr, "CRPPartition") ) stop("Only 'CRPPartition' is currently supported.")
  m <- distr$concentration
  d <- distr$discount
  g <- sapply(seq_len(nrow(partition)), function(i) {
    p <- partition[i,]
    t <- table(p)
    nc <- length(t)
    w <- c(t - d, m + nc*d)
    j <- sample(seq_along(w), 1, prob=w)
    s <- if ( j > nc ) sampleParameter() else parameters[[i]][[j]]
    density(grid, s)
  })
  list(x=grid, y=apply(g, 1, mean))
}
