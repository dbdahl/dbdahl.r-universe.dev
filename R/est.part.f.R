#' Fernando's Splintered Clustering Idea
#'
#' @param post.mean Vector containing estimated posterior means for those subjects under consideration.
#'
#' @return A list giving the number of clusters, cluster labels, and cluster sizes using [mclust::Mclust()].
#' @export
#' @importFrom mclust Mclust
#'
est.part.f <- function(post.mean) {
  out <- mclust::Mclust(post.mean)
  maxclust <- out$G
  clusters <- out$classification
  clust.sizes <- table(clusters)
  output <- list("num.clusters"=maxclust, "clusters"=clusters, "clust.sizes"=clust.sizes)
  return(output)
}
