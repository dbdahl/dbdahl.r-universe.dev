#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances using the CaviarPD method.
#'
#' @param distance An object of class 'dist' or a pairwise distance matrix.
#' @param nClusters A numeric vector that specifies the range for the number of clusters to consider in the search for a clustering estimate.
#' @param mass The mass value to use for sampling. If \code{NULL}, the mass value is found by inverting values from \code{nClusters}.
#' @param nSamples The number of samples drawn per candidate estimate.
#' @param gridLength The number of candidate estimates to consider. The final estimate is obtained from \code{nSamples} \eqn{\times} \code{gridLength} total samples.
#' @param loss The SALSO method (Dahl, Johnson, Müller, 2021) tries to minimize this expected loss when searching the partition space for an optimal estimate. This must be either "binder" or "VI".
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param similarity Either \code{"exponential"} or \code{"reciprocal"} to indicate the desired similarity function.
#' @param maxNClusters The maximum number of clusters that can be considered by the SALSO method.
#' @param nRuns The number of runs of the SALSO algorithm.
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @details
#' A range for the number of clusters to be considered is supplied using the
#' \code{nClusters} argument.
#'
#' @return A object of class \code{salso.estimate}, which provides a clustering estimate (a vector of cluster labels) that can be displayed and plotted.
#'
#' @references
#'
#' D. B. Dahl, J. Andros, J. B. Carter (2023), Cluster Analysis via Random Partition
#' Distributions, \emph{Statistical Analysis and Data Mining}, \doi{10.1002/sam.11602}.
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2022), Search Algorithms and Loss
#' Functions for Bayesian Clustering, \emph{Journal of Computational and
#' Graphical Statistics}, 31(4), 1189-1201, \doi{10.1080/10618600.2022.2069779}.
#''
#' @examples
#' # To reduce load on CRAN servers, limit the number of samples, grid length, and CPU cores.
#' set.seed(34)
#' iris.dis <- dist(iris[,-5])
#' est <- caviarpd(distance=iris.dis, nClusters=c(2,4), nSamples=20, nCores=1)
#' if ( require("salso") ) {
#'   summ <- summary(est, orderingMethod=2)
#'   plot(summ, type="heatmap")
#'   plot(summ, type="mds")
#' }
#'
#' @importFrom stats median
#' @export
#'
caviarpd <- function(distance, nClusters, mass=NULL, nSamples=200, gridLength=5,
                     loss="binder", temperature=100, similarity=c("exponential","reciprocal")[1],
                     maxNClusters=0, nRuns=4, nCores=nRuns) {
  if ( is.matrix(distance) ) {
    if ( !isSymmetric(distance) || !is.numeric(distance) ) stop("'distance' is not a symmetric numerical matrix.")
  } else if ( inherits(distance,'dist') ) {
    distance <- as.matrix(distance)
  } else stop("'distance' argument must be an object of class 'dist' or a symmetric numerical matrix.")
  if ( !is.numeric(nClusters) || !all(is.finite(nClusters)) || any(nClusters<1) ) stop("'nClusters' must a numeric vector of finite values not less than 1")
  if ( !is.null(mass) && ( !is.numeric(mass) || !all(is.finite(mass)) || any(mass<=0.0) ) ) stop("'mass', if non-null, must be a numeric vector of finite values greater than 0")
  if ( !is.numeric(nSamples) || ! length(nSamples) %in% c(1,2) || any(nSamples <= 0) || any(nSamples %% 1 != 0) ) stop("'nSamples' must be a strictly positive and length 1 or 2")
  if ( !is.numeric(gridLength) || length(gridLength) != 1 || gridLength < 2 || gridLength %% 1 != 0 ) stop("'gridLength' must be a strictly positive integer not less than 2")
  if ( !is.character(loss) || length(loss) != 1 || ! loss %in% c("binder","VI") ) stop("'loss' must be either 'binder' or 'VI'")
  if ( !is.numeric(temperature) || !is.vector(temperature) || length(temperature) != 1 || temperature < 0 ) stop("'temperature' must be nonnegative and length 1")
  if ( !is.character(similarity) || length(similarity) != 1 || ! similarity %in% c("exponential","reciprocal") ) stop("'similarity' must be either 'exponential' or 'reciprocal'")
  if ( !is.numeric(maxNClusters) || length(maxNClusters) != 1 || maxNClusters < 0 || maxNClusters %% 1 != 0 ) stop("'maxNClusters' must be 0 or a positive integer")
  if ( maxNClusters == 0 ) maxNClusters <- max(nClusters) + 1
  if ( !is.numeric(nRuns) || length(nRuns) != 1 || nRuns < 1 || nRuns %% 1 != 0 ) stop("'nRuns' must be a strictly positive integer")
  if ( !is.numeric(nCores) || length(nCores) != 1 || nCores < 0 || nCores %% 1 != 0 ) stop("'nCores' must be 0 or a positive integer")
  distance <- distance / median(as.vector(distance))
  similarity <- if ( similarity == "exponential" ) {
    exp( -temperature * distance )
  } else if ( similarity == "reciprocal" ) {
    if ( any(distance == 0.0 ) ) distance <- distance + 0.01
    1/distance^temperature
  } else stop("Unsupported similarity")
  if ( ! all(is.finite(similarity)) ) stop("'distance', 'temperature', and/or 'similarity' yield similarity with nonfinite values")
  result <- .Call(.caviarpd_algorithm2, similarity, min(nClusters), max(nClusters), mass, nSamples, gridLength, getOption("caviarpd.n0",100), getOption("caviarpd.tol",0.01), loss=="VI", maxNClusters, nRuns, nCores)
  structure(result$estimate, class="salso.estimate", draws=result$samples, info=list(loss=loss))
}

mass <- function(expected_number_of_clusters, n_items) {
  .Call(.caviarpd_mass, expected_number_of_clusters, n_items)
}

sampleEPA <- function(similarity, mass, nSamples=500, nCores=0) {
  .Call(.sample_epa, nSamples, similarity, mass, nCores)
}
